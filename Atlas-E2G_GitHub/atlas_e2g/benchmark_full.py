"""Full-input fitting and coordinate-consistent evaluation."""
from pathlib import Path
import time
import numpy as np,pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score as ap
from .features import DATA,FACTORS,load,native
from .model import splits,fit_logistic,score
from .evaluation import partial,valid_coordinates

METHODS=['full','matched','published','without_neighborhood','without_identity','shared','without_profiles','without_ratios']
REPORTER_METHODS=['full','without_neighborhood','without_identity','shared','without_ratios']


def design(raw,b,method,sc,cs):
    if method=='published':return raw
    x=sc.transform(raw)
    if method=='matched':return x
    a=b['pool'];nf=len(FACTORS);n=cs.transform(b['neighborhood'])
    if method=='without_neighborhood':return np.c_[x,a]
    if method=='without_identity':a=np.c_[np.sort(a[:,:nf],axis=1),np.sort(a[:,nf:],axis=1)]
    elif method=='shared':a=(a[:,:nf]+a[:,nf:])/np.sqrt(2)
    elif method=='without_profiles':a=a[:,:0]
    elif method=='without_ratios':n=n[:,2:]
    return np.c_[x,a,n]

def refit(output,methods=None):
    output=Path(output);methods=methods or METHODS
    if set(methods)-set(METHODS):raise ValueError('Unknown methods: '+str(set(methods)-set(METHODS)))
    d,b=load('source');q,qb=load('reporter');raw=native(d);mask=d.mapped.to_numpy();y=d.truth.to_numpy();nf=len(FACTORS)
    # Full-input fitting uses 9604 covered source pairs; valid coordinates are
    # required at evaluation.
    folds=list(splits(d))
    scores={};fitted={};start=time.time()
    for pr,fold,test,train in folds+[('fit_all','all',mask,mask)]:
        sc=StandardScaler().fit(raw[train]);cs=StandardScaler().fit(b['neighborhood'][train])
        for method in methods:
            x=design(raw,b,method,sc,cs);m=fit_logistic(x[train],y[train],method=='published')
            if pr=='fit_all':
                fitted[method]=(m,sc,cs)
                np.savez_compressed(output/(method+'.npz'),**m,native_mean=sc.mean_,native_scale=sc.scale_,context_mean=cs.mean_,context_scale=cs.scale_)
            else:scores.setdefault(pr+'/'+method,np.full(len(d),np.nan))[test]=score(m,x[test])
            if pr=='chromosome' and method in REPORTER_METHODS:
                hit=q.chrom.eq(fold).to_numpy();profile=qb['profile'][hit]
                if method=='shared':coef=m['coef'][8:8+nf]/np.sqrt(2)
                else:coef=m['coef'][8+nf:8+2*nf]
                if method=='without_identity':profile=np.sort(profile,axis=1)
                scores.setdefault('reporter/'+method,np.full(len(q),np.nan))[hit]=profile@coef
        print('Full-input fold',pr,fold,len(methods),'fits',round(time.time()-start,1),'s',flush=True)
    for kind in ['fivecell','astrocyte']:
        td,tb=load(kind)
        for second in ([False,True] if kind=='astrocyte' else [False]):
            label=kind+('_second' if second else '');rr=native(td,second)
            for method,(m,sc,cs) in fitted.items():
                if method=='published':p=td['published_second' if second else 'published'].to_numpy()
                else:p=score(m,design(rr,tb,method,sc,cs))
                scores[label+'/'+method]=p
    np.savez_compressed(output/'scores.npz',**scores);return scores

def summarize(scores,output,bootstrap=0):
    output=Path(output);rows=[]
    def add(dataset,subset,method,metric,value,n):rows.append(dict(dataset=dataset,subset=subset,method=method,metric=metric,value=float(value),n=int(n)))
    for kind,protocol in [('source','chromosome'),('source','study_gene_element_disjoint'),('fivecell','fivecell'),('astrocyte','astrocyte'),('astrocyte','astrocyte_second')]:
        d=pd.read_csv(DATA/f'{kind}.csv.gz');v=valid_coordinates(d);mapped=d['mapped_second' if protocol.endswith('second') else 'mapped'].to_numpy();covered=v&mapped
        subsets={'covered':covered}
        if kind!='astrocyte':subsets['all']=v
        if kind=='fivecell':
            subsets.update(non_K562=covered&d.cell.ne('K562').to_numpy(),non_K562_distal=covered&d.cell.ne('K562').to_numpy()&d.benchmark_distance.gt(100000).to_numpy())
            subsets.update({c:covered&d.cell.eq(c).to_numpy() for c in sorted(d.cell.unique())})
        for method in METHODS:
            key=protocol+'/'+method
            if key not in scores:continue
            p=np.where(mapped,scores[key],0.)
            if not np.isfinite(p[covered]).all():continue
            for subset,ix in subsets.items():
                add(protocol,subset,method,'AP',ap(d.truth[ix],p[ix],sample_weight=d.weight[ix] if kind=='fivecell' else None),ix.sum())
            if not protocol.endswith('second'):
                ii=np.flatnonzero(v if kind!='astrocyte' else covered);selected=ii[np.argsort(-p[ii],kind='stable')[:100]]
                add(protocol,'all' if kind!='astrocyte' else 'covered',method,'top100',d.truth.iloc[selected].sum(),len(selected))
            if kind in ['fivecell','astrocyte'] and not protocol.endswith('second'):
                for subset,ix in subsets.items():
                    if subset not in ['covered','non_K562','non_K562_distal']:continue
                    hit=slots=0
                    for _,g in d[ix].groupby(['cell','study','gene'],sort=True):
                        if len(g)<4:continue
                        ii=g.index.to_numpy();chosen=ii[np.argsort(-p[ii],kind='stable')[:3]];hit+=int(d.truth.iloc[chosen].sum());slots+=3
                    add(protocol,subset,method,'top3_positive',hit,slots)
    q,qb=load('reporter')
    for method in REPORTER_METHODS:
        key='reporter/'+method
        if key not in scores or not np.isfinite(scores[key]).all():continue
        p=scores[key]
        for subset,ix in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
            add('reporter',subset,method,'Spearman',spearmanr(q.response[ix],p[ix]).statistic,ix.sum())
            add('reporter',subset,method,'adjusted_Spearman',partial(q[ix],p[ix],qb['amount'][ix]),ix.sum())
    if bootstrap and 'reporter/full' in scores and np.isfinite(scores['reporter/full']).all():
        ci=[]
        for subset,ix in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
            groups=[g.index.to_numpy() for _,g in q[ix].groupby('gene',sort=True)];rng=np.random.default_rng(20261030);values=[]
            for _ in range(bootstrap):
                ii=np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))]);values.append(spearmanr(q.response.iloc[ii],scores['reporter/full'][ii]).statistic)
            lo,hi=np.quantile(values,[.025,.975]);ci.append(dict(subset=subset,low=lo,high=hi,resamples=bootstrap))
        pd.DataFrame(ci).to_csv(output/'reporter_intervals.csv',index=False)
    result=pd.DataFrame(rows);result.to_csv(output/'results.csv',index=False);return result

def compare_results(actual,expected):
    keys=['dataset','subset','method','metric'];m=actual.merge(expected,on=keys,suffixes=('_actual','_expected'),validate='one_to_one')
    if len(m)!=len(actual):raise AssertionError('Unrecognized full-input result')
    delta=abs(m.value_actual-m.value_expected);limits=np.where(m.metric.isin(['top100','top3_positive']),0,2e-5)
    if (delta>limits).any() or (m.n_actual!=m.n_expected).any():raise AssertionError('Full-input reproduction mismatch:\n'+m[(delta>limits)|(m.n_actual!=m.n_expected)].to_string(index=False))
    return dict(checked=len(m),max_metric_difference=float(delta.max()),status='PASS')
