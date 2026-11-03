"""Portable reproduction of public-input, component-transfer and query analyses."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score as ap
from .features import ROOT,DATA,FACTORS,load,native,public_design,profiles,normalize
from .evaluation import valid_coordinates
from .model import splits
from .rebuild_public import optional_identity,fit_views as public_fit
from .evaluation import partial

NF=len(FACTORS)

def data(kind):
    d,b=load(kind);z=np.load(DATA/'public_features.npz',allow_pickle=False)
    _,g=public_design(b)
    c=optional_identity(d,sparse.load_npz(DATA/'source_overlap.npz'),g.astype(float)) if kind=='source' else g
    return dict(d=d,b=b,lo=z[kind+'/lo'],g=g,c=c,activity=z[kind+'/activity'],available=z[kind+'/available'])

def fit_offset(views,y,offset):
    def objective(w):
        loss=5*np.dot(w[:-1],w[:-1]);grad=10*w[:-1];intercept=0.
        for x in views:
            v=x@w[:-1]+w[-1]+offset;r=(expit(v)-y)/len(views)
            loss+=np.sum(np.logaddexp(0,v)-y*v)/len(views)
            grad=grad+x.T@r;intercept+=r.sum()
        return loss,np.r_[grad,intercept]
    fit=minimize(objective,np.zeros(views[0].shape[1]+1),jac=True,method='L-BFGS-B',options=dict(maxiter=1500,ftol=1e-9))
    if not fit.success:raise RuntimeError(fit.message)
    return dict(coef=fit.x[:-1],intercept=float(fit.x[-1]))

def teacher(v,train):
    raw=np.c_[native(v['d']),v['b']['neighborhood']];sc=StandardScaler().fit(raw[train]);x=np.c_[sc.transform(raw),v['b']['pool']]
    m=fit_offset([x[train]],v['d'].truth.to_numpy()[train],np.zeros(train.sum()))
    m.update(mean=sc.mean_,scale=sc.scale_)
    return m,x

def design(v,m,view,method,profile='g'):
    low=v['lo'] if view=='gene' else np.c_[v['lo'],v['activity']]
    g=v[profile]
    if method=='anchored_promoter':g=g[:,:NF]
    elif method in ['independent_pooled_promoter','pooled_distilled']:g=np.c_[g[:,:NF],v['b']['pool'][:,NF:]]
    return np.c_[(low-m['mean'])/m['scale'],g]

def fit_components(v,train,methods):
    t,x=teacher(v,train);y=v['d'].truth.to_numpy();soft=.5*y+.5*expit(x@t['coef']+t['intercept']);models={}
    for view in ['gene','cell']:
        low=v['lo'] if view=='gene' else np.c_[v['lo'],v['activity']];sc=StandardScaler().fit(low[train]);scale=dict(mean=sc.mean_,scale=sc.scale_)
        for method in methods:
            pc=t['coef'][14+NF:14+2*NF] if method=='anchored_promoter' else np.zeros(NF)
            offset=v['b']['pool'][:,NF:]@pc
            xx=[design(v,scale,view,method,k)[train] for k in ['g','c']]
            m=fit_offset(xx,(soft if method=='pooled_distilled' else y)[train],offset[train]);m.update(**scale,promoter_coef=pc);models[view+'/'+method]=m
    return models,t

def predict(v,models):
    pred={}
    for key,m in models.items():
        view,method=key.split('/');x=design(v,m,view,method);offset=v['b']['pool'][:,NF:]@m['promoter_coef']
        pred[key]=expit(x@m['coef']+m['intercept']+offset)
    for method in {k.split('/')[1] for k in pred}:pred['cell/'+method][~v['available']]=pred['gene/'+method][~v['available']]
    return pred

def fit_transfer(out):
    v=data('source');d=v['d'];valid=d.mapped.to_numpy()&valid_coordinates(d)
    methods=['independent','independent_pooled_promoter','anchored_promoter','pooled_distilled'];scores={}
    models,t=fit_components(v,valid,methods)
    np.savez_compressed(out/'component_models.npz',**{k+'/'+n:a for k,m in models.items() for n,a in m.items()})
    np.savez_compressed(out/'public_components.npz',factors=np.asarray(FACTORS),**{view+'_'+n:a for view in ['gene','cell'] for n,a in models[view+'/anchored_promoter'].items()})
    for kind in ['fivecell','astrocyte']:
        pred=predict(data(kind),models)
        for key,value in pred.items():
            view,method=key.split('/');prefix='controls' if method=='pooled_distilled' else 'transfer'
            scores[prefix+'/'+kind+'/'+key]=value
        print('Component transfer:',kind,flush=True)
    q,qb=load('reporter');_,balanced=profiles(sparse.load_npz(DATA/'reporter_overlap.npz'));public_profile=normalize(balanced[q.uid.to_numpy()])
    for pr,chrom,_,train in splits(d):
        if pr!='chromosome':continue
        hit=q.chrom.eq(chrom).to_numpy()
        if not hit.any():continue
        models,t=fit_components(v,train&valid_coordinates(d),methods)
        values={'transfer/reporter/standard':qb['profile']@t['coef'][14+NF:14+2*NF]}
        for key,m in models.items():
            view,method=key.split('/');n=len(m['mean'])
            if method=='independent':
                values['public/reporter/'+view+'/promoter']=public_profile@m['coef'][n+NF:n+2*NF]
                values['public/reporter/'+view+'/enhancer']=public_profile@m['coef'][n:n+NF]
            elif method in ['independent_pooled_promoter','pooled_distilled']:
                label='pooled_hard' if method=='independent_pooled_promoter' else method
                values['controls/reporter/'+view+'/'+label]=qb['profile']@m['coef'][n+NF:n+2*NF]
        for key,p in values.items():scores.setdefault(key,np.full(len(q),np.nan))[hit]=p[hit]
        print('Reporter chromosome:',chrom,flush=True)
    return scores

def fit_activity(out):
    v=data('source');d=v['d'];train=d.mapped.to_numpy()&valid_coordinates(d);y=d.truth.to_numpy();models={}
    choices={'baseline':[],'activity':[0,1],'activity_nearby':[0,1,3],'complete':[0,1,2,3,4]}
    for key,cols in choices.items():
        x=np.c_[v['lo'],v['activity'][:,cols]];sc=StandardScaler().fit(x[train]);x=sc.transform(x)
        g=v['g'];c=v['c']
        m=public_fit(np.c_[x[train],g[train]],np.c_[x[train],c[train]],y[train]);m.update(mean=sc.mean_,scale=sc.scale_);models[key]=m
    scores={};z=np.load(DATA/'public_features.npz')
    for kind in ['fivecell','astrocyte']:
        v=data(kind);pred={}
        for key,cols in choices.items():
            m=models[key];x=np.c_[v['lo'],v['activity'][:,cols]];g=v['g']
            pred[key]=expit(np.c_[(x-m['mean'])/m['scale'],g]@m['coef']+m['intercept'])
        m=models['complete'];x=np.c_[v['lo'],z[kind+'/k562_activity']]
        pred['K562_for_all_cells']=expit(np.c_[(x-m['mean'])/m['scale'],v['g']]@m['coef']+m['intercept'])
        for key,p in pred.items():
            if key!='baseline':p[~v['available']]=pred['baseline'][~v['available']]
            scores['activity/'+kind+'/'+key]=p
        print('Public activity ablation:',kind,flush=True)
    return scores

def fit_source(out):
    v=data('source');d=v['d'];q,qb=load('reporter');scores={}
    cases=[('removal',s) for s in sorted(d.study.unique())]+[('diagnostic',s) for s in ['Gasperini2019_only','count_matched_3101','count_matched_3102','count_matched_3103']]
    for scope,name in cases:
        for protocol,chrom,_,train in splits(d):
            if protocol!='chromosome':continue
            hit=q.chrom.eq(chrom).to_numpy()
            if not hit.any():continue
            train=train&valid_coordinates(d)
            if scope=='removal':train&=d.study.ne(name).to_numpy()
            elif name.endswith('_only'):train&=d.study.eq('Gasperini2019').to_numpy()
            else:
                n=int((train&d.study.ne('Gasperini2019').to_numpy()).sum());rng=np.random.default_rng(int(name.rsplit('_',1)[-1]));ix=rng.choice(np.flatnonzero(train),n,replace=False);train=np.zeros(len(d),bool);train[ix]=True
            m,_=teacher(v,train)
            for role,k in [('promoter',14+NF),('enhancer',14)]:
                if scope=='diagnostic' and role=='enhancer':continue
                key=scope+'/reporter/'+name+('/'+role if scope=='removal' else '')
                scores.setdefault(key,np.full(len(q),np.nan))[hit]=qb['profile'][hit]@m['coef'][k:k+NF]
        print('Source composition:',name,flush=True)
    return scores

def metrics(scores,out,bootstrap=0):
    q,qb=load('reporter');rows=[];rep={k:v for k,v in scores.items() if '/reporter/' in k}
    for subset,mask in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
        for key,p in rep.items():
            for metric,value in [('Spearman',spearmanr(p[mask],q.response[mask]).statistic),('adjusted_Spearman',partial(q[mask],p[mask],qb['amount'][mask]))]:rows.append(dict(dataset='reporter',subset=subset,method=key,metric=metric,value=float(value),n=int(mask.sum())))
    for kind in ['fivecell','astrocyte']:
        d=pd.read_csv(DATA/f'{kind}.csv.gz');valid=valid_coordinates(d)&d.mapped.to_numpy();subsets={'covered':valid}
        if kind=='fivecell':subsets.update(non_K562=valid&d.cell.ne('K562').to_numpy(),non_K562_distal=valid&d.cell.ne('K562').to_numpy()&d.benchmark_distance.gt(100000).to_numpy())
        subsets.update({cell:valid&d.cell.eq(cell).to_numpy() for cell in sorted(d.cell.unique())})
        for key,p in scores.items():
            if '/'+kind+'/' not in key:continue
            for subset,mask in subsets.items():rows.append(dict(dataset=kind,subset=subset,method=key,metric='AP',value=ap(d.truth[mask],p[mask],sample_weight=d.weight[mask] if kind=='fivecell' else None),n=int(mask.sum())))
    pd.DataFrame(rows).to_csv(out/'metrics.csv',index=False)
    if bootstrap:
        intervals=[]
        for subset,mask in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
            groups=[v.index.to_numpy() for _,v in q[mask].groupby('gene',sort=True)]
            for scope,seed in [('transfer',20261030),('controls',20261030),('public',20261030),('removal',20261031)]:
                chosen={k:v for k,v in rep.items() if k.startswith(scope+'/')};rng=np.random.default_rng(seed);values={k:[] for k in chosen}
                for _ in range(bootstrap):
                    ix=np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))])
                    for k,p in chosen.items():values[k].append([spearmanr(p[ix],q.response.iloc[ix]).statistic,partial(q.iloc[ix],p[ix],qb['amount'][ix])])
                for k,v in values.items():
                    for j,m in enumerate(['Spearman','adjusted_Spearman']):
                        lo,hi=np.quantile(np.asarray(v)[:,j],[.025,.975]);intervals.append(dict(subset=subset,method=k,metric=m,low=lo,high=hi))
        pd.DataFrame(intervals).to_csv(out/'reporter_intervals.csv',index=False)
    return pd.DataFrame(rows)

def check_metrics(result, fresh=False):
    expected=json.loads((ROOT/'reference_results/public_metrics.json').read_text());checks=[]
    def check(dataset,subset,key,metric,value):
        hit=result[(result.dataset==dataset)&(result.subset==subset)&(result.method==key)&(result.metric==metric)]
        if hit.empty:
            if fresh:return
            raise AssertionError('Missing metric: '+key)
        delta=abs(float(hit.iloc[0].value)-value)
        if not np.isfinite(delta) or delta>(2e-4 if fresh else 1e-12):raise AssertionError(f'Metric mismatch {key} {subset} {metric}: {delta}')
        checks.append(delta)
    for table,prefix in [('control_reporter','controls'),('public_reporter','public')]:
        for row in expected[table]:
            if row['method']=='anchored_promoter':key='transfer/reporter/standard'
            elif row['method']=='full/promoter':key='full/reporter/promoter'
            else:key=prefix+'/reporter/'+row['method']
            check('reporter',row['subset'],key,row['metric'],row['value'])
    for row in expected['source_removal']:
        key='removal/reporter/'+row['removed_collection']+'/'+row['role']
        if row['removed_collection']=='all_collections':key='transfer/reporter/standard'
        check('reporter',row['subset'],key,row['metric'],row['value'])
    for row in expected['source_diagnostic']:
        for metric in ['Spearman','adjusted_Spearman']:check('reporter',row['subset'],'diagnostic/reporter/'+row['configuration'],metric,row[metric])
    for table,prefix in [('transfer_external','transfer'),('control_external','controls')]:
        for row in expected[table]:
            if row['method']=='pooled_hard':key='transfer/'+row['dataset']+'/'+row['view']+'/independent_pooled_promoter'
            else:key=prefix+'/'+row['dataset']+'/'+row['view']+'/'+row['method']
            check(row['dataset'],row['subset'],key,'AP',row['AP'])
    for row in expected['activity']:
        check(row['dataset'],row['subset'],'activity/'+row['dataset']+'/'+row['method'],'AP',row['AP'])
    if not checks:raise AssertionError('No expected metrics checked')
    return len(checks),max(checks)


def comparisons(out):
    full=np.load(ROOT/'reference_results/full_scores.npz');public=np.load(ROOT/'reference_results/public_scores.npz');rows=[];ceilings=[]
    for kind in ['fivecell','astrocyte']:
        d=pd.read_csv(DATA/f'{kind}.csv.gz');valid=valid_coordinates(d)&d.mapped.to_numpy();subsets={'covered':valid}
        if kind=='fivecell':subsets.update(non_K562=valid&d.cell.ne('K562').to_numpy(),non_K562_distal=valid&d.cell.ne('K562').to_numpy()&d.benchmark_distance.gt(100000).to_numpy())
        methods={'published_rE2G':full[kind+'/published'],'rE2G_L2':full[kind+'/matched'],'Atlas_full':full[kind+'/full'],'Atlas_gene':public['activity/'+kind+'/baseline'],'Atlas_public_cell':public['activity/'+kind+'/complete'],'ABC_score':d['ABC.Score'].to_numpy()}
        for subset,mask in subsets.items():
            groups=[list(ii) for ii in d.loc[mask].groupby(['cell','study','gene']).groups.values() if len(ii)>=4]
            for k in [1,3,5]:ceilings.append(dict(dataset=kind,subset=subset,k=k,groups=len(groups),positive_groups=sum(d.truth.iloc[ii].any() for ii in groups),all_eligible_positive_pairs=sum(int(d.truth.iloc[ii].sum()) for ii in groups),oracle_hits=sum(min(k,int(d.truth.iloc[ii].sum())) for ii in groups),selections=sum(min(k,len(ii)) for ii in groups)))
            for model,p in methods.items():
                for k in [1,3,5]:
                    hit=slots=ng=0
                    for ii in groups:
                        chosen=sorted(ii,key=lambda i:(-p[i],i))[:k];h=int(d.truth.iloc[chosen].sum());hit+=h;slots+=len(chosen);ng+=h>0
                    rows.append(dict(dataset=kind,subset=subset,model=model,k=k,AP=ap(d.truth[mask],p[mask],sample_weight=d.weight[mask] if kind=='fivecell' else None),pairs=int(mask.sum()),groups=len(groups),selections=slots,hits=hit,groups_with_hit=ng))
    result=pd.DataFrame(rows);result.to_csv(out/'existing_comparisons.csv',index=False);pd.DataFrame(ceilings).to_csv(out/'topk_label_ceiling.csv',index=False)
    c=pd.read_csv(DATA/'candidate_screen.csv.gz');rows=[]
    methods={'Atlas_gene':'score','reference_logistic':'reference_logistic','distance':'distance_score','public_contact':'public_contact','peak_support':'reference_peak_support','public_accessibility':'public_activity','public_accessibility_contact':'public_activity_contact'}
    for gene,d in c[c.regulatory_eligible].groupby('gene',sort=False):
        for subset in ['all','distal']:
            x=d if subset=='all' else d[d.author_tss_distance_bp>100000]
            for model,col in methods.items():
                ranked=x.sort_values([col,'enhancer_start','candidate_id'],ascending=[False,True,True])
                for k in [1,3,5]:
                    y=ranked.head(k);ids={a for v in y.hit_ids.fillna('') for a in v.split(';') if a};over={a for v in y.overlap_hit_ids.fillna('') for a in v.split(';') if a}
                    rows.append(dict(gene=gene,subset=subset,model=model,k=k,selections=len(y),hits=int(y.label.eq('hit').sum()),unique_hits=len(ids),overlap_hits=int(y.overlap_hit_ids.fillna('').ne('').sum()),unique_overlap_hits=len(over),control_overlaps=int(y.control_overlap.sum()),unassayed=int((~y.assayed).sum())))
    four=pd.DataFrame(rows).groupby(['subset','model','k'],sort=False)[['selections','hits','unique_hits','overlap_hits','unique_overlap_hits','control_overlaps','unassayed']].sum().reset_index();four.to_csv(out/'four_gene_topk.csv',index=False)
    expected=json.loads((ROOT/'reference_results/public_metrics.json').read_text())
    for key,value in [('existing_comparisons',result),('topk_label_ceiling',pd.DataFrame(ceilings)),('four_gene_topk',four)]:pd.testing.assert_frame_equal(value,pd.DataFrame(expected[key]),check_dtype=False,atol=1e-12,rtol=1e-12)
    r=pd.read_csv(DATA/'candidate_retrieval.csv.gz');summary=[]
    for kind,v in r.groupby('dataset'):
        subsets={'all_valid':np.ones(len(v),bool),'covered':v.mapped.to_numpy()}
        if kind=='fivecell':subsets.update(non_K562=v.cell.ne('K562').to_numpy(),non_K562_distal=(v.cell.ne('K562')&v.distance.gt(100000)).to_numpy())
        for label,mask in subsets.items():
            s=v.iloc[np.flatnonzero(mask)];s=s[s.truth.eq(1)]
            for col in ['default_center','inclusive_center']:summary.append(dict(dataset=kind,population=label,metric=col,n=len(s),retrieved=int(s[col].sum())))
    pd.DataFrame(summary).to_csv(out/'retrieval_summary.csv',index=False)

def run(out,refit=False,analyses='transfer,activity,source',bootstrap=0):
    out=Path(out);expected=dict(np.load(ROOT/'reference_results/public_scores.npz',allow_pickle=False));scores={};errors={}
    if refit:
        functions={'transfer':fit_transfer,'activity':fit_activity,'source':fit_source}
        for name in analyses.split(','):
            if name not in functions:raise ValueError('Unknown analysis: '+name)
            scores.update(functions[name](out));np.savez_compressed(out/'public_scores.npz',**scores)
        for key,p in scores.items():
            target=expected[key];mask=np.isfinite(target)
            for kind in ['fivecell','astrocyte']:
                if '/'+kind+'/' in key:
                    d=pd.read_csv(DATA/f'{kind}.csv.gz');mask&=valid_coordinates(d)&d.mapped.to_numpy()
            if not np.isfinite(p[mask]).all():raise AssertionError('Non-finite predictions: '+key)
            errors[key]=float(abs(p[mask]-target[mask]).max())
        if max(errors.values(),default=0)>1e-4:raise AssertionError('Fresh-fit score mismatch: '+str(sorted(errors.items(),key=lambda x:-x[1])[:5]))
    else:
        scores=expected
        scores['full/reporter/promoter']=np.load(ROOT/'reference_results/full_scores.npz')['reporter/full']
    result=metrics(scores,out,bootstrap);count,delta=check_metrics(result,refit);comparisons(out)
    print(f'Public-input analyses reproduced: {count} metrics; largest difference {delta:.3g}.')
