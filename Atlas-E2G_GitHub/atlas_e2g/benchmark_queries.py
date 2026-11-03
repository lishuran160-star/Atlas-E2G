"""Reproduce matched public-input ablations and transcription-supported target queries."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from scipy import sparse
from sklearn.metrics import average_precision_score
from .features import ROOT,DATA,pair_features,public_design
from .predict import make_regions
from .benchmark_public import data
from .benchmark_targets import cases
from .transcription import annotate
from .evaluation import valid_coordinates

OUT=None

def load(kind):
    if kind!='targets':
        v=data(kind)
        return {k:(x.astype(np.float32) if isinstance(x,np.ndarray) and x.dtype.kind=='f' else x) for k,x in v.items()}
    c=pd.read_csv(ROOT/'reference_results/target_promoter_scores.csv.gz')
    regions,ei,pi=make_regions(c)
    b=sparse.load_npz(DATA/'target_overlap.npz')
    if b.shape[0]!=len(regions):raise ValueError('Candidate overlap shape differs.')
    _,g=public_design(pair_features(c,b,ei,pi))
    z=np.load(DATA/'target_public_features.npz')
    d=c.rename(columns={'enhancer_chrom':'chrom','enhancer_start':'start','enhancer_end':'end','promoter_chrom':'pchrom'})
    return dict(d=d,g=g,c=g,lo=z['lo'],activity=z['activity'],available=z['available'])

def build_gene_scores(v,z):
    d=v['d'];records=[]
    for view in ['gene','cell']:
        for qid,qq in d.groupby('query_id',sort=True):
            f=annotate(qq,qq.cell.iloc[0] if view=='cell' else None)
            for m in ['Atlas-E2G','without_identity','without_profiles']:
                f[m]=z['targets/'+view+'/'+m][qq.index]
            for gid,q in f.groupby('gene_id',sort=True):
                w=q.public_promoter_cage_cpm.to_numpy()+.1;w/=w.sum()
                support=float(q.public_detection_support.iloc[0])
                r=dict(view=view,query_id=qid,gene_id=gid,support=support,promoters=len(q),minimum_distance=q.distance_bp.min(),unweighted_Atlas=float(q['Atlas-E2G'].max()))
                for m in ['Atlas-E2G','without_identity','without_profiles']:r[m]=float(q[m].to_numpy()@w)*support
                r['distance']=float((1/np.maximum(q.distance_bp.to_numpy(),2500))@w)*support
                r['support_only']=support;records.append(r)
    return pd.DataFrame(records)

from scipy.special import expit
from sklearn.preprocessing import StandardScaler
from .rebuild_public import fit_views
NF=1598

def profile(g,method):
    if method=='without_identity':return np.c_[np.sort(g[:,:NF],axis=1),np.sort(g[:,NF:],axis=1)]
    if method=='without_profiles':return g[:,:0]
    return g

def fit(v,train,methods=('Atlas-E2G','without_identity','without_profiles')):
    models={};y=v['d'].truth.to_numpy()
    for view in ['gene','cell']:
        x=v['lo'] if view=='gene' else np.c_[v['lo'],v['activity']]
        sc=StandardScaler().fit(x[train]);xx=sc.transform(x)
        for method in methods:
            g=profile(v['g'],method);c=profile(v['c'],method)
            m=fit_views(np.c_[xx[train],g[train]],np.c_[xx[train],c[train]],y[train]);m.update(mean=sc.mean_,scale=sc.scale_)
            models[view+'/'+method]=m
    return models

def predict(v,models):
    out={}
    for key,m in models.items():
        view,method=key.split('/');x=v['lo'] if view=='gene' else np.c_[v['lo'],v['activity']]
        xx=np.c_[(x-m['mean'])/m['scale'],profile(v['g'],method)]
        out[key]=expit(xx@m['coef']+m['intercept'])
    methods={k.split('/')[1] for k in models}
    for method in methods:out['cell/'+method][~v['available']]=out['gene/'+method][~v['available']]
    return out

def frozen_models():
    z=np.load(ROOT/'models/public_reference.npz',allow_pickle=False)
    return {view+'/Atlas-E2G':{k:z[view+'_'+k] for k in ['coef','intercept','mean','scale']} for view in ['gene','cell']}

def paired_ap(d, scores, methods, groups, nboot=2000):
    y=d.truth.to_numpy().astype(bool)
    weight=d.eval_weight.to_numpy()
    rng=np.random.default_rng(20261102); estimates=[]
    gs=list(d.groupby(groups,sort=True).indices.values())
    for _ in range(nboot):
        ix=np.concatenate([gs[j] for j in rng.integers(0,len(gs),len(gs))])
        if y[ix].sum()==0 or y[ix].sum()==len(ix):continue
        estimates.append([average_precision_score(y[ix],scores[m][ix],sample_weight=weight[ix]) for m in methods])
    b=np.asarray(estimates)
    return {m:(float(average_precision_score(y,scores['Atlas-E2G'],sample_weight=weight)-average_precision_score(y,scores[m],sample_weight=weight)),*np.quantile(b[:,0]-b[:,j],[.025,.975]).tolist()) for j,m in enumerate(methods) if j}

def evaluate():
    z=np.load(OUT/'raw_scores.npz');methods=['Atlas-E2G','without_identity','without_profiles']
    rows=[];intervals=[];qa={}
    for kind in ['fivecell','astrocyte']:
        v=load(kind);d=v['d'];base=valid_coordinates(d)&d.mapped.to_numpy()
        d['eval_weight']=d.weight if kind=='fivecell' else 1.
        groups={'all':base,'non_K562':base&d.cell.ne('K562').to_numpy(),'non_K562_distal':base&d.cell.ne('K562').to_numpy()&d.benchmark_distance.gt(100000).to_numpy()} if kind=='fivecell' else {'all':base}
        for view in ['gene','cell']:
            scores={m:z[kind+'/'+view+'/'+m] for m in methods}
            for subset,mask in groups.items():
                for m in methods:
                    rows.append(dict(kind=kind,subset=subset,view=view,method=m,n=int(mask.sum()),positives=int(d.truth[mask].sum()),AP=average_precision_score(d.truth[mask],scores[m][mask],sample_weight=d.eval_weight[mask])))
                if view=='cell':
                    delta=paired_ap(d[mask].reset_index(drop=True),{m:s[mask] for m,s in scores.items()},methods,'gene')
                    for m,(diff,lo,hi) in delta.items():intervals.append(dict(kind=kind,subset=subset,view=view,control=m,delta=diff,ci_low=lo,ci_high=hi))
    pd.DataFrame(rows).to_csv(OUT/'pair_summary.csv',index=False)
    pd.DataFrame(intervals).to_csv(OUT/'pair_intervals.csv',index=False)

    v=load('targets');d=v['d'];frame=build_gene_scores(v,z)
    frame.to_csv(OUT/'query_gene_scores.csv.gz',index=False)
    queries,_,labels=cases()
    positive=labels[labels.truth].groupby(['query_id','gene_id']).benchmark_distance.min().reset_index()
    comp=pd.read_csv(DATA/'target_comparators.csv.gz')
    shared=frame.merge(comp,on=['query_id','gene_id'],validate='many_to_one')
    shared['ABC_supported']=shared.ABC_score*shared.support
    shared['rE2G_supported']=shared.rE2G_score*shared.support
    outcomes=[];ranks=[];coverage=[]
    for univ,full in [('all_annotation',frame),('shared_published',shared)]:
        for view,f in full.groupby('view'):
            names=methods+['distance','unweighted_Atlas','support_only']+(['ABC_score','rE2G_score','ABC_supported','rE2G_supported'] if univ=='shared_published' else [])
            for m in names:
                q=f.sort_values(['query_id',m,'gene_id'],ascending=[True,False,True],kind='stable').copy()
                q['rank']=q.groupby('query_id').cumcount()+1
                r=positive.merge(q[['query_id','gene_id','rank']],on=['query_id','gene_id'],how='left')
                r['rank']=r['rank'].fillna(np.inf)
                ranks.append(r.assign(universe=univ,view=view,method=m))
                for endpoint,rr in [('all_targets',r),('distal_targets',r[r.benchmark_distance.gt(100000)])]:
                    for qid,g in rr.groupby('query_id'):
                        rank=float(g['rank'].min());outcomes.append(dict(universe=univ,view=view,method=m,outcome=endpoint,query_id=qid,rank=rank,hit1=int(rank<=1),hit3=int(rank<=3),rr=1/rank))
            r=positive.merge(f[['query_id','gene_id']],on=['query_id','gene_id'],how='left',indicator=True)
            coverage.append(dict(universe=univ,view=view,candidate_genes=len(f),known_targets=len(r),covered=int(r['_merge'].eq('both').sum())))
    out=pd.DataFrame(outcomes).merge(queries,on='query_id',validate='many_to_one')
    out.to_csv(OUT/'query_outcomes.csv.gz',index=False)
    pd.concat(ranks).to_csv(OUT/'query_target_ranks.csv.gz',index=False)
    results=[];ci=[]
    for scope,mask in [('all_queries',np.ones(len(out),bool)),('multiple_assayed_genes',out.assayed_genes.ge(2))]:
        for subset,sm in [('all_cells',np.ones(len(out),bool)),('non_K562',out.cell.ne('K562'))]:
            for key,q in out[mask&sm].groupby(['universe','view','method','outcome']):
                results.append(dict(scope=scope,subset=subset,universe=key[0],view=key[1],method=key[2],outcome=key[3],queries=len(q),hits1=int(q.hit1.sum()),hits3=int(q.hit3.sum()),top1=q.hit1.mean(),top3=q.hit3.mean(),MRR=q.rr.mean()))
            for key,q in out[mask&sm].groupby(['universe','view','outcome']):
                meta=q.drop_duplicates('query_id').set_index('query_id')
                for endpoint in ['hit1','hit3']:
                    pv=q.pivot(index='query_id',columns='method',values=endpoint)
                    clusters=list(meta.loc[pv.index].groupby('enhancer',sort=True).indices.values())
                    rng=np.random.default_rng(20261102)
                    b=np.array([pv.to_numpy()[np.concatenate([clusters[j] for j in rng.integers(0,len(clusters),len(clusters))])].mean(axis=0) for _ in range(2000)])
                    ai=pv.columns.get_loc('Atlas-E2G')
                    for j,m in enumerate(pv.columns):
                        if m=='Atlas-E2G':continue
                        lo,hi=np.quantile(b[:,ai]-b[:,j],[.025,.975])
                        ci.append(dict(scope=scope,subset=subset,universe=key[0],view=key[1],outcome=key[2],endpoint=endpoint,control=m,delta=float(pv['Atlas-E2G'].mean()-pv[m].mean()),ci_low=float(lo),ci_high=float(hi),queries=len(pv)))
    result=pd.DataFrame(results);result.to_csv(OUT/'query_summary.csv',index=False)
    pd.DataFrame(ci).to_csv(OUT/'query_intervals.csv',index=False)
    pd.DataFrame(coverage).to_csv(OUT/'query_coverage.csv',index=False)
    qa.update(coverage=coverage,queries=len(queries),promoter_pairs=len(d),same_candidates_all_controls=True,identical_support_and_usage_weights=True,distance_aggregated_across_same_promoters=True,unassayed_as_negative=False)
    (OUT/'evaluation_summary.json').write_text(json.dumps(qa,indent=2))
    print(pd.DataFrame(rows).to_string(index=False))
    print(result[(result.scope=='multiple_assayed_genes')&(result.subset=='all_cells')&(result.view=='cell')&(result.outcome=='all_targets')].to_string(index=False))
    print(json.dumps(qa,indent=2))

def refit():
    from .model import splits
    v=load('source');d=v['d'];valid=d.mapped.to_numpy()&valid_coordinates(d)
    chrom=d.chrom.map(lambda s:23 if s=='chrX' else 24 if s=='chrY' else int(s[3:])).to_numpy()
    folds=[('chromosome_blocks',str(i),valid&(chrom%5==i),valid&(chrom%5!=i)) for i in range(5)]
    folds.extend((pr,f,te&valid,tr&valid) for pr,f,te,tr in splits(d) if pr=='study_gene_element_disjoint')
    oof={}
    for protocol,fold,test,train in folds:
        models=fit(v,train);pred=predict(v,models)
        for key,value in pred.items():oof.setdefault(protocol+'/'+key,np.full(len(d),np.nan))[test]=value[test]
        print('Fitted',protocol,fold,flush=True)
    np.savez_compressed(OUT/'source_oof_scores.npz',**oof)
    models=fit(v,valid,methods=('without_identity','without_profiles'));models.update(frozen_models())
    np.savez_compressed(OUT/'control_models.npz',**{k+'/'+n:a for k,m in models.items() for n,a in m.items()})
    scores={}
    for kind in ['fivecell','astrocyte','targets']:
        pred=predict(load(kind),models);scores.update({kind+'/'+k:a for k,a in pred.items()})
    np.savez_compressed(OUT/'raw_scores.npz',**scores)

def run(output,from_saved=False):
    import shutil
    global OUT
    OUT=Path(output)
    if from_saved:
        for name in ['raw_scores.npz','source_oof_scores.npz']:
            shutil.copy2(ROOT/'reference_results'/('query_'+name),OUT/name)
    else:refit()
    z=np.load(OUT/'source_oof_scores.npz');d=pd.read_csv(DATA/'source.csv.gz');rows=[]
    for key,p in z.items():
        protocol,view,method=key.split('/');mask=np.isfinite(p)
        rows.append(dict(protocol=protocol,view=view,method=method,n=int(mask.sum()),AP=average_precision_score(d.truth[mask],p[mask])))
    pd.DataFrame(rows).to_csv(OUT/'source_summary.csv',index=False)
    evaluate()
    if from_saved:
        for filename in ['source_summary.csv','pair_summary.csv','query_summary.csv']:
            name=filename if filename.startswith('query_') else 'query_'+filename
            actual=pd.read_csv(OUT/filename);expected=pd.read_csv(ROOT/'reference_results'/name)
            pd.testing.assert_frame_equal(actual,expected,check_exact=False,atol=1e-12,rtol=1e-10)
        print('Matched public-input and query results reproduced.')
