"""Rebuild public features and refit models from public reference files."""
from pathlib import Path
import re
import numpy as np,pandas as pd
from scipy import sparse
from scipy.special import expit
from scipy.optimize import minimize
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import average_precision_score
from .features import DATA,ROOT,FACTORS,load,reference,profiles,normalize
from .model import splits
from .public_atlas import PublicAtlas,match_cell,public_design_matrix

def valid_pairs(d):
    from .evaluation import valid_coordinates
    return d.mapped.to_numpy()&valid_coordinates(d)

def optional_identity(d,B,global_identity):
    meta=reference();ix=meta.factor.map({v:i for i,v in enumerate(FACTORS)}).to_numpy();_,prior=profiles(B);out=global_identity.copy()
    for cell,gg in d.groupby('cell'):
        canonical=lambda s:re.sub('[^a-z0-9]','',str(s).casefold())
        label=canonical(cell);label={'astrocyte':'astrocytes'}.get(label,label);sel=meta.cell.map(canonical).eq(label).to_numpy()
        if not sel.any():continue
        n=np.bincount(ix[sel],minlength=len(FACTORS));g=sparse.csr_matrix((np.ones(sel.sum()),(np.flatnonzero(sel),ix[sel])),shape=(len(meta),len(FACTORS)));local=((B@g).toarray()+5*prior)/(n[None,:]+5)
        out[gg.index]=np.c_[normalize(local[gg.enhancer_uid]),normalize(local[gg.promoter_uid])]
    return out

def fit_views(x1,x2,y):
    # Convert once so mixed-precision inputs do not slow every solver step.
    x1=np.asarray(x1,dtype=np.float64,order='C')
    x2=np.asarray(x2,dtype=np.float64,order='C')
    n=x1.shape[1]
    def objective(w):
        loss=5*np.dot(w[:-1],w[:-1]);grad=10*w[:-1];intercept=0.
        for x in (x1,x2):
            z=x@w[:-1]+w[-1];res=.5*(expit(z)-y);loss+=.5*(np.logaddexp(0,z).sum()-y@z);grad+=x.T@res;intercept+=res.sum()
        return loss,np.r_[grad,intercept]
    fit=minimize(objective,np.zeros(n+1),jac=True,method='L-BFGS-B',options=dict(maxiter=1500,ftol=1e-9))
    if not fit.success:raise RuntimeError(fit.message)
    return dict(coef=fit.x[:-1],intercept=float(fit.x[-1]))

def fit_heads(low,global_identity,cell_identity,activity,mask,y):
    models={}
    for name,x in [('gene',low),('cell',np.c_[low,activity])]:
        sc=StandardScaler().fit(x[mask]);z=sc.transform(x);m=fit_views(np.c_[z[mask],global_identity[mask]],np.c_[z[mask],cell_identity[mask]],y[mask]);m.update(mean=sc.mean_,scale=sc.scale_);models[name]=m
    return models

def predict_heads(models,low,identity,activity,present):
    out={}
    for name,x in [('gene',low),('cell',np.c_[low,activity])]:
        m=models[name];out[name]=expit(np.c_[(x-m['mean'])/m['scale'],identity]@m['coef']+m['intercept'])
    out['cell'][~present]=out['gene'][~present];return out

def data(kind,atlas):
    d,b=load(kind);valid=valid_pairs(d);ix=np.flatnonzero(valid);q=d.iloc[ix].reset_index(drop=True);co=atlas.pair_contact(q);lo,g=public_design_matrix(q,{k:v[ix] for k,v in b.items()},co)
    low=np.zeros((len(d),lo.shape[1]));identity=np.zeros((len(d),g.shape[1]));low[ix]=lo;identity[ix]=g;act=np.zeros((len(d),5));present=np.zeros(len(d),bool)
    for cell,gg in q.groupby('cell'):
        row=match_cell(cell)
        if row is None:continue
        ids=gg.index.to_numpy();act[ix[ids]]=atlas.activity_features(gg.reset_index(drop=True),co[ids],row);present[ix[ids]]=True
    return d,low,identity,act,present,valid

def run(out,cache_dir,offline=False,source_folds=False):
    out=Path(out);atlas=PublicAtlas(cache_dir,offline);d,low,g,activity,present,valid=data('source',atlas);B=sparse.load_npz(DATA/'source_overlap.npz');cell=optional_identity(d,B,g.astype(float));y=d.truth.to_numpy()
    assert present[valid].all(),'Source public accessibility missing'
    if source_folds:
        oof={};chrom=d.chrom.map(lambda s:23 if s=='chrX' else 24 if s=='chrY' else int(s[3:])).to_numpy();folds=[('chromosome_blocks',str(i),valid&(chrom%5==i),valid&(chrom%5!=i)) for i in range(5)]
        folds.extend((pr,f,test&valid,train&valid) for pr,f,test,train in splits(d) if pr=='study_gene_element_disjoint')
        for protocol,fold,test,train in folds:
            models=fit_heads(low,g,cell,activity,train,y);pred=predict_heads(models,low[test],g[test],activity[test],present[test])
            for name,p in pred.items():oof.setdefault(protocol+'/'+name,np.full(len(d),np.nan))[test]=p
            print('Public source fold',protocol,fold,flush=True)
        pd.DataFrame([dict(protocol=k.split('/')[0],input=k.split('/')[1],AP=average_precision_score(y[valid],v[valid])) for k,v in oof.items()]).to_csv(out/'source_folds.csv',index=False)
    models=fit_heads(low,g,cell,activity,valid,y);np.savez_compressed(out/'public_reference.npz',**{name+'_'+key:value for name,m in models.items() for key,value in m.items()});rows=[];saved={}
    for kind in ['fivecell','astrocyte']:
        d,low,g,activity,present,valid=data(kind,atlas);pred=predict_heads(models,low,g,activity,present);expected=np.load(ROOT/'reference_results/public_scores.npz');subsets={'covered':valid}
        if kind=='fivecell':subsets.update(non_K562=valid&d.cell.ne('K562').to_numpy(),non_K562_distal=valid&d.cell.ne('K562').to_numpy()&(d.benchmark_distance>100000).to_numpy())
        for name,p in pred.items():
            saved[kind+'/'+name]=p;error=float(abs(p[valid]-expected['activity/'+kind+'/'+{'gene':'baseline','cell':'complete'}[name]][valid]).max())
            for subset,ix in subsets.items():rows.append(dict(dataset=kind,subset=subset,input=name,n=int(ix.sum()),AP=average_precision_score(d.truth[ix],p[ix],sample_weight=d.weight[ix] if kind=='fivecell' else None),max_score_difference=error))
        print('Public benchmark',kind,int(valid.sum()),flush=True)
    result=pd.DataFrame(rows);result.to_csv(out/'public_results.csv',index=False);np.savez_compressed(out/'public_scores.npz',**saved)
    difference=float(result.max_score_difference.max())
    if not np.isfinite(difference) or difference>=1e-4:raise AssertionError('Public-input refit did not reproduce stored scores')
    print(f'Public-input results reproduced; largest score difference {difference:.3g}.')
