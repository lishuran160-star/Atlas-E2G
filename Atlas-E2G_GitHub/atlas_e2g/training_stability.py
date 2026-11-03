"""Gene-bootstrap refitting of chromosome-excluded promoter predictions."""
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path
import json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits
from .features import load, native, FACTORS, ROOT
from .model import fit_logistic
from .benchmark_full import design
from .evaluation import partial

METHODS = ('full', 'shared', 'without_identity')
_CACHE = None

def inputs():
    global _CACHE
    if _CACHE is None:
        d,b = load('source'); q,qb = load('reporter')
        covered = np.flatnonzero(d.mapped.to_numpy())
        genes, gene_index = np.unique(d.gene.to_numpy()[covered], return_inverse=True)
        _CACHE = d,b,q,qb,native(d),covered,genes,gene_index
    return _CACHE

def refit_one(replicate, seed, output, baseline=False):
    with threadpool_limits(limits=2):
        d,b,q,qb,raw,covered,genes,gene_index = inputs()
        if baseline: counts = np.ones(len(genes), dtype=int)
        else:
            rng = np.random.default_rng(np.random.SeedSequence([seed, replicate]))
            counts = np.bincount(rng.integers(len(genes), size=len(genes)), minlength=len(genes))
        rows = np.repeat(covered, counts[gene_index])
        scores = {m:np.full(len(q), np.nan) for m in METHODS}
        nf = len(FACTORS)
        for chrom in sorted(q.chrom.unique()):
            train = rows[d.chrom.to_numpy()[rows] != chrom]
            assert np.unique(d.truth.to_numpy()[train]).size == 2
            sc = StandardScaler().fit(raw[train])
            cs = StandardScaler().fit(b['neighborhood'][train])
            hit = q.chrom.eq(chrom).to_numpy()
            for method in METHODS:
                x = design(raw,b,method,sc,cs)
                fit = fit_logistic(x[train], d.truth.to_numpy()[train])
                coef = fit['coef'][8:8+nf]/np.sqrt(2) if method=='shared' else fit['coef'][8+nf:8+2*nf]
                profile = qb['profile'][hit]
                if method=='without_identity': profile = np.sort(profile,axis=1)
                scores[method][hit] = profile@coef
        metrics=[]
        for subset,mask in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
            for method,p in scores.items():
                assert np.isfinite(p).all()
                metrics.append(dict(replicate=replicate,subset=subset,method=method,
                    Spearman=float(spearmanr(q.response[mask],p[mask]).statistic),
                    adjusted_Spearman=float(partial(q[mask],p[mask],qb['amount'][mask]))))
        target=Path(output)/('baseline.npz' if baseline else f'refit_{replicate:03d}.npz')
        np.savez_compressed(target,**scores)
        return metrics

def summarize(metrics, output):
    rows=[]
    for (subset,method),grp in metrics.groupby(['subset','method'],sort=True):
        for metric in ['Spearman','adjusted_Spearman']:
            vals=grp[metric].to_numpy()
            rows.append(dict(subset=subset,method=method,metric=metric,n_refits=len(vals),
                median=float(np.median(vals)),p025=float(np.quantile(vals,.025)),p975=float(np.quantile(vals,.975)),
                minimum=float(vals.min()),maximum=float(vals.max()),fraction_positive=float((vals>0).mean())))
    for subset,grp in metrics.groupby('subset',sort=True):
        for metric in ['Spearman','adjusted_Spearman']:
            wide=grp.pivot(index='replicate',columns='method',values=metric)
            for control in METHODS[1:]:
                vals=(wide.full-wide[control]).to_numpy()
                rows.append(dict(subset=subset,method='full-minus-'+control,metric=metric,n_refits=len(vals),
                    median=float(np.median(vals)),p025=float(np.quantile(vals,.025)),p975=float(np.quantile(vals,.975)),
                    minimum=float(vals.min()),maximum=float(vals.max()),fraction_positive=float((vals>0).mean())))
    result=pd.DataFrame(rows);result.to_csv(Path(output)/'summary.csv',index=False)
    return result

def run(output,replicates=100,workers=4,seed=73421):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    protocol={'replicates':replicates,'seed':seed,'sampling':'source genes with replacement; all pairs of each sampled gene; shared resamples across methods',
        'fitting':'repeat chromosome-excluded fitting and fit preprocessing within each bootstrap training partition',
        'fixed':'model architecture, feature definitions, C=0.1, public references and reporter evaluation rows',
        'interpretation':'empirical training-refit distributions, not new independent validation or combined training-and-evaluation confidence intervals'}
    (output/'protocol.json').write_text(json.dumps(protocol,indent=2)+'\n',encoding='utf-8')
    d,b,q,qb,*_=inputs()
    baseline=refit_one(-1,seed,output,baseline=True)
    saved=np.load(ROOT/'reference_results/full_scores.npz',allow_pickle=False)
    actual=np.load(output/'baseline.npz',allow_pickle=False)
    for method in METHODS: np.testing.assert_allclose(actual[method],saved['reporter/'+method],rtol=1e-8,atol=1e-8)
    pd.DataFrame(baseline).to_csv(output/'baseline_metrics.csv',index=False)
    print('Baseline reproduced; starting',replicates,'paired gene-bootstrap refits.',flush=True)
    rows=[];start=time.time()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        tasks={pool.submit(refit_one,i,seed,output):i for i in range(replicates)}
        for j,future in enumerate(as_completed(tasks),1):
            rows.extend(future.result())
            if j%5==0 or j==replicates:
                pd.DataFrame(rows).sort_values(['replicate','subset','method']).to_csv(output/'metrics.csv',index=False)
                print(j,'/',replicates,'complete;',round(time.time()-start,1),'s',flush=True)
    metrics=pd.DataFrame(rows).sort_values(['replicate','subset','method']);metrics.to_csv(output/'metrics.csv',index=False)
    all_scores={m:np.stack([np.load(output/f'refit_{i:03d}.npz',allow_pickle=False)[m] for i in range(replicates)]) for m in METHODS}
    np.savez_compressed(output/'scores.npz',**all_scores)
    result=summarize(metrics,output)
    print(result[result.subset.eq('all')].to_string(index=False),flush=True)

def from_saved(output):
    """Recalculate training-refit metrics from supplied chromosome-excluded scores."""
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    q,qb=load('reporter');rows=[]
    with np.load(ROOT/'reference_results/training_stability_scores.npz',allow_pickle=False) as scores:
        for method in METHODS:
            for replicate,p in enumerate(scores[method]):
                for subset,mask in [('all',np.ones(len(q),bool)),('new_genes',q.novel.to_numpy())]:
                    rows.append(dict(replicate=replicate,subset=subset,method=method,
                        Spearman=float(spearmanr(q.response[mask],p[mask]).statistic),
                        adjusted_Spearman=float(partial(q[mask],p[mask],qb['amount'][mask]))))
    metrics=pd.DataFrame(rows).sort_values(['replicate','subset','method'])
    metrics.to_csv(output/'metrics.csv',index=False)
    result=summarize(metrics,output)
    expected=pd.read_csv(ROOT/'reference_results/training_stability_summary.csv')
    keys=['subset','method','metric'];cols=['n_refits','median','p025','p975','minimum','maximum','fraction_positive']
    a=result.set_index(keys).sort_index();b=expected.set_index(keys).sort_index()
    pd.testing.assert_index_equal(a.index,b.index)
    np.testing.assert_allclose(a[cols].to_numpy(),b[cols].to_numpy(),rtol=1e-10,atol=1e-12)
    print('Training stability reproduced:',len(a),'summaries from',int(expected.n_refits.max()),'gene-bootstrap refits.')
