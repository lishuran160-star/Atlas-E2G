"""Enhancer-centred target recovery with annotation-defined candidates and frozen models."""
from pathlib import Path
import hashlib, json
import numpy as np
import pandas as pd
from scipy import sparse
from .features import ROOT, DATA, annotations
from .evaluation import valid_coordinates
from .target_query import candidates_for_enhancer

ACCESSIONS = {'K562':'ENCFF970QAX','GM12878':'ENCFF260CIO','WTC11':'ENCFF071ZZU',
              'HCT116':'ENCFF407YWT','Jurkat':'ENCFF681NXF'}


def cases():
    """Labels select evaluation queries, never the candidate genes within a query."""
    d = pd.read_csv(DATA/'fivecell.csv.gz')
    d = d[valid_coordinates(d)].copy()
    d['truth'] = d.truth.astype(bool)
    d['source_row'] = d.index
    genes, promoters = annotations()
    queries, pairs, labels = [], [], []
    for key, g in d.groupby(['cell','study','chrom','start','end'], sort=True):
        if not g.truth.any():
            continue
        cell, study, chrom, start, end = key
        region = f'{chrom}:{start}-{end}'
        qid = f'{cell}|{study}|{region}'
        c, info = candidates_for_enhancer(region, 1000000, cell)
        c.insert(0, 'query_id', qid)
        pairs.append(c)
        queries.append(dict(query_id=qid, cell=cell, study=study, chrom=chrom,
                            start=start, end=end, enhancer=region,
                            assayed_genes=g.gene.nunique(),
                            candidate_genes=info['candidate_genes']))
        for r in g.itertuples():
            m = genes[genes.symbol.eq(r.gene) & genes.chrom.eq(chrom)]
            gid = m.gene_id.iloc[0] if len(m)==1 else 'unresolved:'+r.gene
            if len(m)>1:
                t = promoters[promoters.gene_id.isin(m.gene_id)].assign(delta=lambda x:(x.tss-r.tss).abs())
                best = t[t.delta.eq(t.delta.min())].gene_id.unique()
                if len(best)==1 and t.delta.min()<=1000:
                    gid = best[0]
            labels.append(dict(query_id=qid, gene_id=gid, gene=r.gene,
                               truth=bool(r.truth), benchmark_distance=r.benchmark_distance,
                               source_row=r.source_row, mapped=bool(r.mapped)))
    return pd.DataFrame(queries), pd.concat(pairs, ignore_index=True), pd.DataFrame(labels)


def extract_comparators(queries, directory):
    """Use released ABC/rE2G predictions, with no refitting or score rescaling."""
    columns = ['#chr','start','end','TargetGeneEnsemblID','ABC.Score.Feature','Score']
    parts = []
    for cell, accession in ACCESSIONS.items():
        q = queries[queries.cell.eq(cell)]
        for frame in pd.read_csv(Path(directory)/(accession+'.bed.gz'), sep='\t',
                                 usecols=columns, chunksize=200000):
            for chrom, qq in q.groupby('chrom'):
                z = frame[frame['#chr'].eq(chrom)]
                if z.empty:
                    continue
                for r in qq.itertuples():
                    hit = z[z.start.lt(r.end) & z.end.gt(r.start)].copy()
                    if len(hit):
                        hit['query_id'] = r.query_id
                        parts.append(hit)
        print('Read published scores:', cell, flush=True)
    raw = pd.concat(parts, ignore_index=True)
    if raw.TargetGeneEnsemblID.isna().any():
        raise ValueError('Published predictions contain missing gene identifiers.')
    raw['gene_id'] = raw.TargetGeneEnsemblID.str.split('.').str[0]
    if not np.isfinite(raw[['ABC.Score.Feature','Score']]).all().all():
        raise ValueError('Published predictions contain non-finite scores.')
    return raw.groupby(['query_id','gene_id'], sort=True).agg(
        ABC_score=('ABC.Score.Feature','mean'), rE2G_score=('Score','mean'),
        matched_intervals=('Score','size')).reset_index()


def score_candidates(candidates, bed_dir, cache_dir, workers=8, offline=False):
    """Score all annotated candidates with the distributed frozen public model."""
    from .predict import make_regions, overlaps
    from .public_atlas import predict_public
    regions, _, _ = make_regions(candidates)
    folder = Path(cache_dir)/'overlaps'
    folder.mkdir(parents=True, exist_ok=True)
    signature = (DATA/'reference.tsv.gz').read_bytes()
    def cache_path(r):
        key = hashlib.blake2b(r.to_csv(index=False).encode()+signature,digest_size=16).hexdigest()
        return folder/(key+'.npz')
    whole = cache_path(regions)
    if whole.exists():
        matrix = sparse.load_npz(whole)
    else:
        matrix = overlaps(regions, bed_dir, workers)
        sparse.save_npz(whole, matrix)
    nref = len(pd.read_csv(DATA/'reference.tsv.gz', sep='\t'))
    if matrix.shape!=(len(regions),nref) or not np.isin(matrix.data,[0,1]).all():
        raise ValueError('Invalid cached overlap matrix.')
    lookup = {r:i for i,r in enumerate(regions.itertuples(index=False,name=None))}
    parts = []
    for cell, c in candidates.groupby('cell', sort=True):
        r, _, _ = make_regions(c)
        ix = [lookup[x] for x in r.itertuples(index=False,name=None)]
        sparse.save_npz(cache_path(r),matrix[ix])
        result, _ = predict_public(c,bed_dir,cache_dir,cell,workers,offline)
        parts.append(result)
        print('Scored target candidates:',cell,len(c),flush=True)
    return pd.concat(parts,ignore_index=True)


def gene_scores(scored):
    """Pooling rule is fixed before evaluation; every promoter remains in the saved input."""
    if not np.isfinite(scored[['score','gene_only_score','distance_bp']]).all().all():
        raise ValueError('Candidate scores must be finite.')
    return scored.groupby(['query_id','gene_id'],sort=True).agg(
        Atlas_cell=('score','max'), Atlas_gene=('gene_only_score','max'),
        distance_bp=('distance_bp','min')).reset_index().assign(
        nearest_TSS=lambda x:-x.distance_bp)


def universes(scored, comparators, labels):
    full = gene_scores(scored)
    common = full.merge(comparators, on=['query_id','gene_id'], validate='one_to_one')
    common = common[np.isfinite(common[['ABC_score','rE2G_score']]).all(axis=1)]
    d = pd.read_csv(DATA/'fivecell.csv.gz')
    public = np.load(ROOT/'reference_results/public_scores.npz',allow_pickle=False)
    full_saved = np.load(ROOT/'reference_results/full_scores.npz',allow_pickle=False)
    ix = labels.source_row.to_numpy()
    assayed = labels[['query_id','gene_id']].copy()
    assayed['Atlas_cell'] = public['activity/fivecell/complete'][ix]
    assayed['Atlas_gene'] = public['activity/fivecell/baseline'][ix]
    assayed['Atlas_full'] = full_saved['fivecell/full'][ix]
    assayed['ABC_score'] = d['ABC.Score'].to_numpy()[ix]
    assayed['rE2G_score'] = d.published.to_numpy()[ix]
    assayed['nearest_TSS'] = -d.benchmark_distance.to_numpy()[ix]
    cols = ['Atlas_cell','Atlas_gene','Atlas_full','ABC_score','rE2G_score','nearest_TSS']
    assayed = assayed[labels.mapped.to_numpy() & np.isfinite(assayed[cols]).all(axis=1)]
    assayed = assayed.groupby(['query_id','gene_id'],sort=True)[cols].max().reset_index()
    return {
        'all_annotation':(full,['Atlas_cell','Atlas_gene','nearest_TSS']),
        'shared_published':(common,['Atlas_cell','Atlas_gene','ABC_score','rE2G_score','nearest_TSS']),
        'assayed_pairs':(assayed,cols)
    }


def evaluate(scored, comparators, queries, labels, out, bootstrap=2000):
    out = Path(out)
    positives = labels[labels.truth].groupby(['query_id','gene_id'],sort=True).agg(
        benchmark_distance=('benchmark_distance','min')).reset_index()
    outcomes, target_ranks, coverage = [], [], []
    for universe, (frame, methods) in universes(scored,comparators,labels).items():
        for method in methods:
            order = frame.sort_values(['query_id',method,'gene_id'],
                                      ascending=[True,False,True],kind='stable').copy()
            order['rank'] = order.groupby('query_id').cumcount()+1
            ranked = positives.merge(order[['query_id','gene_id','rank']],
                                     on=['query_id','gene_id'],how='left',validate='one_to_one')
            ranked.insert(0,'method',method)
            ranked.insert(0,'universe',universe)
            target_ranks.append(ranked)
            for outcome, eligible in [('all_targets',ranked),('distal_targets',ranked[ranked.benchmark_distance.gt(100000)])]:
                for qid, g in eligible.groupby('query_id',sort=True):
                    ranks = g['rank'].fillna(np.inf).to_numpy()
                    outcomes.append(dict(universe=universe,method=method,outcome=outcome,
                        query_id=qid,known_targets=len(g),covered_targets=int(np.isfinite(ranks).sum()),
                        best_rank=float(ranks.min()),hit1=int((ranks<=1).any()),
                        hit3=int((ranks<=3).any()),targets1=int((ranks<=1).sum()),
                        targets3=int((ranks<=3).sum()),reciprocal_rank=1/float(ranks.min())))
        for r in queries.itertuples():
            f = frame[frame.query_id.eq(r.query_id)]
            p = positives[positives.query_id.eq(r.query_id)]
            coverage.append(dict(universe=universe,query_id=r.query_id,
                candidates=len(f),annotation_candidates=r.candidate_genes,known_targets=len(p),
                covered_targets=int(p.gene_id.isin(f.gene_id).sum())))
    outcomes = pd.DataFrame(outcomes).merge(queries,on='query_id',validate='many_to_one')
    coverage = pd.DataFrame(coverage).merge(queries[['query_id','cell','assayed_genes']],on='query_id')
    results = []
    for scope, mask in [('all_queries',np.ones(len(outcomes),bool)),('multiple_assayed_genes',outcomes.assayed_genes.ge(2))]:
        for subset, smask in [('all_cells',np.ones(len(outcomes),bool)),('non_K562',outcomes.cell.ne('K562'))]:
            for key, g in outcomes[mask & smask].groupby(['universe','method','outcome'],sort=True):
                finite = g.best_rank.replace(np.inf,np.nan)
                results.append(dict(scope=scope,subset=subset,universe=key[0],method=key[1],
                    outcome=key[2],queries=len(g),known_targets=int(g.known_targets.sum()),
                    covered_targets=int(g.covered_targets.sum()),hits1=int(g.hit1.sum()),
                    hits3=int(g.hit3.sum()),hit_rate1=g.hit1.mean(),hit_rate3=g.hit3.mean(),
                    target_recall1=g.targets1.sum()/g.known_targets.sum(),
                    target_recall3=g.targets3.sum()/g.known_targets.sum(),
                    mean_reciprocal_rank=g.reciprocal_rank.mean(),
                    median_best_rank_covered=finite.median()))
    metrics = pd.DataFrame(results)
    metrics.to_csv(out/'metrics.csv',index=False)
    outcomes.to_csv(out/'query_outcomes.csv',index=False)
    coverage.to_csv(out/'coverage.csv',index=False)
    pd.concat(target_ranks,ignore_index=True).to_csv(out/'known_target_ranks.csv',index=False)
    gene_scores(scored).to_csv(out/'candidate_gene_scores.csv.gz',index=False)
    if bootstrap:
        intervals=[]
        # Resample enhancer coordinates, keeping all associated cells/studies together.
        for scope in ['all_queries','multiple_assayed_genes']:
            base=outcomes if scope=='all_queries' else outcomes[outcomes.assayed_genes.ge(2)]
            for subset in ['all_cells','non_K562']:
                chosen=base if subset=='all_cells' else base[base.cell.ne('K562')]
                for key,g in chosen.groupby(['universe','outcome'],sort=True):
                    meta=g.drop_duplicates('query_id').set_index('query_id')
                    for endpoint in ['hit1','hit3']:
                        pivot=g.pivot(index='query_id',columns='method',values=endpoint)
                        meta=meta.loc[pivot.index]
                        groups=list(meta.groupby('enhancer',sort=True).indices.values())
                        rng=np.random.default_rng(84627)
                        estimates=[]
                        for _ in range(bootstrap):
                            ix=np.concatenate([groups[j] for j in rng.integers(len(groups),size=len(groups))])
                            estimates.append(pivot.iloc[ix].mean().to_numpy())
                        estimates=np.asarray(estimates)
                        ai=pivot.columns.get_loc('Atlas_cell')
                        for j,method in enumerate(pivot.columns):
                            delta=estimates[:,ai]-estimates[:,j]
                            lo,hi=np.quantile(delta,[.025,.975])
                            intervals.append(dict(scope=scope,subset=subset,universe=key[0],
                                outcome=key[1],metric=endpoint,reference='Atlas_cell',comparator=method,
                                difference=float(pivot.Atlas_cell.mean()-pivot[method].mean()),
                                low=float(lo),high=float(hi),clusters=len(groups),replicates=bootstrap))
        pd.DataFrame(intervals).to_csv(out/'paired_intervals.csv',index=False)
    return metrics


def run(out, from_saved=False, bed_dir=None, cache_dir=None, workers=8, offline=False,
        published_dir=None, bootstrap=2000):
    out = Path(out)
    queries, candidates, labels = cases()
    if from_saved:
        scored=pd.read_csv(ROOT/'reference_results/target_promoter_scores.csv.gz')
        identity=['query_id','candidate_id','gene_id','promoter_id']
        a=scored[candidates.columns].sort_values(identity).reset_index(drop=True)
        b=candidates.sort_values(identity).reset_index(drop=True)
        pd.testing.assert_frame_equal(a,b,check_dtype=False)
    else:
        scored=score_candidates(candidates,bed_dir,cache_dir,workers,offline)
        scored.to_csv(out/'target_promoter_scores.csv.gz',index=False)
    comparators=extract_comparators(queries,published_dir) if published_dir else pd.read_csv(DATA/'target_comparators.csv.gz')
    if published_dir:
        comparators.to_csv(out/'target_comparators.csv.gz',index=False)
    metrics=evaluate(scored,comparators,queries,labels,out,bootstrap)
    expected=ROOT/'reference_results/target_metrics.csv'
    if expected.exists():
        saved=pd.read_csv(expected)
        pd.testing.assert_frame_equal(metrics,saved,check_exact=False,atol=1e-12,rtol=1e-10)
    print('Target recovery reproduced:',len(queries),'queries;',len(candidates),'promoter pairs;',len(metrics),'summaries.')
    return metrics

