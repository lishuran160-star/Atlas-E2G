"""Public transcription support for annotation-defined target-gene queries."""
import json
from functools import lru_cache
import numpy as np
import pandas as pd
from .features import DATA

@lru_cache(None)
def reference():
    z=np.load(DATA/'transcription.npz',allow_pickle=False)
    meta=json.loads((DATA/'transcription.json').read_text(encoding='utf-8'))
    p=pd.read_csv(DATA/'promoters.tsv.gz',sep='\t')
    if not np.array_equal(p.promoter_id.to_numpy(),z['promoter_ids']):
        raise ValueError('CAGE reference and promoter annotation differ.')
    counts=z['counts'].astype(float)
    cpm=counts/np.array([s['total_tags'] for s in meta['samples']])*1e6
    totals=pd.DataFrame(counts).assign(gene_id=p.gene_id).groupby('gene_id',sort=False).sum()
    mapped=pd.Series(z['mapped'],index=p.index).groupby(p.gene_id).any()
    return p,counts,cpm,totals,mapped,meta

def annotate(scored,cell_type=None):
    """Keep all candidates and attach evidence independently of the pair score."""
    p,counts,cpm,totals,mapped,meta=reference()
    norm=lambda x:''.join(c for c in str(x).casefold() if c.isalnum())
    ix=[i for i,s in enumerate(meta['samples']) if norm(s['cell'])==norm(cell_type)] if cell_type else []
    matched=bool(ix)
    if not ix:ix=list(range(len(meta['samples'])))
    evidence=pd.Series(cpm[:,ix].mean(axis=1),index=p.promoter_id)
    q=scored.copy()
    q['public_promoter_cage_cpm']=q.promoter_id.map(evidence).fillna(0.)
    support=(totals.iloc[:,ix].gt(0).sum(axis=1)+.5)/(len(ix)+1)
    support.loc[~mapped.reindex(support.index).fillna(False)]=1.
    q['public_detection_support']=q.gene_id.map(support).fillna(1.)
    q['public_transcription_reference']='cell_matched' if matched else 'pooled'
    return q

def gene_priorities(scored,cell_type=None):
    """Usage-weighted pair scores and detection support; no expression filtering."""
    q=annotate(scored,cell_type)
    q['public_promoter_weight']=q.public_promoter_cage_cpm+.1
    q['public_promoter_weight']/=q.groupby('gene_id').public_promoter_weight.transform('sum')
    q['priority_contribution']=q.score*q.public_promoter_weight*q.public_detection_support
    ordered=q.sort_values(['priority_contribution','gene_id','promoter_id'],ascending=[False,True,True],kind='stable')
    best=ordered.drop_duplicates('gene_id').copy()
    groups=q.groupby('gene_id',sort=False)
    best['raw_max_pair_score']=best.gene_id.map(groups.score.max())
    best['priority_score']=best.gene_id.map(groups.priority_contribution.sum())
    best['promoters_considered']=best.gene_id.map(groups.size())
    best['closest_tss_distance_bp']=best.gene_id.map(groups.distance_bp.min())
    if 'gene_only_score' in q:
        pooled=annotate(scored,None)
        weight=pooled.public_promoter_cage_cpm+.1
        weight/=weight.groupby(pooled.gene_id).transform('sum')
        gene_priority=(pooled.gene_only_score*weight*pooled.public_detection_support).groupby(pooled.gene_id).sum()
        best['gene_only_priority_score']=best.gene_id.map(gene_priority)
    best=best.rename(columns={'promoter_id':'representative_promoter_id','score':'representative_pair_score'})
    best=best.sort_values(['priority_score','gene_id'],ascending=[False,True],kind='stable')
    best.insert(0,'gene_rank',range(1,len(best)+1))
    promoters=q.sort_values(['score','gene_id','promoter_id'],ascending=[False,True,True],kind='stable')
    promoters.insert(0,'promoter_rank',range(1,len(promoters)+1))
    return promoters.reset_index(drop=True),best.reset_index(drop=True)

