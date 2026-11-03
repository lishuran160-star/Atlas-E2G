"""Find candidate target genes from an hg38 regulatory-region coordinate."""
from pathlib import Path
import json, re
import numpy as np
import pandas as pd
from .features import DATA, annotations


def parse_region(region):
    match = re.fullmatch(r'(chr(?:[1-9]|1[0-9]|2[0-2]|X|Y)):(\d+)-(\d+)', str(region).strip())
    if not match:
        raise ValueError('Use an hg38 interval such as chr8:127735000-127736000 (zero-based, half-open).')
    chrom, start, end = match.group(1), int(match.group(2)), int(match.group(3))
    sizes = dict(line.split() for line in (DATA/'hg38.primary.chrom.sizes').read_text().splitlines())
    if not 0 <= start < end <= int(sizes[chrom]):
        raise ValueError('Enhancer interval is empty or extends beyond the chromosome.')
    return chrom, start, end


def candidates_for_enhancer(region, radius=1000000, cell='unspecified'):
    """Enumerate all nearby annotated promoter clusters without using screen data."""
    if not isinstance(radius, (int, np.integer)) or not 1000 <= radius <= 5000000:
        raise ValueError('Search radius must be between 1 kb and 5 Mb.')
    chrom, start, end = parse_region(region)
    genes, promoters = annotations()
    mid = (start + end)/2
    candidates = promoters[promoters.chrom.eq(chrom) & promoters.tss.sub(mid).abs().le(radius)].copy()
    count = len(candidates)
    sizes = dict(line.split() for line in (DATA/'hg38.primary.chrom.sizes').read_text().splitlines())
    candidates = candidates[candidates.tss.ge(1000) & (candidates.tss+1001).le(int(sizes[chrom]))]
    candidates = candidates.merge(genes[['gene_id','gene_type']], on='gene_id', validate='many_to_one')
    if candidates.empty:
        raise ValueError('No annotated promoters in this interval. Increase --radius within the 5 Mb limit.')
    candidates = candidates.sort_values(['tss','gene_id','promoter_id'], kind='stable').reset_index(drop=True)
    result = pd.DataFrame({
        'candidate_id': chrom+':'+str(start)+'-'+str(end)+'|'+candidates.promoter_id,
        'enhancer_chrom': chrom, 'enhancer_start': start, 'enhancer_end': end,
        'promoter_chrom': chrom, 'tss': candidates.tss, 'gene': candidates.symbol,
        'gene_id': candidates.gene_id, 'gene_type': candidates.gene_type,
        'promoter_id': candidates.promoter_id, 'transcript': candidates.transcript,
        'distance_bp': candidates.tss.sub(mid).abs(), 'cell': cell,
    })
    info = dict(genome='hg38', enhancer=f'{chrom}:{start}-{end}', radius_bp=int(radius),
                candidate_genes=int(result.gene_id.nunique()), candidate_promoters=len(result),
                excluded_invalid_promoter_windows=count-len(result),
                candidate_policy='All annotated promoter clusters within the radius of the enhancer midpoint; no gene-type, expression or screened-gene filter; no candidate cap.',
                gene_aggregation='CAGE-usage-weighted mean pair score multiplied by public detection support; ties ordered by stable identifiers.')
    return result, info


def rank_genes(scored,cell_type=None):
    """Retain every promoter and rank genes using separate public transcription evidence."""
    from .transcription import gene_priorities
    return gene_priorities(scored,cell_type)


def query(region, bed_dir, cache_dir, output, cell_type=None, radius=1000000, workers=8, offline=False, candidates_only=False):
    candidates, info = candidates_for_enhancer(region,radius,cell_type or 'unspecified')
    output = Path(output); output.mkdir(parents=True,exist_ok=False)
    candidates.to_csv(output/'candidates.csv',index=False)
    info.update(model='Atlas-E2G',scored=not candidates_only)
    if not candidates_only:
        from .public_atlas import predict_public
        scored, details = predict_public(candidates,bed_dir,cache_dir,cell_type,workers,offline)
        promoters, genes = rank_genes(scored,cell_type)
        promoters.to_csv(output/'ranked_promoters.csv',index=False)
        genes.to_csv(output/'ranked_genes.csv',index=False)
        genes.head(50).to_csv(output/'top_genes.csv',index=False)
        info.update(details)
        info['transcription_reference']='GSE118075: conventional CAGE, exact cell when present, pooled otherwise'
        info['gene_priority']='Sum of pair score times normalized (promoter CAGE CPM + 0.1), multiplied by gene detection support; no candidate filtering.'
        info['interpretation'] = 'Candidate target-gene priorities, not expression fold changes or calibrated probabilities of regulation. All annotated gene types are retained; unassayed candidates are unknown.'
    (output/'query.json').write_text(json.dumps(info,indent=2)+'\n',encoding='utf-8')
    return info
