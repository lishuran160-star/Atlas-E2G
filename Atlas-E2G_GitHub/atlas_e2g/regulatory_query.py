"""Retrieve and score candidates around a representative annotated transcript."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
import requests
from atlas_e2g.features import DATA,annotations
from atlas_e2g.predict import resolve_gene
from atlas_e2g.public_atlas import predict_public

RADIUS=1000000

def fetch_annotation(gene,cache,offline=False):
    g,_=resolve_gene(gene);cache=Path(cache);cache.mkdir(parents=True,exist_ok=True)
    path=cache/(g.symbol+'.lookup.json')
    url=f'https://rest.ensembl.org/lookup/id/{g.gene_id}?expand=1;mane=1;content-type=application/json'
    if not path.exists():
        if offline:raise FileNotFoundError('Gene annotation not cached: '+g.symbol)
        for attempt in range(2):
            try:
                response=requests.get(url,timeout=45);response.raise_for_status();x=response.json()
                assert x.get('assembly_name')=='GRCh38' and x['id']==g.gene_id
                path.write_text(json.dumps(x,indent=2),encoding='utf8');break
            except (requests.RequestException,ValueError):
                if attempt:raise
    x=json.loads(path.read_text(encoding='utf8'))
    assert x['assembly_name']=='GRCh38' and x['id']==g.gene_id
    mane=[t for t in x['Transcript'] if any(m.get('type')=='MANE_Select' for m in t.get('MANE',[]))]
    canonical=[t for t in x['Transcript'] if t.get('is_canonical')]
    choices=mane or canonical
    if len(choices)!=1:raise ValueError('No unique MANE Select or canonical transcript for '+g.symbol)
    anchor=choices[0]
    if str(anchor['seq_region_name']) not in [str(i) for i in range(1,23)]+['X','Y']:
        raise ValueError('Selected transcript is on a patch or unsupported contig; no silent remapping')
    tss=int(anchor['start']-1 if anchor['strand']==1 else anchor['end']-1)
    cds=[]
    for t in x['Transcript']:
        if 'Translation' not in t:continue
        a,b=t['Translation']['start']-1,t['Translation']['end']
        for e in t.get('Exon',[]):
            s,z=max(a,e['start']-1),min(b,e['end'])
            if s<z:cds.append((int(s),int(z)))
    utr5=[]
    if 'Translation' in anchor:
        boundary=anchor['Translation']['start']-1 if anchor['strand']==1 else anchor['Translation']['end']
        for e in anchor.get('Exon',[]):
            s,z=(e['start']-1,min(e['end'],boundary)) if anchor['strand']==1 else (max(e['start']-1,boundary),e['end'])
            if s<z:utr5.append((int(s),int(z)))
    return g,dict(gene=g.symbol,gene_id=g.gene_id,chrom='chr'+anchor['seq_region_name'],tss=tss,strand=anchor['strand'],transcript=anchor['id']+'.'+str(anchor['version']),selection='MANE_Select' if mane else 'Ensembl_canonical',cds=sorted(set(cds)),utr5=sorted(set(utr5)),annotation_url=url,interpretation='Public representative-transcript prior; not a measurement of cell-specific promoter use.')

def candidates(gene,annotation_cache,radius=RADIUS,offline=False):
    if not 1000<=radius<=1000000:raise ValueError('Annotated candidate search supports up to 1 Mb per side')
    g,a=fetch_annotation(gene,annotation_cache,offline);t=a['tss'];chrom=a['chrom']
    sizes_path=DATA/'hg38.primary.chrom.sizes'
    sizes=dict(line.split() for line in sizes_path.read_text().splitlines())
    chrom_end=int(sizes[chrom])
    starts=np.arange(max(0,(t-radius)//1000)*1000,((t+radius)//1000+1)*1000,1000,dtype=np.int64)
    starts=starts[(np.abs(starts+500-t)<=radius)&(starts+1000<=chrom_end)]
    z=np.load(DATA/'candidates.npz');code=int(chrom[3:])-1 if chrom not in ['chrX','chrY'] else {'chrX':22,'chrY':23}[chrom]
    keys=code*1000000+starts//1000;ix=np.searchsorted(z['ids'],keys);safe=np.minimum(ix,len(z['ids'])-1);present=(ix<len(z['ids']))&(z['ids'][safe]==keys)
    support=np.where(present,z['support'][safe],0)
    _,allprom=annotations();pos=np.sort(allprom.loc[allprom.chrom.eq(chrom),'tss'].unique())
    promoter_overlap=np.searchsorted(pos,starts+2000,'left')>np.searchsorted(pos,starts-1001,'right')
    def ov(intervals):
        out=np.zeros(len(starts),bool)
        for s,e in intervals:out|=(starts<e)&(starts+1000>s)
        return out
    cds=ov(a['cds']);utr=ov(a['utr5']);own=ov([(t-1000,t+1001)])
    c=pd.DataFrame(dict(candidate_id=[f'{chrom}:{s}-{s+1000}|{a["transcript"]}' for s in starts],enhancer_chrom=chrom,enhancer_start=starts,enhancer_end=starts+1000,promoter_chrom=chrom,tss=t,gene=g.symbol,gene_id=g.gene_id,promoter_id=a['transcript'],cell='unspecified',reference_peak_support=support,promoter_overlap=promoter_overlap,target_cds_overlap=cds,target_5utr_overlap=utr,target_promoter_overlap=own,regulatory_eligible=~(cds|utr|own),reference_candidate_supported=present))
    c['selected_tss_distance_bp']=abs(c.enhancer_start+500-t)
    a.update(radius_bp=radius,candidates=len(c),regulatory_candidates=int(c.regulatory_eligible.sum()),separate_coding_or_promoter_candidates=int((~c.regulatory_eligible).sum()),reference_supported_candidates=int(present.sum()),candidate_definition='Exhaustive 1-kb windows, not assumed validated enhancers; reference support is reported separately.')
    return c,a



def rank_regions(frame):
    x=frame.sort_values(['score','enhancer_start','candidate_id'],ascending=[False,True,True]).copy()
    x.insert(0,'rank',range(1,len(x)+1))
    return x


def query(gene,bed_dir,cache_dir,output,cell_type=None,workers=8,offline=False):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    c,a=candidates(gene,Path(cache_dir)/'gene_annotations',offline=offline)
    c.to_csv(output/'candidates.csv',index=False)
    p,r=predict_public(c,bed_dir,cache_dir,cell_type,workers,offline)
    p.to_csv(output/'all_scored_windows.csv',index=False)
    ranked=rank_regions(p[p.regulatory_eligible])
    ranked.to_csv(output/'ranked_pairs.csv',index=False)
    ranked.head(50).to_csv(output/'top_regions.csv',index=False)
    rank_regions(p[~p.regulatory_eligible]).to_csv(output/'coding_and_promoter_regions.csv',index=False)
    r.update(model='Atlas-E2G',annotation=a,ranking='Public-reference score with representative-transcript annotation',interpretation='Candidate ranking for experimental prioritization.')
    (output/'query.json').write_text(json.dumps(r,indent=2)+'\n',encoding='utf8')
    return p,r
