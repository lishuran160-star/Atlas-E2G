"""Prediction from candidate pairs or a gene name using the fixed public reference."""
from pathlib import Path
import concurrent.futures,difflib,urllib.request,os
import numpy as np,pandas as pd
from scipy import sparse
from scipy.special import expit
from .features import ROOT,DATA,NATIVE,FACTORS,reference,annotations,pair_features
COORDINATES=['enhancer_chrom','enhancer_start','enhancer_end','promoter_chrom','tss']



def make_regions(candidates):
    c=candidates
    for name in COORDINATES:
        if name not in c:raise ValueError(f'Missing column: {name}')
    for name in ['enhancer_start','enhancer_end','tss']:
        v=pd.to_numeric(c[name],errors='raise').to_numpy()
        if not np.isfinite(v).all() or not np.equal(v,np.floor(v)).all():raise ValueError(f'Noninteger coordinate: {name}')
    a=c.enhancer_start.to_numpy(dtype=np.int64);b=c.enhancer_end.to_numpy(dtype=np.int64);t=c.tss.to_numpy(dtype=np.int64)
    if np.any(a<0) or np.any(b<=a) or np.any(t<1000):raise ValueError('Invalid interval or promoter window extends below zero')
    valid={f'chr{i}' for i in range(1,23)}|{'chrX','chrY'}
    if not set(c.enhancer_chrom).issubset(valid) or not set(c.promoter_chrom).issubset(valid):raise ValueError('Use canonical hg38 chr1..chr22,chrX,chrY')
    # Benchmarks may annotate a TSS interval rather than a single base. Preserve
    # its exact reference window when supplied; gene queries use point TSSs.
    supplied=['promoter_start' in c, 'promoter_end' in c]
    if any(supplied) and not all(supplied):raise ValueError('Supply both promoter_start and promoter_end')
    ps,pe=t-1000,t+1001
    if all(supplied):
        ps=pd.to_numeric(c.promoter_start,errors='raise').to_numpy()
        pe=pd.to_numeric(c.promoter_end,errors='raise').to_numpy()
        if not np.isfinite(ps).all() or not np.isfinite(pe).all() or not np.equal(ps,np.floor(ps)).all() or not np.equal(pe,np.floor(pe)).all() or np.any(ps<0) or np.any(pe<=ps):raise ValueError('Invalid promoter interval')
        ps,pe=ps.astype(np.int64),pe.astype(np.int64)
    e=list(zip(c.enhancer_chrom,a,b));p=list(zip(c.promoter_chrom,ps,pe))
    regions=pd.DataFrame(sorted(set(e)|set(p)),columns=['chrom','start','end']);lookup={v:i for i,v in enumerate(regions.itertuples(index=False,name=None))}
    return regions,np.array([lookup[v] for v in e]),np.array([lookup[v] for v in p])

def atlas_features(regions,reference,factors,bed_dir,workers=8,verify=False):
    bed_dir=Path(bed_dir);bychr={c:(g.start.to_numpy(),g.end.to_numpy(),g.index.to_numpy()) for c,g in regions.groupby('chrom',sort=False)}
    fi={f:i for i,f in enumerate(factors)};factor_ix=np.array([fi[f] for f in reference.factor]);counts=np.bincount(factor_ix,minlength=len(factors));v=np.zeros((len(regions),len(factors)),dtype=np.float64)
    def intersect(row):
        path=bed_dir/(row.accession+'.05.bed')
        if not path.exists():raise FileNotFoundError(f'Missing fixed reference BED: {path}')
        if verify and path.stat().st_size!=row.size_bytes:raise ValueError(f'Incomplete reference BED: {path}')
        peaks=pd.read_csv(path,sep='\t',header=None,usecols=[0,1,2],names=['chrom','start','end'],comment='#')
        hits=[]
        for chrom,g in peaks.groupby('chrom',sort=False):
            if chrom not in bychr:continue
            s=g.start.to_numpy(dtype=np.int64);e=g.end.to_numpy(dtype=np.int64);order=np.argsort(s,kind='stable');s=s[order];e=np.maximum.accumulate(e[order]);qs,qe,qi=bychr[chrom]
            j=np.searchsorted(s,qe,side='left')-1;keep=(j>=0)&(e[np.maximum(j,0)]>qs);hits.append(qi[keep])
        return np.concatenate(hits) if hits else np.empty(0,dtype=np.int64)
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for i,ids in enumerate(pool.map(intersect,reference.itertuples(index=False))):
            v[ids,factor_ix[i]]+=1.0/counts[factor_ix[i]]
    return v

def resolve_gene(query):
    genes,promoters=annotations();query=query.strip()
    if query.upper().startswith('ENSG'):g=genes[genes.gene_id.eq(query.upper().split('.')[0])]
    else:g=genes[genes.symbol.str.casefold().eq(query.casefold())]
    if len(g)==0:
        near=difflib.get_close_matches(query.upper(),genes.symbol.unique(),n=4)
        raise ValueError('Gene not found in the hg38 GENCODE annotation. Nearby names: '+', '.join(near))
    if len(g)!=1:raise ValueError('Ambiguous symbol; supply one of these Ensembl IDs: '+', '.join(g.gene_id))
    row=g.iloc[0];p=promoters[promoters.gene_id.eq(row.gene_id)].copy()
    if not len(p):raise ValueError('No supported annotated promoter for this gene')
    return row,p

def candidates_for_gene(query,radius=500000,max_candidates=1000,include_promoters=False,cell='unspecified'):
    gene,promoters=resolve_gene(query);allgenes,allprom=annotations();chrom=gene.chrom
    code=int(chrom[3:])-1 if chrom not in ['chrX','chrY'] else {'chrX':22,'chrY':23}[chrom]
    z=np.load(DATA/'candidates.npz');ids=z['ids'];support=z['support'];found=[]
    for p in promoters.itertuples():
        low=max(0,p.tss-radius)//1000;high=(p.tss+radius)//1000
        left=np.searchsorted(ids,code*1000000+low);right=np.searchsorted(ids,code*1000000+high,side='right')
        found.extend(range(left,right))
    ix=np.unique(found).astype(int)
    bins=pd.DataFrame(dict(start=(ids[ix]%1000000)*1000,end=(ids[ix]%1000000)*1000+1000,support=support[ix]))
    positions=np.sort(allprom.loc[allprom.chrom.eq(chrom),'tss'].unique())
    # A distal candidate does not overlap +/-1 kb of any annotated promoter.
    bins['promoter_overlap']=np.searchsorted(positions,bins.end.to_numpy()+1000,'left')>np.searchsorted(positions,bins.start.to_numpy()-1001,'right')
    retrieved=len(bins)
    if not include_promoters:bins=bins[~bins.promoter_overlap]
    eligible=len(bins);bins=bins.sort_values(['support','start'],ascending=[False,True]).head(max_candidates).sort_values('start')
    rows=[]
    for e in bins.itertuples():
        for p in promoters.itertuples():
            if abs((e.start+e.end)/2-p.tss)<=radius:
                rows.append(dict(candidate_id=f'{chrom}:{e.start}-{e.end}|{p.promoter_id}',enhancer_chrom=chrom,enhancer_start=int(e.start),enhancer_end=int(e.end),promoter_chrom=chrom,tss=int(p.tss),gene=gene.symbol,gene_id=gene.gene_id,promoter_id=p.promoter_id,cell=cell,reference_peak_support=int(e.support),promoter_overlap=bool(e.promoter_overlap)))
    if not rows:raise ValueError('No eligible public candidate bins in the specified search interval')
    receipt=dict(gene=gene.symbol,gene_id=gene.gene_id,chromosome=chrom,genome='hg38',radius_bp=radius,annotated_promoter_clusters=len(promoters),public_bins_retrieved=retrieved,distal_or_requested_bins=eligible,retained_bins=len(bins),omitted_by_candidate_cap=max(0,eligible-len(bins)),pairs=len(rows),candidate_definition='1 kb intervals supported by peak midpoints in at least 3 reference experiments; these are putative regulatory regions, not validated enhancers.')
    return pd.DataFrame(rows),receipt

def overlaps(regions,bed_dir,workers=8):
    ref=reference().copy();ref['factor']=ref.accession
    return sparse.csr_matrix(atlas_features(regions,ref,ref.accession.tolist(),bed_dir,workers,verify=True))

def predict(c,bed_dir,workers=8,offline=False):
    c=c.reset_index(drop=True)
    if not (c.enhancer_chrom==c.promoter_chrom).all():raise ValueError('Only cis pairs are supported')
    regions,ei,pi=make_regions(c)
    if any(x not in c for x in NATIVE):raise ValueError('Full inputs require: '+', '.join(NATIVE))
    raw=c[NATIVE].apply(pd.to_numeric,errors='raise').to_numpy(float)
    if not np.isfinite(raw).all():raise ValueError('Full-input predictions require eight finite covariates per pair')
    ensure_reference(bed_dir,workers,offline)
    B=overlaps(regions,bed_dir,workers);blocks=pair_features(c,B,ei,pi)
    out=c.copy()
    w=np.load(ROOT/'models/full.npz',allow_pickle=False)
    x=(np.log(np.abs(raw)+.01)-w['mean'])/w['scale'];a=blocks['pool'];nf=len(FACTORS)
    out['enhancer_component']=a[:,:nf]@w['atlas_coef'][8:8+nf]
    out['promoter_component']=a[:,nf:]@w['atlas_coef'][8+nf:8+2*nf]
    context=(blocks['neighborhood']-w['context_mean'])/w['context_scale']
    out['native_component']=x@w['atlas_coef'][:8]
    out['neighborhood_component']=context@w['atlas_coef'][-6:]
    out['score']=expit(out[['native_component','enhancer_component','promoter_component','neighborhood_component']].sum(axis=1)+w['atlas_intercept'])
    return out

def ensure_reference(bed_dir,workers=8,offline=False):
    """Reuse fixed peaks and acquire missing files before scoring a new query."""
    bed_dir=Path(bed_dir);rows=list(reference().itertuples(index=False))
    missing=[]
    for row in rows:
        path=bed_dir/(row.accession+'.05.bed')
        if not path.exists():missing.append(row)
        elif path.stat().st_size!=row.size_bytes:
            raise ValueError('Existing reference size mismatch; move the incomplete file aside and retry: '+str(path))
    info=dict(experiments=len(rows),downloaded_files=len(missing),
              acquisition='downloaded_missing_peaks' if missing else 'local_peaks')
    if not missing:return info
    if offline:
        raise FileNotFoundError(f'{len(missing)} ChIP-Atlas peak files are not cached in {bed_dir}. Rerun without --offline, or use --bed-dir to select an existing reference.')
    bed_dir.mkdir(parents=True,exist_ok=True)
    size=sum(int(row.size_bytes) for row in missing)
    print(f'Preparing ChIP-Atlas reference: {len(missing)} missing files ({size/1e9:.2f} GB) in {bed_dir}. Completed files are reused.',flush=True)
    def one(row):
        out=bed_dir/(row.accession+'.05.bed')
        tmp=out.with_suffix('.part')
        last=None
        for attempt in range(3):
            try:
                with urllib.request.urlopen(row.url,timeout=120) as src,tmp.open('wb') as dst:
                    while chunk:=src.read(1024*1024):dst.write(chunk)
                if tmp.stat().st_size!=row.size_bytes:raise ValueError('Incomplete reference download')
                os.replace(tmp,out);return
            except Exception as exc:last=exc
        raise RuntimeError(f'{row.accession}: {last}')
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        for i,_ in enumerate(pool.map(one,missing),1):
            if i%100==0 or i==len(missing):print(f'ChIP-Atlas reference {i}/{len(missing)} downloaded',flush=True)
    return info


def download(bed_dir,workers=8):
    return ensure_reference(bed_dir,workers)

def rebuild_overlaps(bed_dir,output,workers=8):
    output=Path(output);output.mkdir(parents=True,exist_ok=False)
    for kind in ['source','fivecell','astrocyte','reporter']:
        regions=pd.read_csv(DATA/f'{kind}_regions.csv.gz')
        actual=overlaps(regions,bed_dir,workers)
        expected=sparse.load_npz(DATA/f'{kind}_overlap.npz')
        if (actual!=expected).nnz:raise AssertionError(f'{kind}: rebuilt overlaps differ')
        sparse.save_npz(output/f'{kind}_overlap.npz',actual)
        print(f'{kind}: exact overlap reproduction',flush=True)
