"""Portable, availability-aware public-reference inference. No benchmark labels are read."""
from pathlib import Path
from functools import lru_cache
import hashlib,json,re,subprocess,shutil,tempfile
import numpy as np,pandas as pd
from scipy import sparse
from scipy.special import expit
from .features import ROOT,DATA,FACTORS,annotations,public_design,pair_features
from .predict import make_regions,overlaps,ensure_reference

EPS=1e-4
CHROMS=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']

def manifest():return json.loads((DATA/'public_atlas.json').read_text(encoding='utf-8'))

def match_cell(cell):
    norm=lambda s:re.sub('[^a-z0-9]','',str(s).casefold())
    name=norm(cell);name={'astrocytes':'astrocyte'}.get(name,name)
    return next((r for r in manifest()['cells'] if norm(r['cell'])==name),None)

def static_features(d):
    gene=pd.read_csv(DATA/'public_gene_class.tsv.gz',sep='\t').drop_duplicates('TargetGene').set_index('TargetGene')
    t=pd.read_csv(DATA/'public_refseq_tss.bed.gz',sep='\t',header=None,names=['chrom','start','end','gene','v','strand'])
    by={k:(np.sort(v.start.to_numpy()),np.sort(v.end.to_numpy())) for k,v in t.groupby('chrom')};out=np.zeros((len(d),4))
    for i,r in enumerate(d.itertuples()):
        a,b=sorted(((r.start+r.end)/2,r.tss));st,en=by.get(r.chrom,(np.array([]),np.array([])))
        out[i,0]=np.log1p(np.searchsorted(st,b,'left')-np.searchsorted(en,a,'right'))
        if r.gene in gene.index:out[i,1]=float(gene.loc[r.gene,'is_ubiquitous_uniform']);out[i,2]=float(gene.loc[r.gene,'P2PromoterClass'])
        else:out[i,3]=1
    return out

def public_design_matrix(d,blocks,contact):
    base,identity=public_design(blocks);co=np.nan_to_num(contact[:,0],nan=0)
    return np.c_[base,static_features(d),np.log1p(np.maximum(co,0)),(np.isfinite(contact[:,0])&(contact[:,1]>0)).astype(float)],identity

class PublicAtlas:
    def __init__(self,cache_dir,offline=False):
        self.cache=Path(cache_dir);self.cache.mkdir(parents=True,exist_ok=True);self.offline=offline;self.validated=set()
        self.contacts=self.cache/'contact_background';self.contacts.mkdir(exist_ok=True)
        self.grids=self.cache/'dnase';self.grids.mkdir(exist_ok=True)
        z=np.load(DATA/'candidates.npz');self.ids=z['ids'];_,promoters=annotations();self.positions={}
        for code,chrom in enumerate(CHROMS):
            ix=np.flatnonzero(self.ids//1000000==code);start=(self.ids[ix]%1000000)*1000;pos=np.sort(promoters.loc[promoters.chrom==chrom,'tss'].unique())
            excluded=np.searchsorted(pos,start+2000,'left')>np.searchsorted(pos,start-1001,'right');ix=ix[~excluded]
            self.positions[chrom]=((self.ids[ix]%1000000)*1000+500,ix)

    def ensure_contact(self,d):
        req=[dict(chr=c,tss=int(t)) for c,t in d[['chrom','tss']].drop_duplicates().itertuples(index=False,name=None)]
        missing=[r for r in req if not all((self.contacts/(r['chr']+'_'+str(r['tss']//5000)+ext)).exists() for ext in ['.bin','.json'])]
        if not missing:return
        if self.offline:raise FileNotFoundError('Public contacts are not cached. Rerun without --offline to retrieve the required ranges.')
        node=shutil.which('node')
        if not node:raise RuntimeError('Public contact retrieval requires Node.js 20.19+ or 22.12+, and npm ci in the Atlas-E2G directory.')
        with tempfile.TemporaryDirectory(prefix='atlas_e2g_',dir=self.cache) as tmp:
            p=Path(tmp)/'requests.json';p.write_text(json.dumps(missing),encoding='utf-8')
            subprocess.run([node,str(ROOT/'atlas_e2g/contact_background.mjs'),str(p),str(self.cache)],check=True,cwd=ROOT)

    @lru_cache(None)
    def contact(self,chrom,tss_bin):
        key=f'{chrom}_{tss_bin}';meta=json.loads((self.contacts/(key+'.json')).read_text());p=self.contacts/(key+'.bin')
        a=np.fromfile(p,dtype='<f4')
        if len(a)!=meta['end_bin']-meta['start_bin']+1 or not np.isfinite(a).all() or (a<0).any():raise ValueError('Invalid contact cache: '+key)
        return meta['start_bin'],a

    def pair_contact(self,d):
        self.ensure_contact(d);out=np.zeros((len(d),2))
        for i,r in enumerate(d.itertuples()):
            start,a=self.contact(r.chrom,int(r.tss)//5000);j=int((r.start+r.end)//2)//5000-start
            if j<0 or j>=len(a):raise ValueError('Public-input pairs must be within 5 Mb of the promoter')
            v=float(a[j]);out[i]=v,int(v>0)
        return out

    def grid(self,row):
        path=self.grids/(row['file']+'.npz')
        if not path.exists():
            if self.offline:raise FileNotFoundError('DNase track '+row['file']+' is not cached; rerun without --offline.')
            from .public_download import prepare_grid
            prepare_grid(row,self.grids)
        if row['file'] not in self.validated:
            sizes=dict(line.split() for line in (DATA/'hg38.primary.chrom.sizes').read_text().splitlines())
            with np.load(path,allow_pickle=False) as grid:
                for chrom in CHROMS:
                    v=grid[chrom]
                    if v.ndim!=1 or len(v)!=(int(sizes[chrom])+999)//1000 or not np.isfinite(v).all() or (v<0).any():
                        raise ValueError('Invalid public DNase grid: '+row['file'])
            self.validated.add(row['file'])
        return np.load(path,allow_pickle=False)

    def activity_features(self,d,contact,row):
        grid=self.grid(row);regions,ei,pi=make_regions(d.rename(columns={'chrom':'enhancer_chrom','start':'enhancer_start','end':'enhancer_end','pchrom':'promoter_chrom'}));raw=np.zeros(len(regions))
        for i,r in enumerate(regions.itertuples(index=False)):
            arr=grid[r.chrom];bins=np.arange(int(r.start)//1000,(int(r.end)-1)//1000+1);valid=(bins>=0)&(bins<len(arr));length=np.maximum(0,np.minimum((bins+1)*1000,r.end)-np.maximum(bins*1000,r.start))
            raw[i]=arr[bins[valid]]@length[valid]/(r.end-r.start)
        obs=np.log1p(raw/row['scale']);vals=np.zeros(len(self.ids),np.float32)
        for code,chrom in enumerate(CHROMS):
            ii=np.flatnonzero(self.ids//1000000==code);bins=self.ids[ii]%1000000;arr=grid[chrom];valid=bins<len(arr);vals[ii[valid]]=arr[bins[valid]]
        vals=np.log1p(vals/row['scale']);memo={};result=np.zeros((len(d),5))
        for i,r in enumerate(d.itertuples()):
            key=(r.chrom,int(r.tss));pos,idx=self.positions[r.chrom]
            if key not in memo:
                left=np.searchsorted(pos,r.tss-5000000,'left');right=np.searchsorted(pos,r.tss+5000000,'right');start,h=self.contact(r.chrom,int(r.tss)//5000)
                memo[key]=max(float((vals[idx[left:right]]+EPS)@h[(pos[left:right]//5000).astype(int)-start]),1e-12)
            den=memo[key];mid=(r.start+r.end)/2;left,right=np.searchsorted(pos,[mid-5000,mid+5000]);near=vals[idx[left:right]].sum();ae=obs[ei[i]];ap=obs[pi[i]]
            result[i]=[np.log1p(ae),np.log1p(ap),np.log((ae+EPS)*contact[i,0]/den+1e-12),np.log1p(near),np.log1p(den)]
        grid.close()
        if not np.isfinite(result).all():raise ValueError('Invalid public activity features')
        return result

def score_arrays(low,identity,activity=None):
    z=np.load(ROOT/'models/public_reference.npz',allow_pickle=False);name='cell' if activity is not None else 'gene'
    x=np.c_[low,activity] if activity is not None else low;x=(x-z[name+'_mean'])/z[name+'_scale'];coef=z[name+'_coef'];n=len(x[0]);nf=len(FACTORS)
    context=x@coef[:n];enhancer=identity[:,:nf]@coef[n:n+nf];promoter=identity[:,nf:]@coef[n+nf:]
    value=expit(context+enhancer+promoter+float(z[name+'_intercept']))
    return value,dict(public_context_component=context,enhancer_component=enhancer,promoter_component=promoter)

def score_retained_component(low,identity,pooled_promoter,activity=None):
    """Evaluate the fixed source promoter term and its public-input adaptation."""
    name='cell' if activity is not None else 'gene';nf=len(FACTORS)
    with np.load(ROOT/'models/public_components.npz',allow_pickle=False) as z:
        if not np.array_equal(z['factors'],np.asarray(FACTORS)):raise ValueError('Component factor order does not match the reference')
        x=np.c_[low,activity] if activity is not None else low
        x=(x-z[name+'_mean'])/z[name+'_scale'];coef=z[name+'_coef'];n=x.shape[1]
        promoter=pooled_promoter@z[name+'_promoter_coef']
        score=expit(x@coef[:n]+identity[:,:nf]@coef[n:]+promoter+float(z[name+'_intercept']))
    return dict(retained_promoter_component=promoter,component_preserving_score=score)


def predict_public(c,bed_dir,cache_dir,cell_type=None,workers=8,offline=False):
    c=c.reset_index(drop=True)
    if 'gene' not in c:
        from .features import annotation_features
        c=c.copy();c['gene']=annotation_features(c)[1]
    if not c.enhancer_chrom.eq(c.promoter_chrom).all():raise ValueError('Only same-chromosome candidate pairs are supported')
    regions,ei,pi=make_regions(c)
    midpoint=(c.enhancer_start.to_numpy(float)+c.enhancer_end.to_numpy(float))/2
    if np.any(np.abs(midpoint-c.tss.to_numpy(float))>5000000):
        raise ValueError('Public-input pairs must be within 5 Mb of the promoter. Interchromosomal and more distant interactions are outside the supported search.')
    cache=Path(cache_dir);folder=cache/'overlaps';folder.mkdir(parents=True,exist_ok=True)
    key=hashlib.blake2b(regions.to_csv(index=False).encode()+(DATA/'reference.tsv.gz').read_bytes(),digest_size=16).hexdigest()
    path=folder/(key+'.npz')
    if path.exists():
        B=sparse.load_npz(path)
        if B.shape!=(len(regions),len(pd.read_csv(DATA/'reference.tsv.gz',sep='\t'))) or not np.isin(B.data,[0,1]).all():
            raise ValueError('Invalid cached ChIP-Atlas overlaps')
        peak_info=dict(experiments=B.shape[1],downloaded_files=0,acquisition='cached_query_overlaps')
    else:
        peak_info=ensure_reference(bed_dir,workers,offline)
        B=overlaps(regions,bed_dir,workers);sparse.save_npz(path,B)
    blocks=pair_features(c,B,ei,pi);d=c.rename(columns={'enhancer_chrom':'chrom','enhancer_start':'start','enhancer_end':'end','promoter_chrom':'pchrom'})
    atlas=PublicAtlas(cache,offline);contact=atlas.pair_contact(d);low,identity=public_design_matrix(d,blocks,contact);gene_score,gparts=score_arrays(low,identity)
    row=match_cell(cell_type) if cell_type else None;out=c.copy();out['gene_only_score']=gene_score
    if row:
        activity=atlas.activity_features(d,contact,row);out['score'],parts=score_arrays(low,identity,activity)
        out['public_enhancer_activity']=np.expm1(activity[:,0]);out['public_promoter_activity']=np.expm1(activity[:,1]);out['public_activity_contact_score']=np.exp(activity[:,2])
        mode='public_cell_reference';reason=None
    else:out['score']=gene_score;parts=gparts;mode='public_gene_reference';reason='cell_not_specified' if not cell_type else 'no_exact_cell_reference'
    for k,v in parts.items():out[k]=v
    retained=score_retained_component(low,identity,blocks['pool'][:,len(FACTORS):],activity if row else None)
    for k,v in retained.items():out[k]=v
    out['public_contact']=contact[:,0];out['prediction_mode']=mode
    meta=dict(input_mode=mode,requested_cell_type=cell_type,matched_cell_type=row['cell'] if row else None,fallback_reason=reason,public_dnase_file=row['file'] if row else None,public_dnase_url=row['url'] if row else None,query_cell_measurement_used=bool(row),user_assays_required=False,imputed_activity_used=False,contact_source=manifest()['contact'],competition_radius_bp=5000000,activity_grid_bp=1000,reference_cells_available=len(manifest()['cells']))
    meta['component_preserving_output']='Separate fixed-source promoter contribution and adapted pair score; default score and ranking use the unconstrained public-input fit.'
    meta['chip_atlas_reference']=peak_info
    meta['candidate_scope']='Contact evidence scores supplied candidates; it does not retrieve pairs outside the query radius. Only same-chromosome pairs within 5 Mb are supported.'
    return out,meta
