"""Fixed-reference features shared by reproduction and inference."""
from pathlib import Path
from functools import lru_cache
import json
import numpy as np
import pandas as pd
from scipy import sparse

ROOT=Path(__file__).resolve().parent.parent
DATA=ROOT/'data'
CONFIG=json.loads((ROOT/'models/config.json').read_text())
NATIVE=CONFIG['native_features']
FACTORS=CONFIG['factors']

def normalize(v,scale=2):
    return v/np.maximum(np.linalg.norm(v,axis=1)[:,None],1e-12)*scale

@lru_cache(None)
def reference():
    return pd.read_csv(DATA/'reference.tsv.gz',sep='\t',keep_default_na=False)

@lru_cache(None)
def annotations():
    return (pd.read_csv(DATA/'genes.tsv.gz',sep='\t',keep_default_na=False),
            pd.read_csv(DATA/'promoters.tsv.gz',sep='\t',keep_default_na=False))

def profiles(B):
    meta=reference();fi={x:i for i,x in enumerate(FACTORS)}
    ix=meta.factor.map(fi).to_numpy();counts=np.bincount(ix,minlength=len(FACTORS))
    mean=sparse.csr_matrix((1/counts[ix],(np.arange(len(meta)),ix)),shape=(len(meta),len(FACTORS)))
    pooled=(B@mean).toarray()
    groups=list(meta.groupby(['factor','cell'],sort=True).indices.items())
    gix=np.empty(len(meta),int);gf=[]
    for j,((factor,cell),ids) in enumerate(groups):gix[ids]=j;gf.append(fi[factor])
    gc=np.bincount(gix);gf=np.array(gf);nc=np.bincount(gf,minlength=len(FACTORS))
    G=sparse.csr_matrix((1/gc[gix],(np.arange(len(meta)),gix)),shape=(len(meta),len(groups)))
    H=sparse.csr_matrix((1/nc[gf],(np.arange(len(groups)),gf)),shape=(len(groups),len(FACTORS)))
    return pooled,(B@G@H).toarray()


def annotation_features(candidates):
    genes,tss=annotations();sym=genes.drop_duplicates('symbol').set_index('symbol');byid=genes.set_index('gene_id');ambiguous=set(genes.loc[genes.symbol.duplicated(False),'symbol']);bychr={c:np.sort(g.tss.unique()) for c,g in tss.groupby('chrom')}
    near={c:g.sort_values('tss') for c,g in tss.groupby('chrom')}
    c=candidates;mid=(c.enhancer_start.to_numpy()+c.enhancer_end.to_numpy())/2;t=c.tss.to_numpy();dist=np.abs(mid-t)
    out=np.zeros((len(c),6),float);out[:,0]=np.log(dist+.01);out[:,1]=np.log(c.enhancer_end.to_numpy()-c.enhancer_start.to_numpy()+.01)
    resolved=[]
    for i,row in enumerate(c.itertuples()):
        chrom=row.enhancer_chrom;positions=bychr.get(chrom,np.empty(0));a,b=sorted([mid[i],t[i]])
        out[i,2]=np.log1p(max(0,np.searchsorted(positions,b,'left')-np.searchsorted(positions,a,'right')))
        symbol=getattr(row,'gene','')
        if not symbol or symbol not in sym.index:
            tt=near.get(row.promoter_chrom)
            if tt is not None:
                j=np.argmin(np.abs(tt.tss.to_numpy()-row.tss))
                if abs(tt.tss.iloc[j]-row.tss)<=2000:symbol=tt.symbol.iloc[j]
        resolved.append(symbol)
        gid=str(getattr(row,'gene_id','')).split('.')[0]
        if gid in byid.index:g=byid.loc[gid]
        elif symbol in sym.index and symbol not in ambiguous:g=sym.loc[symbol]
        elif symbol in sym.index:
            # Resolve duplicate symbols by the supplied promoter coordinate, not
            # an arbitrary first annotation record. No outcome labels are used.
            tt=tss[(tss.symbol==symbol)&(tss.chrom==row.promoter_chrom)]
            if len(tt):g=byid.loc[tt.iloc[np.argmin(np.abs(tt.tss.to_numpy()-row.tss))].gene_id]
            else:g=sym.loc[symbol]
        else:continue
        out[i,3]=np.log1p(g.promoter_count);out[i,4]=np.log1p(g.end-g.start);out[i,5]=g.gene_type=='protein_coding'
    return out,resolved

def competition(c,assay_support):
    z=np.load(DATA/'candidates.npz');ids=z['ids'];support=z['support'].astype(float);bychr={}
    for code in range(24):
        ix=np.flatnonzero(ids//1000000==code);chrom='chr'+(str(code+1) if code<22 else ['X','Y'][code-22]);bychr[chrom]=((ids[ix]%1000000)*1000+500,support[ix])
    memo={};out=[]
    for i,row in enumerate(c.itertuples()):
        pos,w=bychr[row.promoter_chrom];t=row.tss;mid=(row.enhancer_start+row.enhancer_end)/2;dist=abs(mid-t);key=(row.promoter_chrom,t)
        if key not in memo:
            masses=[]
            for radius in [500000,5000000]:
                a,b=np.searchsorted(pos,[t-radius,t+radius]);dd=np.abs(pos[a:b]-t);masses.append(float(np.sum(w[a:b]/np.maximum(dd,2500))))
            memo[key]=masses
        m1,m5=memo[key];a,b=np.searchsorted(pos,[mid-5000,mid+5000]);near=float(w[a:b].sum());u,v=np.searchsorted(pos,sorted([mid,t]));between=v-u
        activity=assay_support[i]/max(dist,2500)
        out.append([np.log(activity/max(m1,1e-8)+1e-8),np.log(activity/max(m5,1e-8)+1e-8),np.log1p(m1),np.log1p(m5),np.log1p(near),np.log1p(between)])
    return np.array(out)

def coordinates(d):
    return d.rename(columns={'chrom':'enhancer_chrom','start':'enhancer_start','end':'enhancer_end','pchrom':'promoter_chrom'})

def pair_features(c,B,ei,pi):
    pooled,balanced=profiles(B)
    base,_=annotation_features(c)
    amounts=np.column_stack([np.log1p(pooled[ids].sum(axis=1)) for ids in [ei,pi]]+
                            [np.log1p((pooled[ids]>0).sum(axis=1)) for ids in [ei,pi]])
    return {'base':np.column_stack([base,amounts]),
            'pool':np.column_stack([normalize(pooled[ei]),normalize(pooled[pi])]),
            'balanced':np.column_stack([normalize(balanced[ei]),normalize(balanced[pi])]),
            'neighborhood':competition(c,np.asarray(B[ei].sum(axis=1)).ravel())}

def load(kind):
    d=pd.read_csv(DATA/f'{kind}.csv.gz')
    B=sparse.load_npz(DATA/f'{kind}_overlap.npz')
    if kind=='reporter':
        pooled,_=profiles(B);raw=pooled[d.uid]
        return d,{'profile':normalize(raw),'amount':np.c_[np.log1p(raw.sum(axis=1)),np.log1p(np.linalg.norm(raw,axis=1))]}
    return d,pair_features(coordinates(d),B,d.enhancer_uid.to_numpy(),d.promoter_uid.to_numpy())

def native(d,second=False):
    cols=[x+'_second' for x in NATIVE] if second else NATIVE
    return np.log(d[cols].fillna(0).abs().to_numpy()+.01)

def public_design(blocks):
    base=blocks['base'].astype(np.float32)
    identity=blocks['balanced'].astype(np.float32)
    neighborhood=blocks['neighborhood'].astype(np.float32)
    return np.column_stack([base,neighborhood]),identity
