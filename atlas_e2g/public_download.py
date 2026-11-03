"""Retrieve a public DNase track and integrate its signal into 1 kb windows."""
from pathlib import Path
from urllib.parse import urlsplit,urlunsplit
import concurrent.futures,json,threading,time,os
import numpy as np
CHROMS=['chr'+str(i) for i in range(1,23)]+['chrX','chrY']

def download(row,path):
    import requests
    tmp=path.with_suffix('.part');state=path.with_suffix('.download.json');block=8*1024*1024
    with requests.get(row['url'],headers={'Range':'bytes=0-15'},stream=True,timeout=45) as r:
        r.raise_for_status();u=urlsplit(r.url);url=urlunsplit((u.scheme,u.netloc,u.path,'',''))
    done=set(json.loads(state.read_text())['done']) if state.exists() and tmp.exists() else set()
    with tmp.open('ab') as f:f.truncate(row['file_size'])
    lock=threading.Lock()
    def checkpoint():
        target=state.with_suffix('.tmp');target.write_text(json.dumps(dict(done=sorted(done),block_bytes=block,file_size=row['file_size'])));os.replace(target,state)
    def part(i):
        a=i*block;b=min(row['file_size'],a+block)-1
        for retry in range(4):
            try:
                with requests.get(url,headers={'Range':f'bytes={a}-{b}'},stream=True,timeout=(30,75)) as r:
                    if r.status_code!=206 or r.headers.get('Content-Range','').split('/')[0]!=f'bytes {a}-{b}':raise ValueError('Public DNase range response mismatch')
                    n=0
                    with tmp.open('r+b') as f:
                        f.seek(a)
                        for chunk in r.iter_content(256*1024):f.write(chunk);n+=len(chunk)
                    if n!=b-a+1:raise ValueError('Truncated DNase download')
                with lock:done.add(i);checkpoint()
                return
            except Exception:
                if retry==3:raise
                time.sleep(retry+1)
    checkpoint();todo=[i for i in range((row['file_size']+block-1)//block) if i not in done]
    print('Retrieving public DNase:',row['cell'],row['file'],flush=True)
    with concurrent.futures.ThreadPoolExecutor(max_workers=6) as pool:list(pool.map(part,todo))
    if tmp.stat().st_size!=row['file_size']:raise ValueError('Incomplete DNase download: '+row['file'])
    os.replace(tmp,path)

def prepare_grid(row,directory):
    import pybigtools
    directory=Path(directory);directory.mkdir(parents=True,exist_ok=True);raw=directory/(row['file']+'.bigWig');dest=directory/(row['file']+'.npz')
    if not raw.exists():download(row,raw)
    if raw.stat().st_size!=row['file_size']:raise ValueError('Incomplete DNase file')
    bw=pybigtools.open(str(raw));sizes=bw.chroms()
    if sizes.get('chr1')!=248956422:raise ValueError('Public DNase track is not the expected hg38 assembly')
    arrays={}
    for chrom in CHROMS:
        size=sizes[chrom];n=size//1000;v=bw.values(chrom,0,n*1000,bins=n,summary='mean',exact=True,uncovered=0.,fillna=0.).astype(np.float32)
        if size%1000:
            tail=bw.values(chrom,n*1000,size,bins=1,summary='mean',exact=True,uncovered=0.,fillna=0.)[0];v=np.r_[v,np.float32(tail*(size%1000)/1000)]
        if not np.isfinite(v).all() or (v<0).any():raise ValueError('Invalid DNase values')
        ix=min(n-1,12345);exact=sum((min(e,(ix+1)*1000)-max(s,ix*1000))*value for s,e,value in bw.records(chrom,ix*1000,(ix+1)*1000))/1000
        if not np.isclose(v[ix],exact,rtol=3e-5,atol=1e-5):raise ValueError('DNase integration check failed')
        arrays[chrom]=v
    bw.close();tmp=dest.with_suffix('.tmp')
    with tmp.open('wb') as f:np.savez_compressed(f,**arrays)
    os.replace(tmp,dest)
    # Remove only this verified, regenerable raw download after the grid is durable.
    raw.unlink();print('Public DNase ready:',row['cell'],flush=True)
