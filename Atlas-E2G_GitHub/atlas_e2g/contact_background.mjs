// Fixed +/-5 Mb contact background for each requested promoter; no outcomes read.
import fs from 'node:fs';import path from 'node:path';
import Straw from 'hic-straw';
const base=process.argv[3];if(!base)throw Error('Provide requests JSON and cache directory');
const cache=path.join(base,'hic_ranges'),out=path.join(base,'contact_background');fs.mkdirSync(out,{recursive:true});fs.mkdirSync(cache,{recursive:true});
const url='https://s3.us-central-1.wasabisys.com/aiden-encode-hic-mirror/bifocals_iter2/tissues.hic';
let transferred=0;const inflight=new Map();
class RangeFile{
 async read(position,length){
  if(!Number.isSafeInteger(position)||!Number.isSafeInteger(length)||position<0||length<0||length>128*1024*1024)throw Error('Invalid range');
  const key=`${position}_${length}`,dest=path.join(cache,key+'.bin');let b;
  if(fs.existsSync(dest))b=fs.readFileSync(dest);
  else{
   if(!inflight.has(key))inflight.set(key,(async()=>{
    for(let retry=0;retry<4;retry++)try{
     const r=await fetch(url,{headers:{Range:`bytes=${position}-${position+length-1}`},signal:AbortSignal.timeout(120000)});
     if(r.status!==206||r.headers.get('content-range')?.split('/')[0]!==`bytes ${position}-${position+length-1}`){await r.body?.cancel();throw Error('Invalid range response');}
     if(transferred+length>5e9){await r.body?.cancel();throw Error('Transfer limit');}
     const data=Buffer.from(await r.arrayBuffer());if(data.length!==length)throw Error('Truncated range');transferred+=length;
     const temp=dest+'.tmp'+process.pid;fs.writeFileSync(temp,data);fs.renameSync(temp,dest);return data;
    }catch(e){if(retry===3)throw e;await new Promise(resolve=>setTimeout(resolve,1000*(retry+1)));}
   })());
   b=await inflight.get(key);
  }
  if(b.length!==length)throw Error('Cached range length mismatch');return b.buffer.slice(b.byteOffset,b.byteOffset+b.byteLength);
 }
}
const straw=new Straw({file:new RangeFile()});const meta=await straw.getMetaData();
const sizes=Object.fromEntries(meta.chromosomes.map(c=>[c.name,c.size]));
const input=process.argv[2];
const requests=JSON.parse(fs.readFileSync(input,'utf8'));const unique=new Map();
for(const r of requests){const tb=Math.floor(r.tss/5000);unique.set(r.chr+'_'+tb,{chr:r.chr,tb});}
let cursor=0,count=0;const tasks=[...unique.entries()].sort((a,b)=>a[0].localeCompare(b[0]));
async function worker(){while(cursor<tasks.length){
 const [key,{chr,tb}]=tasks[cursor++],file=path.join(out,key+'.bin'),idx=path.join(out,key+'.json');
 if(fs.existsSync(file)&&fs.existsSync(idx)){count++;continue;}
 const lo=Math.max(0,tb-1001),hi=Math.min(Math.ceil(sizes[chr]/5000)-1,tb+1001);
 const p={chr,start:tb*5000,end:(tb+1)*5000};
 if(!await straw.hicFile.getNormalizationVector('SCALE',chr,'BP',5000))throw Error('Missing SCALE');
 const up=await straw.getContactRecords('SCALE',{chr,start:lo*5000,end:(tb+1)*5000},p,'BP',5000);
 const down=await straw.getContactRecords('SCALE',p,{chr,start:tb*5000,end:(hi+1)*5000},'BP',5000);
 const a=new Float32Array(hi-lo+1);
 for(const r of [...up,...down]){const j=r.bin1===tb?r.bin2:r.bin1;if(j>=lo&&j<=hi&&Number.isFinite(r.counts))a[j-lo]=r.counts;}
 fs.writeFileSync(file,Buffer.from(a.buffer));fs.writeFileSync(idx,JSON.stringify({chr,tss_bin:tb,start_bin:lo,end_bin:hi,resolution:5000,normalization:'SCALE',radius_bp:5000000}));
 if(++count%50===0)console.log('CONTACT_BACKGROUND',count,'/',tasks.length,'newMB',Math.round(transferred/1e6));
}}
await Promise.all(Array.from({length:3},()=>worker()));console.log('FINISHED_BACKGROUND',count,'newMB',Math.round(transferred/1e6));
