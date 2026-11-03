"""Atlas-E2G: candidate ranking and reproduction of the paper's analyses."""
import os
for name in ['OMP_NUM_THREADS','MKL_NUM_THREADS','OPENBLAS_NUM_THREADS']:os.environ.setdefault(name,'4')
import argparse,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
BED=Path(os.environ.get('ATLAS_E2G_BED_DIR',ROOT/'reference/chip_atlas'))
CACHE=Path(os.environ.get('ATLAS_E2G_CACHE_DIR',ROOT/'reference/public'))

def new_output(path):
    path=Path(path);path.mkdir(parents=True,exist_ok=False);return path

def main():
    p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
    a=sub.add_parser('reproduce',help='Recompute paper results from processed data')
    a.add_argument('--from-saved',action='store_true',help='Calculate metrics from supplied predictions without refitting')
    a.add_argument('--analyses',choices=['all','full','activity','transfer','source'],default='all')
    a.add_argument('--bootstrap',type=int,default=0);a.add_argument('--output',type=Path,default=Path('output/reproduction'))
    a=sub.add_parser('stability',help='Evaluate promoter stability under training-gene bootstrap refitting')
    a.add_argument('--from-saved',action='store_true');a.add_argument('--replicates',type=int,default=100)
    a.add_argument('--workers',type=int,default=4);a.add_argument('--seed',type=int,default=73421)
    a.add_argument('--output',type=Path,default=Path('output/stability'))
    a=sub.add_parser('figures',help='Draw manuscript figures as PNG, PDF and editable SVG')
    a.add_argument('--output',type=Path,default=Path('output/figures'))
    a=sub.add_parser('query-analysis',help='Matched binding ablations and transcription-supported target recovery')
    a.add_argument('--from-saved',action='store_true')
    a.add_argument('--output',type=Path,default=Path('output/query_analysis'))
    a=sub.add_parser('target-benchmark',help='Evaluate enhancer-to-gene recovery with frozen models')
    a.add_argument('--from-saved',action='store_true');a.add_argument('--bootstrap',type=int,default=2000)
    a.add_argument('--bed-dir',type=Path,default=BED);a.add_argument('--cache-dir',type=Path,default=CACHE)
    a.add_argument('--workers',type=int,default=8);a.add_argument('--offline',action='store_true')
    a.add_argument('--published-dir',type=Path);a.add_argument('--output',type=Path,default=Path('output/target_benchmark'))
    a=sub.add_parser('targets',help='Find candidate target genes from an hg38 enhancer coordinate')
    a.add_argument('--enhancer',required=True,help='chr:start-end; zero-based, half-open hg38')
    a.add_argument('--cell-type');a.add_argument('--radius',type=int,default=1000000)
    a.add_argument('--bed-dir',type=Path,default=BED);a.add_argument('--cache-dir',type=Path,default=CACHE)
    a.add_argument('--offline',action='store_true');a.add_argument('--workers',type=int,default=8)
    a.add_argument('--candidates-only',action='store_true');a.add_argument('--output',type=Path,default=Path('output/targets'))
    for command in ['query','predict']:
        a=sub.add_parser(command,help='Find candidates for a gene' if command=='query' else 'Score a CSV of candidate pairs')
        a.add_argument('--cell-type');a.add_argument('--bed-dir',type=Path,default=BED);a.add_argument('--cache-dir',type=Path,default=CACHE)
        a.add_argument('--offline',action='store_true');a.add_argument('--workers',type=int,default=8)
        a.add_argument('--output',type=Path,default=Path('output')/command)
        if command=='query':
            a.add_argument('--gene',required=True)
            a.add_argument('--candidate-set',choices=['reference','annotated'],default='reference')
            a.add_argument('--radius',type=int,default=500000);a.add_argument('--max-candidates',type=int,default=1000)
            a.add_argument('--include-promoters',action='store_true');a.add_argument('--candidates-only',action='store_true')
        else:
            a.add_argument('--input',type=Path,required=True);a.add_argument('--mode',choices=['full','public'],default='full')
    a=sub.add_parser('download',help='Download the ChIP-Atlas reference peaks')
    a.add_argument('--bed-dir',type=Path,default=BED);a.add_argument('--workers',type=int,default=8)
    sub.add_parser('cell-types',help='List available public cell references')
    a=sub.add_parser('prepare-public',help='Download and process one accessibility reference')
    a.add_argument('--cell-type',required=True);a.add_argument('--cache-dir',type=Path,default=CACHE);a.add_argument('--offline',action='store_true')
    a=sub.add_parser('reproduce-public',help='Rebuild public features and refit from public signals')
    a.add_argument('--cache-dir',type=Path,default=CACHE);a.add_argument('--offline',action='store_true');a.add_argument('--source-folds',action='store_true')
    a.add_argument('--output',type=Path,default=Path('output/public_reproduction'))
    a=sub.add_parser('rebuild-overlaps',help='Rebuild benchmark overlaps from downloaded peaks')
    a.add_argument('--bed-dir',type=Path,default=BED);a.add_argument('--workers',type=int,default=8);a.add_argument('--output',type=Path,default=Path('output/overlaps'))
    args=p.parse_args()
    if getattr(args,'workers',1)<1:p.error('--workers must be positive')
    if getattr(args,'bootstrap',0)<0:p.error('--bootstrap must be nonnegative')
    if args.command=='reproduce':
        import numpy as np,pandas as pd
        out=new_output(args.output)
        if args.analyses in ['all','full']:
            from atlas_e2g.benchmark_full import refit,summarize,compare_results
            target=new_output(out/'full')
            scores=dict(np.load(ROOT/'reference_results/full_scores.npz',allow_pickle=False)) if args.from_saved else refit(target)
            result=summarize(scores,target,args.bootstrap);expected=pd.read_csv(ROOT/'reference_results/full_metrics.csv')
            check=compare_results(result,expected)
            if len(result)!=len(expected):raise AssertionError('Incomplete full-input results')
            print(f"Full-input results reproduced: {check['checked']} metrics; largest difference {check['max_metric_difference']:.3g}.")
        if args.analyses!='full':
            from atlas_e2g.benchmark_public import run
            run(new_output(out/'public'),not args.from_saved,'transfer,activity,source' if args.analyses=='all' else args.analyses,args.bootstrap)
        print('Results:',out.resolve())
    elif args.command=='stability':
        from atlas_e2g.training_stability import run,from_saved
        if args.replicates<1:p.error('--replicates must be positive')
        out=new_output(args.output)
        if args.from_saved:from_saved(out)
        else:run(out,args.replicates,args.workers,args.seed)
    elif args.command=='figures':
        from atlas_e2g.figures import run
        run(new_output(args.output));print('Figures:',args.output.resolve())
    elif args.command=='query-analysis':
        from atlas_e2g.benchmark_queries import run
        run(new_output(args.output),args.from_saved)
    elif args.command=='target-benchmark':
        from atlas_e2g.benchmark_targets import run
        run(new_output(args.output),args.from_saved,args.bed_dir,args.cache_dir,args.workers,args.offline,args.published_dir,args.bootstrap)
    elif args.command=='targets':
        from atlas_e2g.target_query import query
        info=query(args.enhancer,args.bed_dir,args.cache_dir,args.output,args.cell_type,args.radius,args.workers,args.offline,args.candidates_only)
        print(f"Retrieved {info['candidate_genes']} candidate genes and {info['candidate_promoters']} promoters: {args.output.resolve()}")
    elif args.command in ['query','predict']:
        import pandas as pd
        from atlas_e2g.predict import predict,candidates_for_gene
        if args.command=='query':
            if args.candidate_set=='annotated':
                if args.include_promoters or args.radius!=500000 or args.max_candidates!=1000:p.error('Annotated candidates use a fixed 1 Mb radius and annotation filtering; omit reference-candidate options.')
                if args.candidates_only:
                    from atlas_e2g.regulatory_query import candidates
                    pairs,receipt=candidates(args.gene,args.cache_dir/'gene_annotations',offline=args.offline)
                else:
                    from atlas_e2g.regulatory_query import query
                    _,receipt=query(args.gene,args.bed_dir,args.cache_dir,args.output,args.cell_type,args.workers,args.offline)
                    print(f"Ranked annotated candidates for {receipt['annotation']['gene']}: {args.output.resolve()}");return
            else:
                if not 1000<=args.radius<=5000000 or not 1<=args.max_candidates<=10000:p.error('Radius: 1 kb–5 Mb; candidate cap: 1–10000')
                pairs,receipt=candidates_for_gene(args.gene,args.radius,args.max_candidates,include_promoters=args.include_promoters,cell=args.cell_type or 'unspecified')
                receipt['include_promoter_overlaps']=args.include_promoters
            mode='public'
        else:pairs=pd.read_csv(args.input);receipt={};mode=args.mode
        if args.cell_type and mode!='public':p.error('--cell-type requires --mode public')
        out=new_output(args.output);pairs.to_csv(out/'candidates.csv',index=False)
        if args.command=='query' and args.candidates_only:
            receipt.update(model='Atlas-E2G',scored=False)
            (out/'query.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8')
            print(f'Retrieved {len(pairs)} candidate pairs: {out.resolve()}');return
        if mode=='public':
            from atlas_e2g.public_atlas import predict_public
            result,info=predict_public(pairs,args.bed_dir,args.cache_dir,args.cell_type,args.workers,args.offline);receipt.update(info)
        else:result=predict(pairs,args.bed_dir,args.workers,args.offline)
        keys=['score','candidate_id'] if 'candidate_id' in result else ['score']
        result=result.sort_values(keys,ascending=[False,True] if len(keys)==2 else False,kind='stable')
        result.insert(0,'rank',range(1,len(result)+1));result.to_csv(out/'ranked_pairs.csv',index=False)
        result.drop_duplicates(['enhancer_chrom','enhancer_start','enhancer_end']).head(50).to_csv(out/'top_regions.csv',index=False)
        receipt.update(model='Atlas-E2G',pairs=len(result),genome='hg38',interpretation='Candidate ranking for experimental prioritization.')
        (out/'query.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf8');print(f'Ranked {len(result)} pairs: {out.resolve()}')
    elif args.command=='reproduce-public':
        from atlas_e2g.rebuild_public import run
        run(new_output(args.output),args.cache_dir,args.offline,args.source_folds)
    elif args.command in ['cell-types','prepare-public']:
        from atlas_e2g.public_atlas import manifest,match_cell,PublicAtlas
        if args.command=='cell-types':
            for row in manifest()['cells']:print(row['cell']+'\t'+row['file'])
        else:
            row=match_cell(args.cell_type)
            if row is None:p.error('No matching public reference; run cell-types to see available cells.')
            PublicAtlas(args.cache_dir,args.offline).grid(row).close();print('Public cell reference ready:',row['cell'])
    else:
        from atlas_e2g.predict import download,rebuild_overlaps
        if args.command=='download':download(args.bed_dir,args.workers)
        else:rebuild_overlaps(args.bed_dir,args.output,args.workers)

if __name__=='__main__':main()

