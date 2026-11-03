# Using Atlas-E2G

## Start with a gene

After installing the dependencies, specify the target. Missing [public references](reference_download.md) are prepared automatically:

```sh
python run.py query --gene MYC
python run.py query --gene MYC --cell-type K562 --output output/MYC_K562
python run.py cell-types
```

The gene can be a human symbol or an Ensembl gene ID. The default search finds reference-supported 1 kb regions within 500 kb of annotated transcription start sites (TSSs), excludes promoter-overlapping regions and retains up to 1,000 regions. Each annotated promoter is scored separately. `--radius`, `--max-candidates` and `--include-promoters` adjust this search.

A recognized cell name selects measured public accessibility for that cell. If no matching reference exists, the model returns the gene-only ranking and records the reason in `query.json`. It does not infer a new cell state. Use `--offline` when all needed public files are cached.

To inspect coordinates before scoring:

```sh
python run.py query --gene MYC --candidates-only --output output/MYC_candidates
```

## Choose candidates around a representative transcript

```sh
python run.py query --gene SIN3A --candidate-set annotated --output output/SIN3A
```

This uses MANE Select, or the Ensembl canonical transcript, and scores every 1 kb window within 1 Mb of its TSS. Windows overlapping the target gene's coding exons, selected transcript's 5′ UTR or selected promoter are placed in `coding_and_promoter_regions.csv`. Other genes' promoter overlaps remain eligible and are flagged.

The scoring model is the same. This option changes the candidate set and transcript annotation, not the fitted parameters. Its representative transcript is an annotation prior rather than a measurement of promoter use in the requested cell.

## Start with an enhancer coordinate

~~~sh
python run.py targets --enhancer chr8:127735000-127736000
python run.py targets --enhancer chr8:127735000-127736000 --cell-type K562 --output output/enhancer_K562
python run.py targets --enhancer chr8:127735000-127736000 --candidates-only --output output/enhancer_candidates
~~~

Coordinates are hg38, zero-based and half-open: the first base is included and the end is excluded. The interval is the region to score; it need not match a stored reference bin. The search finds all GENCODE promoter clusters on the same chromosome within 1 Mb of the interval midpoint. `--radius` accepts 1 kb to 5 Mb. It does not filter by experimental gene lists, gene type or expression, and has no candidate cap. Promoters whose feature windows exceed chromosome bounds are excluded.

A cell name selects the same measured public accessibility reference used by gene queries. Leaving it out uses the frozen public model without cell-specific activity. The software acquires missing public references and constructs the features automatically.

- `candidates.csv`: every candidate enhancer–promoter pair, with gene ID, gene type, transcript and distance.
- `ranked_promoters.csv`: all promoter scores and their component contributions.
- `ranked_genes.csv`: one row per gene, ordered by `priority_score`, which combines promoter pair scores with public CAGE usage and gene detection support. It includes the representative promoter, number of promoters considered and distance to the nearest candidate TSS.
- `top_genes.csv`: the first 50 rows for convenient inspection; the complete ranking remains in `ranked_genes.csv`.
- `query.json`: search settings and the public references used.

Ties use stable gene/promoter identifiers. `representative_promoter_id` identifies the promoter contributing most to the gene priority; component columns describe that promoter's raw pair score. `raw_max_pair_score` preserves the highest unweighted pair score for comparison. `gene_only_priority_score` uses the no-cell model with pooled CAGE support. None of these gene summaries is a calibrated probability of regulation.

Unmeasured candidate genes remain unknown. A high score is a priority for follow-up experiments, not evidence that a regulatory effect has already been established. The [target-recovery evaluation](reproducibility.md#public-input-ablations-and-enhancer-to-gene-recovery) reports how ranks change when candidate genes extend beyond those tested in the screens.

## Longer-range interactions

The search radius is a computational boundary, not a maximum biological reach. Genomic distance is measured along the DNA sequence; chromatin folding can bring distant regions together. Atlas-E2G uses public Hi-C contact values when scoring a candidate, but does not add a candidate outside the search window because it has a contact.

For a broader search with the same fitted model:

```sh
python run.py query --gene MYC --radius 5000000 --max-candidates 10000 --output output/MYC_5Mb
python run.py targets --enhancer chr8:127735000-127736000 --radius 5000000 --output output/enhancer_5Mb
```

Add `--cell-type` if an appropriate DNase reference exists. Gene-query `query.json` reports `omitted_by_candidate_cap`; increasing the radius without increasing the cap can still omit reference-supported regions. The annotation-aware `--candidate-set annotated` mode retains its fixed 1 Mb search.

The ±5 Mb contact background and the candidate search serve different purposes: the background defines contact and activity normalization, while the search determines which pairs receive scores. The enhancer-query evaluation in the paper covered 188/189 known-positive links at the default 1 Mb radius; the remaining link was 1.97 Mb away. This illustrates retrieval loss, not evidence against the distant interaction.

Broader searches have not been shown to improve Top-k recovery in this study and introduce more untested competitors. The current public-input model supports cis pairs within 5 Mb, not whole-genome or interchromosomal discovery. Public tissue Hi-C is shared across queries; only DNase is matched to the optional cell. Contact evidence is not by itself proof of a regulatory effect. A future contact-guided candidate search would require its own matched evaluation.

## Supply candidate pairs

```sh
python run.py predict --input examples/pairs.csv
python run.py predict --input examples/pairs.csv --mode public --cell-type K562 --output output/public_pairs
```

Use zero-based, half-open hg38 intervals on chr1–chr22, chrX or chrY. Required columns are:

| Column | Meaning |
| --- | --- |
| `enhancer_chrom`, `enhancer_start`, `enhancer_end` | Candidate regulatory interval |
| `promoter_chrom`, `tss` | Target promoter chromosome and TSS |
| `gene` | Recommended gene symbol; a nearby annotated TSS is used if absent |

Only pairs on the same chromosome are supported. Optional `promoter_start` and `promoter_end` specify an exact promoter interval; otherwise it is TSS − 1,000 to TSS + 1,001 bp. This feature window is separate from the enhancer search radius. Public-input pairs must be within 5 Mb of their promoter.

Full-input prediction also needs eight **raw**, finite covariates:

| Feature | Information |
| --- | --- |
| `numTSSEnhGene` | Intervening TSS count |
| `distanceToTSS` | Enhancer–TSS distance |
| `normalizedDNase_prom` | Promoter accessibility |
| `sumNearbyEnhancers` | Nearby enhancer activity |
| `ubiquitousExpressedGene` | Broad expression annotation |
| `numCandidateEnhGene` | Candidate count |
| `contactFrequency` | Chromatin contact |
| `ABC.Score` | Activity-by-contact score |

Use the compatible ENCODE-rE2G/ABC feature definitions; the model applies its saved transformation. Both prediction modes automatically acquire missing ChIP-Atlas peaks. Full mode uses the eight values in the input table and does not need Node.js or public Hi-C/DNase downloads. Public mode obtains its own reference features and does not require these eight values. `--offline` requires the necessary references to be present in either mode.

## Read the results

`ranked_pairs.csv` contains coordinates, scores, ranks and regional contributions. `top_regions.csv` lists the top 50 distinct regions. Annotation-aware queries also save all scored windows and the separately listed coding/promoter windows.

A higher score gives a candidate higher priority. Enhancer and promoter components contribute to the model's log-odds score; they are not causal effects of individual proteins. Public-input parameters are fitted separately from full-input parameters. Public predictions expose both the default fit and the component-preserving adaptation:

| Field | Interpretation |
| --- | --- |
| `score` | Default public-input pair score; used for ranking |
| `enhancer_component`, `promoter_component` | Regional log-odds contributions in the default fit |
| `retained_promoter_component` | The fixed source-model promoter contribution, computed from the pooled public binding profile |
| `component_preserving_score` | An alternative pair score that retains the source promoter term and refits the public context and enhancer terms |

The retained component is identical for the same promoter profile with or without a cell reference. It is a model contribution associated with reporter responsiveness, not the experimentally measured responsiveness parameter or an expression change. The paper evaluates its biological association separately from pair ranking. The additional score does not change the default ordering. In `ranked_genes.csv`, these component fields describe the representative promoter, rather than the CAGE-weighted gene priority.

Choose a new `--output` directory when rerunning a command. Reference locations can be set with `--bed-dir` and `--cache-dir`, or with the environment variables `ATLAS_E2G_BED_DIR` and `ATLAS_E2G_CACHE_DIR`.


## Transcription-supported target priorities

For enhancer-coordinate queries, `ranked_genes.csv` uses `priority_score`: the CAGE-usage-weighted mean pair score multiplied by public detection support. The processed CAGE reference is included. `raw_max_pair_score` preserves the unweighted maximum; `representative_pair_score` describes the promoter contributing most to the priority. `ranked_promoters.csv` retains all raw scores, CAGE counts per million and normalized weights.

Support records transcription in reference cells, not proof of enhancer regulation. CAGE and DNase availability are independent: a cell can match DNase while using pooled CAGE. `gene_only_priority_score` uses the no-cell model and pooled CAGE. Gene-name queries retain pair-based region ranking; CAGE weighting is applied when comparing target genes for an enhancer.
