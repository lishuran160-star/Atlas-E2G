# Reproducing the paper

The included processed data are sufficient to refit the paper's models and calculate the reported metrics. Run commands from the repository root after installing `requirements.txt`. Python 3.11 and a CPU are sufficient.

## Recalculate or refit

```sh
# Recalculate metrics from the supplied prediction arrays:
python run.py reproduce --from-saved

# Refit from the processed training data:
python run.py reproduce --output output/refitted

# Export the method overview and two result figures:
python run.py figures
```

Recalculation is the quickest way to explore the results. Refitting additionally trains the models and repeats the holdouts. Both compare their results with `reference_results/`. They write numerical tables to a new output directory.

For one analysis at a time:

| Paper analysis | Command |
| --- | --- |
| Full-input holdouts, external CRISPR evaluation, reporter transfer and model controls | `python run.py reproduce --analyses full --output output/full` |
| Public accessibility ablations, including the fixed-K562 comparison | `python run.py reproduce --analyses activity --output output/activity` |
| Public promoter components, component-preserving transfer and distillation | `python run.py reproduce --analyses transfer --output output/transfer` |
| Training-collection removal and matched-size controls | `python run.py reproduce --analyses source --output output/source` |

Add `--bootstrap 2000` to calculate gene-cluster reporter intervals. The result figures use the supplied original interval estimates. The plotting script exports editable SVG text and vector PDF, as well as a 600 dpi PNG.

The benchmark analyses are implemented in `benchmark_full.py` and `benchmark_public.py`; `rebuild_public.py` obtains the public reference features.

The full-input results are in `full/results.csv`. Public and component analyses produce `public/metrics.csv`, with additional tables for matched-candidate top-k results, label-based ceilings and candidate retrieval. Fresh fitting also saves prediction arrays and fitted parameters.

To find the analyses behind the main text:

| Manuscript result | Command | Output |
| --- | --- | --- |
| Sections 3.1–3.2, Table 1 and Figure 2: full-input pair prediction and reporter association | `python run.py reproduce --from-saved` | `output/reproduction/full/results.csv` |
| Section 3.2: sensitivity to training-gene resampling | `python run.py stability --from-saved` | `output/stability/` |
| Sections 3.3–3.4 and Figure 3: public-input prediction and promoter-component transfer | `python run.py reproduce --from-saved` | `output/reproduction/public/metrics.csv` |
| Sections 3.3 and 3.5, Table 2: matched binding ablations and transcription-supported target queries | `python run.py query-analysis --from-saved` | `output/query_analysis/pair_summary.csv` and `query_summary.csv` |
| Section 3.5: MYC gene-name query | `python run.py query --gene MYC` | `output/query/ranked_pairs.csv` |
| Section 3.5: annotation-wide enhancer query | `python run.py targets --enhancer chr8:127735000-127736000` | `output/targets/ranked_genes.csv` and `ranked_promoters.csv` |

Pair AP evaluates predefined CRISPR-tested pairs. Reporter correlation evaluates a regional component against a different assay. Query Top-k measures recovery of known targets among retrieved candidates. These endpoints answer different questions; pair AP is not the success rate of a new gene-name query. Default public-input rankings and the retained promoter component are returned separately.

`reference_results/full_scores.npz` and `public_scores.npz` contain the supplied predictions. Their expected metrics are in `full_metrics.csv` and `public_metrics.json`. `figure_intervals.csv` and `source_removal_contrasts.csv` retain the reported uncertainty estimates.

## How evaluation works

Training labels come from endogenous CRISPR perturbations. Reporter responses are evaluation measurements and are never used as fitting labels. Source evaluation holds out chromosomes or studies; study holdouts also remove shared genes and enhancer elements.

The coordinate filter retains 10,275 source pairs, including 9,527 with full feature coverage, and 4,354 five-cell pairs, including 3,940 covered pairs. The full-input source fit uses the published covered cohort of 9,604 pairs; public-input and component-transfer fits use its coordinate-valid subset. Reporter evaluation has 285 measurements from 278 intervals and 271 genes.

Five-cell average precision uses the benchmark weights; astrocyte average precision is unweighted. Long-range subsets use the original `benchmark_distance > 100000`. Top-k comparisons share the same candidate set, group by cell/study/gene, and require at least four candidates. The label-based ceiling is the maximum number of labeled positives recoverable within the same selection budget.

Evaluation is retrospective: external results informed method selection. The four-gene CREST-seq analysis is a candidate-search diagnostic rather than an independent test. Its scores and annotations are in `data/candidate_screen.csv.gz`; seven ranking methods use identical eligible windows. `data/candidate_retrieval.csv.gz` records pair-level retrieval results. The reproduction command recalculates these summaries from the supplied candidate results.

## Training stability

```sh
# Recalculate the reported stability summaries offline:
python run.py stability --from-saved

# Repeat the 100 training-gene bootstrap refits:
python run.py stability --workers 4 --output output/stability_refitted
```

Each resample draws source genes with replacement and retains all their pairs. The analysis repeats chromosome-excluded fitting and preprocessing for the separate-role, shared-role and identity-removed models, using identical resamples across models. It holds architecture, regularization, reference features and reporter observations fixed. `training_stability_scores.npz` contains the 100 sets of reporter predictions; `training_stability_summary.csv` contains the reported summaries. These empirical refit distributions describe training-sample sensitivity within the source studies, separately from the evaluation bootstrap intervals. Prediction commands use the fitted parameters supplied in `models/`.

## Figures and query outputs

| Figure | Contents | Source |
| --- | --- | --- |
| 1 | Peak-to-profile construction, CRISPR training, reporter evaluation and a MYC query | Editable `docs/overview.svg`; MYC output from `python run.py query --gene MYC` |
| 2 | Five-cell precision–recall curves, top-100 positives, reporter measurements and role/identity controls | `full_scores.npz`, `figure_intervals.csv` and processed outcomes |
| 3 | Public-input AP and component-preserving reporter transfer | `public_scores.npz`, `public_metrics.json` and processed outcomes |

Figure 1a uses schematic peak intervals to show binary overlap, target-wise pooling and normalization. The expression bars in panel b are also schematic: CRISPR-tested functional links supply fitting labels, whereas reporter measurements evaluate the promoter term and never enter training. Panel c shows the first three rows of the default gene-only MYC query, with no cell type supplied. The displayed pairs share promoter `ENSG00000136997:P3` (hg38 chr8 TSS 127736593); their candidate intervals and scores are read directly from `ranked_pairs.csv`. They are ranked genomic candidates, not experimentally confirmed enhancers. An exact available cell match adds measured public DNase; otherwise the gene-only model is used. Default pair ranking and the retained source promoter contribution remain separate outputs.

The plotting command writes the values behind each panel, including the precision–recall curves. `figure3_values.csv` also preserves within-gene Top-3 counts and label-based ceilings discussed in the text.

Public prediction uses `models/public_reference.npz` for default scores and `models/public_components.npz` for the retained source promoter contribution and its adapted score. The transfer analysis saves a reproducible `public_components.npz` alongside its full `component_models.npz` output. These parameters come from the source CRISPR fit; reporter outcomes are never fitting labels. Reporter performance in the figures uses chromosome-excluded models, while query outputs use models fitted to the complete eligible source cohort.

## Rebuild the public reference features

The compact analysis inputs avoid repeating large downloads. To rebuild from public files instead:

```sh
python run.py download
python run.py rebuild-overlaps
npm ci --ignore-scripts
python run.py reproduce-public --source-folds
```

The last command acquires public contact and accessibility data, refits the public-input parameters and repeats source holdouts. Add `--offline` to reuse a complete cache. See the [download guide](reference_download.md) for details.

Reproduction starts from processed perturbation and reporter measurements; it does not repeat raw sequencing alignment or primary assay calling. Publications, accessions and processing definitions are in [SOURCES.json](../data/SOURCES.json).

## Public-input ablations and enhancer-to-gene recovery

~~~sh
python run.py query-analysis --from-saved
# Refit matched controls and source holdouts from processed public inputs:
python run.py query-analysis --output output/query_refitted
~~~

The first command recomputes source AP, external pair AP, target ranks and paired intervals. The second also refits the matched controls; deployed Atlas-E2G coefficients remain fixed. Both run offline from processed data. The expected summaries are `query_source_summary.csv`, `query_pair_summary.csv` and `query_summary.csv` in `reference_results/`.

The five-cell benchmark defines 174 enhancer/cell/study queries with a known positive. Candidate construction reads hg38 annotations, never the assayed-gene list. It includes all promoters within 1 Mb: 26,287 promoter pairs and 11,338 query–gene candidates. Coverage is 188/189 known-positive links; the remaining MYC link is approximately 1.97 Mb away. The principal target-choice subset contains 60 queries with multiple assayed genes. Unassayed genes are unknown. Distal recovery retains proximal competitors.

| Method | All annotations Top-1 / Top-3 | Shared published candidates Top-1 / Top-3 |
| --- | --- | --- |
| Atlas-E2G | 29/60 / 39/60 | 31/60 / 42/60 |
| Without target identity | 28/60 / 38/60 | 30/60 / 42/60 |
| Without binding profiles | 28/60 / 39/60 | 30/60 / 42/60 |
| Distance | 27/60 / 43/60 | 31/60 / 46/60 |
| ABC | — | 33/60 / 42/60 |
| rE2G | — | 28/60 / 41/60 |

Every row uses identical gene-level CAGE support. Atlas controls and distance also share promoter-use weights. ABC/rE2G retain their published gene-score aggregation and use the identical intersection of annotation and released-score coverage: 3,555 query–gene candidates, covering 171/189 positives. Results do not establish superiority in open-query Top-k recovery. Binding identity improves overall screened-pair AP (0.710 versus 0.607); this is a different endpoint.

`target-benchmark --from-saved` reproduces maximum-pair-score pooling without CAGE support, an explicit comparator (11/60 and 32/60 on all annotations), together with assayed-pair and nearest-TSS controls. Without `--from-saved` it reconstructs raw scores from public tracks; `--published-dir` rebuilds comparator scores from the five ENCODE files in `data/SOURCES.json`.

The transcription reference contains eight unperturbed conventional CAGE libraries from GSE118075, not NET-CAGE or perturbation outcomes. Processed counts, sample accessions, source URLs and preprocessing definitions are in `data/transcription.npz` and `data/transcription.json`. Exact cell matches use two libraries; unmatched or unspecified cells use all eight. Target-recovery intervals use 2,000 enhancer-coordinate resamples. All analyses reuse existing cohorts retrospectively.
