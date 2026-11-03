# Obtaining public reference data

You can reproduce the paper's analyses with the small processed files already included. New public-input queries automatically acquire missing reference data. The commands below also let you prepare references before a query or rebuild peak-overlap matrices.

## What happens after specifying a target

After installing the Python and Node.js dependencies, `python run.py query --gene MYC` performs the complete public-input prediction. `python run.py targets --enhancer chr8:127735000-127736000` runs the reverse search. Add `--cell-type K562` to use a matching public accessibility track.

| Reference | How the query obtains it |
| --- | --- |
| GENCODE genes and promoters, candidate bins, gene classes | Processed tables included in `data/` |
| ChIP-Atlas binding profiles | Downloads missing files from the fixed peak list, then computes overlaps for the requested regions; subsequent queries reuse peaks and cached overlaps |
| Public tissue Hi-C | Retrieves and caches byte ranges needed for the requested promoters, without downloading the entire archive |
| Public DNase | Downloads only the matched cell's track, then retains a compact signal grid |
| Promoter use and transcription support | Processed conventional CAGE summaries included in `data/` |

The first query may download the full 6.3 GB peak reference, which is shared by all targets. The BED service supplies complete experiment files, so peak acquisition is not restricted to the query interval. `--candidates-only` lists candidates without downloading peak files or scoring. `--offline` prohibits new downloads and reports which required reference is absent. `query.json` records the cell match, reference source and whether peaks or cached overlaps were used.

## ChIP-Atlas peaks

```sh
python run.py download --bed-dir reference/chip_atlas --workers 4
```

This optional preparation command retrieves the same 7,070 human hg38 peak files used automatically by queries, covering 1,598 assay targets, about 6.3 GB in total. No login, API key, FASTQ or BAM files are needed. `data/reference.tsv.gz` lists each experiment's accession, target, cell label, public URL and expected file size.

Files come directly from ChIP-Atlas:

```text
https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/ACCESSION.05.bed
```

The `bed05` collection uses MACS2 Q < 10^-5. The [ChIP-Atlas download documentation](https://github.com/inutano/chip-atlas/wiki#10-downloads) describes this convention. For example, [SRX471323.05.bed](https://chip-atlas.dbcls.jp/data/hg38/eachData/bed05/SRX471323.05.bed) contains BCL11A peaks in CAL-1 cells.

If a download stops, run the same command again. Completed files with the expected size are skipped; incomplete downloads are retried. If an existing file has a different size, move that file aside and download it again. Use the listed experiments and peak threshold, since the fitted parameters depend on this reference.

For an interrupted `query`, `targets` or `predict` command, rerun with a new `--output` directory. Downloaded references remain in the shared cache and are reused; existing result directories are not overwritten.

To keep files on another drive:

```powershell
python run.py query --gene MYC --bed-dir "D:/Atlas-reference" --output output/MYC
```

You can also set `ATLAS_E2G_BED_DIR` once in your environment. Do not rename or combine the individual BED files.

### What the model reads

The first three BED columns give chromosome, start and end. Atlas-E2G records whether each candidate enhancer or promoter overlaps a peak in each experiment, then groups these overlaps by assay target. The full-input model averages within targets; the public-input model averages within target–cell groups and across cells.

These features describe binding evidence at a genomic region. Functional enhancer–gene labels come from the CRISPR experiments. The reference experiment list and its selection criteria are provided in [SOURCES.json](../data/SOURCES.json).

## Public accessibility and contact data

Install Node.js 20.19+ or 22.12+ and the contact-file reader:

```sh
npm ci --ignore-scripts
python run.py cell-types
python run.py query --gene MYC --cell-type K562 --output output/MYC_K562
```

Only the requested cell's ENCODE DNase track is downloaded. The code checks file length and genome assembly, integrates signal into 1 kb windows and retains the compact grid. The temporary raw bigWig is removed after processing. Accessions, source URLs and normalization scales are listed in `data/public_atlas.json`.

The public tissue Hi-C file is read through byte-range requests, so the complete archive is not downloaded. Contacts use SCALE normalization at 5 kb resolution. The tissue contact reference is shared across queries; the optional cell name selects DNase, not a new cell-specific Hi-C experiment. Both gene-only and cell-matched predictions use a ±5 Mb contact background. Candidate selection is separate: the defaults are ±500 kb for a gene query and ±1 Mb for an enhancer query.

To prepare a cell reference in advance:

```sh
python run.py prepare-public --cell-type WTC11
```

The 12-cell catalog is not downloaded all at once. A matching cell name selects existing measured accessibility; an unmatched name uses the gene-only score. Donor, culture and disease conditions may differ from those of the public reference.

## Reuse downloaded files

Defaults are `reference/chip_atlas/` for peaks and `reference/public/` for accessibility, contact and annotation caches. Override them with `--bed-dir` and `--cache-dir`, or the environment variables `ATLAS_E2G_BED_DIR` and `ATLAS_E2G_CACHE_DIR`. Use `--offline` to require cached public data.

Reference downloads, caches and generated outputs are excluded by `.gitignore`. They do not belong in the GitHub source folder. Public data retain their original terms; the MIT license covers the original Atlas-E2G software.

## 说中文的看这里

安装依赖后直接运行 `python run.py query --gene MYC`，程序会自动补齐缺失公共峰文件并生成候选评分；添加 `--cell-type K562` 可使用匹配的公共可及性数据。首次峰文件下载总计约 6.3 GB，后续复用；已有参考时用 `--bed-dir` 指定所在目录。也可以提前运行 `python run.py download --workers 4`。`--candidates-only` 只列候选，`--offline` 禁止新增下载。仅复现论文指标或重新拟合时，不需要这些原始公共文件。
