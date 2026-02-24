# SimpleSGEViz

A command-line pipeline for generating standard visualization figures from Saturation Genome Editing (SGE) fitness score data.

---

## Setup

### 1. Clone the repository

```bash
git clone https://github.com/ivanw314/SimpleSGEViz.git
cd SimpleSGEViz
```

### 2. Create and activate a virtual environment (recommended)

```bash
python -m venv .venv
source .venv/bin/activate      # macOS / Linux
.venv\Scripts\activate         # Windows
```

Or with conda:

```bash
conda create -n sgeviz -c conda-forge python=3.11 tzdata
conda activate sgeviz
```

> If you see a `tzdata` unsatisfiable error, your conda index is likely stale. Run `conda update -n base conda && conda clean --all`, then retry the create command.

### 3. Install the package

**Python 3.9 or higher is required** (the conda example above uses 3.11; any 3.9+ version works).

```bash
pip install -e .
```

This installs all dependencies, including PNG/SVG export support (`vl-convert-python`) and Excel output support (`openpyxl`).

After a `git pull`, changes are picked up automatically — no reinstall needed. Re-run `pip install -e .` only if `pyproject.toml` changes (e.g. new dependencies were added).

### 4. Verify the setup

```bash
sgeviz --help
```

You should see the usage message. If you get a "command not found" error, make sure your virtual environment is active. If you get an import error, check that the install step completed without errors.

---

## Usage

```
sgeviz <input_dir> <output_dir> [--format html|png|svg] [--excel]
```

> If you haven't installed via `pip install -e .`, you can also run `python pipeline.py` directly with the same arguments.

### Arguments

| Argument | Description |
|---|---|
| `input_dir` | Directory containing the required input files (see below) |
| `output_dir` | Directory where output figures will be saved (created if it does not exist) |
| `--format` | Output format for figures: `html`, `png`, or `svg` (default: `png`) |
| `--excel` | Also write a multi-sheet Excel workbook for each gene (requires `openpyxl`) |

### Example

```bash
sgeviz ./data/BARD1/ ./output/BARD1/ --format png --excel
```

The pipeline will scan `./data/BARD1/` for gene datasets, generate figures for each gene found, and write everything to `./output/BARD1/`.

---

## Input Files

The pipeline detects genes automatically by scanning for `*allscores.tsv` files and extracting the gene name from each filename (e.g. `20260129_RAD51Dallscores.tsv` → `RAD51D`). Multiple genes can be processed in a single run by placing all their files in the same input directory.

### Required

| Pattern | Description |
|---|---|
| `*{gene}allscores.tsv` | Combined SNV and 3bp deletion fitness scores |
| `*{gene}modelparams.tsv` | SGE model thresholds (non-functional and functional) |
| `*{gene}snvcounts.tsv` | Per-replicate SNV read counts |
| `*{gene}delcounts.tsv` | Per-replicate 3bp deletion read counts |

### Optional (auto-detected; figures generated only when present)

| Pattern | Description |
|---|---|
| `*{gene}*ClinVar*SNV*` | ClinVar germline classification file (tab-delimited `.txt` from ClinVar download). File name matching is case-insensitive. |
| `*{gene}*gnomAD*` | gnomAD allele frequencies (CSV or Excel) |
| `*{gene}*Regeneron*` | Regeneron allele frequencies (CSV or Excel) |
| `*{gene}*domain*` | Protein domain annotations (CSV or Excel). Adds a domain cartoon strip above the AA heatmap. File name matching is case-insensitive. See [Domain annotation file format](#domain-annotation-file-format) below. |

### Required columns in `*allscores.tsv`

| Column | Description |
|---|---|
| `pos` | Genomic position (GRCh38, 1-based) |
| `ref` | Reference allele |
| `alt` | Alternate allele |
| `score` | SGE fitness score |
| `functional_consequence` | GMM classification (`functionally_normal`, `functionally_abnormal`, `intermediate`) |
| `variant_qc_flag` | Quality flag; rows with `WARN` are excluded |
| `exon` | Exon identifier (e.g. `GENE_X2A`) |
| `chrom` | Chromosome |
| `consequence` | VEP consequence term (e.g. `missense_variant`) |

**Optional columns** (enable additional figure outputs):

| Column | Description |
|---|---|
| `amino_acid_change` | Amino acid substitution in `A123G` format (enables AA heatmap) |
| `max_SpliceAI` | Maximum SpliceAI delta score (variants > 0.2 are excluded from AA heatmap) |
| `am_score` | AlphaMissense pathogenicity score (shown as VEP sub-panel in AA heatmap) |
| `revel_score` | REVEL score (shown as VEP sub-panel in AA heatmap) |
| `MutPred2` | MutPred2 score (shown as VEP sub-panel in AA heatmap) |

---

## Outputs

All files are written to `output_dir` with the gene name as a prefix.

### Figures

| File | Description | Condition |
|---|---|---|
| `{gene}_histogram_stripplot` | Score distribution histogram (top) and per-consequence strip plot (bottom) | Always |
| `{gene}_correlation_heatmap` | Pairwise Pearson r heatmap across replicates | Always |
| `{gene}_scores_across_gene` | Per-exon scatter plot of fitness scores vs genomic position | Always |
| `{gene}_aa_heatmap` | Amino acid substitution heatmap (AA position × substitution), optionally with a domain cartoon strip and 3bp deletion scatter panel | If `amino_acid_change` column present |
| `{gene}_clinvar_strip` | Strip plot of SGE scores by ClinVar germline classification | If ClinVar file detected |
| `{gene}_clinvar_roc` | ROC curve for B/LB vs P/LP classification using SGE score | If ClinVar file detected and both classes present |
| `{gene}_maf_vs_score` | Binned heatmap of log10(allele frequency) vs fitness score | If gnomAD or Regeneron file detected |

### Excel workbook (with `--excel`)

| Sheet | Contents |
|---|---|
| `scores` | Full scores table with `Germline classification` (if ClinVar) and `gnomad_af` / `regeneron_af` (if AF files) merged in by `pos_id` |
| `thresholds` | Non-functional and functional threshold values |
| `counts` | Per-replicate SNV and deletion read counts |

---

## ClinVar file format

> **Note:** The pipeline currently supports ClinVar annotations for **SNVs only**. 3bp deletion variants are not matched against ClinVar entries.

The ClinVar input file should be the tab-delimited tabular download from [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/). The pipeline requires these columns:

- `GRCh38Location` — 1-based genomic coordinate
- `Canonical SPDI` — Used to extract the + strand alt allele; correctly handles both + and − strand genes without manual strand detection
- `Germline classification` — ClinVar germline pathogenicity classification

Compound labels are normalized automatically:

| Raw label | Normalized to |
|---|---|
| `Benign/Likely benign` | `Likely benign` |
| `Pathogenic/Likely pathogenic` | `Likely pathogenic` |
| `Conflicting classifications of pathogenicity` | `Uncertain significance` |

Variants with any other classification are excluded.

---

## Domain annotation file format

A domain annotation file adds a colored protein domain cartoon strip above the amino acid heatmap. It can be a CSV or Excel file and must contain these columns:

| Column | Description |
|---|---|
| `region_name` | Label shown in the cartoon (e.g. `RecA`, `Walker A`) |
| `aa_residues` | Residue range in `start-end` format (e.g. `114-209`) |

**Optional columns:**

| Column | Description |
|---|---|
| `color` | Hex color for the segment (e.g. `#B9DBF4`). Auto-assigned from a built-in palette if omitted. |
| `tier` | `0` = main domain (default), `1` = sub-feature nested inside a tier-0 domain. |

### Tier-0 vs tier-1 domains

By default all rows are tier-0. Use `tier = 1` to mark sub-features (e.g. Walker A/B motifs inside a RecA domain). Tier-1 sub-features are rendered inline within the parent domain rather than in a separate row: the parent segment is split around the sub-feature, so the cartoon stays as a single strip.

For example, if RecA spans residues 114–209 and Walker A spans 135–140, the cartoon renders as:

```
[RecA 114–135][Walker A 135–140][RecA 140–209]
```

Gaps between tier-0 domains are automatically filled with gray.

### Label display rules

- Labels are drawn above the colored strip.
- Tier-1 labels are always shown (they take priority over tier-0 labels when space is tight).
- A tier-0 label is hidden if it would overlap a tier-1 label or an adjacent tier-0 label. The segment is still visible and its name appears in the tooltip on hover.

### Example file

| region_name | aa_residues | tier | color |
|---|---|---|---|
| N-terminal domain | 1-113 | 0 | #B9DBF4 |
| RecA | 114-209 | 0 | #C8DBC8 |
| Walker A | 135-140 | 1 | #FF9A00 |
| Walker B | 180-185 | 1 | #ffc976 |
| C-terminal domain | 210-328 | 0 | #D5A6BD |
