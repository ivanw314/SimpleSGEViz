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
sgeviz <input_dir> <output_dir> [--format html|png|svg] [--excel] [--protein-length N] [--px-per-aa N] [--gene-name NAME]
```

> If you haven't installed via `pip install -e .`, you can also run `python pipeline.py` directly with the same arguments.

### Arguments

| Argument | Description |
|---|---|
| `input_dir` | Directory containing the required input files (see below) |
| `output_dir` | Directory where output figures will be saved (created if it does not exist) |
| `--format` | Output format for figures: `html`, `png`, or `svg` (default: `png`) |
| `--excel` | Also write a multi-sheet Excel workbook for each gene (requires `openpyxl`) |
| `--protein-length N` | Known full protein length in amino acids. If the data covers fewer residues the x-axis is extended to this length. If omitted, you are prompted interactively for each gene. |
| `--px-per-aa N` | Pixels allocated per amino acid column in the AA heatmap (default: `4`). Reduce to produce a narrower figure, e.g. `--px-per-aa 2`. |
| `--gene-name NAME` | Override the gene name used in figure titles and output filenames. Useful when the name auto-detected from the filename differs from the preferred display name. Cannot be used when multiple gene datasets are detected in the same input directory. |

### Example

```bash
sgeviz ./data/BARD1/ ./output/BARD1/ --format png --excel
```

The pipeline will scan `./data/BARD1/` for gene datasets, generate figures for each gene found, and write everything to `./output/BARD1/`.

```bash
sgeviz ./data/BRCA1/ ./output/BRCA1/ --px-per-aa 2 --gene-name "BRCA1 (exon 11)"
```

Produces a narrower AA heatmap and overrides the auto-detected gene name in all figure titles and output filenames.

---

## Input Files

The pipeline detects genes automatically by scanning for `*allscores.tsv` files and extracting the gene name from each filename. Two naming conventions are supported:

- **Run-together:** `20260129_RAD51Dallscores.tsv` → `RAD51D`
- **Dot-separated:** `CTCF.allscores.tsv` → `CTCF`

Multiple genes can be processed in a single run by placing all their files in the same input directory.

### Required

| Pattern | Description |
|---|---|
| `*{gene}allscores.tsv` or `*{gene}.allscores.tsv` | Combined SNV and 3bp deletion fitness scores |
| `*{gene}modelparams.tsv` or `*{gene}.modelparams.tsv` | SGE model thresholds (non-functional and functional) |
| `*{gene}snvcounts.tsv` or `*{gene}.snvcounts.tsv` | Per-replicate SNV read counts |
| `*{gene}delcounts.tsv` or `*{gene}.delcounts.tsv` | Per-replicate 3bp deletion read counts |

### Optional (auto-detected; figures generated only when present)

| Pattern | Description |
|---|---|
| `*{gene}*ClinVar*SNV*` | ClinVar germline classification file (tab-delimited `.txt` from ClinVar download). File name matching is case-insensitive. |
| `*{gene}*gnomAD*` | gnomAD allele frequencies (CSV or Excel) |
| `*{gene}*Regeneron*` | Regeneron allele frequencies (CSV or Excel) |
| `*{gene}*domain*` | Protein domain annotations (CSV or Excel). Adds a domain cartoon strip above the AA heatmap. File name matching is case-insensitive. See [Domain annotation file format](#domain-annotation-file-format) below. |
| `*{gene}*editrates*` | Library edit rates (tab-delimited `.tsv`). File name matching is case-insensitive. See [Edit rates file format](#edit-rates-file-format) below. |
| `*{gene}*cartoon*` | Gene structure cartoon (Excel `.xlsx`). Generates a scalable exon cartoon with UTR/CDS distinction, ATG/stop markers, and compressed introns. If a `lib_coords` sheet is present a library amplicon track is added below. File name matching is case-insensitive. See [Gene cartoon file format](#gene-cartoon-file-format) below. |
| `*{gene}*vep*` | VEP output — Excel (`.xlsx`) or tab-delimited text (`.txt`). AlphaMissense, REVEL, CADD, and SpliceAI scores are extracted and merged into the variant data, enabling the VEP predictor sub-panels in the AA heatmap. File name matching is case-insensitive. See [VEP file format](#vep-file-format) below. |

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
| `max_SpliceAI` | Maximum SpliceAI delta score (variants > 0.2 are excluded from missense-specific predictor panels and the AA heatmap) |
| `am_score` | AlphaMissense pathogenicity score (shown in predictor scatter and AA heatmap VEP sub-panel). Can be supplied directly in this file or loaded automatically from a `*{gene}*vep*` file. |
| `revel_score` | REVEL score (shown in predictor scatter and AA heatmap VEP sub-panel). Can be supplied directly in this file or loaded automatically from a `*{gene}*vep*` file. |
| `cadd_score` | CADD phred score (shown in predictor scatter). Can be supplied directly in this file or loaded automatically from a `*{gene}*vep*` file. |
| `MutPred2` | MutPred2 score (shown in predictor scatter and AA heatmap VEP sub-panel). Must be supplied directly in this file. |

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
| `{gene}_predictor_scatter` | Grid of scatter plots: each predictor score (AlphaMissense, REVEL, CADD, MutPred2) vs. fitness score, colored by consequence, with SGE and published predictor threshold lines | If any predictor column present (from VEP file or `*allscores.tsv`) |
| `{gene}_clinvar_strip` | Strip plot of SGE scores by ClinVar germline classification | If ClinVar file detected |
| `{gene}_clinvar_roc` | ROC curve for B/LB vs P/LP classification using SGE score | If ClinVar file detected and both classes present |
| `{gene}_maf_vs_score` | Binned heatmap of log10(allele frequency) vs fitness score | If gnomAD or Regeneron file detected |
| `{gene}_edit_rate_barplot` | Bar chart of library edit rates per SGE target, grouped by replicate | If edit rates file detected |
| `{gene}_exon_cartoon` | Scalable exon structure cartoon with UTR/CDS regions, ATG/stop markers, and `//` break marks at compressed introns | If cartoon file detected (no `lib_coords` sheet) |
| `{gene}_library_cartoon` | Exon structure cartoon stacked above a library amplicon coverage track | If cartoon file detected with `lib_coords` sheet |

### Excel workbook (with `--excel`)

| Sheet | Contents |
|---|---|
| `scores` | Full scores table with `Germline classification` (if ClinVar) and `gnomad_af` / `regeneron_af` (if AF files) merged in by `pos_id` |
| `thresholds` | Non-functional and functional threshold values |
| `counts` | Per-replicate SNV and deletion read counts |
| `edit_rates` | Raw edit rates table (if edit rates file detected) |

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

---

## Edit rates file format

The edit rates file is a tab-delimited `.tsv` with two columns:

| Column | Description |
|---|---|
| `target_rep` | Combined SGE target and replicate identifier (e.g. `CTCF_X10A_R1R4_D05`) |
| `edit_rate` | Fractional library edit rate for that target/replicate |

The `target_rep` field is parsed as `{GENE}_X{target}_{rep}_{day}`. Replicate identifiers are mapped to display labels as follows:

| Raw replicate string | Display label |
|---|---|
| `R1R4`, `R1R2R3`, `R1`, `R4` | Rep. 1 |
| `R2R5`, `R4R5R6`, `R2`, `R5` | Rep. 2 |
| `R3R6`, `R7R8R9`, `R3`, `R6` | Rep. 3 |

### Example file (excerpt)

```
target_rep	edit_rate
CTCF_X10A_R1R4_D05	0.118527
CTCF_X10A_R2R5_D05	0.117841
CTCF_X10A_R3R6_D05	0.121773
CTCF_X3H_R1R4_D05	0.000138443
```

---

## Gene cartoon file format

The gene cartoon file is an Excel workbook (`.xlsx`) whose filename must contain `cartoon` (case-insensitive). Recommended naming: `{gene}.cartoon.xlsx` (e.g. `DARS2.cartoon.xlsx`). Avoid including `domain` in the filename — files matching both `*domain*` and `*cartoon*` are routed to the cartoon figure and will not be treated as domain annotation files. It generates one of two figures:

- **`{gene}_exon_cartoon`** — an exon structure cartoon with UTR/CDS height distinction, ATG/stop codon markers, `//` break marks at compressed introns, and exon number labels. Both introns and UTR regions are compressed using a min/max pixel-width scale, so very long UTRs or introns do not dominate the figure.
- **`{gene}_library_cartoon`** — the above exon track stacked above a library amplicon coverage track (generated when the `lib_coords` sheet is present).

### Sheet: `exon_coords` (required)

| Column | Description |
|---|---|
| `exon` | Exon label (e.g. `GENE_X1`) |
| `start` | Genomic start coordinate |
| `end` | Genomic end coordinate |

### Sheet: `metadata` (required)

A two-column table with `type` and `info` rows:

| `type` value | `info` description |
|---|---|
| `atg` | Genomic position of the start codon (used to delimit 5′ UTR from CDS) |
| `stop` | Genomic position of the stop codon (used to delimit CDS from 3′ UTR) |
| `strand` | Gene strand: `plus` (default) or `minus`. For minus-strand genes all coordinates in all sheets should still be in standard genomic convention (start < end); the figure is automatically drawn 5′→3′ left-to-right. |
| `exon_color` | Hex fill color for exon rectangles (default: `#2E86C1`) |
| `lib_color` | Hex fill color for library amplicons (default: `#888888`) |

### Sheet: `lib_coords` (optional)

| Column | Description |
|---|---|
| `start` | Genomic start of the library amplicon |
| `end` | Genomic end of the library amplicon |

When present, a library amplicon track is drawn below the exon cartoon. Regions covered by multiple overlapping amplicons are progressively darkened. A bracket above the band labels the estimated number of variants covered (`covered_bases × 3` SNVs + `covered_bases ÷ 3` 3-bp deletions), and the left margin shows the total number of amplicons.

---

## VEP file format

A VEP file provides pathogenicity predictor scores for the AA heatmap sub-panels. Place a VEP output file matching `*{gene}*vep*` (e.g. `BRCA1.vep.xlsx` or `BRCA1.VEP.txt`) in the input directory. Both Excel (`.xlsx`) and tab-delimited text (`.txt`) output formats from the Ensembl VEP web tool are supported.

### Generating the VEP file

**Step 1 — Create a VCF from your allscores file:**

```bash
python make_vcf.py path/to/GENE.allscores.tsv
```

This writes a `GENE.allscores.vcf` in the same directory and prints the suggested output filename and next steps.

**Step 2 — Run Ensembl VEP:**

Upload the VCF to [https://www.ensembl.org/Tools/VEP](https://www.ensembl.org/Tools/VEP) with the following options enabled:

- *Identifiers and frequency data*: Gene symbol, MANE
- *Pathogenicity predictions*: AlphaMissense, CADD
- *Plugins*: REVEL, SpliceAI

Download the results as Excel (`.xlsx`) or plain text (`.txt`) — both formats are supported.

**Step 3 — Add to input directory:**

Rename the downloaded file to `{gene}.vep.xlsx` or `{gene}.vep.txt` (e.g. `BRCA1.vep.xlsx`) and place it in your pipeline input directory. The pipeline will detect and load it automatically — both Excel and plain text formats are supported.

---

The pipeline expects columns produced by a standard VEP run with the following plugins:

| Plugin | Column(s) extracted | Renamed to |
|---|---|---|
| AlphaMissense | `am_pathogenicity` | `am_score` |
| REVEL | `REVEL` | `revel_score` |
| CADD | `CADD_PHRED` | `cadd_score` |
| SpliceAI | `SpliceAI_pred_DS_AG/AL/DG/DL` | `max_SpliceAI` (row-wise max) |

All score columns are optional — the pipeline uses whichever are present.

### How it works

The `Location` column (`chr:pos-pos` format, e.g. `2:214728633-214728633`) and `Allele` column are used to construct a `pos_id` matching the `pos:alt` format in `*allscores.tsv`. Scores are merged by left join, so variants without VEP scores receive NaN.

When multiple transcript rows exist for the same variant, the first row is kept (Ensembl VEP typically outputs the MANE_SELECT transcript first).

If a `*{gene}*vep*` file is not present but predictor score columns (`am_score`, `revel_score`, `cadd_score`, `MutPred2`) already exist in `*allscores.tsv`, those columns are used directly. If no predictor scores are available from either source, the predictor scatter figure and AA heatmap VEP sub-panels are omitted.
