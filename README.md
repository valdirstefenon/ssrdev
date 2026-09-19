# SSRdev 4.5

**A Platform for SSR Marker Development and Validation**

[![Python](https://img.shields.io/badge/python-3.8%2B-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](#license)
[![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-lightgrey.svg)](#system-requirements)

**SSRdev** (Simple Sequence Repeat Developer) is a comprehensive desktop application for the discovery, characterization, and validation of SSR molecular markers from genomic sequence data. Built with Python and Tkinter, it provides a user-friendly graphical interface accessible to researchers without programming experience.

---

## Table of Contents

- [Features](#-features)
- [Installation](#-installation)
- [System Requirements](#-system-requirements)
- [Getting Started](#-getting-started)
- [Module 1: Prospecting](#-module-1-prospecting-ssr-discovery--primer-design)
- [Module 2: Primer Settings](#-module-2-primer-settings)
- [Module 3: Characterization](#-module-3-characterization-blast-analysis)
- [Module 4: Validation](#-module-4-validation-in-silico-pcr)
- [Project Management](#-project-management)
- [Output Files](#-output-files)
- [Troubleshooting](#-troubleshooting)
- [Appendix A: Default Parameters](#-appendix-a-default-parameters)
- [Appendix B: File Format Specifications](#-appendix-b-file-format-specifications)
- [Citation](#-citation)
- [License](#-license)
- [Support](#-support)
- [Acknowledgments](#-acknowledgments)

---

## Features

SSRdev 4.5 integrates five complementary analysis capabilities:

- **SSR Mining** — detection of mono- to hexa-nucleotide repeats.
- **Primer Design** — automated `primer3`-based primer design.
- **BLAST Characterization** — NCBI and local BLAST for marker validation.
- **In Silico PCR** — virtual amplification testing across multiple genomes.
- **Phylogenetic Analysis** — dendrogram generation from presence/absence data.

---

## Installation

### Prerequisites

- **Python 3.8 or higher** (Python 3.10+ recommended)
- `pip` package manager
- Internet connection (for NCBI BLAST functionality)

### Step 1 — Install Python

Download and install Python from [python.org](https://www.python.org/). During installation, ensure you check **"Add Python to PATH"**.

### Step 2 — Install required dependencies

Open a terminal (Command Prompt on Windows, Terminal on macOS/Linux) and run:

```bash
pip install biopython pandas openpyxl numpy matplotlib scikit-learn scipy gffutils gtfparse primer3-py
```

### Step 3 — Verify installation

```bash
python -c "import Bio, pandas, primer3, openpyxl; print('All dependencies installed successfully!')"
```

### Step 4 — Run SSRdev

Save the script as `SSRdev.py` and run:

```bash
python SSRdev.py
```

### Dependencies Summary

| Package | Version | Purpose |
|---|---|---|
| **Biopython** | ≥ 1.79 | Sequence parsing, BLAST, GenBank handling |
| **pandas** | ≥ 1.3.0 | Data manipulation and Excel output |
| **openpyxl** | ≥ 3.0.0 | Excel file reading/writing with formatting |
| **numpy** | ≥ 1.21.0 | Numerical operations |
| **matplotlib** | ≥ 3.4.0 | Gel visualization, SSR maps, dendrograms |
| **scikit-learn** | ≥ 1.0.0 | Distance calculations for dendrograms |
| **scipy** | ≥ 1.7.0 | Hierarchical clustering |
| **gffutils** | ≥ 0.9.0 | GFF3 file parsing |
| **gtfparse** | ≥ 1.2.0 | GTF file parsing |
| **primer3-py** | ≥ 2.0.0 | Primer design engine |

> **Note:** If `primer3-py` fails to install, try:
> ```bash
> pip install primer3-py --no-cache-dir
> ```

---

## System Requirements

| Component | Minimum | Recommended |
|---|---|---|
| **Operating System** | Windows 10 / macOS 10.15 / Ubuntu 20.04 | Windows 11 / macOS 12+ / Ubuntu 22.04 |
| **RAM** | 4 GB | 8 GB or more |
| **Storage** | 500 MB free space | 2 GB (for BLAST databases) |
| **Processor** | 2 cores @ 2.0 GHz | 4+ cores @ 2.5+ GHz |
| **Display** | 1280 × 720 | 1920 × 1080 or higher |
| **Python** | 3.8 | 3.10+ |

---

## Getting Started

Launch the application:

```bash
python SSRdev.py
```

The main window will appear with **four tabs**:

1. **Prospecting** — SSR discovery and primer design.
2. **Primer Settings** — configure primer design parameters.
3. **Characterization** — BLAST-based marker validation.
4. **Validation** — in silico PCR testing.

### Interface Overview

- **Left Sidebar** — navigation, SSR search settings, project management.
- **Main Content Area** — active module interface.
- **Status Bar** — progress indicators and status messages.

---

## Module 1: Prospecting (SSR Discovery & Primer Design)

### Purpose

Identify SSR motifs within input sequences and design flanking PCR primers.

### Input Files

#### FASTA File (required)

- **Format:** standard FASTA (`.fasta`, `.fa`, `.fna`)
- **Content:** nucleotide sequences (DNA)
- **Example:**

```text
>Scaffold_1
ATCGATCGATCGATCGATCGATCGATCGATCGATCGATCG
>Scaffold_2
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
```

#### Or Paste Sequence

Directly paste FASTA-formatted sequences into the text area. Multiple sequences are supported.

### SSR Detection Parameters

Configure in the left sidebar:

| SSR Type | Motif Length | Default Min Repeats | Typical Usage |
|---|---|---|---|
| **Mono** (mononucleotide) | 1 bp | 10 | Less common, prone to stutter |
| **Di** (dinucleotide) | 2 bp | 5 | Most common in plants |
| **Tri** (trinucleotide) | 3 bp | 4 | Common in coding regions |
| **Tetra** (tetranucleotide) | 4 bp | 3 | High resolution |
| **Penta** (pentanucleotide) | 5 bp | 3 | Very stable |
| **Hexa** (hexanucleotide) | 6 bp | 3 | Excellent resolution |

> **Pro Tip:** Enable/disable specific SSR types using the checkboxes.

### Running Analysis

1. Enter a **Project Name** (optional).
2. Select **FASTA file** *or* paste sequences.
3. Adjust **SSR detection parameters** as needed.
4. Click **"Run Analysis"**.

### Progress Monitoring

- The progress bar shows completion status.
- The status label updates with the current sequence being processed.
- A log file is created: `ssrdev_YYYYMMDD_HHMMSS.log`.

### Results Display

| Column | Description |
|---|---|
| **Type** | SSR type (mono, di, tri, tetra, penta, hexa) |
| **Start** | Start position in sequence (1-indexed) |
| **End** | End position in sequence |
| **Motif** | Repeat unit sequence (canonical form) |
| **Repeats** | Number of complete repeats |
| **Length** | Total SSR length in base pairs |
| **Fwd Primer** | Designed forward primer sequence |
| **Fwd GC%** | GC content of forward primer |
| **Fwd Tm** | Melting temperature (°C) |
| **Fwd Size** | Primer length (bp) |
| **Rev Primer** | Designed reverse primer sequence |
| **Rev GC%** | GC content of reverse primer |
| **Rev Tm** | Melting temperature (°C) |
| **Rev Size** | Reverse primer length (bp) |
| **Product Size** | Expected amplicon size (bp) |
| **Tm Diff** | Difference between primer Tm values (°C) |

### Filtering Results

- **Type filter** — show only specific SSR types.
- **Search box** — text search across all fields.

### SSR Map Visualization

Each sequence includes a map tab showing:

- Linear representation of the sequence.
- Colored vertical bars at SSR positions.
- Color legend for SSR types.
- Position scale in base pairs.

### Export Options

**Flanking FASTA Export** — right-click or use the **"Export Flanking FASTA"** button to save regions surrounding SSRs (default: 200 bp each side).

---

## Module 2: Primer Settings

### Purpose

Configure and save `primer3` parameters for SSR flanking primer design.

### Primer Size

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Min Size** | 15 bp | 12–30 | Minimum primer length |
| **Max Size** | 30 bp | 15–40 | Maximum primer length |
| **Opt Size** | 20 bp | 15–25 | Optimal primer length |

### GC Content

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Min GC** | 30% | 20–60% | Minimum GC percentage |
| **Max GC** | 70% | 40–80% | Maximum GC percentage |
| **Opt GC** | 50% | 40–60% | Optimal GC percentage |

### Melting Temperature (Tm)

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Min Tm** | 50 °C | 45–60 °C | Minimum melting temperature |
| **Max Tm** | 70 °C | 60–80 °C | Maximum melting temperature |
| **Opt Tm** | 60 °C | 55–65 °C | Optimal melting temperature |
| **Max Diff Tm** | 2 °C | 1–5 °C | Maximum Tm difference between primers |

### Product Size

| Parameter | Default | Range | Description |
|---|---|---|---|
| **Min Product Size** | 100 bp | 50–500 bp | Minimum amplicon length |
| **Max Product Size** | 600 bp | 100–1500 bp | Maximum amplicon length |

### Primer Quality Filters

| Parameter | Default | Description |
|---|---|---|
| **Max Poly-X** | 4 | Maximum length of mononucleotide runs |
| **Max 3' GC Clamp** | 2 | Maximum G/C bases in last 5 positions |
| **Max Self Any** | 8.0 | Maximum self-complementarity score |
| **Max Self End** | 3.0 | Maximum 3' self-complementarity |
| **Max Pair Compl Any** | 8.0 | Maximum pair complementarity |
| **Max Pair Compl End** | 3.0 | Maximum 3' pair complementarity |

### Salt and Concentration

| Parameter | Default | Units | Description |
|---|---|---|---|
| **Monovalent Salt** | 50.0 | mM | KCl or NaCl concentration |
| **Divalent Salt** | 1.5 | mM | MgCl₂ concentration |
| **dNTP Concentration** | 0.8 | mM | dNTP concentration |
| **DNA Concentration** | 50.0 | nM | Primer concentration |

### Parameter Validation

- Invalid entries (non-positive numbers) appear in **red text**.
- Validation errors prevent analysis execution.

### Saving / Loading Settings

- **Save Settings** — export parameters as JSON (`.json`).
- **Reset to Defaults** — restore factory settings.

---

## Module 3: Characterization (BLAST Analysis)

### Purpose

Validate SSR markers by BLAST searching against reference databases to verify uniqueness and identify gene associations.

### Input Requirements

#### Reference FASTA File (required)

- The same FASTA file used in the prospecting module.
- Contains the original reference sequences.

#### SSR Results File (required)

- Excel file (`.xlsx`) from the Prospecting module.
- Must contain `Marker`, `Start`, and `End` columns.

#### Annotation File (optional — one of the following)

| Format | Extension | Description |
|---|---|---|
| **GenBank** | `.gb`, `.gbff`, `.gbk` | Full annotated genome |
| **GFF3 + FASTA** | `.gff`, `.gff3` + `.fasta` | Gene annotation file + sequence |
| **GTF + FASTA** | `.gtf` + `.fasta` | Gene Transfer Format + sequence |

### BLAST Options

#### NCBI BLAST (Online)

- **Database:**
  - `DNA_refSeq` — NCBI nucleotide reference sequences
  - `mRNA_refSeq` — NCBI mRNA reference sequences
- **Organism filter:** select taxonomic family from built-in list.
- **Hit limit:** 5 best hits.
- **Expect threshold:** 10.

#### Local BLAST (Offline)

- Uses loaded GenBank/GFF/GTF annotations.
- Configurable mismatch tolerance (0–100%).
- Gene flank expansion for linkage analysis.

### BLAST Parameters

| Parameter | Default | Description |
|---|---|---|
| **Flanking region size** | 2000 bp | Sequence extracted around SSR for BLAST |
| **Gene flanking region** | 2000 bp | Region expanded around genes for overlap detection |
| **Min overlap** | 50 bp | Minimum overlap for gene linkage |
| **Max mismatches** | 10% | Maximum allowed mismatches in match |

### Running BLAST

1. Select **Family** using the family browser (NCBI BLAST only).
2. Load **Reference FASTA** file.
3. Load **SSR Results** file.
4. (Optional) Load **Annotation** file for gene linkage.
5. Adjust **BLAST parameters**.
6. Click **"Run NCBI BLAST"** or **"Run Local GenBank BLAST"**.

### Results Display

| Column | Description |
|---|---|
| **Marker** | SSR marker identifier |
| **Status** | FOUND / NOT_FOUND / ERROR / LINKED |
| **Accession** | GenBank accession or sequence ID |
| **Gene Link Info** | Associated gene(s) with overlap details |
| **Description** | Match description or position |
| **E-value** | BLAST expectation value |
| **Identity %** | Percent sequence identity |
| **Original Name** | Original sequence identifier |

### Gene Linkage Information

When annotation is provided, results include:

- **Linked** — marker overlaps or is near a gene (green highlight).
- **Intergenic region** — marker in non-coding region (yellow highlight).
- **Gene details** — gene name, type, strand, overlap length.

**Example:**

```text
GeneX (gene, strand +1, overlap 150bp) | GeneY (CDS, strand -1, overlap 25bp)
```

### Output Files

- `[base]_[NCBI/Local]_BLAST_[timestamp].xlsx` — full BLAST results.
- Gene linkage information included for local BLAST.

---

## Module 4: Validation (In Silico PCR)

### Purpose

Simulate PCR amplification across multiple genome sequences to:

- Confirm primer specificity.
- Predict amplicon sizes.
- Generate virtual gel images.
- Compute genetic relationships.

### Input Files

#### Multi-FASTA File (required)

- Contains all genome sequences to test.
- Format: standard FASTA with unique identifiers.
- **Example:**

```text
>SpeciesA_genome
ATCGATCGATCGATCGATCGATCGATCGATCG...
>SpeciesB_genome
GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT...
```

#### Primer Excel File (required)

- Generated from the Prospecting module (`*_Primers.xlsx`).
- Required columns:
  - **Marker** — SSR marker name.
  - **Forward primer** — forward primer sequence.
  - **Reverse primer** — reverse primer sequence.

### Options

#### Virtual Gel Generation

| Option | Description |
|---|---|
| **Generate virtual gel images** | Master toggle for gel generation |
| **Gel for each scaffold** | One gel per species showing all markers |
| **Gel for each marker** | One gel per marker showing all species |

**Gel details:**

- 100 bp DNA ladder included.
- Band colors: **red** (scaffold gels), **blue** (marker gels).
- Amplicons > 1000 bp flagged with asterisk.
- Output: JPEG images (300 DPI).

#### Relationship Analysis

| Option | Description |
|---|---|
| **Compute relationship** | Generate presence/absence matrix |
| **Bootstrap dendrogram** | Perform bootstrap resampling |
| **Replicates** | Number of bootstrap iterations (default: 100) |

### Output Files

#### Excel Workbook (`*_validation_results.xlsx`)

| Sheet Name | Content |
|---|---|
| **Amplicon_Summary** | Consolidated amplicon sizes by marker |
| **[Species Name]** | Per-species amplification results |
| **Presence_Absence** | Binary matrix (1 = present, 0 = absent) |
| **Similarity_Matrix** | Jaccard similarity coefficients |
| **No_Amplicons** | Species with no successful amplifications |

#### Gel Images

- **Scaffold gels:** `[species]_scaffold_gel_[n].jpg`
- **Marker gels:** `marker_[marker_name]_gel_[n].jpg`

#### Dendrogram Outputs

| File | Format | Description |
|---|---|---|
| `dendrogram.jpg` | JPEG | UPGMA dendrogram image |
| `dendrogram.pdf` | Vector | Publication-quality figure |
| `dendrogram.newick` | Newick | Tree file for external software |

### Understanding Results

#### Amplicon Sizes

- Multiple sizes indicate possible alternative binding sites.
- "Not found" means no amplification in that species.
- Sizes listed as comma-separated values.

#### Presence/Absence Matrix

- Used for phylogenetic distance calculation.
- `1` = successful amplification (any size).
- `0` = no amplification.

#### Genetic Distance

- Based on Jaccard similarity: `J(A, B) = |A ∩ B| / |A ∪ B|`.
- Converted to distance: `1 − J(A, B)`.
- UPGMA clustering algorithm.

---

## Project Management

### Saving a Project

- **File extension:** `.ssrdev` (JSON format).
- **Saved data includes:**
  - Project name
  - FASTA file path
  - SSR detection parameters
  - Primer design parameters
  - Analysis results (without full sequences)
  - Marker name mapping
  - BLAST reference file path

### Loading a Project

- Restores all parameters and previously computed results.
- Re-loads reference sequences automatically.
- Analysis results are displayed immediately.

### Log Files

- Created automatically in the working directory.
- **Format:** `ssrdev_YYYYMMDD_HHMMSS.log`.
- **Contains:** timestamps, analysis progress, errors, warnings.

---

## Output Files Summary

### Prospecting Module

| File | Format | Content |
|---|---|---|
| `[project]_SSR_results.txt` | Text | Detailed SSR and primer information |
| `[project]_SSR_results.xlsx` | Excel | Tabular results with summary sheets |
| `[project]_Primers.xlsx` | Excel | Primer sequences only (deduplicated) |

**Excel workbook sheets:**

| Sheet Name | Description |
|---|---|
| **prospecting and primer design** | Main results table |
| **Total summary** | Counts by SSR type with bar chart |
| **Summary by scaffold** | Per-sequence breakdown |

### Characterization Module

| File | Format | Content |
|---|---|---|
| `[base]_[type]_BLAST_[timestamp].xlsx` | Excel | Complete BLAST results |

### Validation Module

| File / Folder | Description |
|---|---|
| `[base]_validation_results.xlsx` | Main validation workbook |
| `gel_images/` | Virtual gel JPEG images |
| `dendrograms/` | Dendrogram images and Newick file |

---

## Troubleshooting

### Common Issues and Solutions

#### "primer3 module not available"

Install `primer3-py`:

```bash
pip install primer3-py
```

#### "Biopython version 1.79 or higher required"

Upgrade Biopython:

```bash
pip install --upgrade biopython
```

#### "Error reading FASTA file"

**Possible causes:**

- File contains non-standard characters.
- Sequence contains spaces or line breaks.
- File encoding issues.

**Solution:**

- Ensure FASTA format with `>` headers.
- Remove spaces from sequence lines.
- Save file as UTF-8.

#### NCBI BLAST fails

**Possible causes:**

- No internet connection.
- Organism family not properly selected.
- NCBI server timeout.

**Solutions:**

- Check internet connection.
- Verify family selection.
- Try local BLAST with a GenBank file.
- Wait and retry (NCBI server limits).

#### "No gene annotations" warning

Download a complete GenBank file:

- From NCBI, select **"GenBank (full)"** format.
- Look for files with `.gbff` extension.
- Ensure the file includes feature annotations.

#### Permission error saving files

**Solution:**

- Close the target Excel file if open.
- Check write permissions for the output directory.
- Save to a different location.

#### Memory errors with large genomes

**Solutions:**

- Split large FASTA files into smaller chunks.
- Increase system RAM.
- Use 64-bit Python.
- Process fewer sequences at once.

### Performance Tips

1. **For large genomes (> 500 MB):**
   - Use local BLAST instead of NCBI.
   - Reduce flanking region size.
   - Process sequences in batches.

2. **For many sequences (> 1000):**
   - Disable primer design if not needed.
   - Disable gel image generation.
   - Reduce bootstrap replicates.

3. **For NCBI BLAST:**
   - Add delays between queries (automatic in SSRdev).
   - Use local GenBank for large datasets.

---

## Appendix A: Default Parameters

### SSR Detection Defaults

```json
{
  "mono":  {"motif_length": 1, "min_repeats": 10, "enabled": true},
  "di":    {"motif_length": 2, "min_repeats": 5,  "enabled": true},
  "tri":   {"motif_length": 3, "min_repeats": 4,  "enabled": true},
  "tetra": {"motif_length": 4, "min_repeats": 3,  "enabled": true},
  "penta": {"motif_length": 5, "min_repeats": 3,  "enabled": true},
  "hexa":  {"motif_length": 6, "min_repeats": 3,  "enabled": true}
}
```

### Primer3 Defaults

```json
{
  "PRIMER_OPT_SIZE": 20,
  "PRIMER_MIN_SIZE": 15,
  "PRIMER_MAX_SIZE": 30,
  "PRIMER_OPT_TM": 60.0,
  "PRIMER_MIN_TM": 50.0,
  "PRIMER_MAX_TM": 70.0,
  "PRIMER_MIN_GC": 30.0,
  "PRIMER_MAX_GC": 70.0,
  "PRIMER_OPT_GC_PERCENT": 50.0,
  "PRIMER_MAX_POLY_X": 4,
  "PRIMER_MAX_END_GC": 2,
  "PRIMER_PRODUCT_SIZE_RANGE": [[100, 600]],
  "PRIMER_MAX_SELF_ANY": 8.0,
  "PRIMER_MAX_SELF_END": 3.0,
  "PRIMER_PAIR_MAX_COMPL_ANY": 8.0,
  "PRIMER_PAIR_MAX_COMPL_END": 3.0,
  "PRIMER_MAX_DIFF_TM": 2.0,
  "PRIMER_SALT_MONOVALENT": 50.0,
  "PRIMER_SALT_DIVALENT": 1.5,
  "PRIMER_DNTP_CONC": 0.8,
  "PRIMER_DNA_CONC": 50.0
}
```

---

## Appendix B: File Format Specifications

### FASTA Format

```text
>[header] [optional description]
SEQUENCE_LINE_1
SEQUENCE_LINE_2
...
```

**Rules:**

- Header line starts with `>`.
- No spaces in sequence lines.
- Nucleotides: A, T, C, G only (case insensitive).
- IUPAC ambiguity codes allowed: R, Y, S, W, K, M, B, D, H, V, N.

### GenBank Format (for annotation)

- Standard NCBI GenBank format.
- Must include `FEATURES` section with gene/CDS entries.
- `.gbff` (full) files preferred over `.gb` (partial).

### Excel Input for Characterization

**Required columns:**

| Column | Type | Description |
|---|---|---|
| **Marker** | String | Unique marker identifier |
| **Start** | Integer | SSR start position (1-indexed) |
| **End** | Integer | SSR end position |

### Excel Input for Validation

**Required columns:**

| Column | Type | Description |
|---|---|---|
| **Marker** | String | Unique marker identifier |
| **Forward primer** | String | Forward primer sequence (5' → 3') |
| **Reverse primer** | String | Reverse primer sequence (5' → 3') |

### Project File (`.ssrdev`)

JSON format containing:

- Analysis parameters
- Results metadata
- File references
- Marker mappings

---

## Citation

If you use SSRdev in your research, please cite:

> STEFENON, V. M. (2026). **SSRdev: A Comprehensive Platform for SSR Marker Development and Validation.** https://github.com/valdirstefenon/ssrdev

---

## License

SSRdev is provided under the **MIT License** for academic and commercial use.

---

## Support

For bug reports, feature requests, or general assistance:

- **GitHub Issues:** [Project Repository URL]
- **Email:** [support@example.com]

When reporting issues, please include:

1. SSRdev version
2. Python version (`python --version`)
3. Operating system
4. Input file samples (if possible)
5. Error message or log file excerpt

---

## Acknowledgments

SSRdev incorporates:

- **Primer3** for primer design (Untergasser et al., 2012).
- **Biopython** for sequence handling (Cock et al., 2009).
- **NCBI BLAST** for sequence similarity searches (Altschul et al., 1990).

---

*README version 1.0 — Compatible with SSRdev 4.5*
