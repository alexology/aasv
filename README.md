# aasv

**aasv** is an R package for Windows that provides a complete pipeline for
assessing whether Amplicon Sequence Variants (ASVs) from COI metabarcoding
studies represent genuine biological sequences or sequencing/PCR artefacts.

Given a set of query ASVs with their taxonomy, the package:

1. Downloads COI reference sequences from BOLD up to the family level
2. Aligns them at species, genus, and family level using amino acid profiles
3. Checks and improves alignments iteratively
4. Compares each ASV's amino acid sequence against the reference database to
   flag unusual substitutions
5. Optionally maps raw reads back to each ASV with VSEARCH to assess codon
   quality
6. Classifies each ASV as likely real or as a potential artefact (lab error
   or NUMT) based on Grantham distance, hydrophobicity, and read quality

---

## Installation

**aasv** depends on several Bioconductor packages and is not available on
CRAN. Install it from GitHub using `remotes`:

```r
# Install Bioconductor dependencies first
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("Biostrings", "DECIPHER", "ShortRead"))

# Then install aasv from GitHub
remotes::install_github("alexology/aasv")
```

### Optional external dependency: VSEARCH

Read-back mapping to assess per-codon sequencing quality requires
[VSEARCH](https://github.com/torognes/vsearch).  
Download a Windows binary from the
[VSEARCH releases page](https://github.com/torognes/vsearch/releases) and
note the path to `vsearch.exe`. The rest of the pipeline runs fully without
it.

---

## Pipeline overview

```
Taxonomy + ASV table
        │
        ▼
align_species_seq()   ← downloads sequences from BOLD, aligns at species level
        │
        ▼
align_genus_seq()     ← merges species alignments into genus-level profiles
        │
        ▼
align_family_seq()    ← merges genus alignments into a family-level profile
        │
        ▼
check_alignment()     ← flags alignments with internal gaps
        │
        ▼
improve_alignment()   ← removes offending sequences and re-aligns iteratively
        │
        ▼
asv_functional_structure()   ← compares each ASV against the reference database
        │
        ▼
classify_asv()        ← returns a per-ASV verdict (IS_REAL TRUE / FALSE)
```

---

## Quick start

### Step 1 — Build the reference database

```r
library(aasv)

# A data frame with at least columns: family, genus, species
taxonomy <- read.csv("my_taxonomy.csv")

# Download and align sequences from BOLD
align_species_seq(
  taxonomy      = taxonomy,
  alignment_dir = "results/alignments",
  min_length    = 640,
  max_length    = 700
)

# Merge into genus-level profiles
align_genus_seq(alignment_dir = "results/alignments")

# Merge into family-level profiles
align_family_seq(alignment_dir = "results/alignments")
```

### Step 2 — Check and fix alignments

```r
# Inspect for internal gaps at the species level
bad <- check_alignment(
  taxon       = "Baetidae",
  folder_path = "results/alignments",
  tax_lev     = "species"
)

# Iteratively remove gap-causing sequences and re-align
if (length(bad) > 0) improve_alignment(bad)
```

### Step 3 — Functional analysis

```r
# asv_taxonomy must include columns: family, genus, species, ASV
# numeric sample-abundance columns are also expected

asv_taxonomy <- read.csv("my_asv_taxonomy.csv")

# Without VSEARCH (quality columns will be NA)
asv_functional_structure(
  asv_taxonomy  = asv_taxonomy,
  output_dir    = "results",
  alignment_dir = "results/alignments"
)

# With VSEARCH read-back mapping
asv_functional_structure(
  asv_taxonomy   = asv_taxonomy,
  output_dir     = "results",
  alignment_dir  = "results/alignments",
  vsearch_path   = "C:/tools/vsearch.exe",
  trimmed_folder = "data/trimmed_reads"
)
```

### Step 4 — Classify ASVs

```r
report <- classify_asv(
  output_dir          = "results",
  grantham_threshold  = 50,
  quality_threshold   = 0.999,
  include_hydro       = TRUE
)

# report has columns: ASV_id, tax_lev, IS_REAL, REASON, ...
head(report)
```

---

## Function reference

| Function | Description |
|---|---|
| `align_species_seq()` | Download COI sequences from BOLD and align at species level |
| `align_genus_seq()` | Merge species alignments into genus-level profiles |
| `align_family_seq()` | Merge genus alignments into a family-level profile |
| `check_alignment()` | Detect alignments with internal gaps |
| `improve_alignment()` | Iteratively remove gap-causing sequences and re-align |
| `re_align()` | Re-align an existing alignment with new parameters |
| `asv_functional_structure()` | Run functional analysis for all ASVs |
| `classify_asv()` | Classify ASVs as real or artefact |
| `aa_similarity_matrix()` | Compute pairwise codon similarity matrix |
| `aa_similarity_score()` | Score ASV substitutions against reference codons |

---

## How classification works

For each ASV, `classify_asv()` applies a hierarchical check:

1. **Alignment gap check** — if the ASV introduces internal gaps when aligned
   to the reference, it is flagged immediately (`IS_REAL = FALSE`)
2. **Substitution check** — for each amino acid position where the ASV
   carries a substitution not seen in the reference database:
   - *Codon position*: mutations at wobble positions (position 3) are less
     penalised than those at positions 1–2 (Forbidden Zone)
   - *Grantham distance*: radical amino acid changes (> threshold) are
     penalised
   - *Sequencing quality* (optional, requires VSEARCH): low-quality codons
     are penalised; `NA` quality (no VSEARCH) passes this check
3. **Hydrophobicity check** (optional) — the global hydrophobicity score and
   the local profile of the ASV are compared against the reference envelope

The worst result across all mutations and all available taxonomic levels
determines the final verdict.
---

## Acknowledgments
This work was funded by the European Union under the NextGeneration EU 
Programme within the Plan “PNRR - Missione 4 “Istruzione e Ricerca” - Componente
C2 Investimento 1.1 “Fondo per il Programma Nazionale di Ricerca e Progetti 
di Rilevante Interesse Nazionale (PRIN)” by the Italian Ministry of University 
and Research (MUR), Project title: “METAbarcoding for METAcommunities: towards 
a genetic approach to community ecology (META2) ”, 
Project code: 2022PA3BS2 (CUP D53D23008270006), MUR D.D. 
financing decree n. 1015 of 07/07/2023

---

## Citation

If you use **aasv** in your research, please cite:

> Laini A., Voyron S., Gruppuso L. (2024). *aasv: an R package for functional
> assessment of COI metabarcoding ASVs*. (preprint / in preparation)

---

## License

GPL (>= 2)
