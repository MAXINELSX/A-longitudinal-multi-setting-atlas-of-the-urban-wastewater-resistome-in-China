# A longitudinal multi-setting atlas of the urban wastewater resistome in China

## Overview

This repository contains analysis and plotting scripts, selected figure summaries, and hospital mixed-effects model results accompanying the manuscript **A longitudinal multi-setting atlas of the urban wastewater resistome in China**.

The study analysed 2,614 wastewater metagenomes collected over 14 months from 130 sites in 12 Chinese cities, covering hospitals, communities, wet markets, and wastewater treatment plant (WWTP) influent. Analyses address ARG abundance and recurrence, putative bacterial hosts and predicted genetic contexts, wastewater–clinical genomic comparisons, candidate surveillance targets, and associations with medication use and urban covariates.

The repository contains R scripts for the main, supplementary, and extended-data figures, together with Python scripts for Tier I sensitivity and internal-consistency analyses. Detailed methods and interpretation are provided in the manuscript and Supplementary Information.

## Repository contents

### Analysis and plotting scripts

| File | Coverage |
| --- | --- |
| [fig1_btoe.R](fig1_btoe.R) | Fig. 1b–e |
| [fig2_a-h.R](fig2_a-h.R) | Fig. 2a–h |
| [Fig3_abcde.R](Fig3_abcde.R) | Fig. 3a–e |
| [Fig4_cde.R](Fig4_cde.R) | Fig. 4c–e |
| [Fig5_b_to_g.R](Fig5_b_to_g.R) | Fig. 5b–g |
| [Fig6_a_to_h.R](Fig6_a_to_h.R) | Fig. 6a–h |
| [SI_S2_to_S14.R](SI_S2_to_S14.R) | Supplementary Figs. S2–S14 |
| [Extended_Data_Figures_1_to_3.R](Extended_Data_Figures_1_to_3.R) | Extended Data Figs. 1–3 |
| [TierI_threshold_sensitivity.py](TierI_threshold_sensitivity.py) | Supplementary Methods S3-5; Supplementary Data S5-2 and S5-3 |
| [TierI_prepare_inputs.py](TierI_prepare_inputs.py) | Input preparation and optional workbook reconciliation for S3-6 |
| [TierI_internal_consistency.py](TierI_internal_consistency.py) | Supplementary Methods S3-6; Supplementary Data S5-4–S5-6 |
| [TierI_export_tables.py](TierI_export_tables.py) | Table-numbered CSV exports of saved S5-2–S5-6 results |

Maps, study-design diagrams, externally rendered trees, and final page assembly are prepared separately.

### Public result files

| File | Contents |
| --- | --- |
| [Fig3b.txt](Fig3b.txt) | Proportions for ARG-type, predicted genetic-compartment, and wastewater-setting combinations underlying Fig. 3b |
| [Fig3c.txt](Fig3c.txt) | Descriptive summaries of sample-level mobile-associated ARG fractions, including sample number, mean, median, quartiles, range, and standard deviation |
| [Fig3d_species_burden_summary.tsv](Fig3d_species_burden_summary.tsv) | Species-by-setting summaries of ARG burden: mean, median, and maximum |
| [Supplementary_Table_ALL_LMM_results_part1of2.xlsx](Supplementary_Table_ALL_LMM_results_part1of2.xlsx) | Hospital mixed-effects model results, part 1 of 2 |
| [Supplementary_Table_ALL_LMM_results_part2of2.xlsx](Supplementary_Table_ALL_LMM_results_part2of2.xlsx) | Hospital mixed-effects model results, part 2 of 2 |

The two hospital result workbooks should be consulted together. They provide the complete hospital model results referenced in the manuscript, including coefficients, confidence intervals, nominal and Benjamini–Hochberg-adjusted P values, model status, and diagnostic codes. Supplementary Data S7-5 contains the selected associations displayed in Fig. 6g,h.

### Supporting documentation

| File | Purpose |
| --- | --- |
|Required input schemas, analysis entry points, and module-specific usage |
| [who_species.txt](who_species.txt) | The predefined clinically relevant host set used in this study |
| [requirements-tieri.txt](requirements-tieri.txt) | Pinned Python dependencies for the Tier I workflows |
| [TIERI_VERSION.md](TIERI_VERSION.md) | Description of the Tier I workflow snapshot |
| [TIERI_SHA256SUMS.txt](TIERI_SHA256SUMS.txt) | File checksums associated with that snapshot |

## Reading the figure summaries

Fractions in `Fig3b.txt` and `Fig3c.txt` are reported on a 0–1 scale; multiply by 100 to express them as percentages.

- **Fig. 3b:** Each fraction uses the pooled number of classified ARG annotations across study samples as its denominator. Individual annotations are counted separately, including those on the same contig. `Display_ARG_type` maps the original ARG categories to the displayed categories. The values describe pooled composition, rather than an average of sample-level fractions.
- **Fig. 3c:** The mobile-associated fraction is the number of ARG annotations assigned to predicted plasmid or phage contexts divided by all classified ARG annotations within a sample. Samples with an undefined denominator are excluded; valid fractions of zero and one are retained. The summaries include 2,607 samples. The `All settings` row pools these samples rather than giving each setting equal weight.
- **Fig. 3d:** In the public summary, `site` denotes the wastewater setting, not an individual sampling location. `mean_burden`, `median_burden`, and `max_burden` summarise ARG burden across the corresponding species–setting MAG collection.

These files are descriptive summaries. Re-fitting statistical models requires the underlying observations and relevant metadata.

## Software requirements

### R

Install the packages required by the selected analysis. Dependencies across the R scripts include:

```text
readr, readxl, openxlsx, data.table, tidyverse, dplyr, tidyr,
tibble, purrr, stringr, forcats, ggplot2, scales, cowplot,
patchwork, ggrepel, ggnewscale, ggalluvial, VennDiagram,
ComplexHeatmap, circlize, RColorBrewer, networkD3, htmlwidgets,
jsonlite, vegan, permute, lme4, lmerTest, multcompView,
rstatix, FSA, performance, randomForest, rfPermute, VIM, nnls
```

`ComplexHeatmap` is distributed through Bioconductor. Some graphics-export options additionally require `svglite`, `ragg`, or `rsvg`. Module-specific requirements are described in `INPUTS.md` and the scripts. Use the software versions reported for the relevant analysis in the manuscript and Supplementary Information, and retain session information when running an analysis.

### Python

The Tier I workflows use Python 3.12. From the repository root, install the pinned dependencies with:

```sh
python -m pip install -r requirements-tieri.txt
```

## Running the analyses

Run commands from the repository root and supply the required input files locally. Paths beginning with `data/`, `local_inputs/`, or `results/` in the examples are local paths, not bundled study-data directories.

Use a separate R session for each figure script. Sourcing a script loads its functions; call the relevant entry point to start the workflow. Check input-validation messages, model diagnostics, and comparison tables when reviewing generated outputs.

### Fig. 1

```sh
Rscript fig1_btoe.R --metadata data/fig1/sample_metadata.csv --subtype-dir data/subtypes --output results/Fig1_b-e.pdf
```

Required metadata columns are `Sample`, `Clean_Reads_GBp`, `City`, `Sample_Type`, and `Sample_Date`. Dates may use `YYYYMMDD` or `YYYY-MM-DD`.


### Fig. 2

```r
Sys.setenv(
  FIG2_DATA_ROOT = "data/fig2",
  FIG2_METADATA_FILE = "data/fig2/sample_metadata.xlsx",
  ARG_EXTERNAL_DATA_ROOT = "data",
  ARG_RESULTS_DIR = "results/fig2"
)
source("fig2_a-h.R")
run_fig2()
```

Required metadata columns are `Sample`, `City`, `Sample_Date` (`YYYYMMDD`), and `Sample_Rename`. The terminal date suffix in `Sample_Rename` is removed to obtain the physical-site identifier.

Panel d additionally requires the external comparison table at `data/c_merge/final_global.xlsx`. To omit that panel, set `FIG2_PANELS` to `"a,b,c,e,f,g,h"` before sourcing the script. Panels b and c run together.

### Fig. 5

```sh
Rscript Fig5_b_to_g.R data/tier1 data/sample_metadata.xlsx data/prioritization results/fig5 bcdefg
```

The abundance workbooks follow the Fig. 2 layout and contain the selected Tier I subtype roster. Set `FIG5_SUBTYPE_SUFFIX` if needed. Metadata may instead be a CSV containing `Sample`, `Setting`, `City`, `PhysicalSite`, `Sample_Date`, and `SamplingMonth`.

The prioritization directory contains `subtype_risk_4level_core_only_final.csv`, `subtype_level_I_core_only.csv`, `core_subtype_list.csv`, and `panel_subtype_mobility_summary.csv`.

### Other figure workflows

For Fig. 3, Fig. 4, Fig. 6, Supplementary Figs. S2–S14, and Extended Data Figs. 1–3, follow the entry points and input definitions in [INPUTS.md](INPUTS.md).

## Tier I sensitivity and internal consistency

These workflows correspond to Supplementary Methods S3-5 and S3-6. Threshold sensitivity uses the original core membership and assessment eligibility. Internal-consistency analyses retain the original candidate panel. Supplementary Data S5-1 remains the original priority list.

The workflow order is:

1. `TierI_threshold_sensitivity.py`: assess sensitivity to the specified selection rules.
2. `TierI_prepare_inputs.py`: prepare the fixed-panel input tables and, optionally, reconcile them against the source workbooks.
3. `TierI_internal_consistency.py`: calculate cross-city, temporal, common-site, and balanced-subsampling summaries.
4. `TierI_export_tables.py`: export the saved results into table-numbered CSV files.

Each script accepts `--help`. Required arguments and input schemas are documented under **Tier I analyses** in [INPUTS.md](INPUTS.md). If prepared inputs are already available, the preparation step can be skipped. If all analysis results have been generated, the export script can be run without repeating the analyses.


## Data availability

Public files in this repository are listed above. De-identified derived data supporting the study are also provided in the manuscript's Supplementary Data, including sample metadata, ARG abundance summaries, MAG annotations, wastewater–clinical comparison matrices, and statistical outputs. Accession numbers, available metadata, and publication links for the reused public wastewater metagenomes are listed in Supplementary Data S1-3.

Raw wastewater sequencing reads, clinical genome resources, detailed site-resolved metadata, and additional analysis inputs are not bundled with this repository. These include detailed annotation and sample-level count tables, individual-MAG records, complete wastewater–clinical cluster membership and screening edges, hospital-quarter inputs, and the external inputs required by the Tier I workflows.

Access is governed by the institutional and data-custodian requirements described in the manuscript's Data Availability statement, including restrictions concerning sensitive surveillance-site information and clinical resources. Requests should be directed to the corresponding authors listed in the manuscript and require the relevant approvals and compliance with applicable data-use conditions. During peer review, editors and reviewers may request access through the editorial office.


## Contact

For scientific questions and data-access requests, contact the corresponding authors tangsong@nieh.chinacdc.cn listed in the manuscript. 
