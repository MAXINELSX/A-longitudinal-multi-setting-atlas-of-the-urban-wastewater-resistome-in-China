# Input formats and entry points

See README.md for installation, Fig. 1, Fig. 2 and Fig. 5 usage, and data requests.

## Fig. 3 and Fig. 4

### Calls

```r
source("Fig3_abcde.R")
fig3_main(c("--data-dir", "local_inputs/fig3", "--output-dir", "results/fig3",
            "--panels", "a,b,c,d,e"))
source("Fig4_cde.R")
run_fig4_cde("local_inputs/fig4", "results/fig4")
```

Command line: `Rscript Fig3_abcde.R --data-dir local_inputs/fig3 --output-dir results/fig3 --panels a,b,c,d,e` and `Rscript Fig4_cde.R local_inputs/fig4 results/fig4 c,d,e`.

### Fig. 3 inputs

All tables are TSV except optional canonical sample metadata (CSV or XLSX).

- a: `Fig3a_species_total_counts.tsv` (`species`, `total_count`); `Fig3a_species_by_sample_counts.tsv` (first column species, remaining columns sample IDs containing integer contig counts); `Fig3a_sample_species_contig_totals.tsv` (`sample`, `total_species_contigs`, all species-assigned contigs, not just displayed species); `Fig3a_species_by_ARG_type_relative.tsv` (first column species, remaining columns ARG types with relative composition); `Fig3a_sample_metadata.tsv` (`sample`, `Setting` or `Sample_Type`, `City`). Explicit setting columns take priority over `Site`, which is a legacy setting-column alias. The heatmap pools counts and denominators before division within city-setting. Left totals and right compositions are provided as separate external plot inputs.
- b: `02_Fig3b_type_context_setting_counts.tsv` (`arg_type`, `Context`, `Setting`, `n_ARG_ORFs`). Contexts are Chromosome, Plasmid and Phage (Virus accepted). These counts must use one final ARG annotation per ORF. Distinct ORFs are not deduplicated by contig or ARG type. Top ten types plus Others are plotted with full-count denominators for nodes and setting-specific denominators for right-hand percentages.
- c: `08_Fig3c_sample_fractions.tsv` (`Sample`, `chromosome_ARG_ORFs`, `plasmid_ARG_ORFs`, `phage_ARG_ORFs`, `mobile_ARG_ORFs`, `total_ARG_ORFs`, `mobile_fraction`, plus embedded `City`, `Setting`, `Sample_Date`, `PhysicalSite`, `SamplingMonth`). `total_ARG_ORFs` is the classified denominator. A separate `Fig3_sample_metadata.csv` or `--metadata-file` CSV/XLSX may provide canonical metadata instead. Undefined denominators are excluded and logged; valid zero and one fractions remain. Optional `unclassified_ARG_ORFs` and `final_ARG_ORFs` are checked for conservation. The original-scale LMM uses Setting + City + SamplingMonth + (1|PhysicalSite), REML and six two-sided asymptotic Wald contrasts with BH adjustment.
- d: `Fig3d_MAG_total_burden.tsv` (`site`, `species`, `total_burden`, one row per MAG); `Fig3d_species_burden_summary.tsv` (`site`, `species`, `n_mag`, `mean_burden`, `prop_pathogen_mags`). Here `site` denotes wastewater setting. Per-species Kruskal–Wallis and Dunn/BH outputs are exported. The supplied MAG collection is not quality-filtered or treated as a mixed-model dataset by this plotting script.
- e: four external tables named `ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_` plus `bubble_pct.tsv` (`site`, `species`, `arg_type_plot`, `n_type_carrying_mags`, `pct_mags`), `site_species_type_mean_burden.tsv` (`site`, `species`, `arg_type` or `arg_type_plot`, `mean_type_burden`), `species_used_by_site.tsv` (`site`, `species`, `n_mag_site`), `argtypes_used.tsv` (`arg_type_plot` or `arg_type`). Supply `Fig3e_MAG_ARG_presence.tsv` (`mag_uid`, `arg_type`; alternatively `Sample` and `bin_id`) and `Fig3e_MAG_registry.tsv` (`mag_uid`, `site`, `species`) for union carrier counts. Alternatively `--bubble-premerged true` accepts counts already merged as a MAG union. Marginal category counts must not be added to make an Others carrier count. Mean type burdens must use the same full species-setting MAG denominator before summing type means.

### Fig. 4 inputs

- c: `Fig4c_cluster_data.tsv` (`Cluster`, `Clinical`, `Wastewater`; optional `Total`, `WW_percent`). One row per mixed-source cluster in intended display order.
- d: `Fig4d_pairs.tsv` (`Pair_ID`, `Clinical_display_ID`, `Clinical_specimen`, `Site`, `Clinical_city`, `Wastewater_city`, `same_city`, `Wastewater_setting`, `ANI_min`, `fixed_SNP_Mbp`, `callable_fraction`, `Direct_remap`, `inStrain_partial`, `Clinical_ARG_count`, `Wastewater_ARG_count`, `Shared_ARG_count`, `High_priority_ARG`). This is an externally preselected display table containing the highest-ranked clinical candidate for each wastewater MAG. The script does not select candidates or reconstruct screening edges. Its SNP-density sort controls row order only. `High_priority_ARG` uses sections `shared: ...; clinical-only: ...`. Read support fields are zero, one or missing. ANI is percent and callable fraction is 0–1.
- e: `e_carbapenemase_overlap.tsv` (`Evidence`, `Allele`, `Supported`, `Support_score`); `e_mobile_genetic_context.tsv` (`Context`, `Contigs`, `Percent`); `e_selected_taxa.tsv` (`Taxon_group`, `Contigs`, `Percent_of_classified`). Percentages are recalculated from the supplied unique-contig counts, not ORF counts. The input category sets follow the panel labels.

### Dependencies

Fig. 3: data.table, dplyr, tibble, ComplexHeatmap, circlize, RColorBrewer, lme4, multcompView, ggplot2, stringr, FSA, readr, scales and tidyverse. readxl is needed only for external XLSX metadata. networkD3/htmlwidgets/jsonlite enable the interactive Sankey; the vector SVG requires base R only. rsvg optionally creates PNG/PDF from SVG. Fig. 4: ggplot2, patchwork and svglite; Cairo-enabled graphics for PDF export. ComplexHeatmap is a Bioconductor package. No script downloads restricted data or installs packages automatically.

## Fig. 6

```r
source("Fig6_a_to_h.R")
cfg <- fig6_network_config()
network_input <- fig6_network_prepare_lmm(network_data, cfg)
network_lmm <- fig6_network_fit_lmm(network_input, "Total", cfg)
network_rf <- fig6_network_fit_rf(fig6_network_prepare_rf(network_input, cfg), "Total", cfg)
run_fig6("output/Fig6", network_rf = list(Total = network_rf),
         network_lmm = list(Total = network_lmm), panels = c("a", "e"))
```

Repeat with the Tier I outcome and use `TierI` named entries for panels b/f. `network_data` requires numeric `arg_abundance`, `city`, globally unique `site_nm`, `setting`, `season`, `sample_date`, and one numeric column per configured predictor and exposure window (`_lag0`, `_lag1`, `_lag2`, `_lag3`). The four suffixes represent sampling day and preceding 30/60/90-day means, respectively. Optional `fig6_network_lags(samples, daily)` constructs windows from daily columns `city,date,Index,value`. Duplicate city-date-predictor records must be resolved explicitly before this step.

```r
inputs <- fig6_hospital_model_inputs("input/hospital_rf")
rf_plan <- fig6_fit_hospital_rf(inputs, dry_run = TRUE)
hospital_rf <- fig6_fit_hospital_rf(inputs, dry_run = FALSE)
run_fig6("output/Fig6", hospital_rf = hospital_rf, panels = c("c", "d"))
```

Hospital RF input files:

| File | Required content |
|---|---|
| `keys.tsv` | Unique `hospital_code,quarter` keys; quarter is an integer code. |
| `drug_exposures.tsv` | The same keys and numeric, original-scale combined-route DDD intensities in `DRUG__001` through `DRUG__097`. |
| `hospital_context.tsv` | The same keys and `HOSP_mean_patient_days,HOSP_mean_los,HOSP_mean_discharges,HOSP_mean_quarterly_amount,QUARTER,HOSPITAL,CITY,HOSP_type`. Mean descriptors are constant within hospital. Expenditure is in RMB (yuan). |
| `outcomes_log10.tsv` | The same keys and all catalogue outcome columns, already transformed as log10(hospital-quarter mean abundance + 1e-8). Sum constituent subtypes within samples before original-scale hospital-quarter averaging. |
| `predictor_dictionary.tsv` | Ordered `source_id,source_label_English,predictor_family,drug_broad_class`; drug predictors followed by the listed context fields. |
| `outcome_catalog.tsv` | Unique `outcome_id,outcome_level`, where levels are `total,type,clinical` and total is `TOTAL__all_ARG`. |

The alternative `hospital_rf_dir` plotting route reads saved `total_importance.tsv`, `group22_importance.tsv`, `predictor_dictionary.tsv`, `outcome_catalog.tsv`, `total_run_settings.tsv`, and `group22_run_settings.tsv`. Settings require `trees,num.rep,quick_test`; the two importance tables require `source_id,outcome_id,%IncMSE_unscaled,%IncMSE.pval`. They are user-supplied results, not bundled data.

```r
hospital_lmm <- fig6_fit_hospital_lmm(
  outcomes, list(oral = oral, injection = injection, combined = combined),
  hospital_metadata, dataset = "group22", outcome_scale = "linear"
)
run_fig6("output/Fig6", hospital_lmm = hospital_lmm,
         hospital_selection = plot_selection, panels = c("g", "h"))
```

Hospital LMM inputs have unique `hospital_code,quarter` keys. Outcome matrices contain complete planned outcome columns in original copies-per-cell hospital-quarter means (`linear`) or the declared already-transformed values (`log10`). Exposure tables contain all planned antibiotic columns on the original DDDs-per-100-patient-days scale. Metadata also require original `HOSP_patient_days`. The plot selection requires `dataset,outcome_id,predictor_id,profile,panel,section,order_key,outcome,drug`; identifiers match fitted results, labels are supplied explicitly, and the layout follows the original panel arrangement. No rows are selected automatically by significance. The saved-results route `hospital_lmm_dir` expects complete planned result families, not display-only tables.

Dependencies: `ggplot2,lme4,lmerTest,rfPermute,randomForest,cowplot,ragg`; `VIM` for RF-only kNN imputation, `data.table` for saved hospital result families, `patchwork` for a combined figure, and optional `performance,svglite`.

## Extended Data Figs. 1–3

```r
source("Extended_Data_Figures_1_to_3.R")
cfg <- extended_data_config("input", "output/extended_data")
run_extended_data(cfg, figures = c(1, 2, 3))
```

| File | Required columns |
|---|---|
| `global_benchmark.tsv` | `sample_id,country,region,source,abundance`; source is `study` or `reference`. Supply the final curated comparison set. |
| `subtype_abundance.tsv` | `sample_id` plus numeric subtype columns on original copies-per-cell scale. |
| `sample_metadata.tsv` | `sample_id,setting,city,month,physical_site`; one row per sample. |
| `subtype_family_map.tsv` | `subtype,family`; supply all constituent subtypes for each requested family, not only the plotted selection. |
| `selected_subtypes.csv` | `family,subtype,plot_order,display_label`; use the retained initial-screen selection and original plotting order. |

Settings are `Hospital,WWTP,Community,Wet market`; families are `OXA,KPC,VIM,NDM,IMP` (a `bla` prefix is accepted by preparation). The default selection mode is `selected`; it does not reselect subtypes from current model results. The optional `screen` mode expects an independently supplied original-screen table with `setting_r2,city_r2,setting_q`. Do not substitute new repeated-site-model results for the original screen.

Blank abundance cells stop by default. Set `cfg$blank_policy <- "nondetection"` only for an input whose blank cells are confirmed nondetections; genuine missing measurements must not be silently converted. `cfg$abundance_orientation <- "subtype_rows"` accepts a transposed matrix with a `subtype` identifier column.

Dependencies: `ggplot2,lme4,multcompView,patchwork`; `readxl` only for Excel input. No actual sample, hospital or genome data are embedded in either script.

## Supplementary Figures S2–S14

Load `SI_S2_to_S14.R`. It defines functions without reading study files or fitting models automatically. `si_check_packages()` lists dependencies; `si_install_packages()` installs missing packages when called explicitly.

### Shared inputs

- `abundance`: a numeric sample-by-subtype matrix with unique sample IDs as row names and subtype IDs as column names. Entries must be finite and non-negative. Missing measurements must be resolved explicitly, not silently replaced by zero.
- `metadata`: `sample_id`, `city`, `setting`, `sample_date` (`YYYY-MM-DD`) and `physical_site`. Site IDs must be globally unique and belong to one city and setting. Setting values are `Hospital`, `WWTP`, `Community` and `Wet market`.
- `meta <- si_metadata(metadata)` returns `meta$lower` for S2/S4/S5/S7/S8/S9 and `meta$upper` for S3/S6/S10/S14. The lower table includes `month`; the upper table uses `Sample`, `City`, `Setting` and `Date`.
- `tier1_subtypes`: the final subtype IDs used for Tier I analyses.
- `subtype_dictionary`: one row per subtype, with `Subtype` and `ARG_type`.

### Figure entry points

| Figure | Analysis and plotting functions | Additional input |
|---|---|---|
| S2 | `si245_s2(abundance, meta$lower, tier1_subtypes)`; `si245_plot_s2(result)` | None |
| S3 | `si3_monthly(abundance, meta$upper, subtype_dictionary)` | None; returns `$plot` and descriptive tables |
| S4 | `si245_s4(abundance, meta$lower, tier1_subtypes)`; `si245_plot_s4(result)` | None; city-specific negative-binomial mixed models |
| S5 | `si245_s5(group_abundance, meta$lower, output_dir=...)`; `si245_plot_s5(result)` | `group_abundance`: sample-by-ARG-group matrix. Alternatively derive it with `si245_group_abundance(abundance, mapping)`, where `mapping` has unique `Subtype`, `Feature` pairs |
| S6 | `si6_family_bars(abundance, meta$upper, variant_map, selected_variants)` | `variant_map`: unique `Subtype`, `Family`, `Variant`; `selected_variants`: `Family`, `Variant`. Default families: `OXA`, `KPC`, `VIM`, `NDM`, `IMP` |
| S7 | `si789_prepare_s7(annotations, meta$lower)`; `si789_plot_s7(prepared)` | Annotation schema below |
| S8 | `si789_run_s8(traits=traits, output_dir=...)`, or `si789_run_s8(annotations=annotations, metadata=meta$lower, core_membership=core_membership, output_dir=...)` | Pooled subtype metrics or annotation records plus core membership, as defined below |
| S9 | `si789_prepare_s9(annotations, meta$lower, host_sets)`; `si789_fit_s9(prepared$fractions, meta$lower)`; `si789_plot_s9(prepared$fractions, meta$lower, fitted)` | `host_sets`: unique `species`, `host_group`, `source_reference`. Host groups: `clinical_pathogen`, `opportunistic_environment`; other assigned species remain in the denominator |
| S10 | `si10_monthly(abundance, meta$upper, subtype_dictionary, tier1_subtypes)` | None; returns `$plot` and descriptive tables |
| S11 | `si111213_s11_activity(activity)`; `si111213_s11_amu(intensities)`; `si111213_plot_s11(activity_result, amu_result)` | Activity and consumption schemas below |
| S12 | WWTP preparation, RF and LMM calls below; `si111213_plot_s12(rf_total, rf_tier1, lmm_total, lmm_tier1)` | Daily covariates and wastewater physicochemical records |
| S13 | `si111213_plot_s13(network_total, network_tier1, wwtp_total, wwtp_tier1)` | Complete LMM coefficient tables for each outcome and analysis scope |
| S14 | `si14_profile_features(abundance, meta$upper, marker_map, tier1_subtypes)`; `si14_fit(prepared$profiles)`; `si14_plot(fitted$coefficients)` | `marker_map`: unique `Marker`, `Subtype` pairs |

Use `si_save_plots(plots, directory, prefix)` for plot lists. S8's runner exports its combined figure and statistical tables directly. S5 exports primary comparisons and site-sensitivity tables when `output_dir` is supplied; `si245_export_s5(result, directory)` can export an existing result without refitting.

### Annotation and S8 inputs

`annotations` contains character columns `sample_id`, `contig_id`, `arg_occurrence_id`, `arg_subtype`, `arg_type`, `host_species`, `host_rank` and `carrier`. `arg_occurrence_id` identifies one final ARG annotation and is unique within a sample. Species-level assignments use `host_rank="S"`; absent host assignments can be `NA`. Carrier values are `chromosome`, `plasmid`, `phage` or `unknown`. Contig classifications must agree within each sample–contig.

S7 counts unique sample–contig–ARG-type combinations and requires resolved carrier assignments. S8 deduplicates sample–contig–subtype combinations before pooling by subtype. S9 calculates ARG-occurrence and unique-contig fractions separately.

`core_membership` contains `arg_subtype`, `setting` and logical `is_core`, with an explicit row for every subtype in every setting. It can be generated from the supplied abundance matrix with `si789_core_from_abundance(abundance, meta$lower)`.

The S8 `traits` long format contains `arg_subtype`, `core_group`, `metric` and numeric `value`, with one row per subtype–metric. Core groups are `four_setting_shared_core`, `setting_restricted_core` and `non_core`. Metrics are `dominant_host_fraction`, `host_species_richness`, `host_shannon`, `mobile_fraction`, `plasmid_fraction` and `phage_fraction`. Undefined values are `NA`; valid zeros and ones are retained. A wide table can instead use `Subtype`, `core_group`, `dominant_host_fraction`, `host_species_richness`, `host_shannon`, `mobile_fraction_all`, `plasmid_fraction_all` and `virus_fraction_all`.

S8 performs Kruskal–Wallis and two-sided Dunn comparisons on original-scale pooled metrics, with BH adjustment across three comparisons per metric. Only the richness display uses `log10(x+1)`. Sample-by-subtype inputs are rejected by the S8 fitting function.

### Hospital activity and consumption for S11

`activity` contains `city`, `hospital`, `quarter`, `patient_days` and `discharges`, with unique hospital–quarter keys. `intensities` contains `city`, `hospital`, `quarter`, `route`, `category` and numeric `intensity` in DDDs per 100 patient-days. Routes are `oral`, `injection` and `combined`. If supplying raw DDDs, use `si111213_ddd_intensity(ddd, activity)` with the same keys, oral/injection routes and a `ddd` column. Missing route records are not treated as zero; a known absence of use must be supplied as an explicit zero.

### Urban covariates for S12 and S13

`samples` contains `sample_id`, `city`, `site_nm`, `setting`, `sample_date`, `season` and numeric `arg_abundance`. Run the workflow separately for total and Tier I outcomes. `daily` contains `city`, `date`, `Index`, `value`, with one value per city–date–predictor. Predictor names are defined in `fig6_network_config()$classes`. `quality` contains the same columns for `CODCr`, `NH3_N`, `TSS`, `pH`, `flow` and `wastewater_temperature`; multiple records within a city–date are averaged. Only the first four enter the models.

```r
prepared <- si111213_prepare_wwtp(samples, daily, quality)
lmm <- si111213_fit_wwtp_lmm(prepared, outcome_label)
rf_input <- si111213_prepare_wwtp_rf(prepared)
rf <- si111213_fit_wwtp_rf(rf_input, outcome_label)
```

Pass `rf$plot_data` and `lmm$results` to the plotting functions. Saved importance tables can be converted with `si111213_rf_table()`: supply a predictor identifier (`Predictor`, `Index` or `variable`), `%IncMSE` and nominal permutation P (`%IncMSE.pval`, `P` or `p`). LMM result tables require `outcome`, `class`, `group`, `lag`, `coef`, `ci_low`, `ci_high`, `p` and `fdr`; `lag0`, `lag1`, `lag2`, `lag3` denote the sampling-day, preceding-30-day, preceding-60-day and preceding-90-day windows. Keep complete class-by-window results when calculating BH adjustment before selecting display panels.

RFs use lag0 predictors, 1,000 trees and 1,000 response permutations. RF predictor imputation is not used for LMMs. Inner RF labels retain signed relative importance; outer class shares use non-negative values and sum to 100%.

## Tier I analyses: Supplementary Methods S3-5 and S3-6

All CSV files use a header row and UTF-8 encoding; a UTF-8 byte-order mark is accepted. Identifiers must match exactly. Scripts do not install dependencies, download study data or infer sample metadata from identifiers. Run Python without `-O` so that validation assertions remain enabled. Commands are in README.md; each script also accepts `--help`.

### Threshold combinations: S5-2 and S5-3

`TierI_threshold_sensitivity.py` requires:

| Argument | Required fields or content |
|---|---|
| `--core` | `arg_subtype,core_flag,Site.Hospital,Site.WWTP,Site.Community,Site.Wet.market`; one row per subtype in the original core union |
| `--assessed` | `arg_subtype,arg_type,risk_level,species_list,host_breadth,chromosome_contigs,plasmid_contigs,virus_contigs,mobile_contigs,total_contigs,mobility_ratio` and the four `Site.*` columns; one row per originally eligible subtype |
| `--original-tier1` | `arg_subtype`; one row per original Tier I candidate, in original roster order. The original complete trait export is also accepted; overlapping fields must agree with `--assessed` |
| `--hosts` | `who_species.txt`: one eligible species name per line, as used in the original analysis |
| `--output` | Directory for result CSV and JSON files |

Core flags retain the source-file convention: an empty field, `NA`, `NaN`, `FALSE` or `0` means not core; any other nonempty value indicates core membership. The `Site.*` fields describe original **core membership**, not detection alone. An assessed subtype must belong to the core union, and its setting flags must agree with that union. Tier labels use `Level I`, `Level II`, `Level III` and `Level IV` in the source tables. Subtypes not in the originally assessed set remain unassessed, with blank scenario fields rather than zero eligibility.

`species_list` is a semicolon-separated list of distinct eligible host species. Its length must equal `host_breadth`, and every species must occur in `who_species.txt`. Host breadth is not the number of all taxonomically assigned species. This workflow uses supplied host evidence; it does not repeat taxonomic assignment.

The subtype-level context counts use distinct sample–contig–subtype records. `mobile_contigs = plasmid_contigs + virus_contigs`, and `total_contigs = chromosome_contigs + plasmid_contigs + virus_contigs`. Thus `total_contigs` is the classified-context denominator; unknown contexts are excluded. `mobility_ratio` must agree with the count-derived fraction. These subtype-carriage counts are separate from Fig. 3b,c's ORF-counting unit.

The nine combinations cross mobility ≥0.60/0.70/0.80 with host breadth ≥2/3/4. The baseline is mobility ≥0.70 and host breadth ≥3. Exact integer comparisons of the numerator and denominator avoid floating-point ambiguity at the mobility thresholds. Baseline membership must reproduce the supplied original Tier I roster. Core membership and original assessment eligibility remain fixed throughout.

Outputs are `TierI_threshold_summary.csv` (S5-2), `TierI_subtype_membership.csv` (S5-3) and `sensitivity_results.json` (checks and source hashes). Scenario flags are 1=selected and 0=not selected; blank=not originally assessed. Retention uses the original Tier I list as denominator. Jaccard overlap uses the intersection divided by the union of the original and alternative lists. Output fields containing `102` or `10` retain the published panel labels; numerical denominators are derived from the supplied original roster, not hardcoded result counts.

### Fixed-panel preparation

`TierI_prepare_inputs.py` reads:

| Argument | Required fields or content |
|---|---|
| `--metadata` | `Sample,City,Setting,PhysicalSite,Sample_Date,SamplingMonth`; one canonical row per sample |
| `--abundance-long` | `Sample,Subtype,abd,Site,City`; one row per sample–subtype combination in the fixed panel |
| `--original-tier1` | The same original roster used in the threshold workflow; `arg_subtype` is required |
| `--membership` | `TierI_subtype_membership.csv` produced by the threshold workflow |
| `--workbooks` | Optional paths to the original setting-specific XLSX matrices, separated by spaces; quote paths containing spaces |
| `--output` | Directory for the prepared CSV files and input audit |

In metadata, `Setting` is exactly `Hospital`, `WWTP`, `Community` or `Wet market`. `PhysicalSite` is globally unique and maps to one city and setting. `Sample_Date` uses `YYYY-MM-DD`; `SamplingMonth` uses `YYYY-MM` and must agree with the date. In the long abundance table, `Site` is the source-file name for **wastewater setting**, not physical site. It must agree with canonical `Setting`; `City` must also agree. `abd` is original-scale copies per cell.

The long table must explicitly contain every canonical sample–candidate combination, including zero abundance for confirmed nondetections. Missing rows, NaN, infinite values and negative values cause an error; they are not replaced by zero. If `--workbooks` is supplied, every sheet starts with `Subtype`, followed by sample columns. Blank abundance cells in these source workbooks are treated as confirmed nondetections for reconciliation, following the source-data convention. Metadata blanks are not imputed. The script compares workbook cells against the long table, rather than silently replacing either source.

Preparation checks that the full-period setting-specific recurrence reproduces original core flags for the fixed candidates. A mismatch is reported in `matrix_audit.json` and stops the workflow; the script does not redefine original membership.

Prepared files:

- `analysis_metadata.csv`: canonical metadata plus `Period`, with labels `2024-11 to 2025-05` and `2025-06 to 2025-12`.
- `analysis_abundance_matrix.csv`: first column `Sample`, followed by one original-scale numeric column per fixed candidate.
- `fixed_panel_traits.csv`: `arg_subtype,arg_type,hospital_core,wwtp_core,community_core,wet_market_core,nonhospital_additional_core` and the carried-through threshold-analysis fields. Core flags are 0/1. A non-hospital addition is core in a non-hospital setting but not in hospitals; it need not be absent from hospital samples.

The prepared CSV files can be supplied directly to `TierI_internal_consistency.py --input-dir ...`. Preserve the original metadata row order and original candidate-roster order for exact reproduction of seeded subsampling. Preparation preserves both orders. Site IDs must retain their original spelling because common sites are traversed in sorted order. Filenames can be relocated without changing their content.

### Fixed-panel internal consistency: S5-4–S5-6

The original panel is fixed throughout. Detection is abundance >0. Within each setting and city, month or period, recurrence requires prevalence ≥0.70 and arithmetic mean abundance ≥10⁻⁵ copies per cell; zeros contribute to both calculations. Prevalence comparisons allow a numerical tolerance of 10⁻¹². Unsampled strata are omitted, not coded as zero.

Sample-weighted profiles give each sample equal weight. Equal-site profiles calculate prevalence and mean abundance within each physical site and then average these site values with equal weights. A city or month is counted as recurrent if at least one sampled setting meets both thresholds; the non-hospital version excludes hospital settings.

The early period is 1 November 2024–31 May 2025, and the late period is 1 June–31 December 2025. Common sites occur in both periods. Temporal persistence requires the same setting to meet both thresholds in both periods. For non-hospital additions, the corresponding measure requires the same non-hospital setting. The workflow evaluates all available sites and common sites, with sample and equal-site weighting. Each setting must have at least one common site.

Balanced subsampling draws `min(n_early, n_late)` samples without replacement from each period within each common site. Within each draw, site means are weighted equally within setting. Defaults are 1,000 repetitions and seed 20260925, using NumPy `default_rng` (PCG64). Both the draws and per-candidate persistence counts are saved. Use the original row order, site IDs and package versions to reproduce the individual draws. Percentiles and retention frequencies summarize subsampling, not confidence intervals, hypothesis-test P values or external validation.

Coverage compares the original hospital-core subset with the full fixed panel on the same samples; it is not an equal-sized-panel optimization. Period-specific prevalence-rank correlations use average ranks for ties and are descriptive, with no hypothesis test.

| Saved analysis file | Content and table mapping |
|---|---|
| `subtype_consistency_summary.csv` | Per-subtype city/month recurrence, temporal flags and balanced-draw retention; selected columns form S5-4 |
| `subtype_stratum_profiles.csv` | Setting-specific prevalence, abundance, denominators and recurrence flags; S5-5 |
| `balanced_subsampling_summary.csv` | Subsampling minima, percentiles, medians and maxima; S5-6 block 2 |
| `period_concordance.csv` | Setting-specific temporal overlap and descriptive rank correlations; S5-6 block 3 |
| `stratum_coverage.csv` | Matched-sample panel coverage; S5-6 block 4 |
| `balanced_subsampling_repetitions.csv` | Counts in individual draws, underlying the summary block |
| `sample_coverage.csv`, `balanced_site_manifest.csv` | Local sample/site-level intermediate outputs; not public code-package contents |
| `results.json`, `run_settings.json` | Summary, input hashes, seed, repetitions and software versions |

`TierI_export_tables.py` selects and orders the S5-4 columns as in the supplementary table, copies the S5-2/S5-3/S5-5 numerical outputs, and writes four separate CSV blocks for S5-6. Block 1 (`S5-6_temporal_schemes.csv`) combines saved subtype-level persistence flags with sample/site denominators from prepared metadata. The exporter performs no fitting or subsampling. It does not replace S5-1, redraw figures, alter the manuscript or embed local input records in the repository. Keep generated sample/site-level outputs outside public uploads.
