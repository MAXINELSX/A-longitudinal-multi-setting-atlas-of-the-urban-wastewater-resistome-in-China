# Fig. 2
suppressWarnings(Sys.setlocale("LC_ALL", "English_United States.utf8"))
set.seed(20260918)

# ----- fig2_data_models.R -----
# Fig. 2 audited cohort and repeated-site inference, 2026-09-18.
# No rounding of continuous abundance and no automatic model-family substitution.
validate_nb_response <- function(y) {
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(abs(y-round(y)) > 1e-8))
    stop('Negative-binomial responses must be nonnegative integer counts; copies/cell is continuous.')
  TRUE
}

resolve_abundance_missing <- function(x, policy='stop') {
  if (any(!is.finite(x[!is.na(x)])) || any(x < 0, na.rm=TRUE))
    stop('Nonfinite or negative abundance: do not silently truncate invalid inputs.')
  if (!policy %in% c('stop','confirmed_nondetection','legacy_zero_review')) stop('Unknown missing-value policy.')
  if (anyNA(x) && policy=='stop') stop('Blank abundance cells found. Review missingness audit and confirm whether blanks denote nondetection.')
  if (policy %in% c('confirmed_nondetection','legacy_zero_review')) x[is.na(x)] <- 0
  x
}

load_fig2_cohort <- function(data_root, audit_dir, blank_policy='confirmed_nondetection') {
  dir.create(audit_dir, recursive=TRUE, showWarnings=FALSE)
  suppressWarnings(Sys.setlocale('LC_ALL', 'English_United States.utf8'))
  suffix <- '_split_by_12city_removed_target_columns_leftjoin_merged_subtype_cleaned.xlsx'
  chinese <- c(Hospital='医院', WWTP='污水处理厂', Community='社区', 'Wet market'='农集贸市场')
  metadata_path <- file.path(data_root, 'standardized_2614samples_cleaned_0427.xlsx')
  raw_meta <- readxl::read_excel(metadata_path)
  required <- c('Sample','City','Sample_Date','Sample_Rename')
  stopifnot(all(required %in% names(raw_meta)), !anyDuplicated(raw_meta$Sample))
  meta <- data.frame(Sample=as.character(raw_meta$Sample), City=as.character(raw_meta$City),
    Sample_Date=as.Date(as.character(raw_meta$Sample_Date), '%Y%m%d'),
    PhysicalSite=sub('_[0-9]{8}$', '', raw_meta$Sample_Rename), stringsAsFactors=FALSE)
  stopifnot(!anyNA(meta), nrow(meta)==2614L, length(unique(meta$PhysicalSite))==130L,
    length(unique(meta$City))==12L, !'SH' %in% meta$City)
  blocks <- list(); blank_cells <- list(); provenance <- list(); samples_seen <- character(); ref_subtypes <- NULL
  for (setting in names(chinese)) {
    path <- file.path(data_root, paste0('subtype_by_site_', chinese[[setting]], suffix))
    if (!file.exists(path)) stop('Missing workbook: ',path)
    for (sheet in readxl::excel_sheets(path)) {
      z <- readxl::read_excel(path, sheet=sheet, .name_repair='minimal')
      ids <- names(z)[-1]; subtype <- as.character(z[[1]])
      stopifnot(length(ids)>0, !anyNA(subtype), !anyDuplicated(subtype), !anyDuplicated(ids),
        !any(ids %in% samples_seen), all(ids %in% meta$Sample))
      samples_seen <- c(samples_seen,ids)
      if (is.null(ref_subtypes)) ref_subtypes <- subtype
      stopifnot(setequal(subtype,ref_subtypes), length(subtype)==2646L)
      z <- z[match(ref_subtypes,subtype),]; subtype <- ref_subtypes
      raw <- as.matrix(z[-1]); x <- suppressWarnings(matrix(as.numeric(raw), nrow=length(subtype)))
      if (any(is.na(x) & !is.na(raw) & trimws(raw)!='')) stop('Nonnumeric abundance in ',path,' / ',sheet)
      wh <- which(is.na(x), arr.ind=TRUE)
      if(nrow(wh)) blank_cells[[length(blank_cells)+1L]] <- data.frame(
        Workbook=basename(path), Sheet=sheet, Setting=setting,
        Sample=ids[wh[,2]], Subtype=subtype[wh[,1]], stringsAsFactors=FALSE)
      # Audit first; abort only after all missing cells have been enumerated.
      x <- resolve_abundance_missing(x, if(blank_policy=='stop') 'legacy_zero_review' else blank_policy)
      idx <- match(ids,meta$Sample)
      b <- data.frame(Sample=rep(ids,each=length(subtype)), Subtype=rep(subtype,length(ids)),
        Abundance=as.vector(x), Setting=setting, stringsAsFactors=FALSE)
      for (v in c('City','Sample_Date','PhysicalSite')) b[[v]] <- rep(meta[[v]][idx],each=length(subtype))
      blocks[[length(blocks)+1L]] <- b
      provenance[[length(provenance)+1L]] <- data.frame(Workbook=basename(path),Sheet=sheet,
        Setting=setting,Samples=length(ids),Subtypes=length(subtype),Blank_cells=nrow(wh))
    }
  }
  blanks <- dplyr::bind_rows(blank_cells)
  if(nrow(blanks)) blanks <- dplyr::left_join(blanks,meta,by='Sample')
  utils::write.csv(blanks,file.path(audit_dir,'raw_blank_cells.csv'),row.names=FALSE)
  utils::write.csv(dplyr::bind_rows(provenance),file.path(audit_dir,'input_workbooks.csv'),row.names=FALSE)
  utils::write.csv(data.frame(Blank_policy=blank_policy,Blank_cells=nrow(blanks),
    Affected_samples=length(unique(blanks$Sample)),Status=if(blank_policy=='legacy_zero_review')
      'PROVISIONAL: reproduces old blank-to-zero convention, not a claim of biological absence'
      else if(blank_policy=='confirmed_nondetection')
        'Author confirmed 2026-09-18: blank subtype cells and explicit numeric zeros denote nondetection'
      else blank_policy),
    file.path(audit_dir,'missing_value_policy.csv'),row.names=FALSE)
  if(nrow(blanks) && blank_policy=='stop') stop('Found ',nrow(blanks),' raw blank cells. Audit exported; confirmation needed before zero filling.')
  long <- dplyr::bind_rows(blocks)
  stopifnot(setequal(samples_seen,meta$Sample),nrow(long)==2614*2646)
  setting_map <- unique(long[c('Sample','Setting')]); meta <- dplyr::left_join(meta,setting_map,by='Sample')
  meta$SamplingMonth <- format(meta$Sample_Date,'%Y-%m')
  site_map <- unique(meta[c('PhysicalSite','City','Setting')]); stopifnot(!anyDuplicated(site_map$PhysicalSite))
  meta$Raw_blank_count <- tabulate(match(blanks$Sample,meta$Sample),nbins=nrow(meta))
  samples <- long |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(Richness=sum(Abundance>0),Total_abundance=sum(Abundance),.groups='drop') |>
    dplyr::left_join(meta,by='Sample')
  samples$Data_status <- ifelse(samples$Raw_blank_count>0,
    paste0('Blank cells resolved by policy: ',blank_policy),'Complete raw abundance profile')
  means <- samples |> dplyr::group_by(Setting) |> dplyr::summarise(
    n=dplyr::n(), sites=dplyr::n_distinct(PhysicalSite),mean_richness=mean(Richness),
    mean_abundance=mean(Total_abundance),.groups='drop')
  expected <- c(Hospital=411,WWTP=1349,Community=416,'Wet market'=438)
  stopifnot(all(means$n==expected[means$Setting]))
  utils::write.csv(meta,file.path(audit_dir,'canonical_sample_metadata.csv'),row.names=FALSE)
  utils::write.csv(samples,file.path(audit_dir,'sample_richness_total_abundance.csv'),row.names=FALSE)
  utils::write.csv(means,file.path(audit_dir,'setting_descriptive_means.csv'),row.names=FALSE)
  list(long=long,meta=meta,samples=samples,blanks=blanks,blank_policy=blank_policy)
}

fit_fig2_model <- function(samples, outcome, family, out_dir) {
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if(!family %in% c('negative_binomial','gaussian_original')) stop('Unsupported or unconfirmed model family.')
  if(!outcome %in% c('Richness','Total_abundance')) stop('Unknown outcome.')
  if(outcome=='Total_abundance' && family=='negative_binomial')
    stop('2c is continuous copies/cell: negative-binomial fitting/rounding is not a valid implementation.')
  d <- samples
  d$Setting <- factor(d$Setting,levels=c('Hospital','WWTP','Community','Wet market'))
  d$City <- factor(d$City);d$SamplingMonth <- factor(d$SamplingMonth);d$PhysicalSite <- factor(d$PhysicalSite)
  d$value <- d[[outcome]]
  required <- c('value','Setting','City','SamplingMonth','PhysicalSite')
  cc <- complete.cases(d[required]);utils::write.csv(d[!cc,],file.path(out_dir,'excluded_incomplete_rows.csv'),row.names=FALSE)
  d <- d[cc,];warnings_seen <- character()
  formula <- value ~ Setting + City + SamplingMonth + (1|PhysicalSite)
  message('Fitting ',outcome,': ',family,'; n=',nrow(d),'; sites=',nlevels(d$PhysicalSite))
  fit <- withCallingHandlers({
    if(family=='negative_binomial') {
      validate_nb_response(d$value)
      lme4::glmer.nb(formula,data=d,control=lme4::glmerControl(optimizer='bobyqa',optCtrl=list(maxfun=200000)))
    } else {
      lme4::lmer(formula,data=d,REML=TRUE,control=lme4::lmerControl(optimizer='bobyqa'))
    }
  },warning=function(w){warnings_seen <<- c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')})
  beta <- lme4::fixef(fit);vc <- as.matrix(vcov(fit))
  nd <- d[rep(1,4),];nd$Setting <- factor(levels(d$Setting),levels=levels(d$Setting))
  X <- model.matrix(~Setting+City+SamplingMonth,nd)[,names(beta),drop=FALSE]
  comparisons <- combn(1:4,2,simplify=FALSE)
  contrasts <- dplyr::bind_rows(lapply(comparisons,function(j){
    L <- X[j[1],]-X[j[2],];est <- sum(L*beta);se <- sqrt(drop(L%*%vc%*%L));z <- est/se
    data.frame(group1=levels(d$Setting)[j[1]],group2=levels(d$Setting)[j[2]],
      estimate=est,SE=se,Wald_z=z,CI_low=est-qnorm(.975)*se,CI_high=est+qnorm(.975)*se,
      P=2*pnorm(-abs(z)),Scale=if(family=='negative_binomial') 'log mean ratio' else 'copies/cell difference')
  }))
  contrasts$q_BH <- p.adjust(contrasts$P,'BH')
  contrasts$Comparison <- paste(contrasts$group1,contrasts$group2,sep='-')
  if(family=='negative_binomial') {contrasts$mean_ratio<-exp(contrasts$estimate);contrasts$ratio_low<-exp(contrasts$CI_low);contrasts$ratio_high<-exp(contrasts$CI_high)}
  pmat <- matrix(1,4,4,dimnames=list(levels(d$Setting),levels(d$Setting)))
  for(i in seq_len(nrow(contrasts))) pmat[contrasts$group1[i],contrasts$group2[i]] <- pmat[contrasts$group2[i],contrasts$group1[i]] <- contrasts$q_BH[i]
  # Ordering affects the alphabet only, not any significance decision.
  ord <- levels(d$Setting)[order(drop(X%*%beta),decreasing=TRUE)]
  cld <- multcompView::multcompLetters(pmat[ord,ord],threshold=.05)$Letters
  letters <- data.frame(Site=names(cld),Letters=unname(cld))
  for(i in seq_len(nrow(contrasts))) {
    a <- strsplit(cld[contrasts$group1[i]],'')[[1]];b <- strsplit(cld[contrasts$group2[i]],'')[[1]]
    stopifnot((length(intersect(a,b))==0)==(contrasts$q_BH[i]<.05))
  }
  conv <- paste(unlist(fit@optinfo$conv$lme4$messages),collapse='; ')
  diag <- data.frame(Outcome=outcome,Family=family,Formula=paste(deparse(formula),collapse=' '),
    Transform='None',Link=if(family=='negative_binomial') 'log' else 'identity',
    n_samples=nrow(d),n_sites=length(unique(d$PhysicalSite)),n_cities=length(unique(d$City)),
    n_months=length(unique(d$SamplingMonth)),samples_with_raw_blanks=sum(d$Raw_blank_count>0),
    Inference='Two-sided asymptotic Wald z contrasts; BH across exactly six setting pairs',
    Singular=lme4::isSingular(fit),Convergence_messages=conv,Warnings=paste(unique(warnings_seen),collapse='; '))
  utils::write.csv(diag,file.path(out_dir,'model_diagnostics.csv'),row.names=FALSE)
  utils::write.csv(contrasts,file.path(out_dir,'six_setting_contrasts_BH.csv'),row.names=FALSE)
  utils::write.csv(letters,file.path(out_dir,'setting_letters.csv'),row.names=FALSE)
  utils::write.csv(as.data.frame(lme4::VarCorr(fit)),file.path(out_dir,'variance_components.csv'),row.names=FALSE)
  fit_stats <- data.frame(AIC=AIC(fit),logLik=as.numeric(logLik(fit)),
    NB_size=if(family=='negative_binomial') lme4::getME(fit,'glmer.nb.theta') else NA_real_,
    Pearson_residual_sum_squares=sum(residuals(fit,type='pearson')^2),
    Optimizer=paste(fit@optinfo$optimizer,collapse='; '))
  utils::write.csv(fit_stats,file.path(out_dir,'model_fit_statistics.csv'),row.names=FALSE)
  utils::write.csv(data.frame(Term=names(beta),Coefficient=beta,SE=sqrt(diag(vc))),file.path(out_dir,'fixed_effects.csv'),row.names=FALSE)
  saveRDS(fit,file.path(out_dir,'fitted_model.rds'))
  pdf(file.path(out_dir,'residual_diagnostics.pdf'),width=8,height=4)
  par(mfrow=c(1,2));plot(fitted(fit),residuals(fit,type='pearson'),xlab='Fitted',ylab='Pearson residual');abline(h=0,lty=2)
  qqnorm(residuals(fit,type='pearson'));qqline(residuals(fit,type='pearson'));dev.off()
  if(nzchar(conv)) warning('Model convergence messages: check output diagnostics before using letters.')
  list(fit=fit,contrasts=contrasts,letters=letters,diagnostics=diag)
}


# ----- panels_ade_fix.R -----
# Scoped adapters for the audited Fig. 2 a, d and e sections.
# Source this after constructing .canonical_long and .canonical_meta.
# Units throughout: ARG copies per microbial cell.

audited_city_order <- c(
  "BJ", "CQ", "GY", "GZ", "HF", "HK", "HR", "NJ", "SZ", "XA", "XM", "ZB"
)

.audited_check_canonical <- function() {
  required_long <- c("Sample", "Subtype", "Abundance", "Setting", "City", "Sample_Date", "PhysicalSite")
  required_meta <- setdiff(required_long, c("Subtype", "Abundance"))
  if (!exists(".canonical_long", inherits = TRUE) || !exists(".canonical_meta", inherits = TRUE)) {
    stop("Build .canonical_long and .canonical_meta before using the a/d/e adapters.")
  }
  if (!all(required_long %in% names(.canonical_long)) ||
      !all(required_meta %in% names(.canonical_meta))) {
    stop("Canonical data are missing required sample, setting or measurement columns.")
  }
  if (nrow(.canonical_meta) == 0L || anyNA(.canonical_meta$Sample) ||
      anyDuplicated(.canonical_meta$Sample)) {
    stop("The canonical metadata must contain exactly one row per nonmissing sample ID.")
  }
  if (anyNA(.canonical_meta$City) ||
      !all(as.character(.canonical_meta$City) %in% audited_city_order)) {
    stop("The canonical metadata contain Shanghai, an unknown city, or a missing city.")
  }
  if (anyNA(.canonical_meta$Setting) ||
      !all(as.character(.canonical_meta$Setting) %in% c("Community", "Hospital", "Wet market", "WWTP"))) {
    stop("Use canonical setting labels: Community, Hospital, Wet market, WWTP.")
  }
  if (!is.numeric(.canonical_long$Abundance) || anyNA(.canonical_long$Abundance) ||
      any(!is.finite(.canonical_long$Abundance)) || any(.canonical_long$Abundance < 0)) {
    stop("Canonical abundance must be finite, nonmissing, nonnegative copies/cell.")
  }
  if (!setequal(as.character(.canonical_long$Sample), as.character(.canonical_meta$Sample))) {
    stop("Canonical abundance and metadata sample sets differ.")
  }
  invisible(TRUE)
}

# Analysis 1:
audited_comp_all <- function(type_mapper = canonical_type) {
  .audited_check_canonical()
  if (!is.function(type_mapper)) stop("A subtype-to-type mapping function is required.")
  out <- dplyr::transmute(
    .canonical_long,
    Sample = as.character(Sample),
    Subtype = as.character(Subtype),
    Abundance = Abundance,
    City = as.character(City),
    Site = as.character(Setting),
    Type_raw = sub("__.*$", "", as.character(Subtype))
  )
  out$Type <- type_mapper(out$Type_raw)
  out
}

# Analysis 4:
audited_core_long <- function() {
  .audited_check_canonical()
  # This simple size guard rejects a sparse positive-only long table. The shared
  # cohort loader owns uniqueness of the Sample/Subtype key and zero completion.
  expected_rows <- as.double(nrow(.canonical_meta)) * dplyr::n_distinct(.canonical_long$Subtype)
  if (nrow(.canonical_long) != expected_rows) {
    stop("Core analysis requires one explicit abundance, including zero, per sample and subtype.")
  }
  dplyr::transmute(
    .canonical_long,
    Sample = as.character(Sample),
    Subtype = as.character(Subtype),
    Abundance = Abundance,
    Site = as.character(Setting),
    City = as.character(City),
    Sample_Date = Sample_Date,
    PhysicalSite = as.character(PhysicalSite)
  )
}

# Analysis 3: 
audited_global_input <- function(df) {
  .audited_check_canonical()
  df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
  names(df) <- trimws(names(df))
  choose_column <- function(candidates, description) {
    hits <- intersect(candidates, names(df))
    if (length(hits) != 1L) stop("Expected one ", description, " column in the global benchmark.")
    hits[[1L]]
  }
  country_col <- choose_column(c("COUNTRY", "COUNTRY/REGION", "Country/Region", "Country", "Regions"), "country")
  abundance_col <- choose_column(c("Copies_per_cell", "copies/cell", "copies_per_cell", "copies_cell"), "copies/cell")
  if (!all(c("Sample_ID", "Habitats", "Sub_habitat", "DOI", "Continents") %in% names(df))) {
    stop("The global benchmark is missing required source metadata columns.")
  }
  if (!"Asia_subregion" %in% names(df)) df$Asia_subregion <- NA_character_
  is_study <- function(x) !is.na(x) & trimws(as.character(x)) %in% c("This Study", "本研究")
  old_study <- is_study(df[[country_col]]) | is_study(df$Asia_subregion) | is_study(df$Sample_ID)
  source_rows <- which(!old_study) + 1L
  external <- df[!old_study, , drop = FALSE]
  duplicate_extra <- duplicated(external)
  duplicate_any <- duplicate_extra | duplicated(external, fromLast = TRUE)

  audit_dir <- file.path(.figure_results_dir, "panel_d_global")
  dir.create(audit_dir, recursive = TRUE, showWarnings = FALSE)
  duplicate_audit <- external[duplicate_any, , drop = FALSE]
  duplicate_audit$Source_excel_row <- source_rows[duplicate_any]
  duplicate_audit$Duplicate_extra_record <- duplicate_extra[duplicate_any]
  utils::write.csv(duplicate_audit, file.path(audit_dir, "audit_external_exact_duplicates.csv"), row.names = FALSE)

  # Author approved removing the eight exact duplicate extra records on 2026-09-18.
  duplicate_policy <- tolower(trimws(Sys.getenv("ARG_GLOBAL_DUPLICATES", unset = "drop")))
  if (!duplicate_policy %in% c("exploratory", "keep", "drop")) {
    stop("ARG_GLOBAL_DUPLICATES must be drop (default), keep, or exploratory.")
  }
  if (duplicate_policy == "drop") {
    removed <- duplicate_audit[duplicate_audit$Duplicate_extra_record, , drop = FALSE]
    utils::write.csv(removed, file.path(audit_dir, "removed_external_duplicate_rows.csv"), row.names = FALSE)
    external <- external[!duplicate_extra, , drop = FALSE]
  }

  wwtp_meta <- .canonical_meta[as.character(.canonical_meta$Setting) == "WWTP", , drop = FALSE]
  if (nrow(wwtp_meta) == 0L) stop("No canonical WWTP samples are available for the global benchmark.")
  totals <- .canonical_long |>
    dplyr::filter(Sample %in% wwtp_meta$Sample) |>
    dplyr::group_by(Sample) |>
    dplyr::summarise(Copies_per_cell = sum(Abundance), .groups = "drop")
  wwtp_meta <- dplyr::left_join(wwtp_meta, totals, by = "Sample")
  if (anyNA(wwtp_meta$Copies_per_cell)) stop("A canonical WWTP sample is missing its ARG total.")

  # Seed zero-row source columns so their original names/types are respected.
  replacement <- df[rep(NA_integer_, nrow(wwtp_meta)), , drop = FALSE]
  replacement$Sample_ID <- as.character(wwtp_meta$Sample)
  replacement$Habitats <- "Sewage"
  replacement$Sub_habitat <- "Sewage"
  replacement$DOI <- "Labdata"
  replacement[[country_col]] <- "This Study"
  replacement$Continents <- "Asia"
  replacement$Asia_subregion <- "This Study"
  replacement[[abundance_col]] <- wwtp_meta$Copies_per_cell
  # Preserve real identities and repeated-measure identifiers in the export.
  for (field in c("City", "Setting", "Sample_Date", "PhysicalSite")) {
    if (!field %in% names(external)) external[[field]] <- NA_character_
    external[[field]] <- as.character(external[[field]])
    replacement[[field]] <- as.character(wwtp_meta[[field]])
  }
  external$Benchmark_origin <- "External benchmark"
  replacement$Benchmark_origin <- "Canonical WWTP cohort"
  out <- dplyr::bind_rows(external, replacement)

  note <- paste(
    "Exploratory global comparison: the retained Kruskal-Wallis/Dunn tests do not model",
    "repeated observations at physical sites or external study effects."
  )
  if (duplicate_policy == "exploratory" && any(duplicate_extra)) {
    note <- paste(note, "Exact external duplicate rows are retained by explicit exploratory override; see the duplicate audit.")
    warning(note, call. = FALSE)
  } else {
    if (duplicate_policy == "drop") note <- paste(note,
      paste0(sum(duplicate_extra), " exact duplicate extra external records removed with author approval; first occurrences retained."))
    message(note)
  }
  policy_audit <- data.frame(
    Removed_old_study_rows = sum(old_study),
    Inserted_canonical_WWTP_rows = nrow(replacement),
    External_exact_duplicate_extra_rows = sum(duplicate_extra),
    External_duplicate_rows_removed = if(duplicate_policy == "drop") sum(duplicate_extra) else 0L,
    Duplicate_policy = duplicate_policy,
    External_rows_retained = nrow(external),
    Inference_status = note,
    stringsAsFactors = FALSE
  )
  utils::write.csv(policy_audit, file.path(audit_dir, "audit_global_cohort_and_policy.csv"), row.names = FALSE)
  utils::write.csv(wwtp_meta, file.path(audit_dir, "audit_canonical_WWTP_totals.csv"), row.names = FALSE)
  attr(out, "global_inference_status") <- note
  out
}



# Fig. 2 

run_panel_a <- function() {
  .panel_previous_wd <- getwd()
  on.exit(setwd(.panel_previous_wd), add = TRUE)
  if (!exists(".figure_results_dir", inherits = TRUE)) stop("Set the new .figure_results_dir before running panels.")

  # Panel a: mean ARG type abundance in copies/cell; canonical cohort.

  pkgs <- c(
    "readxl", "openxlsx", "ggplot2", "dplyr", "tidyr", "stringr",
    "purrr", "scales", "tibble", "patchwork"
  )

  library(readxl)
  library(openxlsx)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(purrr)
  library(scales)
  library(tibble)
  library(patchwork)

  # Input is the shared audited canonical cohort.

  output_dir <- file.path(.figure_results_dir, "panel_a")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  site_order <- c("Community", "Hospital", "Wet market", "WWTP")
  city_order <- audited_city_order

  city_order_left <- setdiff(city_order, "HK")
  city_order_right <- city_order

  site_name_map <- c(
    "农集贸市场" = "Wet market",
    "社区"       = "Community",
    "医院"       = "Hospital",
    "污水处理厂" = "WWTP",
    "wet_market" = "Wet market",
    "community"  = "Community",
    "hospital"   = "Hospital",
    "WWTP"       = "WWTP"
  )

  my_cols_cns <- c(
    "Aminoglycoside" = "#4E79A7",
    "Antibacterial fatty acid" = "#2F5D62",
    "Bacitracin" = "#C7A76C",
    "Bicyclomycin" = "#7BAE7F",
    "Bleomycin" = "#D97A5C",
    "Chloramphenicol" = "#A23B72",
    "Defensin" = "#B9CF8A",
    "Edeine" = "#6B8FB3",
    "Factumycin" = "#C05A78",
    "Florfenicol" = "#8FA8C9",
    "Fosfomycin" = "#A9D6E5",
    "Fusidic acid" = "#D8A7B1",
    "MLS" = "#6FB1A0",
    "Multidrug" = "#5B3F8C",
    "Mupirocin" = "#3C6E71",
    "Novobiocin" = "#4C9F70",
    "Other peptide antibiotics" = "#7FB069",
    "Pleuromutilin/Tiamulin" = "#D8C24A",
    "Polymyxin" = "#5FA8A9",
    "Puromycin" = "#7FA7C9",
    "Quinolone" = "#8E6BBE",
    "Rifamycin" = "#9FC490",
    "Streptothricin" = "#3FA7D6",
    "Sulfonamide" = "#B88A2E",
    "Tetracenomycin C" = "#C8A97E",
    "Tetracycline" = "#E6D690",
    "Trimethoprim" = "#D97B7B",
    "Tunicamycin" = "#6C757D",
    "Vancomycin" = "#A3B18A",
    "beta_lactam" = "#556B2F",
    "Others" = "#BDBDBD"
  )

  canonical_type <- function(x) {
    x <- as.character(x)

    dplyr::case_when(
      x %in% c("aminoglycoside", "Aminoglycoside") ~ "Aminoglycoside",
      x %in% c("antibacterial_fatty_acid", "Antibacterial fatty acid") ~ "Antibacterial fatty acid",
      x %in% c("bacitracin", "Bacitracin") ~ "Bacitracin",
      x %in% c("bicyclomycin", "Bicyclomycin") ~ "Bicyclomycin",
      x %in% c("bleomycin", "Bleomycin") ~ "Bleomycin",
      x %in% c("chloramphenicol", "Chloramphenicol") ~ "Chloramphenicol",
      x %in% c("defensin", "Defensin") ~ "Defensin",
      x %in% c("edeine", "Edeine") ~ "Edeine",
      x %in% c("factumycin", "Factumycin") ~ "Factumycin",
      x %in% c("florfenicol", "Florfenicol") ~ "Florfenicol",
      x %in% c("fosfomycin", "Fosfomycin") ~ "Fosfomycin",
      x %in% c("fusidic_acid", "fusidic acid", "Fusidic acid") ~ "Fusidic acid",
      x %in% c("macrolide-lincosamide-streptogramin", "MLS", "mls") ~ "MLS",
      x %in% c("multidrug", "Multidrug") ~ "Multidrug",
      x %in% c("mupirocin", "Mupirocin") ~ "Mupirocin",
      x %in% c("novobiocin", "Novobiocin") ~ "Novobiocin",
      x %in% c("other_peptide_antibiotics", "Other peptide antibiotics") ~ "Other peptide antibiotics",
      x %in% c("pleuromutilin", "tiamulin", "pleuromutilin/tiamulin", "Pleuromutilin/Tiamulin") ~ "Pleuromutilin/Tiamulin",
      x %in% c("polymyxin", "Polymyxin") ~ "Polymyxin",
      x %in% c("puromycin", "Puromycin") ~ "Puromycin",
      x %in% c("quinolone", "Quinolone") ~ "Quinolone",
      x %in% c("rifamycin", "Rifamycin") ~ "Rifamycin",
      x %in% c("streptothricin", "Streptothricin") ~ "Streptothricin",
      x %in% c("sulfonamide", "Sulfonamide") ~ "Sulfonamide",
      x %in% c("tetracenomycin_c", "Tetracenomycin C") ~ "Tetracenomycin C",
      x %in% c("tetracycline", "Tetracycline") ~ "Tetracycline",
      x %in% c("trimethoprim", "Trimethoprim") ~ "Trimethoprim",
      x %in% c("tunicamycin", "Tunicamycin") ~ "Tunicamycin",
      x %in% c("vancomycin", "Vancomycin") ~ "Vancomycin",
      x %in% c("beta_lactam") ~ "beta_lactam",
      x %in% c("others", "Others") ~ "Others",
      TRUE ~ x
    )
  }

  theme_nature <- function(base_size = 13) {
    theme_bw(base_size = base_size) +
      theme(
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 0.8, colour = "black"),
        axis.line = element_line(linewidth = 0.6, colour = "black"),
        axis.text.x = element_text(
          angle = 45, hjust = 1, vjust = 1,
          colour = "black", face = "bold"
        ),
        axis.text.y = element_text(colour = "black", face = "bold"),
        axis.title = element_text(colour = "black", face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5, colour = "black"),
        legend.title = element_text(face = "bold", colour = "black"),
        legend.text = element_text(face = "bold", colour = "black"),
        legend.position = "right",
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", colour = "black")
      )
  }

  save_plot_dual <- function(plot_obj, file_prefix, width = 14, height = 9) {
    ggsave(
      filename = paste0(file_prefix, ".png"),
      plot = plot_obj,
      width = width,
      height = height,
      dpi = 300,
      bg = "white"
    )
    tryCatch({
      ggsave(
        filename = paste0(file_prefix, ".pdf"),
        plot = plot_obj,
        width = width,
        height = height,
        device = cairo_pdf,
        bg = "white"
      )
    }, error = function(e) {
      ggsave(
        filename = paste0(file_prefix, ".pdf"),
        plot = plot_obj,
        width = width,
        height = height,
        device = "pdf",
        bg = "white"
      )
    })
  }

  fix_site_english <- function(x) {
    out <- x
    for (nm in names(site_name_map)) {
      out <- str_replace_all(out, fixed(nm), site_name_map[[nm]])
    }
    out
  }

  clean_city_name <- function(sheet_name) {
    city <- sheet_name
    city <- str_replace(city, "^Sheet 1_", "")
    city <- str_replace(city, "^Sheet1_", "")
    city <- str_replace(city, "^Sheet_", "")
    city <- str_replace(city, "^X_", "")
    city <- str_replace(city, "^ws$", "SH")
    city <- str_replace(city, "^WS$", "SH")
    city
  }

  extract_site_name <- function(file_name) {
    nm <- tools::file_path_sans_ext(basename(file_name))
    nm <- str_replace(nm, "_split_by_city$", "")
    nm <- str_replace(nm, "^subtype_by_site_", "")
    nm <- fix_site_english(nm)
    nm
  }

  extract_date_from_sample <- function(x) {
    str_extract(x, "20\\d{6}")
  }

  order_by_time <- function(x) {
    d <- extract_date_from_sample(x)
    ord <- order(d, x, na.last = TRUE)
    x[ord]
  }

  read_one_sheet <- function(file, sheet = 1) {
    df <- read_excel(file, sheet = sheet, col_names = TRUE)
    df <- as.data.frame(df, check.names = FALSE, stringsAsFactors = FALSE)
    df <- df[, colSums(!is.na(df)) > 0, drop = FALSE]

    if (ncol(df) < 2) return(NULL)

    colnames(df)[1] <- "Subtype"
    df$Subtype <- as.character(df$Subtype)
    df <- df[!is.na(df$Subtype) & df$Subtype != "", , drop = FALSE]
    if (nrow(df) == 0) return(NULL)

    cn <- colnames(df)
    if (length(cn) >= 2) {
      cn[-1] <- str_replace(cn[-1], "^ws", "SH")
      cn[-1] <- str_replace(cn[-1], "^WS", "SH")
      colnames(df) <- cn
    }

    df
  }

  to_sample_matrix <- function(df) {
    mat <- df
    rownames(mat) <- mat$Subtype
    mat$Subtype <- NULL

    mat[] <- lapply(mat, function(x) as.numeric(as.character(x)))
    mat[is.na(mat)] <- 0

    keep_rows <- rowSums(mat > 0, na.rm = TRUE) > 0
    mat <- mat[keep_rows, , drop = FALSE]
    if (nrow(mat) == 0) return(NULL)

    mat_t <- as.data.frame(t(as.matrix(mat)), check.names = FALSE)
    mat_t$Sample <- rownames(mat_t)
    rownames(mat_t) <- NULL

    mat_t$Sample <- str_replace(mat_t$Sample, "^ws", "SH")
    mat_t$Sample <- str_replace(mat_t$Sample, "^WS", "SH")

    ordered_samples <- order_by_time(mat_t$Sample)
    mat_t <- mat_t[match(ordered_samples, mat_t$Sample), , drop = FALSE]

    mat_t <- mat_t[, c(ncol(mat_t), 1:(ncol(mat_t) - 1)), drop = FALSE]
    rownames(mat_t) <- NULL
    mat_t
  }

  comp_all <- audited_comp_all(type_mapper = canonical_type)

  if (nrow(comp_all) == 0) {
    stop("没有读到有效数据，请检查 input_dir 路径和文件格式。")
  }

  comp_all$Site <- factor(comp_all$Site, levels = site_order)
  comp_all$City <- factor(comp_all$City, levels = city_order)

  if (identical(Sys.getenv("ARG_EXPORT_FULL_TABLES", unset = "0"), "1")) {
    write.xlsx(
      comp_all,
      file.path(output_dir, "00_subtype_long_all_values.xlsx"),
      rowNames = FALSE
    )
  }

  top_types <- comp_all %>%
    group_by(Type) %>%
    summarise(total_abd = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_abd)) %>%
    slice_head(n = 10) %>%
    pull(Type)

  write.xlsx(
    data.frame(Type = top_types),
    file.path(output_dir, "01_top_types_selected.xlsx"),
    rowNames = FALSE
  )

  sample_type_abd <- comp_all %>%
    mutate(Type2 = ifelse(Type %in% top_types, Type, "Others")) %>%
    group_by(Sample, City, Site, Type2) %>%
    summarise(type_abd = sum(Abundance, na.rm = TRUE), .groups = "drop")

  comp_c <- sample_type_abd %>%
    group_by(Site, City, Type2) %>%
    summarise(
      mean_abd = mean(type_abd, na.rm = TRUE),
      n_samples = n(),
      .groups = "drop"
    )

  type_levels_c <- names(my_cols_cns)[names(my_cols_cns) %in% unique(as.character(comp_c$Type2))]
  comp_c$Type2 <- factor(comp_c$Type2, levels = type_levels_c)

  unmapped_types <- setdiff(unique(as.character(comp_c$Type2)), names(my_cols_cns))
  if (length(unmapped_types) > 0) {
    warning("这些 Type 没有在 my_cols_cns 里找到颜色：", paste(unmapped_types, collapse = ", "))
  }

  sum_check <- comp_c %>%
    group_by(Site, City) %>%
    summarise(
      total_mean_abd = sum(mean_abd, na.rm = TRUE),
      .groups = "drop"
    )

  write.xlsx(
    sample_type_abd,
    file.path(output_dir, "02_sample_level_type_abundance.xlsx"),
    rowNames = FALSE
  )

  write.xlsx(
    comp_c,
    file.path(output_dir, "03_panel_a_type_mean_abundance.xlsx"),
    rowNames = FALSE
  )

  write.xlsx(
    sum_check,
    file.path(output_dir, "04_panel_a_total_mean_abundance_check.xlsx"),
    rowNames = FALSE
  )

  plot_one_site <- function(dat, site_name, city_levels,
                            show_legend = FALSE,
                            show_y_title = TRUE,
                            show_x_text = TRUE) {
    dat2 <- dat %>%
      filter(Site == site_name, City %in% city_levels) %>%
      mutate(City = factor(as.character(City), levels = city_levels))

    p <- ggplot(dat2, aes(x = City, y = mean_abd, fill = Type2)) +
      geom_col(width = 0.78, colour = "white", linewidth = 0.2) +
      scale_x_discrete(drop = FALSE) +
      scale_fill_manual(
        values = my_cols_cns,
        breaks = type_levels_c,
        drop = FALSE,
        name = "Type"
      ) +
      labs(
        title = site_name,
        x = NULL,
        y = if (show_y_title) "Mean abundance (copies/cell)" else NULL
      ) +
      theme_nature(base_size = 13) +
      theme(
        legend.position = if (show_legend) "right" else "none",
        plot.margin = margin(8, 8, 8, 8)
      )

    if (!show_y_title) {
      p <- p + theme(axis.title.y = element_blank())
    }

    if (!show_x_text) {
      p <- p + theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
    }

    p
  }

  p_comm <- plot_one_site(
    comp_c,
    "Community",
    city_levels = city_order_left,
    show_legend = FALSE,
    show_y_title = TRUE,
    show_x_text = FALSE
  )

  p_wet <- plot_one_site(
    comp_c,
    "Wet market",
    city_levels = city_order_right,
    show_legend = FALSE,
    show_y_title = FALSE,
    show_x_text = FALSE
  )

  p_hosp <- plot_one_site(
    comp_c,
    "Hospital",
    city_levels = city_order_left,
    show_legend = FALSE,
    show_y_title = TRUE,
    show_x_text = TRUE
  )

  p_wwtp <- plot_one_site(
    comp_c,
    "WWTP",
    city_levels = city_order_right,
    show_legend = TRUE,
    show_y_title = FALSE,
    show_x_text = TRUE
  )

  panel_c_final <- (p_comm | p_wet) / (p_hosp | p_wwtp) +
    plot_annotation(
      title = "Panel a. Mean abundance of top types by city"
    )

  save_plot_dual(
    panel_c_final,
    file.path(output_dir, "panel_a_type_mean_abundance"),
    width = 14,
    height = 9
  )

  message("Panel a canonical-cohort mean abundance (copies/cell) finished.")

  invisible(panel_c_final)
}

run_panel_d <- function() {
  .panel_previous_wd <- getwd()
  on.exit(setwd(.panel_previous_wd), add = TRUE)
  if (!exists(".figure_results_dir", inherits = TRUE)) stop("Set the new .figure_results_dir before running panels.")
  # Analysis 3: Panel d: global WWTP benchmark

  # Read the benchmark by explicit path; keep the calling working directory.

  suppressPackageStartupMessages({
    library(readxl)
    library(readr)
    library(dplyr)
    library(tidyr)
    library(purrr)
    library(stringr)
    library(tibble)
    library(ggplot2)
    library(rstatix)
    library(multcompView)
    library(grid)
  })

  sample_file <- .external_path("c_merge", "final_global.xlsx")

  output_dir <- file.path(.figure_results_dir, "panel_d_global")

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  read_input_auto <- function(path) {
    ext <- tolower(tools::file_ext(path))

    if (ext %in% c("xlsx", "xls")) {
      df <- readxl::read_excel(path)
    } else if (ext %in% c("tsv", "txt")) {
      df <- readr::read_tsv(path, show_col_types = FALSE)
    } else if (ext == "csv") {
      df <- readr::read_csv(path, show_col_types = FALSE)
    } else {
      stop("Unsupported input file: ", path)
    }

    as.data.frame(df)
  }

  to_num <- function(x) {
    suppressWarnings(as.numeric(as.character(x)))
  }

  fix_blank <- function(x) {
    x <- as.character(x)
    x <- stringr::str_trim(x)
    x[x %in% c("", "NA", "N/A", "NULL", "null", "NaN", "nan")] <- NA_character_
    x
  }

  std_country <- function(x) {
    y <- fix_blank(x)
    y <- stringr::str_replace_all(y, "\\s+", " ")

    dplyr::case_when(
      y %in% c(
        "Taiwan", "Taiwan, China", "China: Taiwan", "Chinese Taipei",
        "Taiwan Province of China", "Taiwan, Province of China"
      ) ~ "China Taiwan",

      y %in% c(
        "Hong Kong SAR", "Hong Kong, China", "China: Hong Kong",
        "Hong Kong Special Administrative Region"
      ) ~ "Hong Kong",

      y %in% c("United States", "United States of America") ~ "USA",
      y == "Brasil" ~ "Brazil",
      y == "Republic of Macedonia" ~ "North Macedonia",
      y == "Slovak Republic" ~ "Slovakia",
      y == "Republic of Moldova" ~ "Moldova",
      y == "Czech Republic" ~ "Czechia",

      y %in% c(
        "Cote d'Ivoire", "Côte d'Ivoire", "C?te d’Ivoire",
        "C?te d  Ivoire", "C么te d鈥橧voire", "COTE DIVOIRE"
      ) ~ "Côte d’Ivoire",

      TRUE ~ y
    )
  }

  std_continent <- function(continent, country) {
    continent <- fix_blank(continent)
    country <- std_country(country)

    dplyr::case_when(
      country == "This Study" ~ "This Study",
      country %in% c("China", "China Taiwan", "Hong Kong") ~ "Asia",
      country %in% c("Australia", "New Zealand") ~ "Oceania",
      continent %in% c("Australia", "Australasia") ~ "Oceania",
      TRUE ~ continent
    )
  }

  std_asia_subregion_one <- function(subregion, country, continent) {
    subregion <- fix_blank(subregion)
    country <- std_country(country)
    continent <- std_continent(continent, country)

    if (!is.na(country) && country == "This Study") {
      return("This Study")
    }
    if (!is.na(subregion) && subregion == "This Study") {
      return("This Study")
    }

    subregion <- dplyr::case_when(
      subregion %in% c("中国内地", "Mainland China") ~ "Mainland China",
      subregion %in% c("中国台湾", "China Taiwan") ~ "China Taiwan",
      subregion %in% c("中国Hong Kong SAR", "China Hong Kong SAR") ~ "China Hong Kong SAR",
      subregion %in% c("东亚", "East Asia") ~ "East Asia",
      subregion %in% c("东南亚", "Southeast Asia") ~ "Southeast Asia",
      subregion %in% c("南亚", "South Asia") ~ "South Asia",
      subregion %in% c("中亚", "Central Asia") ~ "Central Asia",
      subregion %in% c("西亚", "West Asia") ~ "West Asia",
      subregion %in% c("This Study", "本研究") ~ "This Study",
      TRUE ~ subregion
    )

    if (is.na(continent) || continent != "Asia") {
      return(NA_character_)
    }

    if (!is.na(subregion) && subregion != "") {
      return(subregion)
    }

    east_asia <- c("Japan", "South Korea", "North Korea", "Mongolia")

    southeast_asia <- c(
      "Brunei", "Cambodia", "Indonesia", "Laos", "Malaysia", "Myanmar",
      "Philippines", "Singapore", "Thailand", "Timor-Leste", "Vietnam"
    )

    south_asia <- c(
      "Afghanistan", "Bangladesh", "Bhutan", "India", "Maldives",
      "Nepal", "Pakistan", "Sri Lanka"
    )

    central_asia <- c(
      "Kazakhstan", "Kyrgyzstan", "Tajikistan", "Turkmenistan", "Uzbekistan"
    )

    west_asia <- c(
      "Armenia", "Azerbaijan", "Bahrain", "Cyprus", "Georgia", "Iran",
      "Iraq", "Israel", "Jordan", "Kuwait", "Lebanon", "Oman",
      "Palestine", "Qatar", "Saudi Arabia", "Syria", "Turkey",
      "United Arab Emirates", "Yemen"
    )

    dplyr::case_when(
      country == "China" ~ "Mainland China",
      country == "China Taiwan" ~ "China Taiwan",
      country == "Hong Kong" ~ "China Hong Kong SAR",
      country %in% east_asia ~ "East Asia",
      country %in% southeast_asia ~ "Southeast Asia",
      country %in% south_asia ~ "South Asia",
      country %in% central_asia ~ "Central Asia",
      country %in% west_asia ~ "West Asia",
      TRUE ~ NA_character_
    )
  }

  df <- audited_global_input(read_input_auto(sample_file))
  names(df) <- stringr::str_trim(names(df))

  cat("\nInput file used:\n")
  print(normalizePath(sample_file))

  cat("\nRaw rows after reading:\n")
  print(nrow(df))

  cat("\nRaw Continents counts directly after reading:\n")
  if ("Continents" %in% names(df)) {
    print(df %>% count(Continents, sort = TRUE))
  }

  if ("COUNTRY/REGION" %in% names(df)) {
    df <- df %>% rename(COUNTRY = `COUNTRY/REGION`)
  }
  if ("Country/Region" %in% names(df)) {
    df <- df %>% rename(COUNTRY = `Country/Region`)
  }
  if ("Country" %in% names(df)) {
    df <- df %>% rename(COUNTRY = Country)
  }
  if ("Regions" %in% names(df) && !"COUNTRY" %in% names(df)) {
    df <- df %>% rename(COUNTRY = Regions)
  }

  if ("copies/cell" %in% names(df)) {
    df <- df %>% rename(Copies_per_cell = `copies/cell`)
  }
  if ("copies_per_cell" %in% names(df)) {
    df <- df %>% rename(Copies_per_cell = copies_per_cell)
  }
  if ("copies_cell" %in% names(df)) {
    df <- df %>% rename(Copies_per_cell = copies_cell)
  }

  if (!"Asia_subregion" %in% names(df)) {
    df$Asia_subregion <- NA_character_
  }

  required_df <- c(
    "Sample_ID", "Habitats", "Sub_habitat", "DOI",
    "COUNTRY", "Continents", "Asia_subregion", "Copies_per_cell"
  )

  missing_df <- setdiff(required_df, names(df))
  if (length(missing_df) > 0) {
    stop("Input file missing required columns: ", paste(missing_df, collapse = ", "))
  }

  asia_as_num <- to_num(df$Asia_subregion)
  copy_as_num <- to_num(df$Copies_per_cell)

  shift_idx <- !is.na(asia_as_num) & is.na(copy_as_num)

  df$Copies_per_cell[shift_idx] <- asia_as_num[shift_idx]
  df$Asia_subregion[shift_idx] <- NA_character_

  df <- df %>%
    mutate(
      Sample_ID = fix_blank(Sample_ID),
      Habitats = fix_blank(Habitats),
      Sub_habitat = fix_blank(Sub_habitat),
      DOI = fix_blank(DOI),
      COUNTRY = std_country(COUNTRY),
      Continents = std_continent(Continents, COUNTRY),
      Copies_per_cell = to_num(Copies_per_cell),

      Asia_subregion = purrr::pmap_chr(
        list(Asia_subregion, COUNTRY, Continents),
        function(Asia_subregion, COUNTRY, Continents) {
          std_asia_subregion_one(Asia_subregion, COUNTRY, Continents)
        }
      ),

      This_Study = dplyr::coalesce(COUNTRY == "This Study", FALSE) |
        dplyr::coalesce(Asia_subregion == "This Study", FALSE) |
        dplyr::coalesce(Sample_ID == "This Study", FALSE),

      COUNTRY_PLOT = dplyr::if_else(This_Study, "This Study", COUNTRY),
      Continents = dplyr::if_else(This_Study, "This Study", Continents),
      Asia_subregion = dplyr::if_else(This_Study, "This Study", Asia_subregion)
    ) %>%
    filter(!is.na(Copies_per_cell), !is.na(COUNTRY_PLOT), !is.na(Continents))

  region_order_fixed <- c(
    "Africa",
    "Oceania",
    "South America",
    "Asia",
    "Europe",
    "North America",
    "This Study",
    "Hong Kong"
  )

  df <- df %>%
    mutate(
      Region = factor(Continents, levels = region_order_fixed),
      Violin_Group = case_when(
        This_Study ~ "This Study",
        COUNTRY_PLOT == "Hong Kong" ~ "Hong Kong",
        TRUE ~ as.character(Region)
      ),
      Violin_Group = factor(Violin_Group, levels = region_order_fixed)
    )

  asia_order <- c(
    "East Asia",
    "Southeast Asia",
    "South Asia",
    "Central Asia",
    "West Asia",
    "This Study",
    "China Hong Kong SAR",
    "China Taiwan",
    "Mainland China"
  )

  df <- df %>%
    mutate(
      Asia_subregion_plot = case_when(
        This_Study ~ "This Study",
        Continents == "Asia" ~ Asia_subregion,
        TRUE ~ NA_character_
      ),
      Asia_subregion_plot = factor(Asia_subregion_plot, levels = asia_order)
    )

  df_all <- df

  df_asia <- df_all %>%
    filter(!is.na(Asia_subregion_plot))

  cat("\nCheck This Study in final df_all:\n")
  print(
    df_all %>%
      filter(This_Study) %>%
      summarise(
        n = n(),
        mean = mean(Copies_per_cell, na.rm = TRUE),
        label = sprintf("%.2f", mean(Copies_per_cell, na.rm = TRUE))
      )
  )

  cat("\nCheck DOI == Labdata but NOT This Study:\n")
  print(
    df_all %>%
      filter(DOI == "Labdata", !This_Study) %>%
      count(COUNTRY, sort = TRUE)
  )

  cat("\nRegion group counts before Plot 1 from df_all:\n")
  region_check <- df_all %>%
    count(Violin_Group, sort = FALSE)

  print(region_check)

  cat("\nAsia subregion counts before Plot 2 from df_asia:\n")
  print(
    df_asia %>%
      count(Asia_subregion_plot, sort = FALSE)
  )

  must_have_region <- c(
    "Africa", "Oceania", "South America", "Asia",
    "Europe", "North America", "This Study", "Hong Kong"
  )

  region_present <- df_all %>%
    filter(!is.na(Violin_Group)) %>%
    pull(Violin_Group) %>%
    as.character() %>%
    unique()

  missing_region <- setdiff(must_have_region, region_present)

  if (length(missing_region) > 0) {
    stop(
      "df_all is missing groups: ",
      paste(missing_region, collapse = ", "),
      "\nYou are reading the wrong input file or globally filtering df before plotting."
    )
  }

  write_csv(
    df_all,
    file.path(output_dir, "merged_sample_with_meta_cleaned.csv")
  )

  unclassified_asia <- df_all %>%
    filter(Continents == "Asia", !This_Study) %>%
    filter(is.na(Asia_subregion) | Asia_subregion == "") %>%
    distinct(COUNTRY, Continents, Asia_subregion)

  write_csv(
    unclassified_asia,
    file.path(output_dir, "unmatched_asia_subregion_after_fix.csv")
  )

  region_palette <- c(
    "Africa" = "#6CBF84",
    "Oceania" = "#B38BCB",
    "South America" = "#F6A6B2",
    "Asia" = "#F3A64A",
    "Europe" = "#999999",
    "North America" = "#6EC5C1",
    "This Study" = "#E15759",
    "Hong Kong" = "#E9C95B"
  )

  asia_palette <- c(
    "East Asia" = "#999999",
    "Southeast Asia" = "#6EC5C1",
    "South Asia" = "#6CBF84",
    "Central Asia" = "#B38BCB",
    "West Asia" = "#C08F5A",
    "This Study" = "#E15759",
    "China Hong Kong SAR" = "#E9C95B",
    "China Taiwan" = "#E9C95B",
    "Mainland China" = "#E9C95B"
  )

  country_palette <- c(
    "Africa" = "#6CBF84",
    "Oceania" = "#B38BCB",
    "South America" = "#F6A6B2",
    "Asia" = "#F3A64A",
    "Europe" = "#999999",
    "North America" = "#6EC5C1",
    "This Study" = "#E15759",
    "China mainland/HK/Taiwan" = "#E9C95B"
  )

  ggname <- function(prefix, grob) {
    grob$name <- grid::grobName(grob, prefix)
    grob
  }

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  GeomFlatViolin <- ggproto(
    "GeomFlatViolin", Geom,

    setup_data = function(data, params) {
      data$width <- data$width %||%
        params$width %||%
        (resolution(data$x, FALSE) * 0.9)

      data %>%
        dplyr::group_by(group) %>%
        dplyr::mutate(
          ymin = min(y),
          ymax = max(y),
          xmin = x,
          xmax = x + width / 2
        )
    },

    draw_group = function(data, panel_scales, coord) {
      data <- transform(
        data,
        xminv = x,
        xmaxv = x + violinwidth * (xmax - x)
      )

      newdata <- rbind(
        data[order(data$y), c(names(data))] %>%
          transform(x = data$xmaxv[order(data$y)]),
        data[order(-data$y), c(names(data))] %>%
          transform(x = data$xminv[order(-data$y)])
      )

      newdata <- rbind(newdata, newdata[1, ])

      ggname(
        "geom_flat_violin",
        GeomPolygon$draw_panel(newdata, panel_scales, coord)
      )
    },

    draw_key = draw_key_polygon,

    default_aes = aes(
      weight = 1,
      colour = NA,
      fill = "grey80",
      linewidth = 0.3,
      alpha = NA,
      linetype = "solid"
    ),

    required_aes = c("x", "y")
  )

  geom_flat_violin <- function(
      mapping = NULL,
      data = NULL,
      stat = "ydensity",
      position = "identity",
      ...,
      trim = TRUE,
      scale = "area",
      na.rm = FALSE,
      show.legend = NA,
      inherit.aes = TRUE
  ) {
    layer(
      data = data,
      mapping = mapping,
      stat = stat,
      geom = GeomFlatViolin,
      position = position,
      show.legend = show.legend,
      inherit.aes = inherit.aes,
      params = list(
        trim = trim,
        scale = scale,
        na.rm = na.rm,
        ...
      )
    )
  }

  make_sig_letters <- function(pairwise_res, group_order, alpha = 0.05) {

    pw <- pairwise_res %>%
      filter(!is.na(p.adj)) %>%
      mutate(
        group1 = as.character(group1),
        group2 = as.character(group2)
      )

    if (nrow(pw) == 0) {
      return(
        tibble(
          Group = factor(character(), levels = group_order),
          Letters = character(),
          x_id = numeric()
        )
      )
    }

    pvec <- pw$p.adj
    names(pvec) <- paste(pw$group1, pw$group2, sep = "-")

    cld <- multcompView::multcompLetters(pvec, threshold = alpha)

    tibble(
      Group = names(cld$Letters),
      Letters = cld$Letters
    ) %>%
      mutate(
        Group = factor(Group, levels = group_order),
        x_id = as.numeric(Group)
      ) %>%
      filter(!is.na(x_id)) %>%
      arrange(x_id)
  }

  make_index_plot <- function(
      data,
      group_var,
      group_order,
      palette,
      file_prefix,
      width = 6,
      height = 3.5,
      y_cap = 10,
      x_angle = 25
  ) {

    dd <- data %>%
      filter(!is.na(.data[[group_var]]), !is.na(Copies_per_cell)) %>%
      mutate(
        Group = as.character(.data[[group_var]]),
        Group = factor(Group, levels = group_order)
      ) %>%
      filter(!is.na(Group)) %>%
      mutate(x_id = as.numeric(Group))

    n_df <- dd %>%
      count(Group, name = "n") %>%
      mutate(Group = factor(Group, levels = group_order))

    kw_res <- tryCatch(
      rstatix::kruskal_test(dd, Copies_per_cell ~ Group),
      error = function(e) tibble()
    )

    pairwise_res <- tryCatch(
      rstatix::dunn_test(dd, Copies_per_cell ~ Group, p.adjust.method = "BH") %>%
        arrange(p.adj),
      error = function(e) tibble(group1 = character(), group2 = character(), p.adj = numeric())
    )

    write_csv(kw_res, file.path(output_dir, paste0(file_prefix, "_kruskal.csv")))
    write_csv(pairwise_res, file.path(output_dir, paste0(file_prefix, "_dunn.csv")))
    write_csv(n_df, file.path(output_dir, paste0(file_prefix, "_group_n.csv")))

    present_order <- dd %>%
      distinct(Group) %>%
      arrange(as.numeric(Group)) %>%
      pull(Group) %>%
      as.character()

    letter_df <- make_sig_letters(pairwise_res, group_order = present_order, alpha = 0.05)

    y_max <- max(dd$Copies_per_cell, na.rm = TRUE)
    y_min <- min(dd$Copies_per_cell, na.rm = TRUE)
    y_range <- y_max - y_min
    if (is.na(y_range) || y_range == 0) y_range <- 1

    mean_df <- dd %>%
      group_by(Group, x_id) %>%
      summarise(
        mean_value = mean(Copies_per_cell, na.rm = TRUE),
        n = n(),
        .groups = "drop"
      ) %>%
      mutate(
        mean_lab = sprintf("%.2f", mean_value)
      )

    mean_y_fixed   <- min(y_max + 0.10 * y_range, y_cap * 0.925)
    letter_y_fixed <- min(y_max + 0.22 * y_range, y_cap * 0.965)

    mean_df <- mean_df %>%
      mutate(label_y = mean_y_fixed)

    letter_df <- letter_df %>%
      mutate(letter_y = letter_y_fixed)

    write_csv(mean_df, file.path(output_dir, paste0(file_prefix, "_descriptive_means.csv")))
    write_csv(letter_df, file.path(output_dir, paste0(file_prefix, "_letters.csv")))

    p <- ggplot() +
      geom_flat_violin(
        data = dd,
        aes(x = x_id, y = Copies_per_cell, fill = Group),
        width = 0.62,
        trim = FALSE,
        alpha = 0.82,
        colour = NA
      ) +
      geom_jitter(
        data = dd %>% mutate(x_point = x_id - 0.17),
        aes(x = x_point, y = Copies_per_cell, colour = Group),
        width = 0.06,
        height = 0,
        alpha = 0.18,
        size = 0.9,
        shape = 16,
        stroke = 0
      ) +
      geom_boxplot(
        data = dd,
        aes(x = x_id, y = Copies_per_cell, group = Group, colour = Group),
        width = 0.082,
        outlier.shape = NA,
        linewidth = 0.6,
        fill = "white"
      ) +
      stat_summary(
        data = dd,
        aes(x = x_id, y = Copies_per_cell, group = Group, colour = Group),
        fun = median,
        geom = "crossbar",
        width = 0.06,
        linewidth = 0.6
      ) +
      geom_text(
        data = mean_df,
        aes(x = x_id, y = label_y, label = mean_lab, colour = Group),
        fontface = "plain",
        size = 3.05,
        vjust = 0
      ) +
      geom_text(
        data = letter_df,
        aes(x = x_id, y = letter_y, label = Letters, colour = Group),
        fontface = "plain",
        size = 4.2,
        vjust = 0
      ) +
      scale_fill_manual(values = palette, drop = TRUE) +
      scale_colour_manual(values = palette, drop = TRUE) +
      scale_x_continuous(
        breaks = seq_along(group_order),
        labels = group_order,
        expand = expansion(mult = c(0.10, 0.12))
      ) +
      labs(
        x = NULL,
        y = "Copies per cell"
      ) +
      scale_y_continuous(
        breaks = seq(0, y_cap, by = 2),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      coord_cartesian(
        ylim = c(0, y_cap),
        clip = "off"
      ) +
      theme_classic(base_size = 12) +
      theme(
        text = element_text(face = "plain"),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.65
        ),
        axis.line = element_blank(),
        axis.title.x = element_text(
          size = 14,
          face = "plain",
          colour = "black",
          margin = margin(t = 12)
        ),
        axis.title.y = element_text(
          size = 14,
          face = "plain",
          colour = "black",
          margin = margin(r = 10)
        ),
        axis.text.x = element_text(
          angle = x_angle,
          hjust = 1,
          vjust = 1,
          size = 13,
          face = "plain",
          colour = "black",
          margin = margin(t = 10)
        ),
        axis.text.y = element_text(
          size = 13,
          face = "plain",
          colour = "black",
          margin = margin(r = 4)
        ),
        axis.ticks = element_line(
          colour = "black",
          linewidth = 0.42
        ),
        axis.ticks.length = unit(0.13, "cm"),
        panel.grid = element_blank(),
        legend.position = "none",
        plot.margin = margin(16, 20, 22, 16)
      )

    ggsave(file.path(output_dir, paste0(file_prefix, ".pdf")), p, width = width, height = height)
    ggsave(file.path(output_dir, paste0(file_prefix, ".png")), p, width = width, height = height, dpi = 600)

    return(p)
  }

  make_country_barplot_by_region <- function(
      data,
      file_prefix,
      min_n = 1,
      width = 16,
      height = 4.8,
      y_cap = 12
  ) {

    dd <- data %>%
      filter(!is.na(COUNTRY_PLOT), !is.na(Copies_per_cell), !is.na(Continents)) %>%
      mutate(
        Plot_Region = case_when(
          COUNTRY_PLOT %in% c("China", "Hong Kong", "China Taiwan") ~ "China mainland/HK/Taiwan",
          COUNTRY_PLOT == "This Study" ~ "This Study",
          TRUE ~ as.character(Continents)
        )
      ) %>%
      group_by(COUNTRY_PLOT, Plot_Region) %>%
      summarise(
        mean_value = mean(Copies_per_cell, na.rm = TRUE),
        sd_value   = ifelse(dplyr::n() > 1, sd(Copies_per_cell, na.rm = TRUE), 0),
        var_value  = ifelse(dplyr::n() > 1, var(Copies_per_cell, na.rm = TRUE), 0),
        sem_value  = ifelse(dplyr::n() > 1, sd(Copies_per_cell, na.rm = TRUE) / sqrt(dplyr::n()), 0),
        n = dplyr::n(),
        .groups = "drop"
      ) %>%
      filter(n >= min_n) %>%
      mutate(
        Country_display = case_when(
          COUNTRY_PLOT == "China" ~ "Mainland China",
          COUNTRY_PLOT == "Hong Kong" ~ "China Hong Kong SAR",
          COUNTRY_PLOT == "China Taiwan" ~ "China Taiwan",
          TRUE ~ COUNTRY_PLOT
        ),
        Country_label = paste0(Country_display, " (n=", n, ")"),
        Plot_Region = factor(Plot_Region, levels = names(country_palette)),
        plot_value = pmin(mean_value, y_cap),

        error_ymin = pmax(mean_value - sem_value, 0),
        error_ymax = mean_value + sem_value,

        error_ymin = pmin(error_ymin, y_cap),
        error_ymax = pmin(error_ymax, y_cap),

        mean_lab = sprintf("%.2f", mean_value),
        label_y = ifelse(
          error_ymax >= y_cap * 0.95,
          y_cap * 0.965,
          error_ymax + y_cap * 0.025
        )
      ) %>%
      arrange(desc(mean_value), Country_label) %>%
      mutate(
        Country_label = factor(Country_label, levels = Country_label)
      )

    write_csv(dd, file.path(output_dir, paste0(file_prefix, "_summary.csv")))

    p <- ggplot(dd, aes(x = Country_label, y = plot_value, fill = Plot_Region)) +
      geom_col(width = 0.78, colour = NA) +
      geom_errorbar(
        aes(ymin = error_ymin, ymax = error_ymax),
        width = 0.18,
        linewidth = 0.45,
        colour = "black"
      ) +
      geom_text(
        aes(y = label_y, label = mean_lab),
        size = 2.6,
        angle = 90,
        vjust = 0.5,
        hjust = 0,
        colour = "black"
      ) +
      scale_fill_manual(values = country_palette, drop = TRUE, name = "Region") +
      labs(
        x = NULL,
        y = "Mean copies per cell ± SEM"
      ) +
      scale_y_continuous(
        breaks = seq(0, y_cap, by = 2),
        expand = expansion(mult = c(0.02, 0.02))
      ) +
      coord_cartesian(ylim = c(0, y_cap), clip = "on") +
      theme_classic(base_size = 12) +
      theme(
        text = element_text(face = "plain"),
        panel.border = element_rect(
          colour = "black",
          fill = NA,
          linewidth = 0.65
        ),
        axis.line = element_blank(),
        axis.title.y = element_text(
          size = 14,
          face = "plain",
          colour = "black",
          margin = margin(r = 10)
        ),
        axis.text.x = element_text(
          angle = 60,
          hjust = 1,
          vjust = 1,
          size = 8.5,
          face = "plain",
          colour = "black"
        ),
        axis.text.y = element_text(
          size = 12,
          face = "plain",
          colour = "black"
        ),
        axis.ticks = element_line(
          colour = "black",
          linewidth = 0.42
        ),
        panel.grid = element_blank(),
        legend.position = "right",
        plot.margin = margin(16, 18, 22, 16)
      )

    ggsave(file.path(output_dir, paste0(file_prefix, ".pdf")), p, width = width, height = height)
    ggsave(file.path(output_dir, paste0(file_prefix, ".png")), p, width = width, height = height, dpi = 600)

    return(p)
  }

  p_region <- make_index_plot(
    data = df_all,
    group_var = "Violin_Group",
    group_order = region_order_fixed,
    palette = region_palette,
    file_prefix = "copies_per_cell_by_region_with_thisstudy_hk",
    width = 4,
    height = 3.5,
    y_cap = 10,
    x_angle = 25
  )

  p_asia <- make_index_plot(
    data = df_asia,
    group_var = "Asia_subregion_plot",
    group_order = asia_order,
    palette = asia_palette,
    file_prefix = "copies_per_cell_by_asia_subregion_with_thisstudy",
    width = 8.5,
    height = 3.5,
    y_cap = 10,
    x_angle = 30
  )

  p_country_all <- make_country_barplot_by_region(
    data = df_all,
    file_prefix = "copies_per_cell_by_country_barplot_region_fill_all",
    min_n = 1,
    width = 22,
    height = 4.8,
    y_cap = 12
  )

  p_country_n10 <- make_country_barplot_by_region(
    data = df_all,
    file_prefix = "copies_per_cell_by_country_barplot_region_fill_n_ge_10",
    min_n = 10,
    width = 16,
    height = 4.8,
    y_cap = 12
  )

  if (interactive()) print(p_region)
  if (interactive()) print(p_asia)
  if (interactive()) print(p_country_all)
  if (interactive()) print(p_country_n10)

  cat("\nDone.\n")
  cat("Outputs saved in: ", output_dir, "\n")
  cat("Generated:\n")
  cat("- copies_per_cell_by_region_with_thisstudy_hk.pdf/png\n")
  cat("- copies_per_cell_by_asia_subregion_with_thisstudy.pdf/png\n")
  cat("- copies_per_cell_by_country_barplot_region_fill_all.pdf/png, mean ± SEM\n")
  cat("- copies_per_cell_by_country_barplot_region_fill_n_ge_10.pdf/png, mean ± SEM\n")
  cat("- cleaned metadata: merged_sample_with_meta_cleaned.csv\n")
  cat("- unmatched Asia check: unmatched_asia_subregion_after_fix.csv\n")

  invisible(list(region = p_region, asia = p_asia, country_all = p_country_all, country_n10 = p_country_n10))
}

run_panel_e <- function() {
  .panel_previous_wd <- getwd()
  on.exit(setwd(.panel_previous_wd), add = TRUE)
  if (!exists(".figure_results_dir", inherits = TRUE)) stop("Set the new .figure_results_dir before running panels.")
  # Analysis 4: Panel e: core prevalence-abundance

  # Input is the shared audited canonical cohort.



  pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "scales")

  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)

  data_long <- audited_core_long()

  if (nrow(data_long) == 0) {
    stop("❌ data_long 为空，请检查文件内容或 site 识别规则")
  }

  cat("长表行数:", nrow(data_long), "\n")
  cat("Site 分布:\n")
  print(table(data_long$Site, useNA = "ifany"))

  prev_df <- data_long %>%
    group_by(Site, Subtype) %>%
    summarise(
      prevalence = mean(Abundance > 0, na.rm = TRUE),
      mean_abd = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    )

  prev_thr <- 0.7
  mean_thr <- 1e-5

  prev_df <- prev_df %>%
    mutate(
      Core = prevalence >= prev_thr & mean_abd >= mean_thr
    )

  top_n <- 25

  top_features <- prev_df %>%
    group_by(Subtype) %>%
    summarise(max_prev = max(prevalence, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(max_prev), Subtype) %>%
    slice_head(n = top_n) %>%
    pull(Subtype)

  plot_df <- prev_df %>%
    filter(Subtype %in% top_features) %>%
    mutate(
      Subtype = factor(Subtype, levels = rev(top_features)),
      Site = factor(Site, levels = c("Community", "Hospital", "Wet market", "WWTP"))
    )

  plot_df <- plot_df %>%
    mutate(
      Subtype_short = str_replace(as.character(Subtype), ".*__", ""),
      Subtype_short = factor(Subtype_short, levels = rev(str_replace(top_features, ".*__", "")))
    )

  use_short_name <- TRUE

  y_var <- if (use_short_name) "Subtype_short" else "Subtype"

  p <- ggplot(plot_df, aes(x = Site, y = .data[[y_var]], fill = prevalence)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_point(
      data = plot_df %>% filter(Core),
      aes(x = Site, y = .data[[y_var]]),
      inherit.aes = FALSE,
      shape = 8,
      size = 2.5,
      color = "black"
    ) +
    scale_fill_gradientn(
      colours = c("#F7FBFF", "#6BAED6", "#2171B5", "#08306B"),
      limits = c(0, 1),
      name = "Prevalence"
    ) +
    labs(
      title = "Core subtype prevalence across site types",
      x = NULL,
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      axis.text.x = element_text(angle = 30, hjust = 1),
      plot.title = element_text(face = "bold", hjust = 0.5)
    )

  if (interactive()) print(p)

  time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(.figure_results_dir, "panel_e_core_prevalence")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  ggsave(
    filename = file.path(out_dir, "core_subtype.png"),
    plot = p,
    width = 9,
    height = 7,
    dpi = 300,
    bg = "white"
  )

  tryCatch(
    {
      ggsave(
        filename = file.path(out_dir, "core_subtype.pdf"),
        plot = p,
        width = 9,
        height = 7,
        device = cairo_pdf,
        bg = "white"
      )
    },
    error = function(e) {
      message("PDF 导出失败，但 PNG 已保存。错误信息：", e$message)
    }
  )

  write.csv(
    plot_df,
    file.path(out_dir, "core_subtype_plot_data.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    prev_df,
    file.path(out_dir, "core_subtype_prevalence_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("✅ 完成：输出目录 -> ", out_dir, "\n")

  core_prevalence_plot <- p

  # Input is the shared audited canonical cohort.



  pkgs <- c("readxl", "dplyr", "tidyr", "stringr", "ggplot2", "scales")

  library(readxl)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(scales)

  has_ggrepel <- requireNamespace("ggrepel", quietly = TRUE)

  data_long <- audited_core_long()

  if (nrow(data_long) == 0) {
    stop("❌ data_long 为空，请检查文件内容")
  }

  cat("长表行数:", nrow(data_long), "\n")
  cat("Site 分布:\n")
  print(table(data_long$Site, useNA = "ifany"))

  prev_thr <- 0.70
  mean_thr <- 1e-5

  core_df <- data_long %>%
    group_by(Site, Subtype) %>%
    summarise(
      prevalence = mean(Abundance > 0, na.rm = TRUE),
      mean_abd   = mean(Abundance, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Core = prevalence >= prev_thr & mean_abd >= mean_thr,
      Site = factor(Site, levels = c("Community", "Hospital", "Wet market", "WWTP"))
    )

  label_n <- 10

  plot_df <- core_df %>%
    filter(is.finite(mean_abd), mean_abd > 0)

  label_df <- plot_df %>%
    group_by(Site) %>%
    arrange(desc(Core), desc(prevalence), desc(mean_abd), .by_group = TRUE) %>%
    slice_head(n = label_n) %>%
    ungroup() %>%
    mutate(
      label_text = str_replace(Subtype, ".*__", "")
    )

  p <- ggplot(plot_df, aes(x = prevalence, y = mean_abd)) +
    geom_vline(
      xintercept = prev_thr,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    geom_hline(
      yintercept = mean_thr,
      linetype = "dashed",
      linewidth = 0.5,
      colour = "grey40"
    ) +
    geom_point(
      aes(fill = Core, size = Core),
      shape = 21,
      colour = "black",
      stroke = 0.35,
      alpha = 0.9
    ) +
    facet_wrap(~ Site, ncol = 2, scales = "fixed") +
    scale_fill_manual(
      values = c("FALSE" = "#D9D9D9", "TRUE" = "#D95F5F"),
      name = "Core"
    ) +
    scale_size_manual(
      values = c("FALSE" = 2.2, "TRUE" = 3.3),
      name = "Core"
    ) +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, 0.2),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_log10(
      labels = label_scientific(digits = 1)
    ) +
    labs(
      title = "Core subtype distribution across site types",
      x = "Prevalence",
      y = "Mean abundance (copies/cell)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0),
      strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.7),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major = element_line(colour = "#D9D9D9", linewidth = 0.4),
      panel.grid.minor = element_blank(),
      axis.text = element_text(colour = "black", size = 11),
      axis.title = element_text(face = "bold", size = 13),
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11),
      panel.border = element_rect(colour = "black", linewidth = 0.8)
    )

  if (has_ggrepel) {
    p <- p +
      ggrepel::geom_text_repel(
        data = label_df,
        aes(label = label_text),
        size = 3.0,
        box.padding = 0.25,
        point.padding = 0.2,
        segment.size = 0.35,
        max.overlaps = 100,
        seed = 123,
        show.legend = FALSE
      )
  } else {
    p <- p +
      geom_text(
        data = label_df,
        aes(label = label_text),
        size = 3.0,
        hjust = 0,
        vjust = -0.2,
        check_overlap = TRUE,
        show.legend = FALSE
      )
  }

  if (interactive()) print(p)

  time_tag <- format(Sys.time(), "%Y%m%d_%H%M%S")
  out_dir <- file.path(.figure_results_dir, "panel_e_core_scatter")

  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }

  ggsave(
    filename = file.path(out_dir, "core_subtype_distribution_by_site.png"),
    plot = p,
    width = 11,
    height = 8.5,
    dpi = 300,
    bg = "white"
  )

  tryCatch(
    {
      ggsave(
        filename = file.path(out_dir, "core_subtype_distribution_by_site.pdf"),
        plot = p,
        width = 11,
        height = 8.5,
        device = cairo_pdf,
        bg = "white"
      )
    },
    error = function(e) {
      message("PDF 导出失败，但 PNG 已保存。错误信息：", e$message)
    }
  )

  write.csv(
    core_df,
    file.path(out_dir, "core_subtype_distribution_table.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  write.csv(
    label_df,
    file.path(out_dir, "core_subtype_labels_used.csv"),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("✅ 完成：输出目录 -> ", out_dir, "\n")

  invisible(list(prevalence = core_prevalence_plot, scatter = p))
}


# ----- panels_bc_original_style.R -----
# Generated from supplied Fig. 2 code; half-violin style retained.
run_panels_bc <- function(cohort) {
# Analysis 2: Panels b-c: richness and total abundance

pkg_needed <- c(
  "readxl", "dplyr", "tidyr", "stringr", "ggplot2",
  "grid", "FSA", "multcompView", "purrr"
)

library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(grid)
library(multcompView)
library(purrr)

ggname <- function(prefix, grob) {
  grob$name <- grid::grobName(grob, prefix)
  grob
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

GeomFlatViolin <- ggproto(
  "GeomFlatViolin", Geom,

  setup_data = function(data, params) {
    data$width <- data$width %||%
      params$width %||%
      (resolution(data$x, FALSE) * 0.9)

    data %>%
      dplyr::group_by(group) %>%
      dplyr::mutate(
        ymin = min(y),
        ymax = max(y),
        xmin = x,
        xmax = x + width / 2
      )
  },

  draw_group = function(data, panel_scales, coord) {
    data <- transform(
      data,
      xminv = x,
      xmaxv = x + violinwidth * (xmax - x)
    )

    newdata <- rbind(
      transform(data[order(data$y), ], x = data$xmaxv[order(data$y)]),
      transform(data[order(-data$y), ], x = data$xminv[order(-data$y)])
    )

    newdata <- rbind(newdata, newdata[1, ])

    ggname(
      "geom_flat_violin",
      GeomPolygon$draw_panel(newdata, panel_scales, coord)
    )
  },

  draw_key = draw_key_polygon,

  default_aes = aes(
    weight = 1,
    colour = NA,
    fill = "grey80",
    linewidth = 0.3,
    alpha = NA,
    linetype = "solid"
  ),

  required_aes = c("x", "y")
)

geom_flat_violin <- function(
    mapping = NULL,
    data = NULL,
    stat = "ydensity",
    position = "identity",
    ...,
    trim = TRUE,
    scale = "area",
    na.rm = FALSE,
    show.legend = NA,
    inherit.aes = TRUE
) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      na.rm = na.rm,
      ...
    )
  )
}

in_dir <- .external_path("c_merge", "2614")
out_dir <- file.path(.figure_results_dir, "panels_b_c")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_list <- c(
  "Hospital"   = file.path(in_dir, "subtype_by_site_医院_split_by_12city_removed_target_columns_leftjoin_merged_subtype_cleaned.xlsx"),
  "WWTP"       = file.path(in_dir, "subtype_by_site_污水处理厂_split_by_12city_removed_target_columns_leftjoin_merged_subtype_cleaned.xlsx"),
  "Community"  = file.path(in_dir, "subtype_by_site_社区_split_by_12city_removed_target_columns_leftjoin_merged_subtype_cleaned.xlsx"),
  "Wet market" = file.path(in_dir, "subtype_by_site_农集贸市场_split_by_12city_removed_target_columns_leftjoin_merged_subtype_cleaned.xlsx")
)

site_order <- c("Hospital", "WWTP", "Community", "Wet market")

site_palette <- c(
  "Community"  = "#6EC5C1",
  "Hospital"   = "#E06C75",
  "Wet market" = "#E9C95B",
  "WWTP"       = "#4E79A7"
)

plot_narrow_ratio <- 0.60
geom_narrow_ratio <- 0.60

sample_richness <- cohort$samples %>% mutate(Site=factor(Setting, levels=site_order),
  value=Richness, Metric="Richness", x_id=as.numeric(Site))
sample_copies <- cohort$samples %>% mutate(Site=factor(Setting, levels=site_order),
  value=Total_abundance, Metric="Copies per cell", x_id=as.numeric(Site))
write.csv(sample_richness,file.path(out_dir,"02_sample_richness.csv"),row.names=FALSE)
write.csv(sample_copies,file.path(out_dir,"03_sample_copies_per_cell.csv"),row.names=FALSE)
get_letter_df <- function(dat,value_col="value",site_order) {
  model_dir <- if(unique(dat$Metric)=="Richness") "panel_b_NB" else "panel_c_LMM"
  base <- file.path(.figure_results_dir,model_dir)
  pairwise <- read.csv(file.path(base,"six_setting_contrasts_BH.csv"),check.names=FALSE)
  letters <- read.csv(file.path(base,"setting_letters.csv"))
  stopifnot(nrow(pairwise)==6L, nrow(letters)==length(site_order),
    all(c("P", "q_BH") %in% names(pairwise)), all(nzchar(letters$Letters)),
    setequal(letters$Site,site_order))
  letters$Site <- factor(letters$Site,levels=site_order)
  letters$x_id <- match(letters$Site,site_order)
  list(model_pairwise=pairwise,letter_df=letters)
}
make_half_violin_plot <- function(
    dat,
    ylab_text,
    out_prefix,
    mean_digits = 1,
    plot_width = 6 * plot_narrow_ratio,
    plot_height = 4.8,
    y_cap = Inf,
    y_breaks = waiver()
) {

  stat_list <- get_letter_df(
    dat = dat,
    value_col = "value",
    site_order = site_order
  )

  model_pairwise  <- stat_list$model_pairwise
  letter_df <- stat_list$letter_df

  write.csv(
    model_pairwise,
    file.path(out_dir, paste0(out_prefix, "_mixed_model_BH_results.csv")),
    row.names = FALSE
  )

  y_max <- max(dat$value, na.rm = TRUE)
  y_min <- min(dat$value, na.rm = TRUE)
  y_range <- y_max - y_min

  if (!is.finite(y_range) || y_range == 0) {
    y_range <- 1
  }

  mean_df <- dat %>%
    group_by(Site, x_id) %>%
    summarise(
      mean_value = mean(value, na.rm = TRUE),
      q3 = quantile(value, 0.75, na.rm = TRUE),
      top = max(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean_lab = sprintf(paste0("%.", mean_digits, "f"), mean_value)
    )

  is_capped <- is.finite(y_cap)

  if (is_capped) {

    y_upper <- y_cap

    mean_df <- mean_df %>%
      mutate(
        mean_y = y_upper * 0.88
      )

    letter_y_fixed <- y_upper * 0.94
    clip_mode <- "on"

  } else {

    mean_df <- mean_df %>%
      mutate(
        mean_y = q3 + 0.30 * y_range
      )

    letter_y_fixed <- max(mean_df$mean_y, na.rm = TRUE) + 0.10 * y_range

    y_upper <- letter_y_fixed + 0.08 * y_range

    if (!is.finite(y_upper)) {
      y_upper <- y_max * 1.15
    }

    if (y_upper <= y_max) {
      y_upper <- y_max + 1
    }

    clip_mode <- "off"
  }

  letter_df <- letter_df %>%
    mutate(
      letter_y = letter_y_fixed
    )

  write.csv(
    letter_df,
    file.path(out_dir, paste0(out_prefix, "_group_letters.csv")),
    row.names = FALSE
  )

  p <- ggplot() +

    geom_flat_violin(
      data = dat,
      aes(x = x_id, y = value, fill = Site),
      width = 0.54 * geom_narrow_ratio,
      trim = TRUE,
      alpha = 0.82,
      colour = NA
    ) +

    geom_jitter(
      data = dat %>%
        mutate(
          x_point = x_id - 0.15 * geom_narrow_ratio
        ),
      aes(x = x_point, y = value, colour = Site),
      width = 0.04 * geom_narrow_ratio,
      height = 0,
      alpha = 0.15,
      size = 0.9,
      shape = 16,
      stroke = 0
    ) +

    geom_boxplot(
      data = dat,
      aes(x = x_id, y = value, group = Site, colour = Site),
      width = 0.082 * geom_narrow_ratio,
      outlier.shape = NA,
      linewidth = 0.6,
      fill = "white"
    ) +

    stat_summary(
      data = dat,
      aes(x = x_id, y = value, group = Site, colour = Site),
      fun = median,
      geom = "crossbar",
      width = 0.06 * geom_narrow_ratio,
      linewidth = 0.6
    ) +

    geom_text(
      data = mean_df,
      aes(x = x_id, y = mean_y, label = mean_lab, colour = Site),
      size = 4.2,
      vjust = 0,
      fontface = "plain"
    ) +

    geom_text(
      data = letter_df,
      aes(x = x_id, y = letter_y, label = Letters, colour = Site),
      size = 6.2,
      vjust = 0,
      fontface = "plain"
    ) +

    scale_fill_manual(values = site_palette, drop = FALSE) +
    scale_colour_manual(values = site_palette, drop = FALSE) +

    scale_x_continuous(
      breaks = seq_along(site_order),
      labels = site_order,
      expand = expansion(mult = c(0.04, 0.06))
    ) +

    scale_y_continuous(
      breaks = y_breaks,
      expand = expansion(mult = c(0, 0))
    ) +

    coord_cartesian(
      ylim = c(0, y_upper),
      clip = clip_mode
    ) +

    labs(
      x = NULL,
      y = ylab_text,
      title = NULL
    ) +

    theme_classic(base_size = 12) +

    theme(
      text = element_text(face = "plain"),

      plot.title = element_blank(),

      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.65
      ),

      axis.line = element_blank(),

      axis.text.x = element_text(
        angle = 25,
        hjust = 1,
        vjust = 1,
        size = 13,
        colour = "black",
        face = "plain",
        margin = margin(t = 10)
      ),

      axis.text.y = element_text(
        size = 13,
        colour = "black",
        face = "plain",
        margin = margin(r = 4)
      ),

      axis.title.x = element_text(
        size = 14,
        colour = "black",
        face = "plain",
        margin = margin(t = 12)
      ),

      axis.title.y = element_text(
        size = 14,
        colour = "black",
        face = "plain",
        margin = margin(r = 10)
      ),

      axis.ticks = element_line(
        colour = "black",
        linewidth = 0.42
      ),

      axis.ticks.length = unit(0.13, "cm"),

      panel.grid = element_blank(),
      legend.position = "none",

      plot.margin = margin(16, 20, 22, 16)
    )

  if(interactive()) print(p)

  ggsave(
    file.path(out_dir, paste0(out_prefix, ".pdf")),
    p,
    width = plot_width,
    height = plot_height
  )

  ggsave(
    file.path(out_dir, paste0(out_prefix, ".png")),
    p,
    width = plot_width,
    height = plot_height,
    dpi = 600
  )

  return(p)
}

p_richness <- make_half_violin_plot(
  dat = sample_richness,
  ylab_text = "Richness",
  out_prefix = "half_violin_site_richness_refstyle_narrow60",
  mean_digits = 1,
  plot_width = 7.5 * plot_narrow_ratio,
  plot_height = 4.8,
  y_cap = Inf
)

p_copies <- make_half_violin_plot(
  dat = sample_copies,
  ylab_text = "Copies per cell",
  out_prefix = "half_violin_site_copiescell_refstyle_narrow60_fullrange",
  mean_digits = 2,
  plot_width = 7.5 * plot_narrow_ratio,
  plot_height = 4.8,
  y_cap = Inf,
  y_breaks = waiver()
)

cat("\nDone.\n")
cat("Outputs saved in:\n", out_dir, "\n")

cat("\nGenerated figures:\n")
cat("- half_violin_site_richness_refstyle_narrow60.pdf/png\n")
cat("- half_violin_site_copiescell_refstyle_narrow60_fullrange.pdf/png\n")

cat("\nGenerated tables:\n")
if (identical(Sys.getenv("ARG_EXPORT_FULL_TABLES", unset = "0"), "1")) {
  cat("- 01_long_table_all_sites.csv\n")
}
cat("- 02_sample_richness.csv\n")
cat("- 03_sample_copies_per_cell.csv\n")
cat("- *_mixed_model_BH_results.csv\n")
cat("- *_group_letters.csv\n")


}


# ----- panel_f_fix.R -----
# Audited Fig. 2f: one canonical core membership table for every plot/export.
# Source this file, then call audited_panel_f() after the shared cohort loader.
# No files are reread between calculation and plotting; no working directory changes.
# Layout change: all core subtypes and all types are shown. The original manually
# reconstructed central dendrogram is omitted because its leaf identities were wrong.

.f_settings <- c("Hospital", "WWTP", "Community", "Wet market")
.f_setting_colours <- c(Hospital = "#E06C75", WWTP = "#4E79A7",
                        Community = "#6EC5C1", `Wet market` = "#E9C95B")
.f_type_colours <- c(
  Aminoglycoside = "#4E79A7", `Antibacterial fatty acid` = "#7FB8B5",
  Bacitracin = "#D6B97A", Bicyclomycin = "#8BCB88", Bleomycin = "#E39B76",
  Chloramphenicol = "#B5658D", Defensin = "#B7D989", Edeine = "#8FB6E8",
  Factumycin = "#D17C98", Florfenicol = "#AFC7E8", Fosfomycin = "#B7E0E5",
  `Fusidic acid` = "#E6BAC4", MLS = "#74C2B3", Multidrug = "#7E6BB5",
  Mupirocin = "#5F9EA0", Novobiocin = "#72B77E",
  `Other peptide antibiotics` = "#9FD58A", `Pleuromutilin/Tiamulin` = "#E4CF62",
  Polymyxin = "#75C0C1", Puromycin = "#9DBDE0", Quinolone = "#9A84D6",
  Rifamycin = "#A9C98A", Streptothricin = "#6DBDE3", Sulfonamide = "#C8A24B",
  `Tetracenomycin C` = "#D8BB92", Tetracycline = "#E8D98F",
  Trimethoprim = "#E39A9A", Tunicamycin = "#9AA3A8", Vancomycin = "#AFC9A0",
  `beta-lactam` = "#7DA34D", Others = "#D0D0D0"
)

.f_type_label <- function(x) {
  key <- tolower(gsub("_", " ", trimws(x), fixed = TRUE))
  lut <- stats::setNames(names(.f_type_colours), tolower(names(.f_type_colours)))
  lut <- c(lut, "beta lactam" = "beta-lactam",
           "macrolide-lincosamide-streptogramin" = "MLS",
           "pleuromutilin tiamulin" = "Pleuromutilin/Tiamulin")
  out <- unname(lut[key])
  out[is.na(out)] <- x[is.na(out)]
  out
}

.f_require_packages <- function() {
  packages <- c("dplyr", "tidyr", "ggplot2", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages for Fig. 2f: ", paste(missing, collapse = ", "))
}

audited_panel_f_core <- function(long = .canonical_long, meta = .canonical_meta,
                                 prevalence_threshold = 0.70, abundance_threshold = 1e-5) {
  .f_require_packages()
  needed <- c("Sample", "Subtype", "Abundance", "Setting")
  if (!all(needed %in% names(long)) || !all(c("Sample", "Setting") %in% names(meta))) {
    stop("Fig. 2f requires canonical Sample, Subtype, Abundance and Setting columns.")
  }
  if (anyNA(meta$Sample) || anyDuplicated(meta$Sample) || anyNA(meta$Setting) ||
      !setequal(as.character(meta$Setting), .f_settings)) {
    stop("Fig. 2f requires one metadata row per sample and exactly four canonical settings.")
  }
  if (!is.numeric(long$Abundance) || anyNA(long$Abundance) ||
      any(!is.finite(long$Abundance)) || any(long$Abundance < 0) ||
      anyNA(long$Sample) || anyNA(long$Subtype) || anyNA(long$Setting)) {
    stop("Canonical measurements must be finite, nonnegative and have complete identifiers.")
  }
  long_sample_map <- dplyr::distinct(long, Sample, Setting)
  meta_sample_map <- dplyr::distinct(meta, Sample, Setting)
  if (nrow(long_sample_map) != nrow(meta) ||
      nrow(dplyr::anti_join(long_sample_map, meta_sample_map, by = c("Sample", "Setting"))) ||
      nrow(dplyr::anti_join(meta_sample_map, long_sample_map, by = c("Sample", "Setting")))) {
    stop("Canonical abundance and metadata sample/setting assignments differ.")
  }
  n_subtypes <- dplyr::n_distinct(long$Subtype)
  if (nrow(long) != as.double(nrow(meta)) * n_subtypes) {
    stop("Core calculation requires a dense table, including zeros, for every sample and subtype.")
  }
  denominators <- dplyr::count(meta, Setting, name = "expected_samples")
  out <- long |>
    dplyr::group_by(Setting, Subtype) |>
    dplyr::summarise(
      n_rows = dplyr::n(), n_samples = dplyr::n_distinct(Sample),
      n_detected = sum(Abundance > 0), prevalence = mean(Abundance > 0),
      mean_abd = mean(Abundance), median_abd = stats::median(Abundance),
      max_abd = max(Abundance), .groups = "drop"
    ) |>
    dplyr::left_join(denominators, by = "Setting")
  if (nrow(out) != length(.f_settings) * n_subtypes ||
      anyNA(out$expected_samples) || any(out$n_rows != out$n_samples) ||
      any(out$n_samples != out$expected_samples)) {
    stop("Duplicate sample/subtype keys or missing setting/subtype measurements in canonical data.")
  }
  out |>
    dplyr::mutate(Core = prevalence >= prevalence_threshold & mean_abd >= abundance_threshold) |>
    dplyr::select(-n_rows, -expected_samples) |>
    dplyr::arrange(match(Setting, .f_settings), dplyr::desc(Core), Subtype)
}

audited_panel_f_tables <- function(core_statistics) {
  .f_require_packages()
  full <- as.data.frame(core_statistics, stringsAsFactors = FALSE)
  if (!"Setting" %in% names(full) && "Site" %in% names(full)) names(full)[names(full) == "Site"] <- "Setting"
  required <- c("Setting", "Subtype", "n_samples", "n_detected", "prevalence", "mean_abd", "Core")
  if (!all(required %in% names(full)) || anyDuplicated(full[c("Setting", "Subtype")]) ||
      !setequal(as.character(full$Setting), .f_settings)) stop("Invalid core-statistics table.")
  if (!is.logical(full$Core)) stop("Core must be a logical flag derived from the stated thresholds.")
  core <- dplyr::filter(full, Core)
  membership <- core |>
    dplyr::select(Subtype, Setting) |>
    dplyr::mutate(value = 1L) |>
    tidyr::pivot_wider(names_from = Setting, values_from = value, values_fill = 0L)
  for (st in .f_settings) if (!st %in% names(membership)) membership[[st]] <- 0L
  conditional <- core |>
    dplyr::group_by(Subtype) |>
    dplyr::summarise(
      core_setting_mean_abundance = mean(mean_abd),
      core_setting_mean_prevalence = mean(prevalence), .groups = "drop"
    )
  pooled <- full |>
    dplyr::group_by(Subtype) |>
    dplyr::summarise(
      all_sample_mean_abundance = stats::weighted.mean(mean_abd, n_samples),
      equal_setting_mean_abundance = mean(mean_abd), .groups = "drop"
    )
  membership <- membership |>
    dplyr::left_join(conditional, by = "Subtype") |>
    dplyr::left_join(pooled, by = "Subtype") |>
    dplyr::mutate(
      Subtype_full = Subtype,
      Type_raw = sub("__.*$", "", Subtype), Type = .f_type_label(Type_raw),
      Subtype_label = sub("^.*__", "", Subtype),
      n_core_settings = rowSums(dplyr::pick(dplyr::all_of(.f_settings))),
      core_log10_mean_abundance = log10(core_setting_mean_abundance)
    ) |>
    dplyr::arrange(Type, Subtype)
  membership$intersection <- apply(as.matrix(membership[.f_settings]), 1L, function(z) {
    paste(.f_settings[as.logical(z)], collapse = " | ")
  })
  setting_counts <- data.frame(
    Setting = .f_settings,
    n_core_subtypes = vapply(.f_settings, function(st) sum(membership[[st]]), numeric(1)),
    row.names = NULL
  )
  intersections <- membership |>
    dplyr::group_by(dplyr::across(dplyr::all_of(.f_settings)), intersection, n_core_settings) |>
    dplyr::summarise(n_subtypes = dplyr::n(), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(n_core_settings), dplyr::desc(n_subtypes), intersection) |>
    dplyr::mutate(intersection_id = dplyr::row_number())
  types <- membership |>
    dplyr::group_by(Type) |>
    dplyr::summarise(
      n_subtypes = dplyr::n(),
      mean_log10_core_setting_mean_abundance = mean(core_log10_mean_abundance),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(n_subtypes), Type)
  totals <- data.frame(
    Metric = c("Union", "Shared_all_four", "Hospital_core_only",
               "Nonhospital_core_not_hospital_core", "Shared_three_nonhospital_not_hospital"),
    n_subtypes = c(nrow(membership), sum(membership$n_core_settings == 4L),
                   sum(membership$Hospital == 1L & membership$n_core_settings == 1L),
                   sum(membership$Hospital == 0L),
                   sum(membership$Hospital == 0L & membership$n_core_settings == 3L))
  )
  stopifnot(sum(intersections$n_subtypes) == nrow(membership))
  list(full = full, core = core, membership = membership, setting_counts = setting_counts,
       intersections = intersections, types = types, totals = totals)
}

audited_panel_f_plots <- function(tables) {
  .f_require_packages()
  m <- tables$membership
  combos <- tables$intersections
  y_map <- stats::setNames(1:4, .f_settings) # Identical map for bars, dots and connector lines.
  set_sizes <- tables$setting_counts
  set_sizes$y <- unname(y_map[set_sizes$Setting])
  matrix <- combos |>
    tidyr::pivot_longer(dplyr::all_of(.f_settings), names_to = "Setting", values_to = "present") |>
    dplyr::mutate(y = unname(y_map[Setting]))
  connectors <- matrix |>
    dplyr::filter(present == 1L) |>
    dplyr::group_by(intersection_id) |>
    dplyr::summarise(ymin = min(y), ymax = max(y), .groups = "drop")
  limits_x <- c(0.4, nrow(combos) + 0.6)
  top <- ggplot2::ggplot(combos, ggplot2::aes(intersection_id, n_subtypes)) +
    ggplot2::geom_col(fill = "#566573", width = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = n_subtypes), vjust = -0.35, size = 3.6) +
    ggplot2::scale_x_continuous(limits = limits_x, expand = c(0, 0)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.13))) +
    ggplot2::labs(x = NULL, y = "Exclusive intersection size") +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  left <- ggplot2::ggplot(set_sizes, ggplot2::aes(n_core_subtypes, y, fill = Setting)) +
    ggplot2::geom_col(width = 0.62, orientation = "y") +
    ggplot2::geom_text(ggplot2::aes(label = n_core_subtypes), hjust = -0.1, size = 3.7) +
    ggplot2::scale_fill_manual(values = .f_setting_colours, guide = "none") +
    ggplot2::scale_x_reverse(expand = ggplot2::expansion(mult = c(0.18, 0.02))) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5), breaks = 1:4, expand = c(0, 0)) +
    ggplot2::labs(x = "Core set size", y = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  dots <- ggplot2::ggplot(matrix, ggplot2::aes(intersection_id, y)) +
    ggplot2::geom_segment(data = connectors,
                          ggplot2::aes(xend = intersection_id, y = ymin, yend = ymax), linewidth = 0.7) +
    ggplot2::geom_point(ggplot2::aes(colour = factor(present)), size = 3.5) +
    ggplot2::scale_colour_manual(values = c("0" = "#DDDDDD", "1" = "#222222"), guide = "none") +
    ggplot2::scale_x_continuous(limits = limits_x, expand = c(0, 0), breaks = NULL) +
    ggplot2::scale_y_continuous(limits = c(0.5, 4.5), breaks = 1:4,
                              labels = .f_settings, expand = c(0, 0)) +
    ggplot2::labs(x = NULL, y = NULL) + ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(axis.text.y = ggplot2::element_text(colour = "#222222"))
  upset <- patchwork::wrap_plots(patchwork::plot_spacer(), top, left, dots,
                                 ncol = 2, widths = c(1.3, 4), heights = c(2.3, 1)) +
    patchwork::plot_annotation(title = paste0("Core ARG sharing across settings (", nrow(m), " subtypes)"))

  # Circular track: the same full membership table, with no top-N type filter.
  ord <- m |>
    dplyr::mutate(type_order = match(Type, tables$types$Type)) |>
    dplyr::arrange(type_order, dplyr::desc(core_log10_mean_abundance), Subtype) |>
    dplyr::mutate(index = dplyr::row_number() + 3L * (type_order - 1L))
  n_arc <- max(ord$index)
  n_total <- ceiling(n_arc * 4 / 3) # Leave one quadrant for the key.
  sectors <- ord |>
    dplyr::group_by(Type) |>
    dplyr::summarise(xmin = min(index) - 0.5, xmax = max(index) + 0.5, .groups = "drop") |>
    dplyr::left_join(tables$types, by = "Type")
  sectors$colour <- unname(.f_type_colours[sectors$Type])
  if (anyNA(sectors$colour)) stop("Missing circle colour for type: ", paste(sectors$Type[is.na(sectors$colour)], collapse = ", "))
  # Preserve the original abundance statistic. Heights have a fixed log10 baseline
  # of -5 (the core mean-abundance threshold), not an unexplained min/max rescale.
  sectors$bar_top <- 1.15 + 0.052 * (sectors$mean_log10_core_setting_mean_abundance + 5)
  ring_radii <- c(`Wet market` = 0.79, Community = 0.88, WWTP = 0.97, Hospital = 1.06)
  tracks <- ord |>
    tidyr::pivot_longer(dplyr::all_of(.f_settings), names_to = "Setting", values_to = "present") |>
    dplyr::mutate(radius = unname(ring_radii[Setting]),
                  colour = ifelse(present == 1L, unname(.f_setting_colours[Setting]), "#E9E9E9"))
  radius_max <- max(1.48, max(sectors$bar_top) + 0.04)
  circle <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = sectors,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0.63, ymax = 0.72, fill = colour)) +
    ggplot2::geom_tile(data = tracks,
                       ggplot2::aes(index, radius, fill = colour), width = 0.98, height = 0.065) +
    ggplot2::geom_rect(data = sectors,
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 1.15, ymax = bar_top, fill = colour)) +
    ggplot2::scale_fill_identity() +
    ggplot2::coord_polar(theta = "x", start = -pi / 2, clip = "off") +
    ggplot2::scale_x_continuous(limits = c(0, n_total), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, radius_max), expand = c(0, 0)) +
    ggplot2::theme_void() + ggplot2::labs(title = paste0("All ", nrow(m), " core ARG subtypes")) +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 16, face = "bold", hjust = 0.5))
  # Draw a readable key separately; very small type sectors still retain all members.
  legend <- tables$types
  legend$y <- nrow(legend):1
  legend$colour <- unname(.f_type_colours[legend$Type])
  legend$label <- paste0(legend$Type, " (", legend$n_subtypes, ")")
  key <- ggplot2::ggplot(legend, ggplot2::aes(0, y)) +
    ggplot2::geom_point(ggplot2::aes(colour = colour), shape = 15, size = 3.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), hjust = 0, nudge_x = 0.07, size = 3.1) +
    ggplot2::scale_colour_identity() + ggplot2::xlim(-0.03, 1.6) +
    ggplot2::theme_void() + ggplot2::labs(title = "ARG type (subtypes)") +
    ggplot2::theme(plot.title = ggplot2::element_text(size = 11, face = "bold"))
  circle_final <- patchwork::wrap_plots(circle, key, widths = c(3.6, 1.65)) +
    patchwork::plot_annotation(
      subtitle = "Core-membership tracks, inner to outer: Wet market, Community, WWTP, Hospital. Grey = not core.",
      caption = paste(
        "Outer bars: type mean of log10 subtype abundance; each subtype abundance is the equal-weight mean",
        "of setting means only where that subtype is core. Radial baseline = -5; 0.052 radial units per log10 unit.",
        "This conditional statistic is not the pooled mean across all samples. No dendrogram is inferred.", sep = "\n")
    )
  list(upset = upset, circle = circle_final, circle_subtype_order = ord,
       circle_type_sectors = sectors, membership_tracks = tracks)
}

audited_panel_f_export <- function(tables, output_dir, make_plots = TRUE,
                                   prevalence_threshold = 0.70, abundance_threshold = 1e-5) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  filenames <- c(full = "02_core_subtype_by_setting_FULL.csv",
                 core = "03_core_subtype_by_setting_ONLY_CORE.csv",
                 membership = "core_membership_ALL_subtypes.csv",
                 setting_counts = "core_counts_by_setting.csv",
                 intersections = "exclusive_core_intersections.csv",
                 types = "core_type_summary_ALL_types.csv", totals = "core_union_summary.csv")
  for (nm in names(filenames)) utils::write.csv(tables[[nm]], file.path(output_dir, filenames[[nm]]),
                                                row.names = FALSE, fileEncoding = "UTF-8")
  definitions <- data.frame(
    Quantity = c("Core membership", "Detection", "Core-setting mean abundance", "All-sample mean abundance",
                 "Outer bar", "Circle scope", "Tree", "UpSet intersection"),
    Definition = c(
      paste0("Within-setting prevalence >= ", prevalence_threshold, " AND sample mean abundance >= ", abundance_threshold, " copies/cell; zeros included."),
      "Abundance > 0; denominator = all retained samples within the setting.",
      "Equal-weight average of sample means only over settings where this subtype passes both core criteria.",
      "Sample-count-weighted mean across all four settings; exported for clarity, not used for the original conditional circle statistic.",
      "Arithmetic mean across subtype log10(core-setting mean abundance) within an ARG type. Baseline -5; 0.052 radial units per log10 unit.",
      "All core subtypes and all ARG types, with four setting-membership tracks. No top-15 filter.",
      "Original incorrectly labelled reconstructed central dendrogram omitted; no substitute hierarchy inferred.",
      "Exact/exclusive membership combination; not an inclusive overlap. Not-core does not mean undetected."
    )
  )
  utils::write.csv(definitions, file.path(output_dir, "plot_and_statistic_definitions.csv"), row.names = FALSE)
  plots <- NULL
  if (make_plots) {
    plots <- audited_panel_f_plots(tables)
    for (name in c("upset", "circle")) {
      stem <- file.path(output_dir, paste0("Fig2f_", name, "_ALL_core_subtypes"))
      ggplot2::ggsave(paste0(stem, ".pdf"), plots[[name]], width = 13, height = 8.5, bg = "white")
      ggplot2::ggsave(paste0(stem, ".png"), plots[[name]], width = 13, height = 8.5, dpi = 300, bg = "white")
    }
    utils::write.csv(plots$circle_subtype_order, file.path(output_dir, "circle_subtype_plot_order.csv"), row.names = FALSE)
    utils::write.csv(plots$circle_type_sectors, file.path(output_dir, "circle_type_sectors.csv"), row.names = FALSE)
  }
  invisible(list(tables = tables, plots = plots, output_dir = output_dir))
}

audited_panel_f <- function(long = .canonical_long, meta = .canonical_meta,
                            output_dir = file.path(.figure_results_dir, "panel_f_audited"),
                            make_plots = TRUE) {
  stats <- audited_panel_f_core(long, meta)
  tables <- audited_panel_f_tables(stats)
  audited_panel_f_export(tables, output_dir, make_plots = make_plots)
}


# ----- panel_f_original_style_fixed.R -----

run_panel_f_original <- function(long = .canonical_long, meta = .canonical_meta,
                                 output_dir = file.path(.figure_results_dir, "panel_f_original_style"),
                                 make_plots = FALSE) {
  required <- c("dplyr", "tidyr", "stringr", "purrr", "tibble", "ggplot2",
                "scales", "forcats", "patchwork", "openxlsx")
  missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Original Fig. 2f needs: ", paste(missing, collapse = ", "))
  if (!exists("audited_panel_f_core", mode = "function")) stop("Source panel_f_fix.R first.")
  if (make_plots && !requireNamespace("ggnewscale", quietly = TRUE)) {
    stop("Original circular styling requires ggnewscale, which is not installed. Run make_plots = FALSE for verified tables/tree coordinates; no packages are installed automatically.")
  }
  suppressPackageStartupMessages(lapply(required, library, character.only = TRUE))
  previous_wd <- getwd()
  on.exit(setwd(previous_wd), add = TRUE)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  plot_dir <- normalizePath(output_dir, winslash = "/", mustWork = TRUE)
  core_output_dir <- file.path(plot_dir, "core_by_setting")
  dir.create(core_output_dir, showWarnings = FALSE)
  plot_input_file <- file.path(plot_dir, "subtype_site_summary.csv")
  core_prev_thr <- 0.70
  core_mean_thr <- 1e-5
  site_order <- c("Hospital", "WWTP", "Community", "Wet market")
  core_statistics <- audited_panel_f_core(long, meta, core_prev_thr, core_mean_thr)
  canonical_tables <- audited_panel_f_tables(core_statistics)
  core_by_site <- dplyr::rename(core_statistics, Site = Setting)
  core_only_by_site <- dplyr::filter(core_by_site, Core)
  write.csv(core_by_site, file.path(core_output_dir, "02_core_subtype_by_site_FULL.csv"), row.names = FALSE)
  write.csv(core_only_by_site, file.path(core_output_dir, "03_core_subtype_by_site_ONLY_CORE.csv"), row.names = FALSE)
  write.csv(canonical_tables$setting_counts, file.path(plot_dir, "verified_core_counts.csv"), row.names = FALSE)
  write.csv(canonical_tables$totals, file.path(plot_dir, "verified_core_union_summary.csv"), row.names = FALSE)
  write.csv(canonical_tables$intersections, file.path(plot_dir, "verified_exclusive_intersections.csv"), row.names = FALSE)
  # Deterministic writes fail clearly; they never silently save an alternative path.
  safe_write_csv <- function(x, file, ...) {
    write.csv(x, file, row.names = FALSE, fileEncoding = "UTF-8", ...)
    invisible(file)
  }
  safe_write_xlsx <- function(x, file, ...) {
    openxlsx::write.xlsx(x, file, overwrite = TRUE, ...)
    invisible(file)
  }
  # Constant columns contribute zero distance; they must not become NaN columns.
  scale_tree_features <- function(x) {
    spread <- apply(x, 2L, stats::sd)
    spread[!is.finite(spread) | spread == 0] <- 1
    scale(x, center = colMeans(x), scale = spread)
  }
  non_empty <- function(x) !is.na(x) & trimws(x) != ""

  mix_hex <- function(cols) {
    rgb_mat <- grDevices::col2rgb(cols)
    rgb(
      red   = mean(rgb_mat[1, ]) / 255,
      green = mean(rgb_mat[2, ]) / 255,
      blue  = mean(rgb_mat[3, ]) / 255
    )
  }

  clean_label <- function(x) {
    x <- trimws(x)
    x <- gsub("[[:space:]]+", " ", x)
    x
  }
  df <- core_only_by_site

  need_cols <- c("Site", "Subtype", "prevalence", "mean_abd")
  miss_cols <- setdiff(need_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("core_only_by_site 缺少必要列：", paste(miss_cols, collapse = ", "))
  }

  df2 <- df %>%
    mutate(
      Site = as.character(Site),
      Subtype_full = as.character(Subtype),

      Type = ifelse(
        str_detect(Subtype_full, "__"),
        sub("__.*$", "", Subtype_full),
        "Unknown"
      ),

      Subtype_clean = ifelse(
        str_detect(Subtype_full, "__"),
        sub("^.*__", "", Subtype_full),
        Subtype_full
      ),

      mean_abd = as.numeric(mean_abd),
      prevalence = as.numeric(prevalence)
    ) %>%
    filter(Site %in% site_order)

  site_matrix <- df2 %>%
    distinct(Type, Subtype_clean, Subtype_full, Site) %>%
    mutate(value = 1L) %>%
    pivot_wider(
      id_cols = c(Type, Subtype_clean, Subtype_full),
      names_from = Site,
      values_from = value,
      values_fill = 0
    )

  for (s in site_order) {
    if (!s %in% colnames(site_matrix)) {
      site_matrix[[s]] <- 0L
    }
  }

  site_matrix <- site_matrix %>%
    mutate(
      Hospital = as.integer(Hospital),
      WWTP = as.integer(WWTP),
      Community = as.integer(Community),
      `Wet market` = as.integer(`Wet market`),

      `Site Hospital` = ifelse(Hospital == 1, "Hospital", ""),
      `Site WWTP` = ifelse(WWTP == 1, "WWTP", ""),
      `Site Community` = ifelse(Community == 1, "Community", ""),
      `Site Wet market` = ifelse(`Wet market` == 1, "Wet market", ""),

      n_core_sites = Hospital + WWTP + Community + `Wet market`,

      core_site_pattern = case_when(
        n_core_sites == 4 ~ "Shared_by_all_4_sites",
        n_core_sites == 3 ~ "Shared_by_3_sites",
        n_core_sites == 2 ~ "Shared_by_2_sites",
        Hospital == 1 & n_core_sites == 1 ~ "Hospital_specific",
        WWTP == 1 & n_core_sites == 1 ~ "WWTP_specific",
        Community == 1 & n_core_sites == 1 ~ "Community_specific",
        `Wet market` == 1 & n_core_sites == 1 ~ "Wet_market_specific",
        TRUE ~ "Other"
      )
    )

  stat_summary <- df2 %>%
    group_by(Type, Subtype_clean, Subtype_full) %>%
    summarise(
      mean_abd_mean = mean(mean_abd, na.rm = TRUE),
      mean_abd_max = max(mean_abd, na.rm = TRUE),
      mean_abd_min = min(mean_abd, na.rm = TRUE),
      prevalence_mean = mean(prevalence, na.rm = TRUE),
      prevalence_min = min(prevalence, na.rm = TRUE),
      prevalence_max = max(prevalence, na.rm = TRUE),
      .groups = "drop"
    )

  df_final <- site_matrix %>%
    left_join(
      stat_summary,
      by = c("Type", "Subtype_clean", "Subtype_full")
    ) %>%
    select(
      Type,
      Subtype = Subtype_clean,
      Subtype_full,
      `Site Hospital`,
      `Site WWTP`,
      `Site Community`,
      `Site Wet market`,
      Hospital,
      WWTP,
      Community,
      `Wet market`,
      n_core_sites,
      core_site_pattern,
      mean_abd_mean,
      mean_abd_max,
      mean_abd_min,
      prevalence_mean,
      prevalence_min,
      prevalence_max
    ) %>%
    arrange(
      desc(n_core_sites),
      Type,
      Subtype
    )

  safe_write_csv(
    df_final,
    file.path(plot_dir, "core_subtypeitol.csv"),
    quote = FALSE
  )

  safe_write_xlsx(
    df_final,
    file.path(plot_dir, "core_subtypeitol.xlsx"),
    rowNames = FALSE
  )

  plot_input <- df_final %>%
    transmute(
      Type = Type,
      Subtype = Subtype,
      Subtype_full = Subtype_full,
      `Site Hospital` = `Site Hospital`,
      `Site WWTP` = `Site WWTP`,
      `Site Community` = `Site Community`,
      `Site Wet market` = `Site Wet market`,
      mean_abd = mean_abd_mean
    )

  safe_write_csv(
    plot_input,
    plot_input_file,
    quote = FALSE
  )

  safe_write_xlsx(
    plot_input,
    file.path(plot_dir, "subtype_site_summary.xlsx"),
    rowNames = FALSE
  )

  safe_write_csv(
    df_final %>% filter(n_core_sites == 4),
    file.path(plot_dir, "core_subtype_shared_by_all_4_sites.csv"),
    quote = FALSE
  )

  safe_write_csv(
    df_final %>% filter(n_core_sites == 1),
    file.path(plot_dir, "core_subtype_site_specific.csv"),
    quote = FALSE
  )

  pattern_count <- df_final %>%
    count(core_site_pattern, name = "n_subtypes") %>%
    arrange(desc(n_subtypes))

  safe_write_csv(
    pattern_count,
    file.path(plot_dir, "core_subtype_site_pattern_counts.csv"),
    quote = FALSE
  )

  type_pattern_count <- df_final %>%
    count(Type, core_site_pattern, name = "n_subtypes") %>%
    arrange(Type, desc(n_subtypes))

  safe_write_csv(
    type_pattern_count,
    file.path(plot_dir, "core_subtype_type_pattern_counts.csv"),
    quote = FALSE
  )

  cat("\n====================================================\n")
  cat("Merged core subtype table finished.\n")
  cat("Plot input: ", plot_input_file, "\n")
  cat("Total unique core subtypes:", nrow(df_final), "\n")
  cat("Shared by all 4 sites:", sum(df_final$n_core_sites == 4), "\n")
  cat("Site-specific core:", sum(df_final$n_core_sites == 1), "\n")
  cat("====================================================\n\n")
  print(pattern_count)
  setwd(plot_dir)

  infile <- "subtype_site_summary.csv"

  out_pdf_circle <- "ARG_type_3quarter_circle_groupedTree_publication_ready.pdf"
  out_png_circle <- "ARG_type_3quarter_circle_groupedTree_publication_ready.png"
  out_pdf_upset  <- "ARG_subtype_upset_groupedTree_publication_ready.pdf"
  out_png_upset  <- "ARG_subtype_upset_groupedTree_publication_ready.png"
  out_csv_sum    <- "ARG_type_3quarter_summary_groupedTree_publication_ready.csv"

  site_palette <- c(
    "Hospital"   = "#E06C75",
    "WWTP"       = "#4E79A7",
    "Community"  = "#6EC5C1",
    "Wet market" = "#E9C95B"
  )

  my_cols_cns <- c(
    "Aminoglycoside" = "#4E79A7",
    "Antibacterial fatty acid" = "#7FB8B5",
    "Bacitracin" = "#D6B97A",
    "Bicyclomycin" = "#8BCB88",
    "Bleomycin" = "#E39B76",
    "Chloramphenicol" = "#B5658D",
    "Defensin" = "#B7D989",
    "Edeine" = "#8FB6E8",
    "Factumycin" = "#D17C98",
    "Florfenicol" = "#AFC7E8",
    "Fosfomycin" = "#B7E0E5",
    "Fusidic acid" = "#E6BAC4",
    "MLS" = "#74C2B3",
    "Multidrug" = "#7E6BB5",
    "Mupirocin" = "#5F9EA0",
    "Novobiocin" = "#72B77E",
    "Other peptide antibiotics" = "#9FD58A",
    "Pleuromutilin/Tiamulin" = "#E4CF62",
    "Polymyxin" = "#75C0C1",
    "Puromycin" = "#9DBDE0",
    "Quinolone" = "#9A84D6",
    "Rifamycin" = "#A9C98A",
    "Streptothricin" = "#6DBDE3",
    "Sulfonamide" = "#C8A24B",
    "Tetracenomycin C" = "#D8BB92",
    "Tetracycline" = "#E8D98F",
    "Trimethoprim" = "#E39A9A",
    "Tunicamycin" = "#9AA3A8",
    "Vancomycin" = "#AFC9A0",
    "beta_lactam" = "#7DA34D",
    "Others" = "#D0D0D0"
  )

  type_map <- c(
    "aminoglycoside" = "Aminoglycoside",
    "antibacterial fatty acid" = "Antibacterial fatty acid",
    "bacitracin" = "Bacitracin",
    "bicyclomycin" = "Bicyclomycin",
    "bleomycin" = "Bleomycin",
    "chloramphenicol" = "Chloramphenicol",
    "defensin" = "Defensin",
    "edeine" = "Edeine",
    "factumycin" = "Factumycin",
    "florfenicol" = "Florfenicol",
    "fosfomycin" = "Fosfomycin",
    "fusidic acid" = "Fusidic acid",
    "macrolide-lincosamide-streptogramin" = "MLS",
    "mls" = "MLS",
    "multidrug" = "Multidrug",
    "mupirocin" = "Mupirocin",
    "novobiocin" = "Novobiocin",
    "other peptide antibiotics" = "Other peptide antibiotics",
    "pleuromutilin/tiamulin" = "Pleuromutilin/Tiamulin",
    "polymyxin" = "Polymyxin",
    "puromycin" = "Puromycin",
    "quinolone" = "Quinolone",
    "rifamycin" = "Rifamycin",
    "streptothricin" = "Streptothricin",
    "sulfonamide" = "Sulfonamide",
    "tetracenomycin c" = "Tetracenomycin C",
    "tetracycline" = "Tetracycline",
    "trimethoprim" = "Trimethoprim",
    "tunicamycin" = "Tunicamycin",
    "vancomycin" = "Vancomycin",
    "beta_lactam" = "beta_lactam",
    "others" = "Others"
  )

  # In-memory canonical plotting table: never reread a stale CSV.
  type_map[c("antibacterial_fatty_acid", "other_peptide_antibiotics",
             "pleuromutilin_tiamulin", "fusidic_acid", "tetracenomycin_c")] <-
    c("Antibacterial fatty acid", "Other peptide antibiotics",
      "Pleuromutilin/Tiamulin", "Fusidic acid", "Tetracenomycin C")
  df <- plot_input

  required_cols <- c(
    "Type", "Subtype", "Site Hospital", "Site WWTP",
    "Site Community", "Site Wet market", "mean_abd"
  )

  miss_cols <- setdiff(required_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("画图输入文件缺少这些列：", paste(miss_cols, collapse = ", "))
  }

  df2 <- df %>%
    mutate(
      Type_raw = clean_label(Type),
      Type_key = tolower(Type_raw),
      Type_std = ifelse(Type_key %in% names(type_map), type_map[Type_key], Type_raw),
      Type_std = clean_label(Type_std),
      Subtype = clean_label(Subtype),
      mean_abd = as.numeric(mean_abd),
      Hospital = as.integer(non_empty(`Site Hospital`)),
      WWTP = as.integer(non_empty(`Site WWTP`)),
      Community = as.integer(non_empty(`Site Community`)),
      `Wet market` = as.integer(non_empty(`Site Wet market`))
    ) %>%
    group_by(Type_std, Subtype) %>%
    summarise(
      mean_abd = mean(mean_abd, na.rm = TRUE),
      Hospital = max(Hospital, na.rm = TRUE),
      WWTP = max(WWTP, na.rm = TRUE),
      Community = max(Community, na.rm = TRUE),
      `Wet market` = max(`Wet market`, na.rm = TRUE),
      .groups = "drop"
    )

  min_pos <- suppressWarnings(min(df2$mean_abd[df2$mean_abd > 0], na.rm = TRUE))
  if (!is.finite(min_pos)) min_pos <- 1e-6

  df2 <- df2 %>%
    mutate(
      mean_abd_adj = ifelse(is.na(mean_abd) | mean_abd <= 0, min_pos / 10, mean_abd),
      mean_log10 = log10(mean_abd_adj)
    )

  type_summary <- df2 %>%
    count(Type_std, name = "n_subtype") %>%
    left_join(
      df2 %>%
        group_by(Type_std) %>%
        summarise(mean_log10_type = mean(mean_log10, na.rm = TRUE), .groups = "drop"),
      by = "Type_std"
    ) %>%
    arrange(desc(n_subtype), Type_std)

  # Retain every ARG type and every core subtype.
  selected_types <- type_summary %>%
    pull(Type_std)

  df_top <- df2 %>%
    filter(Type_std %in% selected_types)

  type_summary <- df_top %>%
    group_by(Type_std) %>%
    summarise(
      n_subtype = n(),
      mean_log10_type = mean(mean_log10, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(desc(n_subtype), Type_std) %>%
    mutate(
      type_col = ifelse(Type_std %in% names(my_cols_cns), my_cols_cns[Type_std], "#BDBDBD")
    )

  type_levels <- type_summary$Type_std

  safe_write_csv(
    type_summary,
    file.path(plot_dir, out_csv_sum),
    quote = FALSE
  )

  gap_n <- 4

  df_ord <- df_top %>%
    mutate(Type_std = factor(Type_std, levels = type_levels)) %>%
    arrange(Type_std, desc(mean_log10), Subtype) %>%
    mutate(Type_std = as.character(Type_std))

  plot_list <- list()

  for (i in seq_along(type_levels)) {

    tp <- type_levels[i]
    tmp <- df_ord %>% filter(Type_std == tp)
    tmp$is_gap <- FALSE
    plot_list[[length(plot_list) + 1]] <- tmp

    if (i < length(type_levels)) {
      gap_df <- data.frame(
        Type_std = paste0("gap_", i),
        Subtype = paste0("gap_", i, "_", seq_len(gap_n)),
        mean_abd = NA_real_,
        Hospital = NA_integer_,
        WWTP = NA_integer_,
        Community = NA_integer_,
        `Wet market` = NA_integer_,
        mean_abd_adj = NA_real_,
        mean_log10 = NA_real_,
        is_gap = TRUE,
        stringsAsFactors = FALSE
      )
      plot_list[[length(plot_list) + 1]] <- gap_df
    }
  }

  plot_df_base <- bind_rows(plot_list)
  n_real <- nrow(plot_df_base)
  blank_n <- ceiling(n_real / 3)

  blank_df <- data.frame(
    Type_std = "blank",
    Subtype = paste0("blank_", seq_len(blank_n)),
    mean_abd = NA_real_,
    Hospital = NA_integer_,
    WWTP = NA_integer_,
    Community = NA_integer_,
    `Wet market` = NA_integer_,
    mean_abd_adj = NA_real_,
    mean_log10 = NA_real_,
    is_gap = TRUE,
    stringsAsFactors = FALSE
  )

  plot_df <- bind_rows(plot_df_base, blank_df) %>%
    mutate(idx = row_number())

  n_total <- nrow(plot_df)

  type_sector <- plot_df %>%
    filter(!is_gap) %>%
    group_by(Type_std) %>%
    summarise(
      xmin = min(idx) - 0.5,
      xmax = max(idx) + 0.5,
      xmid = mean(c(min(idx), max(idx))),
      n_subtype = n(),
      .groups = "drop"
    ) %>%
    left_join(type_summary, by = c("Type_std", "n_subtype")) %>%
    mutate(
      type_col = ifelse(Type_std %in% names(my_cols_cns), my_cols_cns[Type_std], "#BDBDBD")
    )

  r_tree_root <- 0.22
  r_tree_group <- 0.29
  r_tree_tip  <- 0.37

  r_type_in  <- 0.525
  type_h     <- 0.070
  r_type_out <- r_type_in + type_h

  site_h   <- 0.024
  site_gap <- 0.008

  r_wet  <- r_type_out + 0.028 + site_h / 2
  r_comm <- r_wet  + site_h + site_gap
  r_wwtp <- r_comm + site_h + site_gap
  r_hosp <- r_wwtp + site_h + site_gap

  r_bar_base <- r_hosp + site_h / 2 + 0.032
  r_bar_top  <- r_bar_base + 0.12

  ring_long <- plot_df %>%
    filter(!is_gap) %>%
    select(idx, Type_std, Subtype, Hospital, WWTP, Community, `Wet market`) %>%
    pivot_longer(
      cols = c(Hospital, WWTP, Community, `Wet market`),
      names_to = "Site",
      values_to = "present"
    ) %>%
    mutate(
      Site = factor(Site, levels = site_order),
      ring_y = case_when(
        Site == "Hospital"   ~ r_hosp,
        Site == "WWTP"       ~ r_wwtp,
        Site == "Community"  ~ r_comm,
        Site == "Wet market" ~ r_wet
      ),
      fill_col = case_when(
        Site == "Hospital"   & present == 1 ~ unname(site_palette["Hospital"]),
        Site == "WWTP"       & present == 1 ~ unname(site_palette["WWTP"]),
        Site == "Community"  & present == 1 ~ unname(site_palette["Community"]),
        Site == "Wet market" & present == 1 ~ unname(site_palette["Wet market"]),
        TRUE ~ "#ECECEC"
      )
    )

  bar_range <- range(type_sector$mean_log10_type, na.rm = TRUE)

  if (!all(is.finite(bar_range)) || diff(bar_range) == 0) {
    type_sector$bar_top <- r_bar_base + 0.05
  } else {
    type_sector$bar_top <- scales::rescale(
      type_sector$mean_log10_type,
      to = c(r_bar_base + 0.015, r_bar_top),
      from = bar_range
    )
  }

  tree_mat <- df_top %>%
    group_by(Type_std) %>%
    summarise(
      Hospital = mean(Hospital, na.rm = TRUE),
      WWTP = mean(WWTP, na.rm = TRUE),
      Community = mean(Community, na.rm = TRUE),
      `Wet market` = mean(`Wet market`, na.rm = TRUE),
      mean_log10_type = mean(mean_log10, na.rm = TRUE),
      n_subtype = n(),
      .groups = "drop"
    ) %>%
    distinct(Type_std, .keep_all = TRUE)

  tree_num <- tree_mat %>%
    select(Hospital, WWTP, Community, `Wet market`, mean_log10_type, n_subtype) %>%
    as.matrix()

  rownames(tree_num) <- tree_mat$Type_std

  tree_num_scaled <- scale_tree_features(tree_num)
  rownames(tree_num_scaled) <- rownames(tree_num)

  hc <- hclust(dist(tree_num_scaled), method = "ward.D2")

  leaf_pos <- type_sector %>%
    select(Type_std, xmid, type_col) %>%
    distinct() %>%
    mutate(Type_std = as.character(Type_std))

  k_cluster <- min(4, nrow(leaf_pos))
  cluster_assign <- cutree(hc, k = k_cluster)

  cluster_df <- data.frame(
    Type_std = names(cluster_assign),
    cluster_id = as.integer(cluster_assign),
    stringsAsFactors = FALSE
  ) %>%
    left_join(leaf_pos, by = "Type_std") %>%
    arrange(cluster_id, xmid)

  cluster_sector <- cluster_df %>%
    group_by(cluster_id) %>%
    summarise(
      xmin = min(xmid),
      xmax = max(xmid),
      xmid = mean(xmid),
      n_type = n(),
      .groups = "drop"
    ) %>%
    arrange(xmid) %>%
    mutate(rank_id = row_number())

  cluster_centers <- cluster_sector$xmid
  global_center <- mean(range(cluster_centers))

  cluster_offsets <- rep(0, length(cluster_centers))

  if (length(cluster_centers) > 1) {
    for (i in seq_along(cluster_centers)) {
      dir_sign <- ifelse(cluster_centers[i] >= global_center, 1, -1)
      cluster_offsets[i] <- dir_sign * (0.45 + 0.10 * abs(i - mean(seq_along(cluster_centers))))
    }
  }

  cluster_sector$root_x <- cluster_sector$xmid + cluster_offsets

  leaf_cluster <- leaf_pos %>%
    left_join(cluster_df %>% select(Type_std, cluster_id), by = "Type_std") %>%
    left_join(cluster_sector %>% select(cluster_id, root_x), by = "cluster_id")

  tree_segments_radial <- data.frame()
  tree_segments_arc <- data.frame()
  leaf_points <- data.frame()

  for (cid in sort(unique(cluster_df$cluster_id))) {

    cluster_types <- cluster_df %>%
      filter(cluster_id == cid) %>%
      arrange(xmid) %>%
      pull(Type_std)

    sub_leaf <- leaf_pos %>%
      filter(Type_std %in% cluster_types) %>%
      arrange(match(Type_std, cluster_types))

    sub_root_x <- unique(leaf_cluster$root_x[leaf_cluster$cluster_id == cid])[1]

    leaf_points <- bind_rows(
      leaf_points,
      sub_leaf %>%
        transmute(
          x = xmid,
          y = r_tree_tip,
          type_col = type_col
        )
    )

    if (nrow(sub_leaf) == 1) {
      tree_segments_radial <- bind_rows(
        tree_segments_radial,
        data.frame(
          x = sub_leaf$xmid,
          xend = sub_leaf$xmid,
          y = r_tree_tip,
          yend = r_tree_group
        ),
        data.frame(
          x = sub_leaf$xmid,
          xend = sub_root_x,
          y = r_tree_group,
          yend = r_tree_group
        ),
        data.frame(
          x = sub_root_x,
          xend = sub_root_x,
          y = r_tree_group,
          yend = r_tree_root
        )
      )
      next
    }

    # Align the numeric rows before assigning cluster type names.
    sub_num <- tree_mat %>%
      filter(Type_std %in% cluster_types) %>%
      arrange(match(Type_std, cluster_types)) %>%
      select(Hospital, WWTP, Community, `Wet market`, mean_log10_type, n_subtype) %>%
      as.matrix()

    rownames(sub_num) <- cluster_types

    sub_num_scaled <- scale_tree_features(sub_num)
    rownames(sub_num_scaled) <- cluster_types

    sub_hc <- hclust(dist(sub_num_scaled), method = "ward.D2")
    # Negative merge indices refer to original input rows, not leaf display order.
    sub_labels <- sub_hc$labels

    sub_leaf2 <- sub_leaf[match(sub_labels, sub_leaf$Type_std), , drop = FALSE]
    stopifnot(identical(as.character(sub_leaf2$Type_std), as.character(sub_hc$labels)))
    rownames(sub_leaf2) <- NULL

    sub_merge <- sub_hc$merge
    sub_height <- sub_hc$height

    if (max(sub_height) == min(sub_height)) {
      sub_height_scaled <- rep(0.5, length(sub_height))
    } else {
      sub_height_scaled <- (sub_height - min(sub_height)) / (max(sub_height) - min(sub_height))
    }

    sub_height_scaled <- sub_height_scaled ^ 0.55

    sub_env <- new.env(parent = emptyenv())

    for (i in seq_len(nrow(sub_leaf2))) {
      assign(
        as.character(-i),
        list(xmid = sub_leaf2$xmid[i], r = r_tree_tip),
        envir = sub_env
      )
    }

    for (i in seq_len(nrow(sub_merge))) {

      left_id  <- sub_merge[i, 1]
      right_id <- sub_merge[i, 2]

      left_node  <- get(as.character(left_id), envir = sub_env)
      right_node <- get(as.character(right_id), envir = sub_env)

      new_r <- r_tree_tip - sub_height_scaled[i] * (r_tree_tip - r_tree_group)

      x1 <- min(left_node$xmid, right_node$xmid)
      x2 <- max(left_node$xmid, right_node$xmid)
      mid0 <- mean(c(x1, x2))
      span <- x2 - x1

      local_shift <- ifelse(span < 1.5, 0, 0.08 * span * ifelse(mid0 >= sub_root_x, 1, -1))
      parent_x <- mid0 + local_shift

      parent_x <- max(x1 + 0.18 * span, min(x2 - 0.18 * span, parent_x))

      tree_segments_radial <- bind_rows(
        tree_segments_radial,
        data.frame(
          x = left_node$xmid,
          xend = left_node$xmid,
          y = left_node$r,
          yend = new_r
        ),
        data.frame(
          x = right_node$xmid,
          xend = right_node$xmid,
          y = right_node$r,
          yend = new_r
        )
      )

      tree_segments_arc <- bind_rows(
        tree_segments_arc,
        data.frame(
          x = seq(min(left_node$xmid, parent_x), max(left_node$xmid, parent_x), length.out = 100),
          y = new_r,
          branch_id = paste0("C", cid, "_L_", i)
        ),
        data.frame(
          x = seq(min(right_node$xmid, parent_x), max(right_node$xmid, parent_x), length.out = 100),
          y = new_r,
          branch_id = paste0("C", cid, "_R_", i)
        )
      )

      assign(
        as.character(i),
        list(xmid = parent_x, r = new_r),
        envir = sub_env
      )
    }

    cluster_top <- get(as.character(nrow(sub_merge)), envir = sub_env)

    tree_segments_radial <- bind_rows(
      tree_segments_radial,
      data.frame(
        x = cluster_top$xmid,
        xend = cluster_top$xmid,
        y = cluster_top$r,
        yend = r_tree_group
      )
    )

    tree_segments_arc <- bind_rows(
      tree_segments_arc,
      data.frame(
        x = seq(min(cluster_top$xmid, sub_root_x), max(cluster_top$xmid, sub_root_x), length.out = 120),
        y = r_tree_group,
        branch_id = paste0("C", cid, "_to_root")
      )
    )

    tree_segments_radial <- bind_rows(
      tree_segments_radial,
      data.frame(
        x = sub_root_x,
        xend = sub_root_x,
        y = r_tree_group,
        yend = r_tree_root
      )
    )
  }

  cluster_root_df <- cluster_sector %>%
    transmute(
      x = root_x,
      xend = root_x,
      y = r_tree_root,
      yend = r_tree_root - 0.01
    )

  tree_segments_radial <- bind_rows(tree_segments_radial, cluster_root_df)

  # Export the exact numerical inputs and corrected tree geometry for inspection.
  write.csv(tree_mat, file.path(plot_dir, "original_tree_numeric_features.csv"), row.names = FALSE)
  write.csv(tree_segments_radial, file.path(plot_dir, "original_tree_radial_segments.csv"), row.names = FALSE)
  write.csv(tree_segments_arc, file.path(plot_dir, "original_tree_arc_segments.csv"), row.names = FALSE)
  write.csv(cluster_df, file.path(plot_dir, "original_tree_cluster_membership.csv"), row.names = FALSE)
  write.csv(plot_df, file.path(plot_dir, "original_circle_all_subtype_plot_order.csv"), row.names = FALSE)
  writeLines(c(
    paste0("Circle and UpSet display every core subtype: ", nrow(df2), "; ARG types: ", nrow(type_summary), "."),
    "Core: within-setting prevalence >= 70% and mean abundance >= 1e-5 copies/cell; sample zeros retained.",
    "Original grouped circular tree layout, Ward.D2 clustering, four-cluster cut, and branching geometry retained.",
    "Fixed numeric-row/type-label alignment and negative hclust merge indices; constant feature columns scale to zero.",
    "Tree is a descriptive clustering of setting-core fractions, conditional log10 abundance and subtype counts, not a phylogeny.",
    "Outer bars retain the original statistic: type average of log10 subtype means, with subtype means averaged only across settings where that subtype is core.",
    "All-source sample-weighted abundance is a different statistic and is not represented by these bars.",
    "Only the corrected ALL-core UpSet is rendered. Its left bars and membership rows share one y map.",
    "No external intermediate CSV is used; plot input stays in memory."
  ), file.path(plot_dir, "original_style_scope_and_definitions.txt"))
  stopifnot(nrow(df2) == nrow(canonical_tables$membership))
  if (!make_plots) {
    return(invisible(list(tables = canonical_tables, plot_input = plot_input,
                          tree_features = tree_mat, tree_radial = tree_segments_radial,
                          tree_arcs = tree_segments_arc, tree_clusters = cluster_df,
                          output_dir = plot_dir, rendered = FALSE)))
  }

  p_circle <- ggplot() +
    geom_path(
      data = tree_segments_arc,
      aes(x = x, y = y, group = branch_id),
      color = "#8F8F8F",
      linewidth = 0.62,
      lineend = "round"
    ) +
    geom_segment(
      data = tree_segments_radial,
      aes(x = x, xend = xend, y = y, yend = yend),
      color = "#8F8F8F",
      linewidth = 0.62,
      lineend = "round"
    ) +
    geom_point(
      data = leaf_points,
      aes(x = x, y = y),
      shape = 16,
      size = 3.8,
      color = leaf_points$type_col
    ) +

    ggnewscale::new_scale_fill() +
    geom_rect(
      data = type_sector,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = r_type_in, ymax = r_type_out,
        fill = Type_std
      ),
      color = "white",
      linewidth = 0.45
    ) +
    scale_fill_manual(
      values = setNames(type_sector$type_col, type_sector$Type_std),
      guide = "none"
    ) +

    ggnewscale::new_scale_fill() +
    geom_tile(
      data = ring_long,
      aes(x = idx, y = ring_y, fill = fill_col),
      width = 0.992,
      height = site_h,
      color = NA
    ) +
    scale_fill_identity() +

    ggnewscale::new_scale_fill() +
    geom_rect(
      data = type_sector,
      aes(
        xmin = xmin, xmax = xmax,
        ymin = r_bar_base, ymax = bar_top,
        fill = Type_std
      ),
      color = "white",
      linewidth = 0.50
    ) +
    scale_fill_manual(
      values = setNames(type_sector$type_col, type_sector$Type_std),
      guide = "none"
    ) +

    coord_polar(theta = "x", start = -pi/2, clip = "off") +
    xlim(0.5, n_total + 0.5) +
    ylim(0.10, r_bar_top + 0.05) +
    theme_void() +
    theme(
      plot.margin = margin(4, 4, 4, 4),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  p_circle <- p_circle + labs(caption = paste0("All ", nrow(df2), " core subtypes; ", nrow(type_summary), " ARG types. Corrected descriptive clustering."))
  ggsave(out_pdf_circle, p_circle, width = 10.5, height = 10.5, bg = "white")
  ggsave(out_png_circle, p_circle, width = 10.5, height = 10.5, dpi = 700, bg = "white")
  pkg_needed <- c(
    "dplyr", "tidyr", "stringr", "ggplot2",
    "patchwork", "purrr", "tibble", "scales"
  )

  for (p in pkg_needed) {
  }

  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(patchwork)
  library(purrr)
  library(tibble)
  library(scales)

  root_dir <- plot_dir
  setwd(root_dir)

  infile <- "subtype_site_summary.csv"

  out_pdf_upset <- "ARG_subtype_upset_ALL_core_subtypes_FIXED_label_order.pdf"
  out_png_upset <- "ARG_subtype_upset_ALL_core_subtypes_FIXED_label_order.png"
  out_check_csv <- "ARG_subtype_upset_ALL_core_subtypes_check_table.csv"

  row_order <- c("Wet market", "Community", "WWTP", "Hospital")

  y_map <- c(
    "Hospital"   = 1,
    "WWTP"       = 2,
    "Community"  = 3,
    "Wet market" = 4
  )

  site_palette <- c(
    "Hospital"   = "#E06C75",
    "WWTP"       = "#4E79A7",
    "Community"  = "#6EC5C1",
    "Wet market" = "#E9C95B"
  )

  non_empty <- function(x) {
    !is.na(x) & trimws(as.character(x)) != ""
  }

  mix_hex <- function(cols) {
    cols <- cols[!is.na(cols)]
    if (length(cols) == 1) return(cols)
    rgb_mat <- grDevices::col2rgb(cols)
    rgb(
      red   = mean(rgb_mat[1, ]) / 255,
      green = mean(rgb_mat[2, ]) / 255,
      blue  = mean(rgb_mat[3, ]) / 255
    )
  }

  # Reuse the same canonical in-memory table as the circular plot.
  df <- plot_input

  required_cols <- c(
    "Type",
    "Subtype",
    "Site Hospital",
    "Site WWTP",
    "Site Community",
    "Site Wet market",
    "mean_abd"
  )

  miss_cols <- setdiff(required_cols, colnames(df))
  if (length(miss_cols) > 0) {
    stop("输入文件缺少必要列：", paste(miss_cols, collapse = ", "))
  }

  upset_df <- df %>%
    mutate(
      Type = as.character(Type),
      Subtype = as.character(Subtype),
      mean_abd = as.numeric(mean_abd),

      Hospital = non_empty(`Site Hospital`),
      WWTP = non_empty(`Site WWTP`),
      Community = non_empty(`Site Community`),
      `Wet market` = non_empty(`Site Wet market`)
    ) %>%
    select(Type, Subtype, Hospital, WWTP, Community, `Wet market`) %>%
    distinct() %>%
    filter(Hospital | WWTP | Community | `Wet market`)

  set_size_check <- tibble(
    group = c("Hospital", "WWTP", "Community", "Wet market"),
    set_size = c(
      sum(upset_df$Hospital),
      sum(upset_df$WWTP),
      sum(upset_df$Community),
      sum(upset_df$`Wet market`)
    )
  )

  write.csv(
    set_size_check,
    out_check_csv,
    row.names = FALSE,
    quote = FALSE,
    fileEncoding = "UTF-8"
  )

  cat("\nSet size check:\n")
  print(set_size_check)

  cat("\nTotal unique core subtypes in UpSet:", nrow(upset_df), "\n")

  comb_df <- upset_df %>%
    transmute(
      `Wet market`,
      Community,
      WWTP,
      Hospital
    ) %>%
    count(`Wet market`, Community, WWTP, Hospital, name = "intersection_size") %>%
    filter(intersection_size > 0) %>%
    rowwise() %>%
    mutate(
      present_sets = list(c(
        if (`Wet market`) "Wet market",
        if (Community) "Community",
        if (WWTP) "WWTP",
        if (Hospital) "Hospital"
      )),
      fill_col = if (length(present_sets) == 1) {
        unname(site_palette[present_sets])
      } else {
        mix_hex(unname(site_palette[present_sets]))
      },
      degree = sum(c(`Wet market`, Community, WWTP, Hospital))
    ) %>%
    ungroup() %>%
    arrange(desc(degree), desc(intersection_size)) %>%
    mutate(
      intersection_name = paste(
        ifelse(`Wet market`, "1", NA),
        ifelse(Community, "2", NA),
        ifelse(WWTP, "3", NA),
        ifelse(Hospital, "4", NA),
        sep = "-"
      ),
      intersection_name = gsub("NA-", "", intersection_name),
      intersection_name = gsub("-NA", "", intersection_name),
      intersection_name = gsub("NA", "", intersection_name),
      intersection_name = gsub("--", "-", intersection_name),
      intersection_name = factor(intersection_name, levels = intersection_name)
    )

  cat("\nIntersection size sum:", sum(comb_df$intersection_size), "\n")
  cat("This should equal total unique core subtypes:", nrow(upset_df), "\n")

  set_size_df <- tibble(
    group = row_order,
    set_size = c(
      sum(upset_df$`Wet market`),
      sum(upset_df$Community),
      sum(upset_df$WWTP),
      sum(upset_df$Hospital)
    )
  ) %>%
    mutate(
      y = unname(y_map[group])
    )

  top_bar <- ggplot(
    comb_df,
    aes(
      x = intersection_name,
      y = intersection_size,
      fill = fill_col
    )
  ) +
    geom_col(width = 0.88, color = NA) +
    geom_text(
      aes(label = intersection_size),
      vjust = ifelse(comb_df$intersection_size > 30, 1.12, -0.30),
      color = ifelse(comb_df$intersection_size > 30, "white", "black"),
      size = 6
    ) +
    scale_fill_identity() +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    labs(y = "Intersection size", x = NULL) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 18, face = "bold"),
      axis.text.y = element_text(size = 13, color = "#3A3A3A"),
      plot.margin = margin(6, 8, 0, 0)
    )

  left_bar <- ggplot(
    set_size_df,
    aes(x = set_size, y = y, fill = group)
  ) +
    geom_col(width = 0.62, orientation = "y") +
    geom_text(
      aes(label = set_size),
      hjust = -0.18,
      color = "#3A3A3A",
      size = 5.6,
      fontface = "bold"
    ) +
    scale_fill_manual(values = site_palette, guide = "none") +
    scale_y_continuous(
      breaks = unname(y_map[row_order]),
      labels = row_order,
      expand = c(0, 0)
    ) +
    coord_cartesian(ylim = c(0.5, 4.5), clip = "off") +
    scale_x_reverse(expand = expansion(mult = c(0.18, 0))) +
    labs(x = "Set size", y = NULL) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      axis.title.x = element_text(size = 18, face = "bold"),
      axis.text.x = element_text(size = 13, color = "#3A3A3A"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.margin = margin(0, 0, 8, 0)
    )

  matrix_df <- comb_df %>%
    mutate(
      x = intersection_name
    ) %>%
    select(x, `Wet market`, Community, WWTP, Hospital) %>%
    pivot_longer(
      cols = c(`Wet market`, Community, WWTP, Hospital),
      names_to = "group",
      values_to = "present"
    ) %>%
    mutate(
      y = unname(y_map[group])
    )

  line_df <- comb_df %>%
    rowwise() %>%
    mutate(
      x = intersection_name,
      y_present = list(c(
        if (`Wet market`) unname(y_map["Wet market"]),
        if (Community) unname(y_map["Community"]),
        if (WWTP) unname(y_map["WWTP"]),
        if (Hospital) unname(y_map["Hospital"])
      )),
      ymin = min(unlist(y_present)),
      ymax = max(unlist(y_present))
    ) %>%
    ungroup()

  matrix_plot <- ggplot(matrix_df, aes(x = x, y = y)) +
    geom_rect(
      data = tibble(
        ymin = c(0.5, 2.5),
        ymax = c(1.5, 3.5)
      ),
      inherit.aes = FALSE,
      aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
      fill = "#F4F4F4",
      color = NA
    ) +
    geom_segment(
      data = line_df,
      aes(x = x, xend = x, y = ymin, yend = ymax),
      inherit.aes = FALSE,
      linewidth = 0.8,
      color = "black"
    ) +
    geom_point(
      aes(fill = present),
      shape = 21,
      size = 4.5,
      stroke = 0.6,
      color = "#9E9E9E"
    ) +
    scale_fill_manual(
      values = c("TRUE" = "black", "FALSE" = "#D9D9D9"),
      guide = "none"
    ) +
    scale_x_discrete(expand = c(0.01, 0.01)) +
    scale_y_continuous(
      breaks = unname(y_map[row_order]),
      labels = row_order,
      limits = c(0.5, 4.5),
      expand = c(0, 0)
    ) +
    labs(x = NULL, y = NULL) +
    theme_minimal(base_size = 18) +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold", color = "#3A3A3A"),
      axis.text.y = element_text(size = 16, color = "#3A3A3A"),
      axis.ticks = element_blank(),
      plot.margin = margin(0, 0, 8, 0)
    )

  top_left_blank <- ggplot() + theme_void()

  top_row <- top_left_blank + top_bar +
    plot_layout(widths = c(0.88, 4.32))

  bottom_row <- left_bar + matrix_plot +
    plot_layout(widths = c(0.88, 4.32))

  p_upset_final <- top_row / bottom_row +
    plot_layout(heights = c(2.15, 1))

  ggsave(
    out_pdf_upset,
    p_upset_final,
    width = 16,
    height = 9,
    bg = "white"
  )

  ggsave(
    out_png_upset,
    p_upset_final,
    width = 16,
    height = 9,
    dpi = 700,
    bg = "white"
  )

  print(p_upset_final)

  cat("\nSaved files:\n")
  cat(" - ", file.path(root_dir, out_pdf_upset), "\n", sep = "")
  cat(" - ", file.path(root_dir, out_png_upset), "\n", sep = "")
  cat(" - ", file.path(root_dir, out_check_csv), "\n", sep = "")

  invisible(list(tables = canonical_tables, circle = p_circle, upset = p_upset_final,
                 tree_features = tree_mat, tree_radial = tree_segments_radial,
                 tree_arcs = tree_segments_arc, output_dir = plot_dir, rendered = TRUE))
}


# ----- panels_gh_fix.R -----
# Audited Fig. 2g/h. Source after the shared canonical cohort has been loaded.
# All abundance values are ARG copies per microbial cell. Detection means > 0.
# Outputs are a review version; no source workbook or original figure is modified.

.gh_gene_name <- function(subtype) {
  toupper(trimws(sub("^.*__", "", as.character(subtype))))
}

.gh_marker_specs <- function() {
  data.frame(
    Marker_Group = c("tetX", "mcr", "vim", "van", "linezolid", "imp", "ndm", "kpc", "oxa48"),
    Marker_Label = c("tet(X)", "mcr", "VIM", "van", "Linezolid", "IMP", "NDM", "KPC", "OXA-48"),
    Gene_pattern = c(
      "^TET\\(X[0-9]*\\)$|^TETX[0-9]*$",
      "^MCR(?:[-._]|$)", "^VIM(?:[-._]|$)", "^VAN[A-Z0-9_.-]*$",
      "^CFR(?:$|[-._(]|[A-Z])|^OPTRA(?:$|[-._])|^POXTA(?:$|[-._])",
      "^IMP(?:[-._]|$)", "^NDM(?:[-._]|$)", "^KPC(?:[-._]|$)",
      "^OXA[-_]?48$"
    ), stringsAsFactors = FALSE
  )
}

.gh_caption_targets <- function() {
  # These are the original eight caption-specified alleles, not representatives
  # chosen from the observed abundance or prevalence ranking.
  data.frame(
    Marker_Group = c("vim", "ndm", "imp", "kpc", "van", "mcr", "tetX", "linezolid"),
    Gene = c("VIM-1", "NDM-10", "IMP-1", "KPC-2", "VANYM", "MCR-3.1", "TET(X3)", "OPTRA"),
    Gene_Label = c("blaVIM-1", "blaNDM-10", "blaIMP-1", "blaKPC-2", "vanYM", "mcr-3.1", "tet(X3)", "optrA"),
    stringsAsFactors = FALSE
  )
}

.gh_check_counts <- function(observed, expected, keys, value, description) {
  comparison <- dplyr::full_join(
    dplyr::select(observed, dplyr::all_of(c(keys, value))),
    dplyr::rename(expected, expected_value = dplyr::all_of(value)), by = keys
  )
  comparison$passed <- !is.na(comparison[[value]]) & !is.na(comparison$expected_value) &
    comparison[[value]] == comparison$expected_value
  if (!all(comparison$passed)) {
    print(comparison[!comparison$passed, , drop = FALSE])
    stop(description, " failed; do not export figures until the cohort or matching difference is resolved.")
  }
  comparison$Check <- description
  comparison
}

run_panels_gh <- function(
    canonical_long = .canonical_long,
    canonical_meta = .canonical_meta,
    out_dir = file.path(.figure_results_dir, "panels_gh")) {
  needed <- c("dplyr", "tidyr", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Required installed packages are unavailable: ", paste(missing, collapse = ", "))
  long_required <- c("Sample", "Subtype", "Abundance", "Setting", "City", "Sample_Date", "PhysicalSite")
  meta_required <- setdiff(long_required, c("Subtype", "Abundance"))
  if (!all(long_required %in% names(canonical_long)) || !all(meta_required %in% names(canonical_meta))) {
    stop("The canonical long table or metadata is missing required columns.")
  }
  settings <- c("Hospital", "WWTP", "Community", "Wet market")
  cities <- c("BJ", "CQ", "GY", "GZ", "HF", "HK", "HR", "NJ", "SZ", "XA", "XM", "ZB")
  months <- format(seq(as.Date("2024-11-01"), as.Date("2025-12-01"), by = "month"), "%Y-%m")
  palette <- c(Hospital = "#E06C75", WWTP = "#4E79A7", Community = "#6EC5C1", "Wet market" = "#E9C95B")
  meta <- dplyr::transmute(
    canonical_meta, Sample = as.character(Sample), Setting = as.character(Setting),
    City = as.character(City), Sample_Date = as.Date(Sample_Date), PhysicalSite = as.character(PhysicalSite)
  )
  if (nrow(meta) != 2614L || anyNA(meta) || anyDuplicated(meta$Sample) ||
      !setequal(meta$City, cities) || !setequal(meta$Setting, settings)) {
    stop("Fig. 2g/h require the validated 2614-sample, 12-city, four-setting canonical cohort.")
  }
  meta$YearMonth <- format(meta$Sample_Date, "%Y-%m")
  if (!all(meta$YearMonth %in% months)) stop("Canonical dates lie outside November 2024 to December 2025.")
  if (!is.numeric(canonical_long$Abundance) || anyNA(canonical_long$Abundance) ||
      any(!is.finite(canonical_long$Abundance)) || any(canonical_long$Abundance < 0)) {
    stop("Canonical abundance must be finite, nonmissing and nonnegative.")
  }
  if (!setequal(as.character(canonical_long$Sample), meta$Sample) ||
      nrow(canonical_long) != nrow(meta) * dplyr::n_distinct(canonical_long$Subtype)) {
    stop("The canonical abundance table must retain every sample and explicit zero observations.")
  }

  specs <- .gh_marker_specs()
  genes <- data.frame(Subtype = unique(as.character(canonical_long$Subtype)), stringsAsFactors = FALSE)
  genes$Gene <- .gh_gene_name(genes$Subtype)
  membership <- dplyr::bind_rows(lapply(seq_len(nrow(specs)), function(i) {
    hit <- genes[grepl(specs$Gene_pattern[i], genes$Gene, perl = TRUE), , drop = FALSE]
    hit$Marker_Group <- rep(specs$Marker_Group[i], nrow(hit))
    hit$Marker_Label <- rep(specs$Marker_Label[i], nrow(hit))
    hit$Gene_pattern <- rep(specs$Gene_pattern[i], nrow(hit))
    hit
  }))
  if (!setequal(membership$Marker_Group, specs$Marker_Group) || anyDuplicated(membership$Subtype)) {
    stop("Every marker must have unambiguous subtype membership.")
  }
  membership_all <- dplyr::left_join(genes, membership, by = c("Subtype", "Gene"))
  membership_all$Included_in_clinical9 <- !is.na(membership_all$Marker_Group)
  clinical <- canonical_long |>
    dplyr::select(Sample, Subtype, Abundance) |>
    dplyr::inner_join(membership, by = "Subtype")
  clinical_key_counts <- clinical |> dplyr::count(Sample, Subtype, name = "n_rows")
  if (nrow(clinical) != nrow(meta) * nrow(membership) || any(clinical_key_counts$n_rows != 1L)) {
    stop("Clinical data contain duplicate or missing sample/subtype observations.")
  }
  marker_observed <- clinical |>
    dplyr::group_by(Sample, Marker_Group) |>
    dplyr::summarise(detected = any(Abundance > 0), total_abundance = sum(Abundance), .groups = "drop")
  marker_sample <- tidyr::crossing(meta, specs[c("Marker_Group", "Marker_Label")]) |>
    dplyr::left_join(marker_observed, by = c("Sample", "Marker_Group"))
  if (anyNA(marker_sample$detected) || anyNA(marker_sample$total_abundance)) {
    stop("Marker-level completeness failed; missing measurements must not silently become zeros.")
  }
  setting_counts <- meta |> dplyr::count(Setting, name = "n_samples")
  detection <- marker_sample |>
    dplyr::group_by(Marker_Group, Marker_Label, Setting) |>
    dplyr::summarise(
      positive_n = sum(detected), n_sample_rows = dplyr::n(),
      mean_abundance_copies_per_cell = mean(total_abundance), .groups = "drop"
    ) |>
    dplyr::left_join(setting_counts, by = "Setting") |>
    dplyr::mutate(detection_fraction = positive_n / n_samples, detection_threshold = 0) |>
    dplyr::arrange(match(Marker_Group, specs$Marker_Group), match(Setting, settings))
  if (nrow(detection) != 36L || any(detection$n_sample_rows != detection$n_samples)) {
    stop("Fig. 2g must have 36 marker/setting rows with complete setting denominators.")
  }

  targets <- .gh_caption_targets()
  target_membership <- dplyr::inner_join(genes, targets, by = "Gene")
  if (nrow(target_membership) != 8L || !setequal(target_membership$Gene, targets$Gene)) {
    stop("The original eight caption-specified genes must each identify exactly one subtype.")
  }
  h_sample <- clinical |>
    dplyr::select(Sample, Subtype, Abundance) |>
    dplyr::inner_join(target_membership, by = "Subtype") |>
    dplyr::left_join(meta, by = "Sample")
  if (nrow(h_sample) != 2614L * 8L || anyNA(h_sample)) stop("Fig. 2h sample/gene table is incomplete.")
  monthly_counts <- meta |> dplyr::count(YearMonth, Setting, name = "n_samples") |>
    dplyr::arrange(YearMonth, match(Setting, settings))
  monthly_totals <- meta |> dplyr::count(YearMonth, name = "n_samples") |> dplyr::arrange(YearMonth)
  h_monthly <- h_sample |>
    dplyr::group_by(Marker_Group, Gene, Gene_Label, Subtype, Setting, YearMonth) |>
    dplyr::summarise(
      n_samples = dplyr::n(), positive_n = sum(Abundance > 0),
      mean_abundance_copies_per_cell = mean(Abundance),
      median_abundance_copies_per_cell = stats::median(Abundance), .groups = "drop"
    ) |>
    dplyr::mutate(detection_fraction = positive_n / n_samples) |>
    dplyr::arrange(match(Gene, targets$Gene), match(Setting, settings), YearMonth)
  if (nrow(monthly_counts) != 56L || nrow(h_monthly) != 448L) {
    stop("Fig. 2h requires all 14 months in all four settings for each of eight genes.")
  }

  # Independent expected counts were recomputed directly from the four source
  # workbooks and the intact XLSX metadata; these are not generated by this code.
  expected_g <- data.frame(
    Marker_Group = rep(specs$Marker_Group, each = 4L), Setting = rep(settings, 9L),
    positive_n = c(411,1349,415,435, 408,1346,414,431, 410,1343,415,432,
                   409,1343,411,431, 399,1330,405,424, 384,1099,255,196,
                   378,882,257,241, 359,626,127,100, 65,62,23,73)
  )
  expected_settings <- data.frame(Setting = settings, n_samples = c(411L, 1349L, 416L, 438L))
  expected_months <- data.frame(YearMonth = months,
    n_samples = c(42L,72L,143L,176L,205L,194L,194L,209L,229L,220L,241L,228L,227L,234L))
  check_g <- .gh_check_counts(detection, expected_g, c("Marker_Group", "Setting"), "positive_n", "clinical9_positive_counts")
  check_settings <- .gh_check_counts(setting_counts, expected_settings, "Setting", "n_samples", "setting_denominators")
  check_months <- .gh_check_counts(monthly_totals, expected_months, "YearMonth", "n_samples", "14_month_sample_counts")
  monthly_join_check <- dplyr::left_join(h_monthly, monthly_counts,
    by = c("YearMonth", "Setting"), suffix = c("", "_metadata"))
  if (any(monthly_join_check$n_samples != monthly_join_check$n_samples_metadata)) {
    stop("Fig. 2h means do not use all samples in their corresponding setting and month.")
  }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  table_dir <- file.path(out_dir, "tables")
  figure_dir <- file.path(out_dir, "figures")
  dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  write_table <- function(x, name) utils::write.csv(x, file.path(table_dir, name), row.names = FALSE, na = "", fileEncoding = "UTF-8")
  write_table(specs, "g_marker_definitions.csv")
  write_table(membership_all, "g_all_subtype_marker_membership.csv")
  write_table(membership, "g_marker_matched_subtypes.csv")
  write_table(marker_sample, "g_marker_detection_by_sample.csv")
  write_table(detection, "g_marker_detection_fraction_by_setting.csv")
  write_table(target_membership, "h_original_caption_subtypes.csv")
  write_table(h_sample, "h_caption_subtype_sample_abundance.csv")
  write_table(h_monthly, "h_monthly_mean_abundance_by_setting.csv")
  write_table(monthly_counts, "h_monthly_sample_counts_by_setting.csv")
  write_table(monthly_totals, "h_monthly_sample_counts_total.csv")
  write_table(check_g, "check_g_known_positive_counts.csv")
  write_table(check_settings, "check_setting_denominators.csv")
  write_table(check_months, "check_14_month_sample_counts.csv")

  g_plot_data <- detection
  g_plot_data$Setting <- factor(g_plot_data$Setting, levels = settings)
  g_plot_data$Marker_Label <- factor(g_plot_data$Marker_Label, levels = rev(specs$Marker_Label))
  p_g <- ggplot2::ggplot(g_plot_data, ggplot2::aes(x = Setting, y = Marker_Label)) +
    ggplot2::geom_point(ggplot2::aes(size = detection_fraction, fill = detection_fraction),
      shape = 21, colour = "white", stroke = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", detection_fraction)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#EAF0F7", high = "#E06C75", limits = c(0, 1), name = "Detection fraction") +
    ggplot2::scale_size(range = c(3, 12), limits = c(0, 1), name = "Detection fraction") +
    ggplot2::guides(size = "none", fill = ggplot2::guide_colourbar(
      title.position = "top", barwidth = grid::unit(3, "in"), barheight = grid::unit(0.12, "in"))) +
    ggplot2::labs(x = NULL, y = NULL,
      caption = "Detection: abundance > 0.\nDenominators include all samples in each setting.") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 35, hjust = 1), legend.position = "bottom",
      plot.caption = ggplot2::element_text(size = 8, hjust = 0))

  h_plot_data <- h_monthly
  h_plot_data$Month_Date <- as.Date(paste0(h_plot_data$YearMonth, "-01"))
  h_plot_data$Setting <- factor(h_plot_data$Setting, levels = settings)
  h_plot_data$Gene_Label <- factor(h_plot_data$Gene_Label, levels = targets$Gene_Label)
  p_h <- ggplot2::ggplot(h_plot_data,
    ggplot2::aes(x = Month_Date, y = mean_abundance_copies_per_cell, colour = Setting, group = Setting)) +
    ggplot2::geom_line(linewidth = 0.7) + ggplot2::geom_point(size = 1.2) +
    ggplot2::facet_wrap(~ Gene_Label, ncol = 2, scales = "free_y") +
    ggplot2::scale_colour_manual(values = palette, breaks = settings, drop = FALSE) +
    ggplot2::scale_x_date(breaks = as.Date(paste0(months, "-01")), date_labels = "%Y-%m",
      expand = ggplot2::expansion(mult = c(0.01, 0.01))) +
    ggplot2::scale_y_continuous(limits = c(0, NA), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::labs(x = NULL, y = "Mean abundance (copies per cell)", colour = "Setting",
      caption = "Monthly sample means include zero observations; points are joined by straight lines. Sample counts vary by month.") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 55, hjust = 1, size = 7),
      strip.background = ggplot2::element_blank(), legend.position = "bottom")
  save_both <- function(plot, stem, width, height) {
    ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot = plot, width = width, height = height, bg = "white")
    ggplot2::ggsave(file.path(figure_dir, paste0(stem, ".png")), plot = plot, width = width, height = height, dpi = 300, bg = "white")
  }
  save_both(p_g, "Fig2g_clinical9_detection_corrected", 6.8, 6.2)
  save_both(p_h, "Fig2h_caption8_monthly_abundance_corrected", 13, 12)
  message("Audited Fig. 2g/h: 2614 samples; 36 detection fractions; 448 monthly gene/setting means. Checks passed.")
  invisible(list(
    marker_membership = membership, g_detection = detection, h_monthly = h_monthly,
    monthly_sample_counts = monthly_counts, monthly_sample_totals = monthly_totals,
    checks = list(detection = check_g, settings = check_settings, months = check_months),
    plots = list(g = p_g, h = p_h), output_dir = out_dir
  ))
}


# ----- fig2_entrypoint.R -----
FIG2_DATA_ROOT <- Sys.getenv('FIG2_DATA_ROOT','E:/c_merge/2614')
BLANK_POLICY <- Sys.getenv('FIG2_BLANK_POLICY','confirmed_nondetection')
C_MODEL <- Sys.getenv('FIG2_C_MODEL','gaussian_original')
PANELS <- strsplit(Sys.getenv('FIG2_PANELS','a,b,c,d,e,f,g,h'),',',fixed=TRUE)[[1]]

run_fig2_audited <- function() {
  if(any(c('b','c') %in% PANELS) && !all(c('b','c') %in% PANELS))
    stop('Select both b and c to regenerate paired half-violin plots with fresh model letters.')
  if('c' %in% PANELS && C_MODEL!='gaussian_original')
    stop('Fig. 2c is confirmed as an original-scale Gaussian LMM; use FIG2_C_MODEL=gaussian_original.')
  args <- commandArgs(FALSE); filearg <- grep('^--file=',args,value=TRUE)
  if(length(filearg)) base <- dirname(normalizePath(sub('^--file=','',filearg[1]),winslash='/'))
  else base <- getwd()
  results <- Sys.getenv('ARG_RESULTS_DIR',file.path(base,'Fig2_corrected_results'))
  dir.create(results,recursive=TRUE,showWarnings=FALSE)
  results <- normalizePath(results,winslash='/')
  .figure_results_dir <<- results
  .external_data_root <<- Sys.getenv('ARG_EXTERNAL_DATA_ROOT','E:/')
  .external_path <<- function(...) file.path(.external_data_root,...)
  cohort <- load_fig2_cohort(FIG2_DATA_ROOT,file.path(results,'cohort_audit'),BLANK_POLICY)
  .canonical_long <<- cohort$long; .canonical_meta <<- cohort$meta
  if(BLANK_POLICY=='legacy_zero_review') warning('REVIEW OUTPUTS ONLY: reproducing old blank-to-zero convention pending biological confirmation.')
  if('b' %in% PANELS) fit_fig2_model(cohort$samples,'Richness','negative_binomial',file.path(results,'panel_b_NB'))
  if('c' %in% PANELS) {
    fit_fig2_model(cohort$samples,'Total_abundance','gaussian_original',file.path(results,'panel_c_LMM'))
  }
  if(any(c('b','c') %in% PANELS)) {
    run_panels_bc(cohort)
  }
  if('a' %in% PANELS) run_panel_a()
  if('d' %in% PANELS) run_panel_d()
  if('e' %in% PANELS) run_panel_e()
  if('f' %in% PANELS) {
    run_panel_f_original(make_plots=requireNamespace('ggnewscale',quietly=TRUE))
    if(!requireNamespace('ggnewscale',quietly=TRUE)) warning('2f tables/tree coordinates calculated; original-style rendering not run because ggnewscale is unavailable.')
  }
  if(any(c('g','h') %in% PANELS)) run_panels_gh()
  writeLines(capture.output(sessionInfo()),file.path(results,'R_sessionInfo.txt'))
  writeLines(c(paste('Blank policy:',BLANK_POLICY),paste('2c model:',C_MODEL),
    '2b: negative-binomial mixed model; setting+city+sampling month fixed, physical site random intercept.',
    '2c: original-scale Gaussian linear mixed model fitted by REML; same fixed/random effects; no log transformation or rounding.',
    'Six two-sided asymptotic Wald contrasts; BH adjustment across these six comparisons.',
    'Both b/c letter displays use the newly fitted mixed-model contrasts; displayed means are descriptive arithmetic means.',
    'Global comparisons retain author KW/Dunn approach and are exploratory; exact duplicate extra external records are removed by default and audited.',
    'Monthly profiles are descriptive; sample/site composition is not constant.',
    if(BLANK_POLICY=='confirmed_nondetection')
      'Author-confirmed encoding: blank subtype abundance cells and explicit numeric zeros both denote nondetection; blanks are encoded as zero.'
    else 'Non-default blank policy selected: consult cohort_audit/missing_value_policy.csv.'),file.path(results,'READ_ME_FIRST.txt'))
  invisible(cohort)
}

if(!identical(Sys.getenv('FIG2_FUNCTIONS_ONLY','0'),'1')) run_fig2_audited()
