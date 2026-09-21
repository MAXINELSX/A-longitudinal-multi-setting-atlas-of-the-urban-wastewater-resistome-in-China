# Figure 6 a-h: analysis, result preparation and plotting in one R file.


# =============================================================================
# Fig6_network.R
# =============================================================================
# Figure 6 a,b,e,f: surveillance-network RF and LMM analysis.

fig6_network_config <- function() {
  list(
    classes = list(
      "Meteorological factors" = c("prec_sum", "temp_mean"),
      "Air pollutions" = c("PM10", "PM2.5"),
      "Digital prescriptions" = c("Aminoglycoside_class", "Antihistamine_drugs", "NSAIDs",
        "Other_types_of_antibiotics", "Peptide_class", "Phenolic_compounds", "Quinolone_class",
        "Sulfonamide_class_and_synergists", "Tetracycline_class", "beta_Lactam_antibiotics",
        "non_NSAIDs", "Macrolide_class"),
      "Population migration" = c("baidu.emigration", "baidu.immigration")),
    lags = c("lag0", "lag1", "lag2", "lag3"),
    lag_days = c(lag0 = 0L, lag1 = 30L, lag2 = 60L, lag3 = 90L),
    outcome = "arg_abundance", ntree = 1000L, num.rep = 1000L, num.cores = 4L,
    seed = 123L, knn_k = 5L, zero_replacement = 1e-8,
    panel_predictors = c("baidu.immigration", "baidu.emigration", "Phenolic_compounds", "NSAIDs",
      "Aminoglycoside_class", "beta_Lactam_antibiotics", "Macrolide_class", "Tetracycline_class"),
    labels = c(prec_sum = "Precipitation", temp_mean = "Air temperature", PM10 = "PM10", "PM2.5" = "PM2.5",
      Aminoglycoside_class = "Aminoglycoside", Antihistamine_drugs = "Antihistamine", NSAIDs = "NSAIDs",
      Other_types_of_antibiotics = "Other antibiotics", Peptide_class = "Peptide", Phenolic_compounds = "Phenolic biocides",
      Quinolone_class = "Quinolone", Sulfonamide_class_and_synergists = "Sulfonamide", Tetracycline_class = "Tetracycline",
      beta_Lactam_antibiotics = "Beta-lactam", non_NSAIDs = "Non-NSAIDs", Macrolide_class = "MLS",
      baidu.emigration = "Emigration", baidu.immigration = "Immigration", city = "City", time = "Month", setting = "Setting")
  )
}

fig6_network_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages before this operation: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

fig6_network_columns <- function(data, required) {
  missing <- setdiff(required, names(data))
  if (length(missing)) stop("Required columns missing: ", paste(missing, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

fig6_network_minmax <- function(x) {
  if (!is.numeric(x)) stop("Min-max scaling requires numeric columns.", call. = FALSE)
  if (any(is.infinite(x))) stop("Infinite values must be resolved before min-max scaling.", call. = FALSE)
  if (all(is.na(x))) return(rep(NA_real_, length(x)))
  span <- range(x, na.rm = TRUE)
  if (span[1] == span[2]) return(ifelse(is.na(x), NA_real_, 0))
  (x - span[1]) / diff(span)
}

# daily: city, date, Index, value, already harmonized and aggregated to one row
# per city/date/Index. Missing days are ignored in window means, as in the source.
# Explicit boundaries: lag0 = sampling date; lag30/60/90 exclude sampling date.
fig6_network_lags <- function(samples, daily, config = fig6_network_config()) {
  fig6_network_columns(samples, c("city", "sample_date"))
  fig6_network_columns(daily, c("city", "date", "Index", "value"))
  samples <- as.data.frame(samples)
  daily <- as.data.frame(daily)
  samples$sample_date <- as.Date(samples$sample_date)
  daily$date <- as.Date(daily$date)
  if (anyNA(samples$sample_date) || anyNA(samples$city)) stop("Sample city/date must be complete.", call. = FALSE)
  if (anyNA(daily[c("city", "date", "Index")])) stop("Daily exposure keys must be complete.", call. = FALSE)
  if (anyDuplicated(daily[c("city", "date", "Index")])) stop("Aggregate duplicate city/date/Index exposure rows explicitly first.", call. = FALSE)
  if (!is.numeric(daily$value)) stop("Daily exposure value must be numeric.", call. = FALSE)
  for (index in unique(as.character(daily$Index))) {
    daily_index <- daily[daily$Index == index, , drop = FALSE]
    for (lag in config$lags) {
      days <- unname(config$lag_days[lag])
      samples[[paste0(index, "_", lag)]] <- vapply(seq_len(nrow(samples)), function(i) {
        delta <- as.numeric(samples$sample_date[i] - daily_index$date)
        keep <- daily_index$city == samples$city[i] & if (days == 0L) delta == 0 else delta >= 1 & delta <= days
        vals <- daily_index$value[keep]
        if (!length(vals) || all(is.na(vals))) NA_real_ else mean(vals, na.rm = TRUE)
      }, numeric(1))
    }
  }
  samples
}


fig6_network_prepare_lmm <- function(data, config = fig6_network_config()) {
  fig6_network_columns(data, c(config$outcome, "city", "site_nm", "season"))
  if (isTRUE(attr(data, "fig6_rf_imputed"))) stop("Network LMM must use the data before RF kNN imputation.", call. = FALSE)
  if (isTRUE(attr(data, "fig6_minmax"))) stop("Data have already been min-max scaled.", call. = FALSE)
  out <- as.data.frame(data)
  outcome <- out[[config$outcome]]
  if (!is.numeric(outcome)) stop("Outcome must be numeric.", call. = FALSE)
  out[[config$outcome]][!is.na(outcome) & outcome == 0] <- config$zero_replacement
  exposures <- unlist(lapply(unlist(config$classes, use.names = FALSE), function(v) paste0(v, "_", config$lags)), use.names = FALSE)
  scale_cols <- intersect(c(config$outcome, exposures), names(out))
  out[scale_cols] <- lapply(out[scale_cols], fig6_network_minmax)
  if ("sample_date" %in% names(out)) out$sample_date <- as.Date(out$sample_date)
  attr(out, "fig6_minmax") <- TRUE
  attr(out, "fig6_preprocessing") <- data.frame(
    Rows = nrow(out), Missing_outcomes = sum(is.na(outcome)), Zeros_replaced = sum(outcome == 0, na.rm = TRUE),
    Missing_outcome_rule = "excluded modelwise; never replaced or kNN imputed",
    Scaling = "outcome and included exposures independently min-max 0-1 before modelwise complete cases",
    stringsAsFactors = FALSE)
  out
}

fig6_network_model_spec <- function(data, class, lag, config = fig6_network_config()) {
  if (!class %in% names(config$classes) || !lag %in% config$lags) stop("Unknown network class or lag.", call. = FALSE)
  expected <- paste0(config$classes[[class]], "_", lag)
  predictors <- intersect(expected, names(data))
  if (!length(predictors)) stop("No predictors available for ", class, "/", lag, call. = FALSE)
  columns <- c(config$outcome, predictors, "season", "city", "site_nm")
  fig6_network_columns(data, columns)
  quote_name <- function(x) paste0("`", x, "`")
  formula_text <- paste(quote_name(config$outcome), "~", paste(quote_name(predictors), collapse = " + "),
                        "+ factor(season) + (1 | city) + (1 | site_nm)")
  list(predictors = predictors, missing_predictors = setdiff(expected, predictors),
       rows = which(stats::complete.cases(data[columns])), formula_text = formula_text,
       formula = stats::as.formula(formula_text))
}

# One BH family = one outcome x predictor class, spanning all predictors/lags.
fig6_network_adjust_lmm <- function(results) {
  fig6_network_columns(results, c("outcome", "class", "p", "coef"))
  results$fdr <- NA_real_
  families <- interaction(results$outcome, results$class, drop = TRUE, lex.order = TRUE)
  for (rows in split(seq_len(nrow(results)), families)) results$fdr[rows] <- stats::p.adjust(results$p[rows], method = "BH")
  significant <- !is.na(results$fdr) & results$fdr < .05
  results$P1 <- ifelse(significant, "*", "")
  results$Trend <- ifelse(significant & results$coef > 0, "Positive", ifelse(significant & results$coef < 0, "Negative", "No Significant"))
  results
}


fig6_network_fit_lmm <- function(prepared, outcome_label = "arg_abundance", config = fig6_network_config()) {
  fig6_network_require(c("lme4", "lmerTest"))
  if (!isTRUE(attr(prepared, "fig6_minmax"))) stop("Use fig6_network_prepare_lmm() on unimputed data first.", call. = FALSE)
  if (isTRUE(attr(prepared, "fig6_rf_imputed"))) stop("RF kNN data cannot be used for network LMM.", call. = FALSE)
  results <- diagnostics <- models <- list()
  for (class in names(config$classes)) for (lag in config$lags) {
    id <- paste(class, lag, sep = "__")
    spec <- tryCatch(fig6_network_model_spec(prepared, class, lag, config), error = identity)
    if (inherits(spec, "error")) {
      diagnostics[[id]] <- data.frame(Class = class, Lag = lag, N = 0L, Status = conditionMessage(spec), Singular = NA,
        Formula = NA_character_, Missing_predictors = NA_character_, Warnings = "")
      next
    }
    dat <- prepared[spec$rows, , drop = FALSE]
    warnings <- character()
    fit <- tryCatch(withCallingHandlers(lmerTest::lmer(spec$formula, data = dat, REML = TRUE, na.action = stats::na.fail),
      warning = function(w) { warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning") }), error = identity)
    failure <- inherits(fit, "error")
    diagnostics[[id]] <- data.frame(Class = class, Lag = lag, N = nrow(dat),
      Status = if (failure) conditionMessage(fit) else "ok", Singular = if (failure) NA else lme4::isSingular(fit),
      Formula = spec$formula_text, Missing_predictors = paste(spec$missing_predictors, collapse = ";"), Warnings = paste(unique(warnings), collapse = ";"))
    if (failure) next
    models[[id]] <- fit
    coefficients <- stats::coef(summary(fit))
    exps <- intersect(spec$predictors, rownames(coefficients))
    intervals <- suppressMessages(stats::confint(fit, parm = exps, method = "Wald"))
    r2 <- if (requireNamespace("performance", quietly = TRUE)) tryCatch(suppressWarnings(performance::r2_nakagawa(fit)), error = function(e) NULL) else NULL
    results[[id]] <- data.frame(outcome = outcome_label, class = class, lag = lag, exposure = exps,
      group = sub("_lag[0-3]$", "", exps), coef = coefficients[exps, "Estimate"], SE = coefficients[exps, "Std. Error"],
      ci_low = intervals[exps, 1], ci_high = intervals[exps, 2], p = coefficients[exps, "Pr(>|t|)"],
      R2_marginal = if (is.list(r2)) unname(r2$R2_marginal) else NA_real_,
      R2_conditional = if (is.list(r2)) unname(r2$R2_conditional) else NA_real_, N = nrow(dat), stringsAsFactors = FALSE)
  }
  result <- if (length(results)) fig6_network_adjust_lmm(do.call(rbind, results)) else data.frame()
  rownames(result) <- NULL
  list(results = result, diagnostics = do.call(rbind, diagnostics), models = models,
       preprocessing = attr(prepared, "fig6_preprocessing"))
}

fig6_network_prepare_rf <- function(prepared, config = fig6_network_config(), impute = TRUE) {
  if (!isTRUE(attr(prepared, "fig6_minmax"))) stop("Use pre-kNN, min-max prepared network data.", call. = FALSE)
  fig6_network_columns(prepared, c("city", "setting", "sample_date", config$outcome))
  keep <- !is.na(prepared[[config$outcome]])
  date <- as.Date(prepared$sample_date[keep])
  out <- data.frame(city = prepared$city[keep], setting = prepared$setting[keep],
                    time = format(date, "%Y-%m"), stringsAsFactors = FALSE)
  out[[config$outcome]] <- prepared[[config$outcome]][keep]
  lag_columns <- grep("_lag0$", names(prepared), value = TRUE)
  out[sub("_lag0$", "", lag_columns)] <- prepared[keep, lag_columns, drop = FALSE]
  if (anyNA(out[c("city", "setting", "time")])) stop("RF context fields must be complete; fix metadata before imputation.", call. = FALSE)
  before <- vapply(out, function(x) sum(is.na(x)), integer(1))
  if (impute && anyNA(out)) {
    fig6_network_require("VIM")
    # The original k=5/median rule is retained. Context fields have no missing
    # entries, so no categorical imputation or ambiguous base::mode is needed.
    out <- VIM::kNN(out, k = config$knn_k, numFun = stats::median, imp_var = FALSE)
  }
  out[c("city", "setting", "time")] <- lapply(out[c("city", "setting", "time")], factor)
  attr(out, "fig6_rf_imputed") <- impute && any(before > 0L)
  attr(out, "fig6_rf_preprocessing") <- data.frame(Predictor = names(out), Missing_before = before,
    Missing_after = vapply(out, function(x) sum(is.na(x)), integer(1)), stringsAsFactors = FALSE)
  attr(out, "fig6_missing_outcomes_excluded") <- sum(!keep)
  out
}

fig6_network_fit_rf <- function(prepared_rf, outcome_label = "arg_abundance", config = fig6_network_config()) {
  fig6_network_require(c("rfPermute", "randomForest"))
  if (!"num.rep" %in% names(formals(getS3method("rfPermute", "default", envir = asNamespace("rfPermute")))))
    stop("Installed rfPermute API lacks the expected num.rep argument.", call. = FALSE)
  if (anyNA(prepared_rf)) stop("RF inputs contain missing values; complete RF preprocessing first.", call. = FALSE)
  if (config$ntree != 1000L || config$num.rep != 1000L) stop("Confirmed network analysis uses 1000 trees and 1000 permutations.", call. = FALSE)
  set.seed(config$seed)
  fit <- rfPermute::rfPermute(stats::reformulate(setdiff(names(prepared_rf), config$outcome), response = config$outcome),
    data = prepared_rf, importance = TRUE, ntree = config$ntree, num.rep = config$num.rep, num.cores = config$num.cores)
  importance <- as.data.frame(randomForest::importance(fit, scale = TRUE), check.names = FALSE)
  importance$variable <- rownames(importance)
  list(model = fit, importance = importance, plot_data = fig6_network_rf_table(importance, outcome_label, config),
       metadata = data.frame(Outcome = outcome_label, N = nrow(prepared_rf), Trees = config$ntree,
         Permutations = config$num.rep, Seed = config$seed, Package_version = as.character(utils::packageVersion("rfPermute")),
         Significance = "nominal permutation P", stringsAsFactors = FALSE),
       preprocessing = attr(prepared_rf, "fig6_rf_preprocessing"))
}

fig6_network_label <- function(predictor, config = fig6_network_config()) {
  label <- unname(config$labels[predictor])
  label[is.na(label)] <- predictor[is.na(label)]
  label
}

fig6_network_class <- function(predictor, config = fig6_network_config()) {
  class <- rep("Other factors", length(predictor))
  for (name in names(config$classes)) class[predictor %in% config$classes[[name]]] <- name
  class[class == "Air pollutions"] <- "Air pollution"
  class[predictor %in% c("city", "setting", "time", "month", "year_month")] <- "Context factors"
  class
}

fig6_network_stars <- function(p) {
  ifelse(is.na(p), "", ifelse(p < .001, "***", ifelse(p < .01, "**", ifelse(p < .05, "*", ""))))
}


fig6_network_rf_table <- function(importance, outcome_label, config = fig6_network_config()) {
  importance <- as.data.frame(importance, check.names = FALSE)
  fig6_network_columns(importance, "%IncMSE")
  key_col <- intersect(c("Predictor", "Index", "variable"), names(importance))
  key <- if (length(key_col)) as.character(importance[[key_col[1]]]) else rownames(importance)
  # Older saved tables replace variable with the display label; invert known
  # labels/legacy aliases without changing any stored numerical values.
  aliases <- c(Phenol = "Phenolic_compounds", "Beta lactam antibiotics" = "beta_Lactam_antibiotics",
               "Non NSAIDs" = "non_NSAIDs", Emmigration = "baidu.emigration")
  label_key <- stats::setNames(names(config$labels), config$labels)
  mapping <- c(aliases, label_key)
  replace <- !key %in% names(config$labels) & key %in% names(mapping)
  key[replace] <- unname(mapping[key[replace]])
  if (anyNA(key) || any(!nzchar(key)) || anyDuplicated(key))
    stop("RF predictor keys must be present and unique; duplicate keys would distort class shares.", call. = FALSE)
  score <- importance[["%IncMSE"]]
  if (!is.numeric(score) || any(!is.finite(score))) stop("All RF importance scores must be finite numeric values.", call. = FALSE)
  if (sum(score) <= 0) stop("Cannot apply source RF percentage rule when signed importance sum is nonpositive.", call. = FALSE)
  p_col <- intersect(c("%IncMSE.pval", "P", "p"), names(importance))
  if (!length(p_col)) stop("Nominal permutation P values are required for network RF panels.", call. = FALSE)
  p <- importance[[p_col[1]]]
  if (any(!is.na(p) & (p < 0 | p > 1))) stop("Invalid RF P values.", call. = FALSE)
  raw_relative <- score / sum(score) * 100
  data.frame(Predictor = key, Label = fig6_network_label(key, config), Class = fig6_network_class(key, config),
    Importance = pmax(raw_relative, 0), Raw_importance = score, Relative_unclipped = raw_relative,
    P = p, Q = NA_real_, StarP = p, Stars = fig6_network_stars(p), Outcome = outcome_label, Display = TRUE,
    stringsAsFactors = FALSE, check.names = FALSE)
}

fig6_network_rf_classes <- function(table) {
  fig6_network_columns(table, c("Class", "Importance"))
  # Do not filter Display: class shares always include every model predictor.
  totals <- stats::aggregate(table$Importance, list(Class = table$Class), sum)
  names(totals)[2] <- "Importance"
  totals$Percent <- totals$Importance / sum(totals$Importance) * 100
  totals[order(-totals$Percent), , drop = FALSE]
}

# Preserve fitted coefficients, intervals, nominal P and stored BH values exactly.
# `panel_only=TRUE` selects the eight confirmed main-figure predictors after BH.
fig6_network_lmm_table <- function(results, outcome_label = NULL, config = fig6_network_config(), panel_only = FALSE) {
  fig6_network_columns(results, c("outcome", "class", "lag", "group", "coef", "ci_low", "ci_high", "p", "fdr"))
  if (any(!is.na(results$p) & (results$p < 0 | results$p > 1)) ||
      any(!is.na(results$fdr) & (results$fdr < 0 | results$fdr > 1)))
    stop("Saved LMM P and BH q values must fall within 0-1.", call. = FALSE)
  lag_labels <- c(lag0 = "Lag0", lag1 = "Lag30", lag2 = "Lag60", lag3 = "Lag90")
  lag <- unname(lag_labels[as.character(results$lag)])
  if (anyNA(lag)) stop("Unknown LMM lag labels; expected lag0-lag3.", call. = FALSE)
  out <- data.frame(Predictor = as.character(results$group), Label = fig6_network_label(as.character(results$group), config),
    Class = as.character(results$class), Lag = factor(lag, levels = lag_labels), Estimate = results$coef,
    SE = if ("SE" %in% names(results)) results$SE else NA_real_, CI_low = results$ci_low, CI_high = results$ci_high,
    P = results$p, Q = results$fdr, Outcome = if (is.null(outcome_label)) as.character(results$outcome) else outcome_label,
    StarQ = ifelse(!is.na(results$fdr) & results$fdr < .05, "*", ""),
    Display = results$group %in% config$panel_predictors, stringsAsFactors = FALSE)
  out$Class[out$Class == "Air pollutions"] <- "Air pollution"
  if (panel_only) {
    out <- out[out$Display, , drop = FALSE]
    out$Label <- factor(out$Label, levels = fig6_network_label(config$panel_predictors, config))
  }
  out
}

fig6_network_read_table <- function(path) {
  if (!file.exists(path)) stop("Result file does not exist: ", path, call. = FALSE)
  extension <- tolower(tools::file_ext(path))
  if (extension == "rds") {
    value <- readRDS(path)
    if (!is.data.frame(value)) stop("Result RDS must contain a table; extract model importance explicitly first.", call. = FALSE)
    return(value)
  }
  if (extension %in% c("tsv", "txt")) return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (extension == "csv") return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  stop("Result tables must be CSV, TSV/TXT or data-frame RDS.", call. = FALSE)
}

# Named paths assign panel outcome identities (e.g. c(Total=..., Tier_I=...)).
# No placeholder tables are created when previously fitted outputs are missing.
fig6_network_load_results <- function(rf_paths = character(), lmm_paths = character(), config = fig6_network_config()) {
  read_set <- function(paths, type) {
    if (!length(paths)) return(list(data = data.frame(), provenance = data.frame()))
    if (is.null(names(paths)) || any(!nzchar(names(paths)))) stop("Result paths must be named by outcome.", call. = FALSE)
    data <- lapply(seq_along(paths), function(i) {
      result <- fig6_network_read_table(paths[[i]])
      if (type == "RF") fig6_network_rf_table(result, names(paths)[i], config)
      else fig6_network_lmm_table(result, names(paths)[i], config)
    })
    list(data = do.call(rbind, data), provenance = data.frame(Analysis = type, Outcome = names(paths),
      Path = normalizePath(paths, winslash = "/", mustWork = TRUE), MD5 = unname(tools::md5sum(paths)),
      stringsAsFactors = FALSE))
  }
  rf <- read_set(rf_paths, "RF")
  lm <- read_set(lmm_paths, "LMM")
  list(rf = rf$data, lmm = lm$data, provenance = rbind(rf$provenance, lm$provenance))
}



fig6_hospital_assert <- function(ok, message) {
  if (!isTRUE(ok)) stop(message, call. = FALSE)
}

fig6_hospital_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  fig6_hospital_assert(!length(missing), paste(label, "lacks", paste(missing, collapse = ", ")))
}

fig6_hospital_cauchy <- function(p) {
  p <- as.numeric(p)
  fig6_hospital_assert(length(p) > 0L && all(is.finite(p)) && all(p >= 0 & p <= 1),
                       "Cauchy combination requires finite P values in [0,1] for every outcome")
  # Preserve the archived runner's endpoint clipping exactly.
  p <- pmin(pmax(p, 1e-15), 1 - 1e-15)
  statistic <- mean(tan((0.5 - p) * pi))
  pmin(pmax(0.5 - atan(statistic) / pi, 0), 1)
}

fig6_hospital_summarize <- function(components, predictor_ids, outcome_ids) {
  fig6_hospital_columns(components, c("source_id", "outcome_id", "%IncMSE_unscaled", "%IncMSE.pval"), "RF results")
  fig6_hospital_assert(length(predictor_ids) > 0 && !anyDuplicated(predictor_ids) &&
                         length(outcome_ids) > 0 && !anyDuplicated(outcome_ids), "Predictor and outcome IDs must be unique")
  fig6_hospital_assert(!anyNA(components[c("source_id", "outcome_id")]) &&
                         !anyDuplicated(components[c("source_id", "outcome_id")]), "Duplicate or missing RF predictor/outcome identifiers")
  fig6_hospital_assert(nrow(components) == length(predictor_ids) * length(outcome_ids) &&
                         setequal(components$source_id, predictor_ids) && setequal(components$outcome_id, outcome_ids),
                       "Expected complete predictor-by-outcome grid; no reduced outcome set is allowed")
  raw <- as.numeric(components[["%IncMSE_unscaled"]])
  p <- as.numeric(components[["%IncMSE.pval"]])
  fig6_hospital_assert(all(is.finite(raw)), "Every fitted importance must be finite")
  fig6_hospital_assert(all(is.finite(p)) && all(p >= 0 & p <= 1), "Every fitted P value must be finite and in [0,1]")
  rows <- lapply(predictor_ids, function(id) {
    at <- which(components$source_id == id)
    fig6_hospital_assert(setequal(components$outcome_id[at], outcome_ids), paste("Incomplete outcome grid for", id))
    data.frame(Predictor = id, RawImportanceMean = mean(raw[at]),
               PositiveImportanceMean = mean(pmax(raw[at], 0)),
               P = if (length(outcome_ids) == 1L) p[at] else fig6_hospital_cauchy(p[at]),
               OutcomeN = length(outcome_ids), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  denominator <- sum(out$PositiveImportanceMean)
  fig6_hospital_assert(denominator > 0, "Positive importance sum must exceed zero")
  out$Importance <- 100 * out$PositiveImportanceMean / denominator
  out$Q <- stats::p.adjust(out$P, method = "BH", n = length(predictor_ids))
  out$StarP <- out$P
  out
}

fig6_hospital_class_order <- function() {
  c("Context factors", "Hospital profile", "Antibiotic consumption", "Hospital continuity")
}

fig6_hospital_classes <- function(ids) {
  result <- rep(NA_character_, length(ids))
  result[ids %in% c("CITY", "QUARTER")] <- "Context factors"
  result[ids %in% c("HOSPITAL", "HOSP_type")] <- "Hospital profile"
  result[grepl("^DRUG__", ids) | ids == "HOSP_mean_quarterly_amount"] <- "Antibiotic consumption"
  # Retain the original figure's class name. These three predictors are
  # hospital-level means, not longitudinal within-hospital changes.
  result[ids %in% c("HOSP_mean_patient_days", "HOSP_mean_los", "HOSP_mean_discharges")] <- "Hospital continuity"
  fig6_hospital_assert(!anyNA(result), paste("Unknown predictor class:", paste(ids[is.na(result)], collapse = ", ")))
  result
}

fig6_hospital_add_metadata <- function(x, dictionary) {
  metadata <- dictionary[match(x$Predictor, dictionary$source_id), , drop = FALSE]
  x$Label <- as.character(metadata$source_label_English)
  x$Class <- fig6_hospital_classes(x$Predictor)
  x$PredictorN <- 1L
  x[, c("Predictor", "Label", "Class", "Importance", "P", "Q", "StarP",
        "RawImportanceMean", "PositiveImportanceMean", "OutcomeN", "PredictorN")]
}

fig6_hospital_class_summary <- function(x) {
  do.call(rbind, lapply(fig6_hospital_class_order(), function(class_name) {
    at <- x$Class == class_name
    data.frame(Class = class_name, Importance = sum(x$Importance[at]),
               PredictorN = as.integer(sum(x$PredictorN[at])), stringsAsFactors = FALSE)
  }))
}

fig6_hospital_display_rows <- function(x, display_drug_ids) {
  drug <- grepl("^DRUG__", x$Predictor)
  retained <- !drug | x$Predictor %in% display_drug_ids
  other <- x[drug & !retained, , drop = FALSE]
  out <- x[retained, , drop = FALSE]
  if (nrow(other)) {
    aggregate <- data.frame(Predictor = "OTHER_ANTIBIOTICS", Label = "Other antibiotics",
      Class = "Antibiotic consumption", Importance = sum(other$Importance),
      P = NA_real_, Q = NA_real_, StarP = NA_real_,
      RawImportanceMean = sum(other$RawImportanceMean),
      PositiveImportanceMean = sum(other$PositiveImportanceMean),
      OutcomeN = unique(other$OutcomeN), PredictorN = as.integer(nrow(other)), stringsAsFactors = FALSE)
    out <- rbind(out, aggregate)
  }
  out <- out[order(match(out$Class, fig6_hospital_class_order()), -out$Importance, out$Predictor), , drop = FALSE]
  rownames(out) <- NULL
  fig6_hospital_assert(abs(sum(out$Importance) - 100) < 1e-9, "Display grouping changed the all-predictor denominator")
  out
}

fig6_prepare_hospital_rf <- function(total, group, dictionary, catalog) {
  fig6_hospital_columns(dictionary, c("source_id", "source_label_English"), "Predictor dictionary")
  fig6_hospital_columns(catalog, c("outcome_id", "outcome_level"), "Outcome catalogue")
  predictor_ids <- as.character(dictionary$source_id)
  expected_predictors <- c(sprintf("DRUG__%03d", 1:97), "HOSP_mean_patient_days", "HOSP_mean_los",
                           "HOSP_mean_discharges", "HOSP_mean_quarterly_amount", "QUARTER", "HOSPITAL", "CITY", "HOSP_type")
  fig6_hospital_assert(identical(predictor_ids, expected_predictors), "Hospital analysis requires the archived ordered 97 + 8 predictors")
  outcome_ids <- as.character(catalog$outcome_id)
  counts <- table(factor(catalog$outcome_level, levels = c("type", "clinical", "total")))
  fig6_hospital_assert(length(outcome_ids) == 22L && !anyDuplicated(outcome_ids) &&
                         identical(as.integer(counts), c(12L, 9L, 1L)) &&
                         "TOTAL__all_ARG" %in% outcome_ids, "Expected 12 ARG types, 9 clinical groups and total ARG")
  c_rows <- fig6_hospital_add_metadata(fig6_hospital_summarize(total, predictor_ids, "TOTAL__all_ARG"), dictionary)
  d_rows <- fig6_hospital_add_metadata(fig6_hospital_summarize(group, predictor_ids, outcome_ids), dictionary)
  display_drug_ids <- sort(unique(c(c_rows$Predictor[grepl("^DRUG__", c_rows$Predictor) & c_rows$P < 0.05],
                                    d_rows$Predictor[grepl("^DRUG__", d_rows$Predictor) & d_rows$P < 0.05])))
  components <- group
  components$importance_nonnegative <- pmax(as.numeric(components[["%IncMSE_unscaled"]]), 0)
  # The component file preserves all nominal and BH P values from the archive.
  list(c = c_rows, d = d_rows,
       class_c = fig6_hospital_class_summary(c_rows), class_d = fig6_hospital_class_summary(d_rows),
       display_c = fig6_hospital_display_rows(c_rows, display_drug_ids),
       display_d = fig6_hospital_display_rows(d_rows, display_drug_ids),
       display_drug_ids = display_drug_ids, components = components,
       outcome_catalog = catalog, predictor_dictionary = dictionary)
}

fig6_hospital_read_tsv <- function(path) {
  fig6_hospital_assert(file.exists(path), paste("Missing hospital RF input:", path))
  out <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, encoding = "UTF-8")
  names(out) <- sub("^\ufeff", "", names(out), useBytes = TRUE)
  out
}

fig6_read_hospital_rf <- function(data_dir = "data/hospital_rf") {
  files <- c(total = "total_importance.tsv", group = "group22_importance.tsv",
             dictionary = "predictor_dictionary.tsv", catalog = "outcome_catalog.tsv",
             original_summary = "group22_original_summary.tsv",
             total_settings = "total_run_settings.tsv", group_settings = "group22_run_settings.tsv",
             total_run = "total_run_summary.tsv", group_run = "group22_run_summary.tsv")
  out <- lapply(files, function(filename) fig6_hospital_read_tsv(file.path(data_dir, filename)))
  for (settings in list(out$total_settings, out$group_settings)) {
    fig6_hospital_assert(nrow(settings) == 1L && settings$trees == 800L &&
                           settings$num.rep == 999L && !settings$quick_test,
                         "Only the completed formal 800-tree / 999-permutation archive is accepted")
  }
  fig6_hospital_assert(out$total_run$pending_tasks == 0L && out$total_run$completed_outcomes == 1L &&
                         out$group_run$pending_tasks == 0L && out$group_run$successful_outcomes == 22L,
                       "Formal hospital RF runs are incomplete")
  out
}

fig6_hospital_model_inputs <- function(data_dir = "data/hospital_rf") {
  read <- function(filename) fig6_hospital_read_tsv(file.path(data_dir, filename))
  keys <- read("keys.tsv")
  outcomes <- read("outcomes_log10.tsv")
  drugs <- read("drug_exposures.tsv")
  context <- read("hospital_context.tsv")
  dictionary <- read("predictor_dictionary.tsv")
  catalog <- read("outcome_catalog.tsv")
  make_key <- function(x) paste(trimws(as.character(x$hospital_code)), as.integer(x$quarter), sep = "__")
  key <- make_key(keys)
  fig6_hospital_assert(nrow(keys) == 43L && !anyDuplicated(key) &&
                         length(unique(keys$hospital_code)) == 15L, "Expected 43 unique hospital-quarter rows from 15 hospitals")
  align <- function(x) {
    x_key <- make_key(x)
    fig6_hospital_assert(!anyDuplicated(x_key) && setequal(x_key, key), "Hospital-quarter input keys disagree")
    x[match(key, x_key), , drop = FALSE]
  }
  drugs <- align(drugs); context <- align(context); outcomes <- align(outcomes)
  ids <- as.character(dictionary$source_id)
  numeric_ids <- c(sprintf("DRUG__%03d", 1:97), "HOSP_mean_patient_days", "HOSP_mean_los",
                   "HOSP_mean_discharges", "HOSP_mean_quarterly_amount")
  context_ids <- ids[!grepl("^DRUG__", ids)]
  raw <- cbind(drugs[, sprintf("DRUG__%03d", 1:97), drop = FALSE], context[, context_ids, drop = FALSE])
  fig6_hospital_assert(identical(names(raw), ids) && ncol(raw) == 105L && !anyNA(raw), "Expected complete ordered 43 x 105 predictor matrix")
  X <- raw
  for (id in numeric_ids) {
    values <- as.numeric(X[[id]])
    fig6_hospital_assert(all(is.finite(values) & values >= 0), paste("Invalid continuous predictor:", id))
    X[[id]] <- log1p(values)
  }
  for (id in setdiff(ids, numeric_ids)) X[[id]] <- factor(as.character(X[[id]]), levels = sort(unique(as.character(X[[id]]))))
  for (id in context_ids[grepl("^HOSP_mean_", context_ids)]) {
    fig6_hospital_assert(all(vapply(split(raw[[id]], keys$hospital_code), function(x) length(unique(x)) == 1L, logical(1))),
                         paste(id, "must be a fixed hospital-level mean"))
  }
  Y <- outcomes[, as.character(catalog$outcome_id), drop = FALSE]
  fig6_hospital_assert(ncol(Y) == 22L && all(vapply(Y, function(y) all(is.finite(y)), logical(1))), "Expected 22 finite archived log10 ARG outcomes")
  # Outcomes are already log10(abundance + 1e-8); never transform them twice.
  list(keys = keys, raw_predictors = raw, X = X, Y = Y, catalog = catalog, dictionary = dictionary)
}

fig6_export_hospital_rf <- function(hospital, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(Fig6c_all_105 = hospital$c, Fig6d_all_105 = hospital$d,
                 Fig6c_display_24 = hospital$display_c, Fig6d_display_24 = hospital$display_d,
                 Fig6c_class_shares = hospital$class_c, Fig6d_class_shares = hospital$class_d,
                 Fig6d_all_2310_components = hospital$components)
  for (name in names(tables)) utils::write.table(tables[[name]], file.path(output_dir, paste0(name, ".tsv")),
    sep = "\t", quote = TRUE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(tables)
}

fig6_hospital_fit_seed <- function(dataset, outcome_id) {
  base_seed <- switch(dataset, TotalARG2646 = 20260812L, group22 = 20260813L,
                      stop("Unknown hospital dataset", call. = FALSE))
  raw <- utf8ToInt(paste(dataset, outcome_id, sep = "::"))
  as.integer((as.double(base_seed) + sum(raw * seq_along(raw)) * 1009) %% 2147483000 + 1)
}

fig6_fit_hospital_rf <- function(model_inputs, dry_run = TRUE) {
  if (!exists("fig6_prepare_hospital_rf", mode = "function"))
    stop("Source Fig6_hospital_rf.R before fitting hospital models", call. = FALSE)
  check <- function(ok, message) if (!isTRUE(ok)) stop(message, call. = FALSE)
  check(is.logical(dry_run) && length(dry_run) == 1L && !is.na(dry_run), "dry_run must be TRUE or FALSE")
  check(all(c("X", "Y", "keys", "catalog", "dictionary") %in% names(model_inputs)),
        "Use the complete output of fig6_hospital_model_inputs()")
  X <- model_inputs$X; Y <- model_inputs$Y
  dictionary <- model_inputs$dictionary; catalog <- model_inputs$catalog
  ids <- c(sprintf("DRUG__%03d", 1:97), "HOSP_mean_patient_days", "HOSP_mean_los",
           "HOSP_mean_discharges", "HOSP_mean_quarterly_amount", "QUARTER", "HOSPITAL", "CITY", "HOSP_type")
  check(is.data.frame(X) && nrow(X) == 43L && identical(names(X), ids) &&
          identical(as.character(dictionary$source_id), ids), "Require 43 rows and the ordered 105 hospital predictors")
  check(all(c("source_label_English", "predictor_family", "drug_broad_class") %in% names(dictionary)),
        "Predictor dictionary metadata is incomplete")
  check(nrow(model_inputs$keys) == 43L && length(unique(model_inputs$keys$hospital_code)) == 15L,
        "Require 43 hospital-quarter observations from 15 hospitals")
  outcome_ids <- as.character(catalog$outcome_id)
  counts <- table(factor(catalog$outcome_level, levels = c("type", "clinical", "total")))
  check(is.data.frame(Y) && nrow(Y) == 43L && ncol(Y) == 22L &&
          identical(names(Y), outcome_ids) && !anyDuplicated(outcome_ids) &&
          "TOTAL__all_ARG" %in% outcome_ids && identical(as.integer(counts), c(12L, 9L, 1L)),
        "Require all 12 ARG types, nine clinical groups and total ARG in catalogue order")
  usable <- function(z) !anyNA(z) && length(unique(z)) > 1L &&
    (is.factor(z) || (is.numeric(z) && all(is.finite(z))))
  check(all(vapply(X, usable, logical(1))), "All 105 predictors must be finite, nonmissing and nonconstant; none may be silently dropped")
  check(all(vapply(X[tail(ids, 4)], is.factor, logical(1))), "Quarter, hospital, city and hospital type must be factors")
  check(all(vapply(Y, function(y) is.numeric(y) && all(is.finite(y)) && length(unique(y)) > 1L, logical(1))),
        "All 22 log10 outcomes must be finite and nonconstant")

  plan <- data.frame(dataset = c("TotalARG2646", rep("group22", 22L)),
                     outcome_id = c("TOTAL__all_ARG", outcome_ids),
                     trees = 800L, num.rep = 999L, num.cores = 1L,
                     predictor_n = 105L, observation_n = 43L, stringsAsFactors = FALSE)
  plan$model_seed <- vapply(seq_len(nrow(plan)), function(i)
    fig6_hospital_fit_seed(plan$dataset[i], plan$outcome_id[i]), integer(1))
  if (dry_run) return(plan)
  for (package in c("rfPermute", "randomForest"))
    check(requireNamespace(package, quietly = TRUE), paste("Install required R package:", package))

  fit_one <- function(task) {
    set.seed(task$model_seed)
    fit <- rfPermute::rfPermute(x = X, y = as.numeric(Y[[task$outcome_id]]),
      importance = TRUE, ntree = 800L, num.rep = 999L, num.cores = 1L)
    raw <- randomForest::importance(fit, scale = FALSE)
    scaled <- randomForest::importance(fit, scale = TRUE)
    row_ids <- rownames(raw)
    check(length(row_ids) == 105L && setequal(row_ids, ids), "Fitted model did not retain all 105 predictors")
    get_column <- function(matrix, column) {
      check(column %in% colnames(matrix), paste("Missing rfPermute result column:", column))
      as.numeric(matrix[match(row_ids, rownames(matrix)), column])
    }
    metadata <- dictionary[match(row_ids, ids), , drop = FALSE]
    result <- data.frame(dataset = task$dataset, outcome_id = task$outcome_id,
      source_id = row_ids, source_label_English = metadata$source_label_English,
      predictor_family = metadata$predictor_family, drug_broad_class = metadata$drug_broad_class,
      `%IncMSE_unscaled` = get_column(raw, "%IncMSE"),
      IncNodePurity_unscaled = get_column(raw, "IncNodePurity"),
      `%IncMSE_scaled` = get_column(scaled, "%IncMSE"),
      IncNodePurity_scaled = get_column(scaled, "IncNodePurity"),
      `%IncMSE.pval` = get_column(raw, "%IncMSE.pval"),
      IncNodePurity.pval = get_column(raw, "IncNodePurity.pval"),
      num.rep = 999L, trees = 800L, model_seed = task$model_seed,
      check.names = FALSE, stringsAsFactors = FALSE)
    result$`%IncMSE.q_BH_within_outcome` <- stats::p.adjust(result$`%IncMSE.pval`, method = "BH", n = 105L)
    result$IncNodePurity.q_BH_within_outcome <- stats::p.adjust(result$IncNodePurity.pval, method = "BH", n = 105L)
    result
  }
  fitted <- lapply(seq_len(nrow(plan)), function(i) fit_one(plan[i, , drop = FALSE]))
  total <- fitted[[1L]]
  group <- do.call(rbind, fitted[-1L])
  check(nrow(total) == 105L && nrow(group) == 2310L, "Incomplete fitted hospital output")
  list(total = total, group = group, dictionary = dictionary, catalog = catalog, plan = plan,
       prepared = fig6_prepare_hospital_rf(total, group, dictionary, catalog),
       package_versions = vapply(c("rfPermute", "randomForest"), function(x)
         as.character(utils::packageVersion(x)), character(1)))
}

fig6_hospital_key <- function(x) {
  required <- c("dataset", "outcome_id", "predictor_id", "profile")
  if (!all(required %in% names(x))) stop("Missing hospital LMM key columns", call. = FALSE)
  if (anyNA(x[required])) stop("Hospital LMM keys contain missing values", call. = FALSE)
  do.call(paste, c(x[required], sep = "|"))
}

fig6_hospital_validate_family <- function(d, planned_n) {
  d <- as.data.frame(d, stringsAsFactors = FALSE)
  required <- c("status", "beta_within", "se_within", "p_within", "ci_low", "ci_high",
                "fdr_global_profile_dataset", "n", "hospital_n")
  if (!all(required %in% names(d))) stop("Incomplete hospital LMM result schema", call. = FALSE)
  if (nrow(d) != planned_n) stop("Full planned family is required before BH adjustment", call. = FALSE)
  if (anyDuplicated(fig6_hospital_key(d))) stop("Duplicate hospital LMM result key", call. = FALSE)
  if (length(unique(paste(d$dataset,d$profile))) != 1L) stop("BH input mixes families", call. = FALSE)
  ok <- d$status == "success" & is.finite(d$p_within)
  if (!any(ok)) stop("No estimable models in family", call. = FALSE)
  if (any(d$n[ok] != 43L | d$hospital_n[ok] != 15L)) stop("Expected 43 observations and 15 hospitals", call. = FALSE)
  if (any(!is.finite(d$beta_within[ok]) | !is.finite(d$se_within[ok]))) stop("Non-finite successful estimate", call. = FALSE)
  if (any(d$p_within[ok] < 0 | d$p_within[ok] > 1)) stop("Invalid P value", call. = FALSE)
  q <- rep(NA_real_,nrow(d)); q[ok] <- p.adjust(d$p_within[ok],method="BH",n=planned_n)
  if (any(!is.finite(d$fdr_global_profile_dataset[ok])) ||
      max(abs(q[ok]-d$fdr_global_profile_dataset[ok])) > 1e-10)
    stop("Saved global FDR differs from BH over the full planned family", call. = FALSE)
  if (any(!is.finite(d$ci_low[ok]) | !is.finite(d$ci_high[ok])) ||
      max(abs(d$ci_low[ok]-(d$beta_within[ok]-1.96*d$se_within[ok])),
          abs(d$ci_high[ok]-(d$beta_within[ok]+1.96*d$se_within[ok]))) > 1e-10)
    stop("Saved confidence intervals are not beta +/- 1.96 SE", call. = FALSE)
  d$q <- q
  d
}

fig6_hospital_read_results <- function(full_results_dir) {
  if (!requireNamespace("data.table",quietly=TRUE)) stop("Install data.table",call.=FALSE)
  families <- expand.grid(dataset=c("full2646","group22"),profile=c("oral","injection","combined"),stringsAsFactors=FALSE)
  keep <- c("outcome_id","predictor_id","status","reason","beta_within","se_within",
            "p_within","ci_low","ci_high","fdr_global_profile_dataset","n","hospital_n",
            "singular","convergence_message","drug_name_English","outcome_label_English",
            "aggregation","observed_routes")
  results <- audit <- vector("list",nrow(families))
  for (i in seq_len(nrow(families))) {
    ds <- families$dataset[i]; profile <- families$profile[i]
    file <- file.path(full_results_dir,paste0("LMM_",ds,"_",profile,"_all.csv"))
    if (!file.exists(file)) stop("Missing full-family results: ",file,call.=FALSE)
    d <- as.data.frame(data.table::fread(file,select=keep,encoding="UTF-8",na.strings=c("","NA")))
    d$dataset <- ds; d$profile <- profile
    planned_n <- c(full2646=113778L,group22=946L)[[ds]]
    d <- fig6_hospital_validate_family(d,planned_n)
    if (any(!is.na(d$aggregation) & !grepl("log10\\(x\\+1e-08\\)",d$aggregation)))
      stop("Outcome transformation differs from archived log10(x+1e-08)",call.=FALSE)
    ok <- d$status=="success" & is.finite(d$p_within)
    audit[[i]] <- data.frame(dataset=ds,profile=profile,planned_models=planned_n,
      successful_models=sum(ok),not_estimable=sum(!ok),singular_successful=sum(d$singular[ok],na.rm=TRUE),
      bh_max_abs_error=max(abs(d$q[ok]-d$fdr_global_profile_dataset[ok])),
      ci_max_abs_error=max(abs(d$ci_low[ok]-d$beta_within[ok]+1.96*d$se_within[ok]),
                           abs(d$ci_high[ok]-d$beta_within[ok]-1.96*d$se_within[ok])),
      file_md5=unname(tools::md5sum(file)))
    results[[i]] <- d
  }
  list(results=do.call(rbind,results),audit=do.call(rbind,audit))
}

fig6_hospital_select <- function(full, selection) {
  full <- as.data.frame(full,stringsAsFactors=FALSE)
  selection <- as.data.frame(selection,stringsAsFactors=FALSE)
  key <- fig6_hospital_key(full); wanted <- fig6_hospital_key(selection)
  if (anyDuplicated(key)) stop("Duplicate hospital LMM result key",call.=FALSE)
  if (anyDuplicated(wanted)) stop("Duplicate selection key",call.=FALSE)
  matched <- match(wanted,key)
  if (anyNA(matched)) stop("Selected hospital LMM keys absent from full family",call.=FALSE)
  d <- full[matched,,drop=FALSE]
  if (any(d$status!="success" | !is.finite(d$beta_within))) stop("Selection includes non-estimable model",call.=FALSE)
  extra <- setdiff(names(selection),names(d))
  d <- cbind(selection[extra],d)
  rownames(d) <- NULL
  d
}

fig6_hospital_plots <- function(d, family="Arial") {
  for (pkg in c("ggplot2","cowplot","ragg")) if(!requireNamespace(pkg,quietly=TRUE)) stop("Install ",pkg,call.=FALSE)
  # Use the platform text renderer for measuring labels; restore the caller's
  # cowplot device afterward. The default PostScript device cannot measure Arial.
  previous_device <- cowplot::set_null_device("agg")
  on.exit(cowplot::set_null_device(previous_device),add=TRUE)
  # Retain the final historical layout, category order, route symbols and colors.
  d$beta <- d$beta_within; d$low <- d$ci_low; d$high <- d$ci_high
  d$route <- factor(d$profile,levels=c("oral","injection","combined"),labels=c("Oral","Injection","Combined"))
  offsets <- c(Oral=-.22,Injection=0,Combined=.22)
  d$x <- d$order_key+ifelse(d$section %in% c("Combined","h"),offsets[as.character(d$route)],0)
  d$star_y <- d$high+pmax(pmax(d$high,0)-pmin(d$low,0),.05)*.06
  d$star <- ifelse(d$q<.05,"*","")
  d$point_colour <- ifelse(d$q<.05,ifelse(d$beta<0,"#4E79A7","#B24C3D"),"#7A7A7A")
  panel <- function(z,ncat,limits,breaks,show_y=TRUE) {
    ggplot2::ggplot(z,ggplot2::aes(x,beta))+
      ggplot2::geom_hline(yintercept=0,linetype=2,linewidth=.38,colour="grey35")+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=low,ymax=high,colour=point_colour),width=.075,linewidth=.55)+
      ggplot2::geom_point(ggplot2::aes(shape=route,colour=point_colour),fill="white",size=2.55,stroke=.82)+
      ggplot2::geom_text(ggplot2::aes(y=star_y,label=star),family=family,size=3.75,fontface="bold",colour="black",vjust=0)+
      ggplot2::scale_shape_manual(values=c(Oral=24,Injection=21,Combined=22),drop=FALSE)+
      ggplot2::scale_colour_identity()+
      ggplot2::scale_x_continuous(limits=c(.45,ncat+.55),breaks=seq_len(ncat),expand=ggplot2::expansion(mult=0))+
      ggplot2::scale_y_continuous(limits=limits,breaks=breaks,labels=function(x)sprintf("%.1f",x),expand=ggplot2::expansion(mult=0))+
      ggplot2::labs(x=NULL,y=NULL)+ggplot2::theme_classic(base_family=family,base_size=7)+
      ggplot2::theme(panel.background=ggplot2::element_rect(fill=NA,colour=NA),
        panel.border=ggplot2::element_rect(colour="black",fill=NA,linewidth=.39),axis.line=ggplot2::element_blank(),
        axis.ticks=ggplot2::element_line(linewidth=.34,colour="black"),axis.ticks.length=grid::unit(1.6,"pt"),
        axis.text.x=ggplot2::element_blank(),axis.text.y=if(show_y)ggplot2::element_text(size=7.2,colour="black") else ggplot2::element_blank(),
        axis.ticks.y=if(show_y)ggplot2::element_line(linewidth=.34) else ggplot2::element_blank(),
        legend.position="none",plot.margin=ggplot2::margin(0,0,0,0),plot.background=ggplot2::element_rect(fill=NA,colour=NA))
  }
  labels <- function(z,ncat,font_pt) {
    z <- unique(z[c("order_key","outcome","drug")]);z<-z[order(z$order_key),]
    z$drug_display <- gsub("Amoxicillin/clavulanic acid","Amoxicillin/\nclavulanic acid",z$drug,fixed=TRUE)
    z$drug_display <- gsub("Benzathine benzylpenicillin","Benzathine benzyl-\npenicillin",z$drug_display,fixed=TRUE)
    z$drug_display <- gsub("Cefoperazone/sulbactam","Cefoperazone/\nsulbactam",z$drug_display,fixed=TRUE)
    ggplot2::ggplot(z,ggplot2::aes(order_key))+
      ggplot2::geom_text(ggplot2::aes(y=1.02,label=outcome),family=family,size=font_pt/(72.27/25.4),angle=30,hjust=1,vjust=.5,colour="#1A1A1A",lineheight=.86)+
      ggplot2::geom_text(ggplot2::aes(y=.70,label=drug_display),family=family,size=font_pt/(72.27/25.4),angle=30,hjust=1,vjust=.5,colour="#B044A8",lineheight=.84)+
      ggplot2::scale_x_continuous(limits=c(.45,ncat+.55),expand=ggplot2::expansion(mult=0))+
      ggplot2::coord_cartesian(ylim=c(0,1),clip="off")+ggplot2::theme_void()+
      ggplot2::theme(plot.margin=ggplot2::margin(0,0,0,0),plot.background=ggplot2::element_rect(fill=NA,colour=NA))
  }
  go<-d[d$panel=="g"&d$section=="Oral",];gc<-d[d$panel=="g"&d$section=="Combined",]
  gi<-d[d$panel=="g"&d$section=="Injection",];h<-d[d$panel=="h",]
  if (!identical(c(length(unique(go$order_key)),length(unique(gc$order_key)),length(unique(gi$order_key)),length(unique(h$order_key))),c(4L,16L,19L,12L)))
    stop("Selection no longer matches the archived 6g/6h layout",call.=FALSE)
  # Include the full star glyph above the largest CI; the historical upper
  # limit put the QnrD1 significance marker against the frame.
  gl<-c(-.60,max(1.08,d$star_y[d$panel=="g"&d$star=="*"]+.09));gb<-c(-.5,0,.5,1)
  g <- cowplot::ggdraw()+
    cowplot::draw_label("Oral",x=.130,y=.977,fontfamily=family,fontface="bold",size=9.2)+
    cowplot::draw_label("Oral + Injection",x=.617,y=.977,fontfamily=family,fontface="bold",size=9.2)+
    cowplot::draw_plot(panel(go,4,gl,gb),x=.0495,y=.6335,width=.1613,height=.3134)+
    cowplot::draw_plot(panel(gc,16,gl,gb,FALSE),x=.2359,y=.6335,width=.7626,height=.3134)+
    cowplot::draw_plot(labels(go,4,5.95),x=.0495,y=.4700,width=.1613,height=.1600)+
    cowplot::draw_plot(labels(gc,16,4.15),x=.2359,y=.4700,width=.7626,height=.1600)+
    cowplot::draw_plot(panel(gi,19,gl,gb),x=.0514,y=.1421,width=.9485,height=.3261)+
    cowplot::draw_plot(labels(gi,19,4.55),x=.0514,y=0,width=.9485,height=.1380)+
    cowplot::draw_label("Injection",x=.527,y=.484,fontfamily=family,fontface="bold",size=8.7)
  hp <- cowplot::ggdraw()+
    cowplot::draw_plot(panel(h,12,c(-.15,.85),c(0,.2,.4,.6,.8)),x=.0491,y=.2717,width=.9509,height=.7083)+
    cowplot::draw_plot(labels(h,12,5.25),x=.0491,y=0,width=.9509,height=.2650)
  list(g=g,h=hp,data=d)
}

fig6_hospital_lmm <- function(full_results_dir, selection_file, out_dir,
                              supplementary_selection_file=NULL, family="Arial") {
  loaded<-fig6_hospital_read_results(full_results_dir)
  selection<-read.delim(selection_file,check.names=FALSE,stringsAsFactors=FALSE)
  selected<-fig6_hospital_select(loaded$results,selection)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  write.csv(loaded$audit,file.path(out_dir,"hospital_LMM_full_family_audit.csv"),row.names=FALSE)
  plots<-fig6_hospital_plots(selected,family)
  write.csv(plots$data,file.path(out_dir,"Fig6gh_displayed_estimates.csv"),row.names=FALSE)
  if(!is.null(supplementary_selection_file)){
    supplement<-fig6_hospital_select(loaded$results,read.delim(supplementary_selection_file,check.names=FALSE,stringsAsFactors=FALSE))
    write.csv(supplement,file.path(out_dir,"S7_5_selected_estimates.csv"),row.names=FALSE)
  }
  for (pn in c("g","h")) {
    size<-if(pn=="g")c(491,234.5) else c(490.25,141.7)
    stem<-file.path(out_dir,paste0("Fig6",pn))
    ggplot2::ggsave(paste0(stem,".png"),plots[[pn]],width=size[1]/72,height=size[2]/72,dpi=450,device=ragg::agg_png,bg="white")
    ggplot2::ggsave(paste0(stem,".pdf"),plots[[pn]],width=size[1]/72,height=size[2]/72,device=grDevices::cairo_pdf,bg="white")
    if(requireNamespace("svglite",quietly=TRUE))ggplot2::ggsave(paste0(stem,".svg"),plots[[pn]],width=size[1]/72,height=size[2]/72,device=svglite::svglite,bg="transparent")
  }
  invisible(list(plots=plots[c("g","h")],data=plots$data,audit=loaded$audit))
}

fig6_hospital_fit_keys <- function(d,label) {
  keys<-c("hospital_code","quarter")
  if(!all(keys %in% names(d)))stop(label,": missing hospital_code/quarter",call.=FALSE)
  if(anyNA(d[keys]))stop(label,": missing key",call.=FALSE)
  key<-paste(as.character(d$hospital_code),as.character(d$quarter),sep="|")
  if(anyDuplicated(key))stop(label,": duplicate hospital-quarter key",call.=FALSE)
  key
}

fig6_hospital_fit_z <- function(x) {
  x<-as.numeric(x);s<-stats::sd(x,na.rm=TRUE)
  if(!is.finite(s)||s<=0)return(rep(0,length(x)))
  (x-mean(x,na.rm=TRUE))/s
}

fig6_hospital_model_frame <- function(outcomes,exposure,metadata,outcome_id,predictor_id,
                                      outcome_scale=c("linear","log10"),
                                      expected_n=43L,expected_hospitals=15L) {
  outcome_scale<-match.arg(outcome_scale)
  outcomes<-as.data.frame(outcomes);exposure<-as.data.frame(exposure);metadata<-as.data.frame(metadata)
  km<-fig6_hospital_fit_keys(metadata,"metadata")
  ky<-fig6_hospital_fit_keys(outcomes,"outcomes");kx<-fig6_hospital_fit_keys(exposure,"exposure")
  if(length(km)!=expected_n||length(unique(metadata$hospital_code))!=expected_hospitals)
    stop("Unexpected hospital-quarter or hospital count",call.=FALSE)
  if(!setequal(km,ky)||!setequal(km,kx))stop("Input hospital-quarter keys do not match",call.=FALSE)
  if(!outcome_id %in% names(outcomes)||!predictor_id %in% names(exposure)||!"HOSP_patient_days" %in% names(metadata))
    stop("Missing requested outcome, exposure, or patient-days column",call.=FALSE)
  y<-outcomes[[outcome_id]][match(km,ky)]
  xr<-exposure[[predictor_id]][match(km,kx)]
  patient<-metadata$HOSP_patient_days
  if(!is.numeric(y)||!is.numeric(xr)||!is.numeric(patient)||any(!is.finite(c(y,xr,patient))))
    stop("Outcome, exposure and patient-days values must be complete finite numbers",call.=FALSE)
  if(outcome_scale=="linear"){
    if(any(y<0))stop("Linear outcome abundance must be non-negative",call.=FALSE)
    y<-log10(y+1e-8)
  }
  hospital<-factor(as.character(metadata$hospital_code))
  xl<-log1p(pmax(xr,0))
  xb<-ave(xl,hospital,FUN=function(v)mean(v,na.rm=TRUE));xw<-xl-xb
  data.frame(y=as.numeric(y),hospital=hospital,quarter_f=factor(metadata$quarter),
    patient_z=fig6_hospital_fit_z(log1p(pmax(patient,0))),
    drug_within_z=fig6_hospital_fit_z(xw),drug_between_z=fig6_hospital_fit_z(xb),
    exposure_nonnegative=pmax(xr,0),row.names=NULL)
}

fig6_fit_hospital_lmm <- function(outcomes,exposures,metadata,dataset=c("full2646","group22"),
                                 outcome_ids=setdiff(names(outcomes),c("hospital_code","quarter")),
                                 predictor_ids=setdiff(names(exposures[[1]]),c("hospital_code","quarter")),
                                 outcome_scale=c("linear","log10"),
                                 expected_outcomes=NULL,expected_drugs=43L,
                                 expected_n=43L,expected_hospitals=15L,min_rows=30L,
                                 return_models=FALSE) {
  dataset<-as.character(dataset)[1];outcome_scale<-match.arg(outcome_scale)
  if(is.null(expected_outcomes)){
    expected_outcomes<-c(full2646=2646L,group22=22L)[dataset]
    if(is.na(expected_outcomes))stop("Set expected_outcomes for an explicitly defined alternative family",call.=FALSE)
  }
  if(!is.list(exposures)||!length(exposures)||is.null(names(exposures))||anyDuplicated(names(exposures))||
     any(!names(exposures) %in% c("oral","injection","combined")))
    stop("exposures must be a named oral/injection/combined list",call.=FALSE)
  if(length(outcome_ids)!=expected_outcomes||length(predictor_ids)!=expected_drugs||
     anyDuplicated(outcome_ids)||anyDuplicated(predictor_ids))
    stop("Incomplete planned outcome/drug family; do not fit only the displayed subset",call.=FALSE)
  if(!all(outcome_ids %in% names(outcomes))||any(!vapply(exposures,function(x)all(predictor_ids %in% names(x)),logical(1))))
    stop("Input matrix is missing a planned outcome or drug column",call.=FALSE)
  for(pkg in c("lme4","lmerTest"))if(!requireNamespace(pkg,quietly=TRUE))stop("Install ",pkg,call.=FALSE)
  formula<-stats::as.formula("y ~ drug_within_z + drug_between_z + quarter_f + patient_z + (1 | hospital)")
  planned_n<-length(outcome_ids)*length(predictor_ids)
  results<-vector("list",length(exposures)*planned_n);models<-list();index<-0L
  for(profile in names(exposures)){
    family_start<-index+1L
    for(oid in outcome_ids)for(did in predictor_ids){
      index<-index+1L
      d<-fig6_hospital_model_frame(outcomes,exposures[[profile]],metadata,oid,did,
        outcome_scale=outcome_scale,expected_n=expected_n,expected_hospitals=expected_hospitals)
      row<-data.frame(dataset=dataset,profile=profile,task_id=index-family_start+1L,
        outcome_id=oid,predictor_id=did,n=nrow(d),hospital_n=nlevels(droplevels(d$hospital)),
        status="skipped",reason="",beta_within=NA_real_,se_within=NA_real_,p_within=NA_real_,
        beta_between=NA_real_,p_between=NA_real_,ci_low=NA_real_,ci_high=NA_real_,
        fdr_within_outcome=NA_real_,fdr_global_profile_dataset=NA_real_,singular=NA,
        convergence_message="",aggregation="hospital_quarter_mean->log10(x+1e-08)",stringsAsFactors=FALSE)
      if(!any(d$exposure_nonnegative>0))row$reason<-"structural_zero"
      else if(stats::sd(d$drug_within_z)<=0)row$reason<-"no_within_hospital_variation"
      else if(nrow(d)<min_rows||nlevels(droplevels(d$hospital))<5L)row$reason<-"insufficient_rows"
      else if(stats::sd(d$y)<=0)row$reason<-"constant_outcome"
      else {
        warnings<-character()
        fit<-tryCatch(withCallingHandlers(lmerTest::lmer(formula,data=d,REML=FALSE),
          warning=function(w){warnings<<-c(warnings,conditionMessage(w));invokeRestart("muffleWarning")}),
          error=function(e)e)
        row$convergence_message<-paste(unique(warnings),collapse=" | ")
        if(inherits(fit,"error")){row$status<-"error";row$reason<-conditionMessage(fit)}
        else {
          co<-summary(fit)$coefficients
          value<-function(term,column)if(term %in% rownames(co)&&column %in% colnames(co))co[term,column] else NA_real_
          row$beta_within<-value("drug_within_z","Estimate")
          row$se_within<-value("drug_within_z","Std. Error")
          row$p_within<-value("drug_within_z","Pr(>|t|)")
          row$beta_between<-value("drug_between_z","Estimate")
          row$p_between<-value("drug_between_z","Pr(>|t|)")
          row$singular<-lme4::isSingular(fit,tol=1e-4)
          estimable<-all(is.finite(c(row$beta_within,row$se_within,row$p_within)))
          row$status<-if(estimable)"success" else "skipped"
          row$reason<-if(estimable)"" else "coefficient_not_estimable"
          row$ci_low<-row$beta_within-1.96*row$se_within
          row$ci_high<-row$beta_within+1.96*row$se_within
          if(return_models)models[[paste(profile,oid,did,sep="|")]]<-fit
        }
      }
      results[[index]]<-row
    }
  }
  results<-do.call(rbind,results)
  for(profile in names(exposures)){
    ok<-which(results$profile==profile&results$status=="success"&is.finite(results$p_within))
    results$fdr_global_profile_dataset[ok]<-stats::p.adjust(results$p_within[ok],"BH",n=planned_n)
    for(oid in outcome_ids){
      ix<-ok[results$outcome_id[ok]==oid]
      results$fdr_within_outcome[ix]<-stats::p.adjust(results$p_within[ix],"BH",n=length(predictor_ids))
    }
  }
  # Singular successful fits remain in BH, exactly as in the archived analysis.
  # No profile rows are deduplicated, even when their consumption vectors match.
  results$q<-results$fdr_global_profile_dataset
  audit<-do.call(rbind,lapply(names(exposures),function(profile){
    d<-results[results$profile==profile,,drop=FALSE]
    data.frame(dataset=dataset,profile=profile,planned_models=planned_n,
      successful_models=sum(d$status=="success"),skipped_models=sum(d$status=="skipped"),
      error_models=sum(d$status=="error"),singular_successful=sum(d$singular[d$status=="success"],na.rm=TRUE),
      expected_hospital_quarters=expected_n,expected_hospitals=expected_hospitals,
      stringsAsFactors=FALSE)
  }))
  list(results=results,audit=audit,models=if(return_models)models else NULL,
    specification=list(dataset=dataset,formula=paste(deparse(formula),collapse=" "),REML=FALSE,
      outcome_transform="log10(hospital-quarter mean abundance + 1e-8)",input_outcome_scale=outcome_scale,
      exposure_units="DDDs/100 patient-days",exposure_transform="log1p -> within/between -> standardize",
      patient_days_transform="log1p -> standardize",planned_per_profile=planned_n,
      BH_family="all planned outcome x drug comparisons, separately by dataset and route",
      confidence_interval="beta_within +/- 1.96 * se_within",p_value="lmerTest summary Pr(>|t|)",
      route_note="Coincident injection and combined estimates are not independent evidence."))
}



fig6_class_palette <- function() {
  c('Context factors'='#C97868','Population migration'='#6F95B8',
    'Air pollution'='#8F9EAA','Digital prescriptions'='#B494C8',
    'Meteorological factors'='#AFCB86','Antibiotic consumption'='#D6B4D9',
    'Hospital profile'='#E36774','Hospital continuity'='#FAA842')
}

fig6_p_stars <- function(p) {
  if (!is.numeric(p) || any(!is.na(p) & (!is.finite(p) | p<0 | p>1)))
    stop('StarP must contain numerical P values in [0,1], or NA for aggregates.')
  ifelse(is.na(p),'',ifelse(p<.001,'***',ifelse(p<.01,'**',ifelse(p<.05,'*',''))))
}

fig6_rf_geometry <- function(table, class_table=NULL) {
  need <- c('Predictor','Label','Class','Importance','StarP')
  if (!all(need %in% names(table))) stop('RF table needs ',paste(need,collapse=', '))
  d <- as.data.frame(table)
  d$Class <- as.character(d$Class); d$Label <- as.character(d$Label)
  if (!nrow(d) || anyNA(d[need[c(1,2,3,4)]]) || anyDuplicated(d$Predictor))
    stop('RF display keys must be unique and its labels/classes/scores complete.')
  if (!is.numeric(d$Importance) || any(!is.finite(d$Importance) | d$Importance<0) || sum(d$Importance)<=0)
    stop('RF display scores must be nonnegative with positive total; retain raw signed scores separately.')
  palette <- fig6_class_palette()
  if (any(!d$Class %in% names(palette))) stop('Unmapped predictor class: ',paste(setdiff(d$Class,names(palette)),collapse=', '))
  d$Stars <- fig6_p_stars(d$StarP)
  d <- d[order(match(d$Class,names(palette)),-d$Importance,d$Predictor),,drop=FALSE]
  d$id <- seq_len(nrow(d)); d$xmin <- d$id-.5; d$xmax <- d$id+.5
  d$ymin <- .82; d$ymax <- d$ymin+log10(d$Importance+1)
  d$name_y <- pmax(d$ymax+.12,d$ymin+.30)
  d$angle <- 90-360*(d$id-.5)/nrow(d)
  d$hjust <- ifelse(d$angle< -90,1,0)
  d$text_angle <- ifelse(d$angle< -90,d$angle+180,d$angle)
  d$ValueLabel <- paste0(sprintf('%.1f%%',d$Importance),d$Stars)
  totals <- aggregate(d$Importance,list(Class=d$Class),sum)
  names(totals)[2] <- 'ScoreSum'; totals$Percent <- totals$ScoreSum/sum(totals$ScoreSum)*100
  if (!is.null(class_table)) {
    if (!all(c('Class','Importance') %in% names(class_table)) || anyDuplicated(class_table$Class))
      stop('Class table requires unique Class and Importance columns.')
    expected <- class_table$Importance/sum(class_table$Importance)*100
    names(expected)<-class_table$Class
    if (!setequal(totals$Class,names(expected)) || any(abs(totals$Percent-expected[totals$Class])>1e-7))
      stop('Display wedges do not reproduce the all-predictor class shares; do not renormalize a selected subset.')
  }
  totals$xmin <- vapply(totals$Class,function(s) min(d$xmin[d$Class==s]),numeric(1))
  totals$xmax <- vapply(totals$Class,function(s) max(d$xmax[d$Class==s]),numeric(1))
  totals$xmid <- (totals$xmin+totals$xmax)/2
  totals$ymin <- max(d$ymax)+1.65; totals$ymax <- totals$ymin+.45
  totals$angle <- 90-360*(totals$xmid-.5)/nrow(d)
  totals$text_angle <- ifelse(totals$angle< -90,totals$angle+180,totals$angle)
  list(variables=d,classes=totals,palette=palette,n=nrow(d))
}

fig6_plot_rf <- function(table,title,class_table=NULL) {
  if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Package ggplot2 is required.')
  g <- fig6_rf_geometry(table,class_table); d<-g$variables; cl<-g$classes
  # Equal angular slots: neither sector width nor outer-ring area is an
  # importance measure. Importance is shown by numerical annotations/bar height.
  ggplot2::ggplot()+
    ggplot2::geom_rect(data=d,ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Class),
      colour='white',linewidth=.28)+
    ggplot2::geom_text(data=d,ggplot2::aes(x=id,y=ymin+(ymax-ymin)*.5,label=ValueLabel,angle=text_angle),
      size=2.1,colour='black')+
    ggplot2::geom_text(data=d,ggplot2::aes(x=id,y=name_y,label=Label,angle=text_angle,hjust=hjust),
      size=2.2,colour='black')+
    ggplot2::geom_rect(data=cl,ggplot2::aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax,fill=Class),
      colour='white',linewidth=.4)+
    ggplot2::geom_text(data=cl,ggplot2::aes(x=xmid,y=(ymin+ymax)/2,label=sprintf('%.1f%%',Percent),angle=text_angle),
      fontface='bold',size=3.5)+
    ggplot2::annotate('text',x=.5,y=0,label=title,size=3.1,fontface='bold',lineheight=.95)+
    ggplot2::scale_fill_manual(values=g$palette,name='Predictor class')+
    ggplot2::scale_x_continuous(limits=c(.5,g$n+.5),expand=c(0,0))+
    ggplot2::scale_y_continuous(limits=c(0,max(cl$ymax)+.5),expand=c(0,0))+
    ggplot2::coord_polar(theta='x',clip='off')+
    ggplot2::theme_void(base_size=10)+
    ggplot2::theme(legend.position='bottom',legend.text=ggplot2::element_text(size=8),
      plot.margin=ggplot2::margin(12,20,12,20))+
    ggplot2::guides(fill=ggplot2::guide_legend(nrow=2,byrow=TRUE))
}

fig6_network_plot_data <- function(table) {
  need<-c('Predictor','Label','Class','Lag','Estimate','CI_low','CI_high','P','Q')
  if (!all(need %in% names(table))) stop('Network LMM table needs ',paste(need,collapse=', '))
  d<-as.data.frame(table)
  if (!nrow(d) || anyDuplicated(d[c('Predictor','Lag')])) stop('Each predictor-lag estimate must occur once.')
  if (anyNA(d[need]) || any(!is.finite(d$Estimate)) || any(d$CI_low>d$Estimate | d$CI_high<d$Estimate))
    stop('LMM coefficients and intervals must be complete and ordered.')
  if (any(d$Q<0 | d$Q>1)) stop('Invalid BH-adjusted P value.')
  d$Lag<-factor(as.character(d$Lag),levels=c('Lag0','Lag30','Lag60','Lag90'))
  if (anyNA(d$Lag)) stop('Unknown lag; use Lag0, Lag30, Lag60 or Lag90.')
  d$Direction<-factor(ifelse(d$Q<.05 & d$Estimate>0,'Positive',
    ifelse(d$Q<.05 & d$Estimate<0,'Negative','Not significant')),
    levels=c('Positive','Negative','Not significant'))
  d$Stars<-ifelse(d$Q<.05,'*','')
  preferred<-c('Immigration','Phenolic biocides','Aminoglycoside','MLS',
    'Emigration','NSAIDs','Beta-lactam','Tetracycline')
  present<-unique(as.character(d$Label))
  d$Label<-factor(d$Label,levels=c(intersect(preferred,present),setdiff(present,preferred)))
  d
}

fig6_plot_network_lmm <- function(table,title) {
  if (!requireNamespace('ggplot2',quietly=TRUE)) stop('Package ggplot2 is required.')
  d<-fig6_network_plot_data(table)
  ranges<-tapply(d$CI_high-d$CI_low,d$Label,max)
  d$star_y<-d$CI_high+as.numeric(ranges[as.character(d$Label)])*.10
  ggplot2::ggplot(d,ggplot2::aes(Lag,Estimate,colour=Direction))+
    ggplot2::geom_hline(yintercept=0,linetype='dashed',colour='grey55',linewidth=.4)+
    ggplot2::geom_errorbar(ggplot2::aes(ymin=CI_low,ymax=CI_high),width=.16,linewidth=.6)+
    ggplot2::geom_point(shape=21,fill='white',size=2.6,stroke=.8)+
    ggplot2::geom_text(ggplot2::aes(y=star_y,label=Stars),colour='black',size=5.2)+
    ggplot2::facet_wrap(~Label,ncol=4,scales='free_y')+
    ggplot2::scale_colour_manual(values=c(Positive='#D86B71',Negative='#507EA5','Not significant'='#9B9B9B'),drop=FALSE)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.08,.20)))+
    ggplot2::labs(x=NULL,y='Fixed-effect coefficient (0-1 scaled)',title=title,colour=NULL)+
    ggplot2::theme_bw(base_size=11)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),strip.background=ggplot2::element_blank(),
      strip.text=ggplot2::element_text(face='bold',size=10),
      axis.text.x=ggplot2::element_text(angle=45,hjust=1),legend.position='bottom',
      plot.title=ggplot2::element_text(size=11),panel.spacing=grid::unit(.7,'lines'))
}

fig6_save_plot <- function(plot,stem,width,height) {
  dir.create(dirname(stem),recursive=TRUE,showWarnings=FALSE)
  dev<-if (capabilities('cairo')) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(paste0(stem,'.pdf'),plot,width=width,height=height,device=dev,bg='white')
  ggplot2::ggsave(paste0(stem,'.png'),plot,width=width,height=height,dpi=400,bg='white')
}


fig6_input_table <- function(x, type, outcome, config) {
  if (is.character(x) && length(x)==1L) x <- fig6_network_read_table(x)
  if (is.matrix(x)) x <- as.data.frame(x, check.names=FALSE)
  if (is.list(x) && !is.data.frame(x)) {
    if (!is.null(x$plot_data)) x <- x$plot_data
    else if (type=='rf' && !is.null(x$importance)) x <- x$importance
    else if (type=='lmm' && !is.null(x$results)) x <- x$results
  }
  if (!is.data.frame(x)) stop('Supply a result table, table path, or fitted result list for ',outcome)
  if (!nrow(x)) stop('No estimable results for ',outcome,'; inspect the fitting diagnostics.')
  if (type=='rf') {
    if (all(c('Predictor','Label','Class','Importance','StarP') %in% names(x))) return(x)
    return(fig6_network_rf_table(x,outcome,config))
  }
  if (!all(c('Predictor','Label','Class','Lag','Estimate','CI_low','CI_high','P','Q') %in% names(x)))
    x <- fig6_network_lmm_table(x,outcome,config,panel_only=FALSE)
  # BH values originate from the full class/window family, not this display set.
  x$Display <- x$Predictor %in% config$panel_predictors
  x
}

run_fig6 <- function(output_dir,
                     network_rf=NULL, network_lmm=NULL,
                     hospital_rf_dir=NULL, hospital_lmm_dir=NULL,
                     hospital_selection_file=NULL,
                     panels=letters[1:8], config=fig6_network_config(),
                     combine=TRUE, family='Arial',
                     hospital_rf=NULL, hospital_lmm=NULL,
                     hospital_selection=NULL) {
  if (!length(panels) || anyNA(panels) || any(!panels %in% letters[1:8]) || anyDuplicated(panels))
    stop('panels must contain unique letters a-h.')
  needed <- function(input, name, map) {
    keys <- unname(map[intersect(names(map),panels)])
    if (length(keys) && (is.null(input) || !all(keys %in% names(input))))
      stop(name,' must have named entries: ',paste(keys,collapse=', '),'.')
  }
  needed(network_rf,'network_rf',c(a='Total',b='TierI'))
  needed(network_lmm,'network_lmm',c(e='Total',f='TierI'))
  if (any(c('c','d') %in% panels) && is.null(hospital_rf_dir) && is.null(hospital_rf))
    stop('Supply hospital_rf_dir or hospital_rf for panels c,d.')
  if (any(c('g','h') %in% panels) && is.null(hospital_lmm_dir) && is.null(hospital_lmm))
    stop('Supply hospital_lmm_dir or hospital_lmm for panels g,h.')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  table_dir <- file.path(output_dir,'plot_data')
  dir.create(table_dir,recursive=TRUE,showWarnings=FALSE)
  plots <- tables <- list()
  save_table <- function(x,name) utils::write.csv(x,file.path(table_dir,paste0(name,'.csv')),row.names=FALSE,na='')

  for (pn in intersect(c('a','b'),panels)) {
    key <- c(a='Total',b='TierI')[[pn]]
    d <- fig6_input_table(network_rf[[key]],'rf',key,config)
    tables[[pn]] <- d
    title <- if(pn=='a') 'Total ARG\nNetwork' else 'Tier I ARG\nNetwork'
    plots[[pn]] <- fig6_plot_rf(d,title)
    save_table(d,paste0('Fig6',pn,'_predictors'))
    save_table(fig6_network_rf_classes(d),paste0('Fig6',pn,'_classes'))
  }
  if (any(c('c','d') %in% panels)) {
    raw <- if(is.null(hospital_rf)) fig6_read_hospital_rf(hospital_rf_dir) else hospital_rf
    if (!all(c('total','group','dictionary','catalog') %in% names(raw)))
      stop('hospital_rf needs total, group, dictionary and catalog tables.')
    hospital <- fig6_prepare_hospital_rf(raw$total,raw$group,raw$dictionary,raw$catalog)
    for (pn in intersect(c('c','d'),panels)) {
      d <- hospital[[paste0('display_',pn)]]
      tables[[pn]] <- d
      title <- if(pn=='c') 'Total ARG\nHospital' else '22 ARG\nmodels\nHospital'
      plots[[pn]] <- fig6_plot_rf(d,title,hospital[[paste0('class_',pn)]])
      save_table(d,paste0('Fig6',pn,'_display'))
      save_table(hospital[[pn]],paste0('Fig6',pn,'_all_predictors'))
      save_table(hospital[[paste0('class_',pn)]],paste0('Fig6',pn,'_classes'))
    }
  }
  for (pn in intersect(c('e','f'),panels)) {
    key <- c(e='Total',f='TierI')[[pn]]
    all <- fig6_input_table(network_lmm[[key]],'lmm',key,config)
    d <- all[all$Display,,drop=FALSE]
    tables[[pn]] <- d
    title <- if(pn=='e') 'Total ARG abundance' else 'Tier I ARG abundance'
    plots[[pn]] <- fig6_plot_network_lmm(d,title)
    save_table(all,paste0('Fig6',pn,'_all_coefficients'))
    save_table(d,paste0('Fig6',pn,'_display'))
  }
  if (any(c('g','h') %in% panels)) {
    if (is.null(hospital_selection)) {
      if (is.null(hospital_selection_file) && !is.null(hospital_lmm_dir))
        hospital_selection_file <- file.path(hospital_lmm_dir,'hospital_lmm_plot_selection.tsv')
      if (is.null(hospital_selection_file) || !file.exists(hospital_selection_file))
        stop('Supply hospital_selection or an existing hospital_selection_file.')
      selection <- utils::read.delim(hospital_selection_file,check.names=FALSE,stringsAsFactors=FALSE)
    } else selection <- hospital_selection
    full <- if(is.null(hospital_lmm)) fig6_hospital_read_results(hospital_lmm_dir) else hospital_lmm
    if (!is.list(full) || !is.data.frame(full$results)) stop('hospital_lmm needs a results table and optional audit table.')
    if (!'q' %in% names(full$results) && 'fdr_global_profile_dataset' %in% names(full$results))
      full$results$q <- full$results$fdr_global_profile_dataset
    selected <- fig6_hospital_select(full$results,selection)
    hp <- fig6_hospital_plots(selected,family=family)
    for (pn in intersect(c('g','h'),panels)) plots[[pn]] <- hp[[pn]]
    tables$hospital_lmm <- selected
    save_table(selected,'Fig6gh_display')
    if(is.data.frame(full$audit)) save_table(full$audit,'Fig6gh_family_checks')
  }
  for (pn in panels) {
    size <- if(pn %in% letters[1:4]) c(7.2,7.2) else if(pn %in% c('e','f')) c(10,5.4) else if(pn=='g') c(491/72,234.5/72) else c(490.25/72,141.7/72)
    fig6_save_plot(plots[[pn]],file.path(output_dir,paste0('Fig6',pn)),size[1],size[2])
  }
  # Optional assembly. Individual vector panels remain available for layout edits.
  if (isTRUE(combine) && setequal(panels,letters[1:8]) && requireNamespace('patchwork',quietly=TRUE)) {
    order <- c('a','e','b','f','c','g','d','h')
    tagged <- lapply(order,function(pn) plots[[pn]]+ggplot2::labs(tag=pn))
    combined <- patchwork::wrap_plots(tagged,ncol=2,widths=c(1,1.65))
    fig6_save_plot(combined,file.path(output_dir,'Fig6_a_to_h'),18,23)
  }
  writeLines(capture.output(sessionInfo()),file.path(output_dir,'sessionInfo.txt'))
  invisible(list(plots=plots[panels],tables=tables))
}


