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
  missing <- setdiff(expected, names(data))
  if (length(missing)) stop("Missing predictors for ", class, "/", lag, ": ",
                            paste(missing, collapse = ", "), call. = FALSE)
  columns <- c(config$outcome, predictors, "season", "city", "site_nm")
  fig6_network_columns(data, columns)
  quote_name <- function(x) paste0("`", x, "`")
  formula_text <- paste(quote_name(config$outcome), "~", paste(quote_name(predictors), collapse = " + "),
                        "+ factor(season) + (1 | city) + (1 | site_nm)")
  list(predictors = predictors, missing_predictors = setdiff(expected, predictors),
       rows = which(stats::complete.cases(data[columns])), formula_text = formula_text,
       formula = stats::as.formula(formula_text))
}

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

fig6_rf_permute <- function(..., num.rep = 1000L) {
  fig6_network_require("rfPermute")
  if (length(num.rep) != 1L || !is.finite(num.rep) || num.rep < 1L || num.rep != as.integer(num.rep))
    stop("num.rep must be a positive integer.", call. = FALSE)
  parameters <- names(formals(getS3method("rfPermute", "default", envir = asNamespace("rfPermute"))))
  parameter <- if ("num.rep" %in% parameters) "num.rep" else if ("nrep" %in% parameters) "nrep" else NA_character_
  if (is.na(parameter)) stop("Installed rfPermute has no recognized explicit permutation-count parameter.", call. = FALSE)
  arguments <- list(...)
  arguments[[parameter]] <- as.integer(num.rep)
  do.call(rfPermute::rfPermute, arguments)
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

si245_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Install required packages: ", paste(absent, collapse = ", "))
}

si245_input <- function(abundance, metadata) {
  a <- as.matrix(abundance); m <- as.data.frame(metadata)
  required <- c("sample_id", "city", "setting", "month", "physical_site")
  if (!all(required %in% names(m))) stop("Metadata requires: ", paste(required, collapse = ", "))
  if (!is.numeric(a) || is.null(rownames(a)) || is.null(colnames(a)) ||
      anyDuplicated(rownames(a)) || anyDuplicated(colnames(a)) || anyNA(a) ||
      any(!is.finite(a)) || any(a < 0)) stop("Use a finite nonnegative named abundance matrix; resolve missing measurements explicitly.")
  if (anyNA(m[required]) || anyDuplicated(m$sample_id) || !setequal(m$sample_id, rownames(a)))
    stop("Metadata/sample matrix keys must be complete, unique and identical.")
  m <- m[match(rownames(a), m$sample_id), , drop = FALSE]
  if (any(!m$setting %in% c("Hospital", "WWTP", "Community", "Wet market"))) stop("Unknown setting.")
  sites <- unique(m[c("physical_site", "city", "setting")])
  if (anyDuplicated(sites$physical_site)) stop("Physical site maps to multiple cities/settings.")
  list(a = a, m = m)
}

si245_site_means <- function(abundance, metadata) {
  z <- si245_input(abundance, metadata)
  sites <- unique(as.character(z$m$physical_site))
  a <- vapply(sites, function(id) colMeans(z$a[z$m$physical_site == id, , drop = FALSE]), numeric(ncol(z$a)))
  dim(a) <- c(ncol(z$a), length(sites))
  a <- t(a)
  colnames(a) <- colnames(z$a); rownames(a) <- sites
  m <- z$m[match(sites, z$m$physical_site), , drop = FALSE]
  m$sample_id <- sites
  list(a = a, m = m)
}

si245_ordination <- function(a, m, site_level = FALSE, permutations = 999L, seed = 123L) {
  si245_require(c("vegan", "permute"))
  removed <- m$sample_id[rowSums(a) == 0]
  keep <- rowSums(a) > 0
  a <- a[keep, , drop = FALSE]; m <- droplevels(m[keep, , drop = FALSE])
  m$setting <- factor(m$setting); m$city <- factor(m$city)
  if (nrow(a) < 4L || nlevels(m$setting) < 2L) stop("Too few nonzero profiles or settings for ordination/testing.")
  distance <- vegan::vegdist(a, method = "bray")
  if (all(as.numeric(distance) == 0)) stop("All profiles have zero Bray-Curtis distance.")
  ord <- vegan::wcmdscale(distance, k = 2L, eig = TRUE, add = if (site_level) "lingoes" else FALSE)
  points <- as.matrix(ord$points)
  if (ncol(points) < 2L) points <- cbind(points, 0)
  positive <- ord$eig[ord$eig > 0]
  explained <- 100 * ord$eig[seq_len(min(2L, length(ord$eig)))] / sum(positive)
  if (length(explained) < 2L) explained <- c(explained, 0)
  control <- permute::how(nperm = as.integer(permutations))
  if (site_level) permute::setBlocks(control) <- m$city
  set.seed(seed)
  model <- if (site_level && nlevels(m$city) > 1L)
    vegan::adonis2(distance ~ city + setting, data = m, permutations = control, by = "margin") else
    vegan::adonis2(distance ~ setting, data = m, permutations = control)
  dispersion <- vegan::betadisper(distance, m$setting, type = "median")
  set.seed(seed)
  dp <- vegan::permutest(dispersion, permutations = control)
  row <- match("setting", rownames(model))
  stats <- data.frame(n = nrow(m), n_sites = length(unique(m$physical_site)),
    R2 = model$R2[row], F = model$F[row], P = model$`Pr(>F)`[row],
    Dispersion_F = dp$tab$F[1], Dispersion_P = dp$tab$`Pr(>F)`[1],
    permutations = permutations, Site_level = site_level)
  list(points = cbind(m, PCoA1 = points[, 1], PCoA2 = points[, 2]), explained = explained,
    stats = stats, permanova = model, dispersion = dp, distance = distance,
    excluded_zero_profiles = removed)
}

si245_s2 <- function(abundance, metadata, tier1_subtypes, seed = 123L) {
  z <- si245_input(abundance, metadata)
  if (!length(tier1_subtypes) || any(!tier1_subtypes %in% colnames(z$a))) stop("Tier I list does not match abundance columns.")
  outcomes <- list(Total = z$a, TierI = z$a[, unique(tier1_subtypes), drop = FALSE])
  sample <- site <- list()
  for (outcome in names(outcomes)) {
    sample[[outcome]] <- lapply(split(seq_len(nrow(z$m)), z$m$city), function(i)
      si245_ordination(outcomes[[outcome]][i, , drop = FALSE], z$m[i, ], FALSE, 999L, seed))
    aggregated <- si245_site_means(outcomes[[outcome]], z$m)
    site[[outcome]] <- si245_ordination(aggregated$a, aggregated$m, TRUE, 9999L, seed)
  }
  list(sample = sample, site = site)
}

si245_plot_s2 <- function(result, ellipse_level = .68) {
  si245_require("ggplot2")
  colours <- c(Hospital = "#E06C75", WWTP = "#4E79A7", Community = "#6EC5C1", "Wet market" = "#E4BE54")
  panels <- list()
  for (outcome in names(result$sample)) {
    fits <- result$sample[[outcome]]
    d <- do.call(rbind, lapply(names(fits), function(city) {
      q <- fits[[city]]$points
      q$panel <- sprintf("%s | axes %.1f%%, %.1f%%", city, fits[[city]]$explained[1], fits[[city]]$explained[2])
      q
    }))
    annotations <- do.call(rbind, lapply(names(fits), function(city) {
      f <- fits[[city]]
      data.frame(panel = sprintf("%s | axes %.1f%%, %.1f%%", city, f$explained[1], f$explained[2]),
        label = sprintf("R\u00b2 = %.3f; P = %.3g", f$stats$R2, f$stats$P))
    }))
    panels[[paste0(outcome, "_sample")]] <- ggplot2::ggplot(d, ggplot2::aes(PCoA1, PCoA2, colour = setting)) +
      ggplot2::geom_point(size = 1, alpha = .6) +
      ggplot2::stat_ellipse(type = "norm", level = ellipse_level, linewidth = .4) +
      ggplot2::facet_wrap(~panel, scales = "free") + ggplot2::scale_colour_manual(values = colours) +
      ggplot2::geom_text(data = annotations, ggplot2::aes(x = -Inf, y = Inf, label = label),
        inherit.aes = FALSE, hjust = -.02, vjust = 1.3, size = 3) +
      ggplot2::theme_bw() + ggplot2::labs(title = outcome, colour = NULL)
    s <- result$site[[outcome]]
    panels[[paste0(outcome, "_site")]] <- ggplot2::ggplot(s$points, ggplot2::aes(PCoA1, PCoA2, colour = setting)) +
      ggplot2::geom_point(size = 2) + ggplot2::scale_colour_manual(values = colours) + ggplot2::theme_bw() +
      ggplot2::labs(x = sprintf("PCoA1 (%.1f%%)", s$explained[1]), y = sprintf("PCoA2 (%.1f%%)", s$explained[2]),
        title = paste(outcome, "site means"), subtitle = sprintf("City-adjusted R\u00b2 = %.3f; P = %.3g", s$stats$R2, s$stats$P), colour = NULL)
  }
  panels
}

si245_nb <- function(d) {
  si245_require(c("lme4", "multcompView"))
  d$setting <- droplevels(factor(d$setting, levels = c("Hospital", "WWTP", "Community", "Wet market")))
  d$month <- factor(d$month); d$physical_site <- factor(d$physical_site)
  if (nlevels(d$setting) < 2L || nlevels(d$month) < 2L) stop("S4 requires at least two settings and two sampling months in each analysed city.")
  fit <- lme4::glmer.nb(richness ~ setting + month + (1 | physical_site), data = d,
    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
  beta <- lme4::fixef(fit); V <- as.matrix(stats::vcov(fit))
  grid <- d[rep(1L, nlevels(d$setting)), , drop = FALSE]
  grid$setting <- factor(levels(d$setting), levels = levels(d$setting))
  X <- stats::model.matrix(~setting + month, data = grid)
  X <- X[, names(beta), drop = FALSE]
  pairs <- utils::combn(seq_len(nrow(grid)), 2L)
  comparisons <- do.call(rbind, lapply(seq_len(ncol(pairs)), function(j) {
    i <- pairs[1, j]; k <- pairs[2, j]; L <- X[i, ] - X[k, ]
    estimate <- sum(L * beta); se <- sqrt(drop(L %*% V %*% L))
    data.frame(group1 = as.character(grid$setting[i]), group2 = as.character(grid$setting[k]),
      log_ratio = estimate, SE = se, ratio = exp(estimate),
      ratio_low = exp(estimate - 1.96 * se), ratio_high = exp(estimate + 1.96 * se),
      P = 2 * stats::pnorm(-abs(estimate / se)))
  }))
  comparisons$Q <- stats::p.adjust(comparisons$P, "BH")
  p <- stats::setNames(comparisons$Q, paste(comparisons$group1, comparisons$group2, sep = "-"))
  letters <- multcompView::multcompLetters(p)$Letters
  list(model = fit, comparisons = comparisons, letters = data.frame(setting = names(letters), letter = unname(letters)),
    diagnostics = data.frame(n = nrow(d), n_sites = nlevels(d$physical_site),
      convergence = paste(fit@optinfo$conv$opt, collapse = ";"),
      messages = paste(fit@optinfo$conv$lme4$messages, collapse = ";"),
      singular = lme4::isSingular(fit), theta = lme4::getME(fit, "glmer.nb.theta")), data = d)
}

si245_s4 <- function(abundance, metadata, tier1_subtypes) {
  z <- si245_input(abundance, metadata)
  if (!length(tier1_subtypes) || any(!tier1_subtypes %in% colnames(z$a))) stop("Supply the exact Tier I subtype list.")
  data <- rbind(transform(z$m, outcome = "Total", richness = rowSums(z$a > 0)),
    transform(z$m, outcome = "TierI", richness = rowSums(z$a[, unique(tier1_subtypes), drop = FALSE] > 0)))
  fits <- lapply(split(data, interaction(data$city, data$outcome, drop = TRUE)), si245_nb)
  list(fits = fits, data = data)
}

si245_plot_s4 <- function(result) {
  si245_require("ggplot2")
  cld <- do.call(rbind, lapply(result$fits, function(f) transform(f$letters,
    city = f$data$city[1], outcome = f$data$outcome[1], y = max(f$data$richness) * 1.1)))
  means <- stats::aggregate(richness ~ setting + city + outcome, result$data, mean)
  d <- result$data; d$setting <- factor(d$setting, levels = c("Hospital", "WWTP", "Community", "Wet market"))
  ggplot2::ggplot(d, ggplot2::aes(setting, richness, fill = setting)) +
    ggplot2::geom_violin(trim = TRUE, alpha = .4) + ggplot2::geom_boxplot(width = .14, outlier.shape = NA) +
    ggplot2::geom_jitter(width = .08, size = .5, alpha = .2) +
    ggplot2::geom_text(data = cld, ggplot2::aes(y = y, label = letter), size = 3) +
    ggplot2::geom_text(data = means, ggplot2::aes(label = sprintf("%.1f", richness)), vjust = -1, size = 2.5) +
    ggplot2::facet_grid(outcome ~ city, scales = "free_y") + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1), legend.position = "none") +
    ggplot2::labs(x = NULL, y = "ARG subtype richness")
}

si245_group_abundance <- function(abundance, mapping) {
  a <- as.matrix(abundance)
  if (!all(c("Subtype", "Feature") %in% names(mapping)) || anyNA(mapping[c("Subtype", "Feature")]) ||
      anyDuplicated(mapping[c("Subtype", "Feature")])) stop("mapping needs unique Subtype/Feature pairs.")
  if (any(!mapping$Subtype %in% colnames(a))) stop("Mapped subtypes are absent from abundance input.")
  features <- unique(as.character(mapping$Feature))
  out <- vapply(features, function(f) rowSums(a[, mapping$Subtype[mapping$Feature == f], drop = FALSE]), numeric(nrow(a)))
  dim(out) <- c(nrow(a), length(features))
  rownames(out) <- rownames(a); colnames(out) <- features
  out
}

si245_s5_tests <- function(abundance, metadata, seed = 123L) {
  z <- si245_input(abundance, metadata); a <- z$a; m <- z$m
  set.seed(seed); rows <- overall <- list(); index <- 0L
  for (setting in unique(m$setting)) for (feature in colnames(a)) {
    i <- which(m$setting == setting); x <- a[i, feature]; city <- as.character(m$city[i]); positive <- x[x > 0]
    if (length(positive) < 10L || length(unique(city)) < 2L) next
    threshold <- as.numeric(stats::quantile(positive, .9, type = 7))
    event <- x > 0 & x >= threshold
    counts <- table(factor(city), factor(event, levels = c(FALSE, TRUE)))
    global_p <- if (any(colSums(counts) == 0)) 1 else
      suppressWarnings(stats::chisq.test(counts, simulate.p.value = TRUE, B = 9999L)$p.value)
    index <- index + 1L
    overall[[index]] <- data.frame(Setting = setting, Feature = feature, Overall_P = global_p,
      threshold = threshold, n_positive = length(positive), n_total = length(x))
    rows[[index]] <- do.call(rbind, lapply(unique(city), function(target) {
      inside <- city == target
      aa <- sum(event & inside); bb <- sum(!event & inside)
      cc <- sum(event & !inside); dd <- sum(!event & !inside)
      rate <- aa / (aa + bb); other <- cc / (cc + dd)
      ratio <- if (other == 0) if (rate == 0) 0 else Inf else rate / other
      fisher <- stats::fisher.test(matrix(c(aa, bb, cc, dd), nrow = 2, byrow = TRUE), alternative = "greater")
      data.frame(Setting = setting, Feature = feature, City = target, a = aa, b = bb, c = cc, d = dd,
        N = aa + bb, Events = aa, Rate = rate, Other_rate = other, Ratio = ratio,
        Mean_abundance = mean(x[inside]), OR = unname(fisher$estimate), P = fisher$p.value)
    }))
  }
  if (!length(rows)) {
    empty <- data.frame(Setting = character(), Feature = character(), City = character(),
      a = integer(), b = integer(), c = integer(), d = integer(), N = integer(), Events = integer(),
      Rate = numeric(), Other_rate = numeric(), Ratio = numeric(), Mean_abundance = numeric(), OR = numeric(), P = numeric(),
      p_enriched_FDR_within_feature = numeric(), Overall_P = numeric(), Overall_Q = numeric(),
      threshold = numeric(), Enriched = logical(), Evidence = character())
    return(list(comparisons = empty, overall = data.frame(), winners = empty,
      recurrence = data.frame(City = character(), Feature = character(), Supported_settings = integer(), Max_rate = numeric()),
      status = "No setting-group has >=10 positive observations and at least two cities"))
  }
  tab <- do.call(rbind, rows); global <- do.call(rbind, overall)
  global$Overall_Q <- ave(global$Overall_P, global$Setting, FUN = function(x) stats::p.adjust(x, "BH"))
  tab$p_enriched_FDR_within_feature <- ave(tab$P, interaction(tab$Setting, tab$Feature, drop = TRUE), FUN = function(x) stats::p.adjust(x, "BH"))
  key <- paste(tab$Setting, tab$Feature, sep = "\r")
  j <- match(key, paste(global$Setting, global$Feature, sep = "\r"))
  tab$Overall_P <- global$Overall_P[j]; tab$Overall_Q <- global$Overall_Q[j]
  tab$threshold <- global$threshold[j]
  tab$Enriched <- tab$p_enriched_FDR_within_feature < .05 & tab$Ratio >= 1.5
  high <- tab$Enriched & tab$Overall_Q < .05 & tab$Ratio >= 2 & tab$N >= 10 & tab$Events >= 3
  moderate <- !high & tab$Enriched & tab$Overall_Q < .1 & tab$N >= 5 & tab$Events >= 2
  tab$Evidence <- ifelse(high, "High", ifelse(moderate, "Moderate", "Low"))
  choose <- function(d) {
    supported <- d$Evidence %in% c("High", "Moderate")
    if (any(supported)) {
      d <- d[supported, , drop = FALSE]
      ord <- order(-match(d$Evidence, c("Low", "Moderate", "High")), -d$Ratio, -d$Rate, -d$Events, d$City)
    } else if (any(d$Enriched)) {
      d <- d[d$Enriched, , drop = FALSE]; ord <- order(-d$Ratio, -d$Rate, -d$Events, d$City)
    } else ord <- order(-d$Rate, -d$Events, -d$Mean_abundance, d$City)
    d[ord[1], , drop = FALSE]
  }
  winners <- do.call(rbind, lapply(split(tab, interaction(tab$Setting, tab$Feature, drop = TRUE)), choose))
  supported <- tab[tab$Evidence %in% c("High", "Moderate"), , drop = FALSE]
  recurrence <- if (nrow(supported)) merge(
    stats::aggregate(list(Supported_settings = supported$Setting), supported[c("City", "Feature")], function(x) length(unique(x))),
    stats::aggregate(list(Max_rate = supported$Rate), supported[c("City", "Feature")], max), by = c("City", "Feature")) else
    data.frame(City = character(), Feature = character(), Supported_settings = integer(), Max_rate = numeric())
  list(comparisons = tab, overall = global, winners = winners, recurrence = recurrence)
}

si245_s5_site <- function(group_abundance, metadata, primary) {
  z <- si245_input(group_abundance, metadata)
  site <- si245_site_means(z$a, z$m)
  counts <- table(z$m$physical_site)
  profiles <- comparisons <- list()
  index <- 0L
  for (setting in unique(as.character(site$m$setting))) for (feature in colnames(site$a)) {
    i <- which(site$m$setting == setting)
    d <- data.frame(Setting = setting, Feature = feature, City = as.character(site$m$city[i]),
      PhysicalSite = as.character(site$m$physical_site[i]), Site_mean_burden = site$a[i, feature],
      Site_n_samples = as.integer(counts[site$m$physical_site[i]]), stringsAsFactors = FALSE)
    positive <- d$Site_mean_burden > 0
    eligible <- sum(positive) >= 10L && length(unique(d$City)) >= 2L
    threshold <- if (eligible) as.numeric(stats::quantile(d$Site_mean_burden[positive], .9, type = 7)) else NA_real_
    d$Site_positive <- positive
    d$Site_outlier_threshold <- threshold
    d$Site_event <- if (eligible) positive & d$Site_mean_burden >= threshold else NA
    index <- index + 1L
    profiles[[index]] <- d
    comparisons[[index]] <- do.call(rbind, lapply(sort(unique(d$City)), function(city) {
      target <- d$City == city
      n_target <- sum(target); n_other <- sum(!target)
      aa <- bb <- cc <- dd <- odds <- p <- target_rate <- other_rate <- ratio <- NA_real_
      if (eligible) {
        aa <- sum(d$Site_event[target]); bb <- n_target - aa
        cc <- sum(d$Site_event[!target]); dd <- n_other - cc
        ft <- stats::fisher.test(matrix(c(aa, bb, cc, dd), nrow = 2L, byrow = TRUE), alternative = "greater")
        odds <- unname(ft$estimate); p <- ft$p.value
        target_rate <- aa / n_target; other_rate <- cc / n_other
        ratio <- if (other_rate > 0) target_rate / other_rate else if (target_rate > 0) Inf else NA_real_
      }
      data.frame(Setting = setting, Feature = feature, City = city,
        Site_a = aa, Site_b = bb, Site_c = cc, Site_d = dd,
        Target_sites = n_target, Other_sites = n_other,
        Target_site_event_rate = target_rate, Other_site_event_rate = other_rate,
        Site_rate_ratio = ratio, Site_OR = odds, Site_P_enriched = p,
        Site_outlier_threshold = threshold, N_sites_in_group = nrow(d),
        N_positive_sites_in_group = sum(positive),
        N_site_events_in_group = if (eligible) sum(d$Site_event) else NA_integer_,
        Median_samples_per_site = stats::median(d$Site_n_samples),
        Site_status = if (eligible) "estimated" else if (sum(positive) < 10L) "fewer_than_10_positive_sites" else "fewer_than_two_cities",
        stringsAsFactors = FALSE)
    }))
  }
  tab <- do.call(rbind, comparisons)
  tab$Site_Q_within_setting_ARG_group <- ave(tab$Site_P_enriched,
    interaction(tab$Setting, tab$Feature, drop = TRUE), FUN = function(x) stats::p.adjust(x, "BH"))
  tab$Site_positive_q05 <- tab$Site_OR > 1 & tab$Site_Q_within_setting_ARG_group < .05
  key <- c("Setting", "Feature", "City")
  p <- primary$comparisons[c(key, "a", "b", "c", "d", "Rate", "Other_rate", "Ratio", "OR", "P", "p_enriched_FDR_within_feature", "Evidence")]
  names(p)[-(1:3)] <- c("Sample_a", "Sample_b", "Sample_c", "Sample_d", "Sample_event_rate", "Other_sample_event_rate",
    "Sample_rate_ratio", "Sample_OR", "Sample_P_enriched", "Sample_Q_within_setting_ARG_group", "Sample_evidence")
  tab <- merge(tab, p, by = key, all.x = TRUE, sort = FALSE)
  tab$Direction_concordant <- sign(tab$Target_site_event_rate - tab$Other_site_event_rate) ==
    sign(tab$Sample_event_rate - tab$Other_sample_event_rate)
  tab$Sample_positive_q05 <- tab$Sample_OR > 1 & tab$Sample_Q_within_setting_ARG_group < .05
  winners <- primary$winners[c(key, "Evidence")]
  names(winners)[4] <- "Display_evidence"
  displayed <- merge(winners, tab, by = key, all.x = TRUE, sort = FALSE)
  tab$Selected_in_primary_display <- si789_key(tab, key) %in% si789_key(winners, key)
  tab$Displayed_High_or_Moderate <- si789_key(tab, key) %in%
    si789_key(winners[winners$Display_evidence %in% c("High", "Moderate"), , drop = FALSE], key)
  scopes <- list(All_city_comparisons = seq_len(nrow(tab)),
    Primary_Fisher_positive = which(tab$Sample_positive_q05),
    Displayed_High_or_Moderate = which(tab$Displayed_High_or_Moderate))
  summary <- do.call(rbind, lapply(names(scopes), function(scope) {
    d <- tab[scopes[[scope]], , drop = FALSE]
    direction_n <- sum(!is.na(d$Direction_concordant))
    data.frame(Scope = scope, N = nrow(d), Direction_evaluable = direction_n,
      Same_direction = sum(d$Direction_concordant, na.rm = TRUE),
      Same_direction_percent = if (direction_n > 0) 100 * mean(d$Direction_concordant, na.rm = TRUE) else NA_real_,
      Site_level_q05 = sum(d$Site_positive_q05, na.rm = TRUE))
  }))
  list(profiles = do.call(rbind, profiles), comparisons = tab, displayed_winners = displayed, summary = summary)
}

si245_export_s5 <- function(result, directory) {
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  tables <- list(S5_primary_city_comparisons = result$primary$comparisons,
    S5_primary_displayed_winners = result$primary$winners)
  if (!is.null(result$site)) tables <- c(tables,
    list(S5_site_level_profiles_and_events = result$site$profiles,
      S5_site_level_all_city_comparisons = result$site$comparisons,
      S5_site_level_displayed_winners = result$site$displayed_winners,
      S5_site_level_concordance_summary = result$site$summary))
  paths <- file.path(directory, paste0(names(tables), ".csv"))
  for (i in seq_along(tables)) utils::write.csv(tables[[i]], paths[i], row.names = FALSE, na = "", fileEncoding = "UTF-8")
  invisible(paths)
}

si245_s5 <- function(group_abundance, metadata, seed = 123L, site_sensitivity = TRUE, output_dir = NULL) {
  primary <- si245_s5_tests(group_abundance, metadata, seed)
  result <- list(primary = primary)
  if (site_sensitivity) {
    result$site <- si245_s5_site(group_abundance, metadata, primary)
    result$selected_sensitivity <- result$site$displayed_winners[
      result$site$displayed_winners$Display_evidence %in% c("High", "Moderate"), , drop = FALSE]
  }
  if (!is.null(output_dir)) si245_export_s5(result, output_dir)
  result
}

si245_plot_s5 <- function(result) {
  si245_require("ggplot2")
  d <- result$primary$winners; r <- result$primary$recurrence
  d$Setting <- factor(d$Setting, levels = c("Hospital", "WWTP", "Community", "Wet market"))
  a <- ggplot2::ggplot(d, ggplot2::aes(Setting, Feature, fill = Rate)) + ggplot2::geom_tile(colour = "white") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%s\n%.1f%%, %s", City, 100 * Rate, Evidence)), size = 3) +
    ggplot2::scale_fill_gradient(low = "#F7F2F9", high = "#A47AB8", labels = function(x) paste0(100 * x, "%")) +
    ggplot2::theme_minimal() + ggplot2::labs(x = NULL, y = NULL, fill = "Event rate", tag = "a")
  b <- ggplot2::ggplot(r, ggplot2::aes(City, Feature, size = Supported_settings, fill = Max_rate)) +
    ggplot2::geom_point(shape = 21, colour = "grey25") + ggplot2::geom_text(ggplot2::aes(label = Supported_settings), size = 3) +
    ggplot2::scale_size_area(max_size = 13, breaks = 1:4) +
    ggplot2::scale_fill_gradient(low = "#F7F2F9", high = "#A47AB8", labels = function(x) paste0(100 * x, "%")) +
    ggplot2::theme_bw() + ggplot2::labs(x = NULL, y = NULL, size = "Supported settings", fill = "Maximum event rate", tag = "b")
  list(a = a, b = b)
}

si_root_require <- function(packages) {
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop('Required packages: ',paste(missing,collapse=', '))
}

si_root_input <- function(abundance,metadata) {
  a <- as.matrix(abundance); m <- as.data.frame(metadata,stringsAsFactors=FALSE)
  need <- c('Sample','City','Setting','Date')
  if(!all(need %in% names(m))) stop('metadata requires ',paste(need,collapse=', '))
  if(!is.numeric(a)||is.null(rownames(a))||is.null(colnames(a))||
     anyDuplicated(rownames(a))||anyDuplicated(colnames(a))) stop('Use a numeric matrix with unique sample/subtype names.')
  if(anyNA(a)||any(!is.finite(a))||any(a<0)) stop('Resolve missing abundance values explicitly; require nonnegative finite measurements.')
  if(anyDuplicated(m$Sample)||anyNA(m[need])||!setequal(rownames(a),m$Sample)) stop('Metadata and matrix must have identical unique sample keys.')
  m <- m[match(rownames(a),m$Sample),,drop=FALSE]
  m$Date <- as.Date(m$Date)
  if(anyNA(m$Date)) stop('Date must be Date or YYYY-MM-DD.')
  m$Setting <- as.character(m$Setting)
  m$Setting[m$Setting %in% c('Wet Market','Wet_market')] <- 'Wet market'
  if(any(!m$Setting %in% c('Hospital','WWTP','Community','Wet market'))) stop('Unknown wastewater setting.')
  list(abundance=a,metadata=m)
}

si_root_palette <- function(n) {
  base <- c('#4E79A7','#7FB8B5','#D6B97A','#8BCB88','#E39B76','#B5658D',
    '#B7D989','#8FB6E8','#D17C98','#AFC7E8','#B7E0E5','#E6BAC4','#74C2B3',
    '#7E6BB5','#5F9EA0','#72B77E','#9FD58A','#E4CF62','#75C0C1','#9DBDE0',
    '#9A84D6','#A9C98A','#6DBDE3','#C8A24B','#D8BB92','#E8D98F','#E39A9A',
    '#9AA3A8','#AFC9A0','#7DA34D')
  if(n<=length(base)) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
}

si3_monthly <- function(abundance,metadata,subtype_dictionary,tier1_subtypes=NULL,
                        months=NULL) {
  si_root_require(c('dplyr','tidyr','ggplot2'))
  input <- si_root_input(abundance,metadata); a<-input$abundance;m<-input$metadata
  if (is.null(months)) {
    month_start <- as.Date(format(m$Date, '%Y-%m-01'))
    months <- format(seq(min(month_start), max(month_start), by='month'), '%Y-%m')
  }
  dict<-as.data.frame(subtype_dictionary,stringsAsFactors=FALSE)
  if(!all(c('Subtype','ARG_type') %in% names(dict))||anyDuplicated(dict$Subtype)) stop('Dictionary requires unique Subtype and ARG_type.')
  if(!is.null(tier1_subtypes)) {
    if(!length(tier1_subtypes)||any(!tier1_subtypes %in% colnames(a))) stop('Tier I subtype list must match the abundance columns.')
    a<-a[,colnames(a) %in% tier1_subtypes,drop=FALSE]
  }
  j<-match(colnames(a),dict$Subtype)
  if(anyNA(j)||anyNA(dict$ARG_type[j])) stop('Every analysed subtype needs an ARG type.')
  types<-unique(as.character(dict$ARG_type[j]));m$Month<-format(m$Date,'%Y-%m')
  if(any(!m$Month %in% months)) stop('Samples fall outside requested monthly display interval; subset explicitly.')
  by_type<-sapply(types,function(t) rowSums(a[,dict$ARG_type[j]==t,drop=FALSE]))
  if(is.null(dim(by_type))) by_type<-matrix(by_type,ncol=1,dimnames=list(m$Sample,types))
  long<-data.frame(Sample=rep(m$Sample,times=length(types)),ARG_type=rep(types,each=nrow(a)),Abundance=as.vector(by_type))
  long<-dplyr::left_join(long,m,by='Sample')
  observed<-long |>
    dplyr::group_by(City,Setting,Month,ARG_type) |>
    dplyr::summarise(Mean=mean(Abundance),n=dplyr::n(),.groups='drop')
  grid<-tidyr::expand_grid(City=sort(unique(m$City)),Setting=c('Hospital','WWTP','Community','Wet market'),Month=months,ARG_type=types)
  summary<-dplyr::left_join(grid,observed,by=c('City','Setting','Month','ARG_type'))
  summary$n[is.na(summary$n)]<-0L
  summary$Status<-ifelse(summary$n==0,'unsampled',ifelse(summary$Mean==0,'observed_zero','observed_positive'))
  summary$Month<-factor(summary$Month,levels=months)
  summary$Setting<-factor(summary$Setting,levels=c('Hospital','WWTP','Community','Wet market'))
  colors<-stats::setNames(si_root_palette(length(types)),types)
  plot<-ggplot2::ggplot(summary[summary$n>0,],ggplot2::aes(Month,Mean,fill=ARG_type))+
    ggplot2::geom_col(width=.82)+ggplot2::facet_grid(City~Setting,drop=FALSE)+
    ggplot2::scale_x_discrete(drop=FALSE)+ggplot2::scale_fill_manual(values=colors)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,.06)))+
    ggplot2::labs(x=NULL,y='Mean ARG abundance (copies per cell)',fill='ARG type')+
    ggplot2::theme_bw(base_size=9)+ggplot2::theme(panel.grid=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_text(angle=90,hjust=1,vjust=.5),legend.position='bottom')
  list(summary=as.data.frame(summary),sample_type=long,plot=plot,
       method='Descriptive sample-weighted monthly arithmetic means; absent sampling cells remain missing.')
}

si10_monthly <- function(abundance,metadata,subtype_dictionary,tier1_subtypes,...) {
  si3_monthly(abundance,metadata,subtype_dictionary,tier1_subtypes=tier1_subtypes,...)
}

si6_select_variants <- function(all_setting_stats,top1_setting,top1_city,n_major=8L) {
  si_root_require('dplyr')
  needs<-list(c('Family','Variant','Mean_family_composition'),
              c('Family','Representative_variant'),c('Family','Representative_variant'))
  inputs<-list(all_setting_stats,top1_setting,top1_city)
  for(i in seq_along(inputs)) if(!all(needs[[i]] %in% names(inputs[[i]]))) stop('Selection table is missing required columns.')
  out<-lapply(unique(all_setting_stats$Family),function(family) {
    z<-all_setting_stats[all_setting_stats$Family==family,] |>
      dplyr::group_by(Variant) |>
      dplyr::summarise(Max=max(Mean_family_composition,na.rm=TRUE),.groups='drop') |>
      dplyr::arrange(dplyr::desc(Max),Variant)
    selected<-unique(c(top1_setting$Representative_variant[top1_setting$Family==family],
      top1_city$Representative_variant[top1_city$Family==family],head(z$Variant,n_major)))
    selected<-selected[!is.na(selected)&nzchar(selected)]
    data.frame(Family=family,Variant=selected)
  })
  do.call(rbind,out)
}

si6_family_bars <- function(abundance,metadata,variant_map,selected_variants,
                            families=c('OXA','KPC','VIM','NDM','IMP')) {
  si_root_require(c('dplyr','tidyr','ggplot2'))
  inp<-si_root_input(abundance,metadata);a<-inp$abundance;m<-inp$metadata
  if(!all(c('Subtype','Family','Variant') %in% names(variant_map))||
     anyDuplicated(variant_map$Subtype)||anyNA(variant_map[c('Subtype','Family','Variant')])) stop('Use a complete unique subtype-family-variant map.')
  if(!all(c('Family','Variant') %in% names(selected_variants))) stop('Supply selected variants.')
  map<-variant_map[variant_map$Family %in% families,,drop=FALSE]
  if(any(!map$Subtype %in% colnames(a))) stop('Mapped subtype absent from abundance matrix.')
  data<-plots<-list()
  for(fam in families) {
    fm<-map[map$Family==fam,,drop=FALSE]
    if(!nrow(fm)) stop('No mapping for family ',fam)
    keep<-selected_variants$Variant[selected_variants$Family==fam]
    groups<-ifelse(fm$Variant %in% keep,fm$Variant,'Others')
    display<-unique(groups)
    vals<-lapply(display,function(v) {
      data.frame(m,Variant=v,Abundance=rowSums(a[,fm$Subtype[groups==v],drop=FALSE]))
    })
    summarized<-do.call(rbind,vals) |>
      dplyr::group_by(City,Setting,Variant) |>
      dplyr::summarise(Mean=mean(Abundance),n=dplyr::n(),.groups='drop')
    grid<-tidyr::expand_grid(City=sort(unique(m$City)),Setting=c('Hospital','WWTP','Community','Wet market'),Variant=display)
    z<-dplyr::left_join(grid,summarized,by=c('City','Setting','Variant'))
    z$n[is.na(z$n)]<-0L;z$Family<-fam
    z$Setting<-factor(z$Setting,levels=c('Hospital','WWTP','Community','Wet market'))
    colors<-stats::setNames(si_root_palette(length(display)),display)
    if('Others' %in% display) colors['Others']<-'#D0D0D0'
    missing_cells<-unique(z[z$n==0,c('City','Setting')])
    p<-ggplot2::ggplot(z[z$n>0,],ggplot2::aes(City,Mean,fill=Variant))+
      ggplot2::geom_col(width=.78)+ggplot2::facet_wrap(~Setting,ncol=1,scales='free_y',drop=FALSE)+
      ggplot2::geom_text(data=missing_cells,ggplot2::aes(x=City,y=0,label='NA'),inherit.aes=FALSE,vjust=-.4)+
      ggplot2::scale_fill_manual(values=colors)+ggplot2::scale_x_discrete(limits=sort(unique(m$City)),drop=FALSE)+
      ggplot2::labs(x=NULL,y='Mean abundance (copies per cell)',title=paste0('bla',fam),fill=NULL)+
      ggplot2::theme_bw(base_size=10)+ggplot2::theme(panel.grid=ggplot2::element_blank(),legend.position='bottom')
    data[[fam]]<-as.data.frame(z);plots[[fam]]<-p
  }
  list(data=do.call(rbind,data),plots=plots,selection=selected_variants,
       method='Sample-level variant sums followed by city-setting arithmetic means. No hypothesis test.')
}

si14_profile_features <- function(abundance,metadata,marker_map,tier1_subtypes,bin_days=14L) {
  inp<-si_root_input(abundance,metadata);a<-inp$abundance;m<-inp$metadata
  si_root_require('dplyr')
  if(!all(c('Marker','Subtype') %in% names(marker_map))||anyNA(marker_map)||
     anyDuplicated(marker_map[c('Marker','Subtype')])) stop('marker_map needs unique Marker/Subtype rows.')
  if(any(!marker_map$Subtype %in% colnames(a))||any(!tier1_subtypes %in% colnames(a))||!length(tier1_subtypes)) stop('Mapped NNLS subtypes are absent.')
  settings<-c('Hospital','WWTP','Community','Wet market')
  presence<-table(m$City,factor(m$Setting,levels=settings))
  eligible<-rownames(presence)[apply(presence>0,1,all)]
  excluded<-setdiff(unique(m$City),eligible)
  if(!length(eligible)) stop('No city has samples from all four settings.')
  rows<-m$City %in% eligible;m<-m[rows,,drop=FALSE];a<-a[rows,,drop=FALSE]
  origin<-min(m$Date)
  m$Bin<-as.Date(origin+bin_days*floor(as.numeric(m$Date-origin)/bin_days))
  cells<-unique(m[c('City','Bin')])
  cells$Complete<-vapply(seq_len(nrow(cells)),function(i) all(settings %in% m$Setting[m$City==cells$City[i]&m$Bin==cells$Bin[i]]),logical(1))
  profiles<-list();k<-0L
  add_profiles<-function(indices,cols,marker,city,temporal=FALSE) {
    for(setting in settings) {
      rr<-indices[m$Setting[indices]==setting]
      if(!length(rr)) next
      values<-colMeans(a[rr,cols,drop=FALSE])
      feature<-if(temporal) paste(cols,as.character(m$Bin[rr[1]]),sep='___') else cols
      k<<-k+1L
      profiles[[k]]<<-data.frame(Marker=marker,City=city,Setting=setting,Feature=feature,
        Subtype=cols,Abundance=as.numeric(values),n_samples=length(rr),Kind=if(temporal)'marker' else 'aggregate')
    }
  }
  for(city in eligible) {
    rr<-which(m$City==city)
    add_profiles(rr,colnames(a),'Total ARGs',city)
    add_profiles(rr,tier1_subtypes,'Tier I ARGs',city)
    bins<-cells$Bin[cells$City==city&cells$Complete]
    for(bin in as.character(bins)) for(marker in unique(marker_map$Marker)) {
      ii<-which(m$City==city&as.character(m$Bin)==bin)
      add_profiles(ii,marker_map$Subtype[marker_map$Marker==marker],marker,city,TRUE)
    }
  }
  list(profiles=do.call(rbind,profiles),bin_presence=cells,excluded_cities=excluded,
       bin_origin=origin,marker_map=marker_map)
}

si14_fit_one <- function(profile,min_features=6L) {
  si_root_require('nnls')
  required<-c('Marker','City','Setting','Feature','Abundance','Kind')
  if(!all(required %in% names(profile))||length(unique(profile$Marker))!=1||length(unique(profile$City))!=1) stop('Pass one marker-city profile.')
  if(anyNA(profile[required])||anyDuplicated(profile[c('Feature','Setting')])||any(profile$Abundance<0)) stop('Invalid or duplicate NNLS feature rows.')
  sources<-c('Community','Hospital','Wet market');settings<-c(sources,'WWTP')
  features<-unique(profile$Feature)
  mat<-matrix(0,length(features),4,dimnames=list(features,settings))
  mat[cbind(match(profile$Feature,features),match(profile$Setting,settings))]<-profile$Abundance
  y<-mat[,'WWTP'];X<-mat[,sources,drop=FALSE]
  use<-y>0|rowSums(X)>0;y<-y[use];X<-X[use,,drop=FALSE]
  nonzero<-colSums(X)>0
  result<-data.frame(Marker=profile$Marker[1],City=profile$City[1],Source=sources,
    Coefficient=NA_real_,Percent=NA_real_,fit_r2=NA_real_,fit_bray=NA_real_,
    n_features=length(y),n_nonzero_sources=sum(nonzero),Status='not_fitted')
  if(length(y)<min_features||sum(y)<=0||!any(nonzero)) {
    result$Status<-'insufficient_nonzero_features';return(list(coefficients=result))
  }
  if(sum(nonzero)==1L&&profile$Kind[1]=='marker') {
    result$Coefficient<-as.numeric(nonzero);result$Percent<-100*result$Coefficient
    result$Status<-'single_nonzero_source_no_fit_diagnostics';return(list(coefficients=result))
  }
  if(sum(nonzero)<2L) {result$Status<-'fewer_than_two_nonzero_sources';return(list(coefficients=result))}
  yc<-y/sum(y);xc<-sweep(X[,nonzero,drop=FALSE],2,colSums(X[,nonzero,drop=FALSE]),'/')
  fit<-nnls::nnls(sqrt(xc),sqrt(yc))
  if(sum(fit$x)<=0) {result$Status<-'zero_coefficients';return(list(coefficients=result))}
  weights<-fit$x/sum(fit$x)
  pred<-pmax(as.vector(xc%*%weights),1e-12)
  denom<-sum((yc-mean(yc))^2)
  result$Coefficient<-0;result$Coefficient[nonzero]<-weights
  result$Percent<-100*result$Coefficient
  result$fit_r2<-if(denom>0) 1-sum((yc-pred)^2)/denom else NA_real_
  result$fit_bray<-sum(abs(yc-pred))/sum(yc+pred)
  result$Status<-ifelse(is.finite(result$fit_r2)&result$fit_r2<0,'poor_fit','fitted')
  list(coefficients=result,observed=yc,reconstructed=pred,source_composition=xc,
       raw_coefficients=fit$x,hellinger_RSS=sum((sqrt(yc)-as.vector(sqrt(xc)%*%fit$x))^2))
}

si14_fit <- function(profiles) {
  if(!all(c('Marker','City','Kind') %in% names(profiles))) stop('Use profiles from si14_profile_features().')
  groups<-split(profiles,interaction(profiles$Marker,profiles$City,drop=TRUE,lex.order=TRUE))
  fits<-lapply(groups,function(z) si14_fit_one(z,if(z$Kind[1]=='marker')6L else 10L))
  list(coefficients=do.call(rbind,lapply(fits,`[[`,'coefficients')),fits=fits)
}

si14_plot <- function(coefficients) {
  si_root_require('ggplot2')
  d<-as.data.frame(coefficients)
  if(!all(c('Marker','City','Source','Coefficient','fit_r2') %in% names(d))) stop('Incomplete NNLS coefficient table.')
  if(anyDuplicated(d[c('Marker','City','Source')])) stop('Duplicate source coefficients.')
  keys<-unique(d[c('Marker','City')]);cities<-sort(unique(d$City))
  for(i in seq_len(nrow(keys))) {
    z<-d[d$Marker==keys$Marker[i]&d$City==keys$City[i],]
    if(nrow(z)!=3L||!setequal(z$Source,c('Community','Hospital','Wet market'))) stop('Retain all three sources for every fit.')
    if(!anyNA(z$Coefficient)&&abs(sum(z$Coefficient)-1)>1e-8) stop('NNLS display coefficients must sum to one.')
  }
  d$x<-match(d$City,cities);d$Source<-factor(d$Source,levels=c('Community','Hospital','Wet market'))
  bad<-unique(d[is.finite(d$fit_r2)&d$fit_r2<0,c('Marker','City','x')])
  segments<-list();n<-0L
  for(i in seq_len(nrow(bad))) for(direction in c(-1,1)) for(offset in seq(-.15,1.15,by=.055)) {
    xs<-c(-.4,.4);ys<-offset+direction*.24*xs
    lo<-max(-.4,min((c(0,1)-offset)/(direction*.24)))
    hi<-min(.4,max((c(0,1)-offset)/(direction*.24)))
    if(lo>=hi) next
    n<-n+1L;segments[[n]]<-data.frame(Marker=bad$Marker[i],x=bad$x[i]+lo,xend=bad$x[i]+hi,
      y=offset+direction*.24*lo,yend=offset+direction*.24*hi)
  }
  p<-ggplot2::ggplot(d[!is.na(d$Coefficient),],ggplot2::aes(x,Coefficient,fill=Source))+
    ggplot2::geom_col(width=.8)+ggplot2::facet_wrap(~Marker,ncol=3)+
    ggplot2::scale_fill_manual(values=c(Community='#6EC5C1',Hospital='#E06C75','Wet market'='#E9C95B'),drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=seq_along(cities),labels=cities)+
    ggplot2::scale_y_continuous(limits=c(0,1),labels=function(x)paste0(x*100,'%'),expand=c(0,0))+
    ggplot2::labs(x=NULL,y='Normalized source-profile coefficient',fill='Reference setting')+
    ggplot2::theme_bw(base_size=10)+ggplot2::theme(panel.grid=ggplot2::element_blank(),legend.position='bottom',axis.text.x=ggplot2::element_text(angle=45,hjust=1))
  if(length(segments))p<-p+ggplot2::geom_segment(data=do.call(rbind,segments),
    ggplot2::aes(x=x,xend=xend,y=y,yend=yend),inherit.aes=FALSE,colour='grey35',linewidth=.23)
  p
}

si789_require <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Missing packages: ", paste(absent, collapse = ", "), call. = FALSE)
}

si789_settings <- function() c("Hospital", "WWTP", "Community", "Wet market")
si789_groups <- function() c("four_setting_shared_core", "setting_restricted_core", "non_core")
si789_metrics <- function() c("dominant_host_fraction", "host_species_richness", "host_shannon",
                             "mobile_fraction", "plasmid_fraction", "phage_fraction")

si789_schemas <- function() list(
  annotations = c(sample_id = "character", contig_id = "character", arg_occurrence_id = "character",
                  arg_subtype = "character", arg_type = "character", host_species = "character",
                  host_rank = "character", carrier = "character"),
  metadata = c(sample_id = "character", setting = "character", city = "character",
               month = "character", physical_site = "character"),
  core_membership = c(arg_subtype = "character", setting = "character", is_core = "logical"),
  host_sets = c(species = "character", host_group = "character", source_reference = "character"),
  s8_model = c(arg_subtype = "character", core_group = "character", metric = "character", value = "numeric"),
  s8_plot = c(arg_subtype = "character", core_group = "character", metric = "character", value = "numeric")
)

si789_check <- function(x, schema, name = "table", nullable = character()) {
  if (!is.data.frame(x)) stop(name, " must be a data.frame.")
  missing <- setdiff(names(schema), names(x))
  if (length(missing)) stop(name, " lacks required columns: ", paste(missing, collapse = ", "))
  for (key in names(schema)) {
    ok <- switch(schema[[key]], character = is.character(x[[key]]),
                 numeric = is.numeric(x[[key]]), logical = is.logical(x[[key]]), FALSE)
    if (!ok) stop(name, "$", key, " must be ", schema[[key]], ".")
    if (!key %in% nullable && anyNA(x[[key]])) stop(name, "$", key, " contains missing values.")
    if (is.character(x[[key]]) && any(!is.na(x[[key]]) & !nzchar(trimws(x[[key]]))))
      stop(name, "$", key, " contains empty identifiers.")
  }
  invisible(x)
}

si789_key <- function(x, columns) {
  do.call(paste0, lapply(x[columns], function(z) paste0(nchar(z, type = "chars"), ":", z)))
}

si789_metadata <- function(metadata) {
  si789_check(metadata, si789_schemas()$metadata, "metadata")
  if (anyDuplicated(metadata$sample_id)) stop("Metadata must contain one row per sample.")
  if (any(!metadata$setting %in% si789_settings())) stop("Unknown wastewater setting.")
  if (any(!grepl("^[0-9]{4}-(0[1-9]|1[0-2])$", metadata$month))) stop("month must be YYYY-MM.")
  sites <- unique(metadata[c("physical_site", "setting", "city")])
  if (anyDuplicated(sites$physical_site)) stop("A physical site maps to more than one setting or city.")
  metadata
}

si789_annotations <- function(annotations, metadata = NULL) {
  si789_check(annotations, si789_schemas()$annotations, "annotations", c("host_species", "host_rank"))
  a <- annotations
  if (!nrow(a)) stop("Annotations contain no ARG occurrences.")
  if (any(!a$carrier %in% c("chromosome", "plasmid", "phage", "unknown")))
    stop("carrier must be chromosome, plasmid, phage or unknown; explicitly map source 'virus' to phage.")
  if (anyDuplicated(si789_key(a, c("sample_id", "arg_occurrence_id"))))
    stop("Duplicate ARG occurrence identifiers; do not double-count input annotation rows.")
  a$contig_uid <- si789_key(a, c("sample_id", "contig_id"))
  cm <- unique(a[c("contig_uid", "carrier", "host_species", "host_rank")])
  if (anyDuplicated(cm$contig_uid)) stop("Contradictory carrier or host annotation for the same sample-contig.")
  if (any(a$host_rank == "S" & is.na(a$host_species), na.rm = TRUE))
    stop("Species-rank assignments require a host species name.")
  if (!is.null(metadata)) {
    m <- si789_metadata(metadata)
    i <- match(a$sample_id, m$sample_id)
    if (anyNA(i)) stop("Annotations contain samples outside the supplied metadata cohort.")
    for (key in setdiff(names(si789_schemas()$metadata), "sample_id")) a[[key]] <- m[[key]][i]
  }
  a
}

si789_prepare_s7 <- function(annotations, metadata) {
  a <- si789_annotations(annotations, metadata)
  if (any(a$carrier == "unknown")) stop("S7 requires resolved three-compartment assignments; unknown carriers cannot be silently dropped.")
  units <- unique(a[c("contig_uid", "sample_id", "contig_id", "arg_type", "carrier", "city")])
  count <- function(columns) {
    out <- stats::aggregate(rep(1L, nrow(units)), units[columns], sum)
    names(out)[ncol(out)] <- "n_contigs"
    out
  }
  totals <- count("arg_type")
  carriers <- count(c("arg_type", "carrier"))
  cities <- count(c("arg_type", "city"))
  for (kind in c("carriers", "cities")) {
    z <- get(kind)
    z$total_contigs <- totals$n_contigs[match(z$arg_type, totals$arg_type)]
    z$fraction <- z$n_contigs / z$total_contigs
    assign(kind, z)
  }
  list(units = units, carriers = carriers, totals = totals, cities = cities,
       audit = data.frame(unique_contigs = length(unique(units$contig_uid)),
                          contig_type_units = nrow(units), annotation_occurrences = nrow(a)))
}

si789_plot_s7 <- function(prepared, city_colours = NULL) {
  si789_require("ggplot2")
  order <- prepared$totals$arg_type[order(prepared$totals$n_contigs, decreasing = TRUE)]
  z <- lapply(prepared[c("carriers", "totals", "cities")], function(d) {
    d$arg_type <- factor(d$arg_type, levels = order); d
  })
  common <- ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 60, hjust = 1))
  a <- ggplot2::ggplot(z$carriers, ggplot2::aes(arg_type, fraction, fill = carrier)) +
    ggplot2::geom_col() + ggplot2::scale_fill_manual(values = c(chromosome = "#8C8C8C", plasmid = "#E39A54", phage = "#70A5CB")) +
    ggplot2::labs(x = NULL, y = "Fraction of unique contigs", tag = "a", fill = "Carrier") + common
  b <- ggplot2::ggplot(z$totals, ggplot2::aes(arg_type, n_contigs)) + ggplot2::geom_col(fill = "#4E79A7") +
    ggplot2::scale_y_log10() + ggplot2::labs(x = NULL, y = "Unique contigs (log scale)", tag = "b") + common
  c <- ggplot2::ggplot(z$cities, ggplot2::aes(arg_type, fraction, fill = city)) + ggplot2::geom_col() +
    ggplot2::labs(x = NULL, y = "Fraction of unique contigs", tag = "c", fill = "City") + common
  if (!is.null(city_colours)) c <- c + ggplot2::scale_fill_manual(values = city_colours)
  list(a = a, b = b, c = c)
}

si789_core_groups <- function(core_membership) {
  si789_check(core_membership, si789_schemas()$core_membership, "core_membership")
  d <- core_membership
  if (!nrow(d) || any(!d$setting %in% si789_settings()) ||
      anyDuplicated(si789_key(d, c("arg_subtype", "setting")))) stop("Invalid or duplicate core membership cells.")
  if (any(table(d$arg_subtype) != 4L)) stop("Each subtype needs explicit membership flags for all four settings.")
  z <- stats::aggregate(d$is_core, list(arg_subtype = d$arg_subtype), sum)
  names(z)[2] <- "n_core_settings"
  z$core_group <- ifelse(z$n_core_settings == 4L, si789_groups()[1],
                         ifelse(z$n_core_settings > 0L, si789_groups()[2], si789_groups()[3]))
  z
}

si789_core_from_abundance <- function(abundance, metadata) {
  m <- si789_metadata(metadata)
  if (!is.matrix(abundance) || !is.numeric(abundance) || is.null(rownames(abundance)) ||
      is.null(colnames(abundance)) || anyDuplicated(rownames(abundance)) || anyDuplicated(colnames(abundance)) ||
      !setequal(rownames(abundance), m$sample_id) || any(!is.finite(abundance)) || any(abundance < 0))
    stop("abundance must be a complete nonnegative numeric sample-by-subtype matrix with unique names.")
  if (!setequal(m$setting, si789_settings())) stop("Core classification requires all four settings.")
  abundance <- abundance[match(m$sample_id, rownames(abundance)), , drop = FALSE]
  do.call(rbind, lapply(si789_settings(), function(s) {
    a <- abundance[m$setting == s, , drop = FALSE]
    prevalence <- colMeans(a > 0); mean_abundance <- colMeans(a)
    data.frame(arg_subtype = colnames(a), setting = s,
               is_core = prevalence >= .70 & mean_abundance >= 1e-5,
               prevalence = prevalence, mean_abundance = mean_abundance, row.names = NULL)
  }))
}

si789_metric_values <- function(units) {
  assigned <- units$host_rank == "S" & !is.na(units$host_rank) & !is.na(units$host_species)
  h <- table(units$host_species[assigned])
  p <- if (length(h)) as.numeric(h) / sum(h) else numeric()
  carrier <- units$carrier[units$carrier %in% c("chromosome", "plasmid", "phage")]
  c(dominant_host_fraction = if (length(p)) max(p) else NA_real_,
    host_species_richness = length(h), host_shannon = if (length(p)) -sum(p * log(p)) else NA_real_,
    mobile_fraction = if (length(carrier)) mean(carrier %in% c("plasmid", "phage")) else NA_real_,
    plasmid_fraction = if (length(carrier)) mean(carrier == "plasmid") else NA_real_,
    phage_fraction = if (length(carrier)) mean(carrier == "phage") else NA_real_)
}

si789_prepare_s8 <- function(annotations, metadata, core_membership) {
  a <- si789_annotations(annotations, metadata)
  groups <- si789_core_groups(core_membership)
  if (any(!a$arg_subtype %in% groups$arg_subtype)) stop("Missing core classification for annotated ARG subtype.")
  units <- a[!duplicated(si789_key(a, c("contig_uid", "arg_subtype"))), , drop = FALSE]
  summarize <- function(columns) {
    groups_index <- split(seq_len(nrow(units)), si789_key(units, columns))
    do.call(rbind, lapply(groups_index, function(i) {
      values <- si789_metric_values(units[i, , drop = FALSE])
      ids <- units[rep(i[1], length(values)), columns, drop = FALSE]
      data.frame(ids, metric = names(values), value = unname(values), row.names = NULL)
    }))
  }
  plot <- summarize("arg_subtype")
  plot$core_group <- groups$core_group[match(plot$arg_subtype, groups$arg_subtype)]
  plot <- plot[c("arg_subtype", "core_group", "metric", "value")]
  list(model_input = plot, plot_data = plot, core_groups = groups,
       audit = data.frame(annotation_occurrences = nrow(a), subtype_contig_units = nrow(units),
                          unknown_carrier_units = sum(units$carrier == "unknown"),
                          mobility_denominator = "Subtype-carrying contigs assigned to chromosome, plasmid or phage",
                          host_denominator = "Species-rank assigned subtype-carrying contigs"))
}

si789_s8_source_plot_table <- function(traits) {
  fields <- c(dominant_host_fraction = "dominant_host_fraction", host_species_richness = "host_species_richness",
              host_shannon = "host_shannon", mobile_fraction = "mobile_fraction_all",
              plasmid_fraction = "plasmid_fraction_all", phage_fraction = "virus_fraction_all")
  si789_check(traits, c(Subtype = "character", core_group = "character",
                        stats::setNames(rep("numeric", length(fields)), unname(fields))),
              "source trait table", unname(fields))
  if (anyDuplicated(traits$Subtype)) stop("Source plot table has duplicate subtypes.")
  d <- do.call(rbind, lapply(names(fields), function(metric) data.frame(
    arg_subtype = traits$Subtype, core_group = traits$core_group, metric = metric,
    value = traits[[fields[[metric]]]], stringsAsFactors = FALSE)))
  si789_validate_s8_values(d)
  d
}

si789_validate_s8_values <- function(d) {
  if (any(!d$metric %in% si789_metrics())) stop("Unknown S8 metric.")
  if ("core_group" %in% names(d) && any(!d$core_group %in% si789_groups())) stop("Unknown core group.")
  if (any(!is.na(d$value) & !is.finite(d$value)) || any(d$value < 0, na.rm = TRUE)) stop("Invalid S8 metric value.")
  fraction <- d$metric %in% c("dominant_host_fraction", "mobile_fraction", "plasmid_fraction", "phage_fraction")
  if (any(d$value[fraction] > 1, na.rm = TRUE)) stop("S8 fractions must lie in [0,1].")
  richness <- d$value[d$metric == "host_species_richness"]
  if (any(richness != round(richness), na.rm = TRUE)) stop("Host species richness must be an integer.")
  invisible(d)
}

si789_host_sets <- function(host_sets) {
  si789_check(host_sets, si789_schemas()$host_sets, "host_sets")
  if (anyDuplicated(host_sets$species)) stop("Host species groups overlap or contain duplicate species.")
  if (!setequal(host_sets$host_group, c("clinical_pathogen", "opportunistic_environment")))
    stop("Provide both predefined host groups using clinical_pathogen and opportunistic_environment.")
  host_sets
}

si789_prepare_s9 <- function(annotations, metadata, host_sets) {
  a <- si789_annotations(annotations, metadata)
  m <- si789_metadata(metadata); hs <- si789_host_sets(host_sets)
  eligible <- !is.na(a$host_rank) & a$host_rank == "S" & !is.na(a$host_species)
  species <- a[eligible, , drop = FALSE]
  species$host_group <- hs$host_group[match(species$host_species, hs$species)]
  species$host_group[is.na(species$host_group)] <- "other_species"
  occurrences <- species
  contigs <- species[!duplicated(species$contig_uid), , drop = FALSE]
  calc <- function(d, counting) {
    denominator <- tabulate(match(d$sample_id, m$sample_id), nbins = nrow(m))
    do.call(rbind, lapply(c("clinical_pathogen", "opportunistic_environment"), function(g) {
      numerator <- tabulate(match(d$sample_id[d$host_group == g], m$sample_id), nbins = nrow(m))
      data.frame(sample_id = m$sample_id, counting = counting, host_group = g,
                 numerator = numerator, denominator = denominator,
                 fraction = ifelse(denominator > 0, numerator / denominator, NA_real_))
    }))
  }
  fractions <- rbind(calc(occurrences, "ARG_occurrence"), calc(contigs, "unique_contig"))
  fractions$panel <- c(a = "ARG_occurrence:clinical_pathogen", b = "ARG_occurrence:opportunistic_environment",
                       c = "unique_contig:clinical_pathogen", d = "unique_contig:opportunistic_environment") |>
    (function(map) names(map)[match(paste(fractions$counting, fractions$host_group, sep = ":"), map)])()
  list(fractions = fractions, metadata = m, host_sets = hs,
       audit = data.frame(all_occurrences = nrow(a), species_assigned_occurrences = nrow(species),
                          species_assigned_contigs = nrow(contigs), excluded_non_species_occurrences = sum(!eligible),
                          same_subtype_same_contig_additional_occurrences = sum(duplicated(si789_key(species, c("contig_uid", "arg_subtype"))))))
}

si789_lmm <- function(d, focal, levels, expected_pairs) {
  si789_require(c("lme4", "multcompView"))
  for (key in c(focal, "setting", "city", "month", "physical_site")) d[[key]] <- factor(d[[key]])
  d[[focal]] <- factor(d[[focal]], levels = levels)
  if (!setequal(as.character(d[[focal]]), levels)) stop("All comparison groups are required for the specified BH family.")
  fixed <- unique(c(focal, "setting", "city", "month"))
  if (any(vapply(d[fixed], nlevels, integer(1)) < 2L)) stop("Each fixed effect requires at least two observed levels.")
  if (nlevels(d$physical_site) < 2L || nlevels(d$physical_site) >= nrow(d)) stop("Random-site grouping requires multiple sites and repeated observations.")
  form <- stats::as.formula(paste("value ~", paste(fixed, collapse = " + "), "+ (1 | physical_site)"))
  warnings_seen <- messages_seen <- character()
  fit <- withCallingHandlers(lme4::lmer(form, d, REML = TRUE, na.action = stats::na.fail,
                         control = lme4::lmerControl(optimizer = "bobyqa")),
    warning = function(w) { warnings_seen <<- c(warnings_seen, conditionMessage(w)); invokeRestart("muffleWarning") },
    message = function(m) { messages_seen <<- c(messages_seen, conditionMessage(m)); invokeRestart("muffleMessage") })
  if (length(attr(lme4::getME(fit, "X"), "col.dropped"))) stop("Rank-deficient fixed-effects design; inference needs review.")
  beta <- lme4::fixef(fit); covariance <- as.matrix(stats::vcov(fit))
  nd <- d[rep(1L, length(levels)), , drop = FALSE]; nd[[focal]] <- factor(levels, levels = levels)
  X <- stats::model.matrix(stats::as.formula(paste("~", paste(fixed, collapse = " + "))), nd)[, names(beta), drop = FALSE]
  contrasts <- do.call(rbind, lapply(utils::combn(seq_along(levels), 2L, simplify = FALSE), function(pair) {
    L <- X[pair[1], ] - X[pair[2], ]; estimate <- sum(L * beta)
    se <- sqrt(drop(L %*% covariance %*% L)); z <- estimate / se
    data.frame(group1 = levels[pair[1]], group2 = levels[pair[2]], estimate = estimate, SE = se,
               Wald_z = z, CI_low = estimate - stats::qnorm(.975) * se,
               CI_high = estimate + stats::qnorm(.975) * se, p_value = 2 * stats::pnorm(-abs(z)))
  }))
  if (nrow(contrasts) != expected_pairs || any(!is.finite(contrasts$p_value))) stop("Incomplete or non-estimable contrast family.")
  contrasts$q_BH <- stats::p.adjust(contrasts$p_value, method = "BH")
  pmat <- matrix(1, length(levels), length(levels), dimnames = list(levels, levels))
  for (i in seq_len(nrow(contrasts))) pmat[contrasts$group1[i], contrasts$group2[i]] <-
    pmat[contrasts$group2[i], contrasts$group1[i]] <- contrasts$q_BH[i]
  ordered <- levels[order(drop(X %*% beta), decreasing = TRUE)]
  cld <- multcompView::multcompLetters(pmat[ordered, ordered], threshold = .05)$Letters
  for (i in seq_len(nrow(contrasts))) {
    shared <- length(intersect(strsplit(cld[contrasts$group1[i]], "")[[1]],
                               strsplit(cld[contrasts$group2[i]], "")[[1]])) > 0L
    if (shared != (contrasts$q_BH[i] >= .05)) stop("Compact-letter display contradicts BH-adjusted comparisons.")
  }
  list(fit = fit, contrasts = contrasts, letters = data.frame(group = names(cld), Letters = unname(cld)),
       diagnostics = data.frame(formula = paste(deparse(form), collapse = " "), transform = "identity",
         REML = TRUE, inference = "Two-sided asymptotic Wald z; unadjusted 95% CI; within-outcome BH",
         BH_family_size = expected_pairs, n_observations = nrow(d), n_samples = length(unique(d$sample_id)),
         n_sites = nlevels(d$physical_site), singular = lme4::isSingular(fit),
         min_fitted = min(stats::fitted(fit)), max_fitted = max(stats::fitted(fit)),
         convergence = paste(unlist(fit@optinfo$conv$lme4$messages), collapse = "; "),
         warnings = paste(unique(warnings_seen), collapse = "; "), messages = paste(unique(messages_seen), collapse = "; ")))
}

si789_fit_s8 <- function(model_input) {
  si789_require("multcompView")
  si789_check(model_input, si789_schemas()$s8_model, "subtype-level S8 input", "value")
  si789_validate_s8_values(model_input)
  if ("sample_id" %in% names(model_input)) stop("S8 requires pooled subtype metrics, not sample-by-subtype observations.")
  if (!setequal(model_input$metric, si789_metrics())) stop("S8 input requires all six metrics.")
  if (anyDuplicated(si789_key(model_input, c("arg_subtype", "metric"))))
    stop("S8 requires one aggregate value per subtype and metric.")
  if (anyDuplicated(unique(model_input[c("arg_subtype", "core_group")])$arg_subtype))
    stop("Each subtype must belong to exactly one core group.")
  d <- model_input
  result <- lapply(si789_metrics(), function(metric) {
    q <- d[d$metric == metric & !is.na(d$value), , drop = FALSE]
    levels <- si789_groups()
    if (!setequal(q$core_group, levels)) stop("All three core groups require defined values for metric: ", metric)
    q$core_group <- factor(q$core_group, levels = levels)
    n <- nrow(q)
    ranks <- rank(q$value, ties.method = "average")
    group_n <- table(q$core_group)
    mean_ranks <- tapply(ranks, q$core_group, mean)
    ties <- as.numeric(table(q$value))
    rank_variance <- max(0, n * (n + 1) / 12 - sum(ties^3 - ties) / (12 * (n - 1)))
    global <- if (rank_variance > 0) stats::kruskal.test(value ~ core_group, data = q) else NULL
    contrasts <- do.call(rbind, lapply(utils::combn(levels, 2L, simplify = FALSE), function(pair) {
      difference <- unname(mean_ranks[pair[1]] - mean_ranks[pair[2]])
      se <- sqrt(rank_variance * sum(1 / group_n[pair]))
      z <- if (se > 0) difference / se else 0
      data.frame(group1 = pair[1], group2 = pair[2], n1 = unname(group_n[pair[1]]),
        n2 = unname(group_n[pair[2]]), mean_rank_difference = difference, SE = se,
        Dunn_z = z, p_value = 2 * stats::pnorm(-abs(z)), stringsAsFactors = FALSE)
    }))
    contrasts$q_BH <- stats::p.adjust(contrasts$p_value, "BH")
    pmat <- matrix(1, 3L, 3L, dimnames = list(levels, levels))
    for (i in seq_len(nrow(contrasts))) pmat[contrasts$group1[i], contrasts$group2[i]] <-
      pmat[contrasts$group2[i], contrasts$group1[i]] <- contrasts$q_BH[i]
    ordered <- levels[order(mean_ranks[levels], decreasing = TRUE)]
    cld <- multcompView::multcompLetters(pmat[ordered, ordered], threshold = .05)$Letters
    for (i in seq_len(nrow(contrasts))) {
      shared <- length(intersect(strsplit(cld[contrasts$group1[i]], "")[[1]],
        strsplit(cld[contrasts$group2[i]], "")[[1]])) > 0L
      if (shared != (contrasts$q_BH[i] >= .05)) stop("Compact-letter display contradicts Dunn/BH comparisons.")
    }
    list(global = data.frame(metric = metric, statistic = if (is.null(global)) 0 else unname(global$statistic),
      df = 2L, p_value = if (is.null(global)) 1 else global$p.value),
      contrasts = contrasts, letters = data.frame(group = names(cld), Letters = unname(cld)),
      diagnostics = data.frame(metric = metric, transform = "identity", inference = "Kruskal-Wallis; two-sided tie-corrected Dunn; BH within metric",
        BH_family_size = 3L, n_subtypes = n, missing_subtypes = sum(d$metric == metric & is.na(d$value)),
        random_effects = "none", covariates = "none", stringsAsFactors = FALSE))
  })
  names(result) <- si789_metrics()
  list(results = result, model_input = d,
       excluded_missing = d[is.na(d$value), c("arg_subtype", "metric"), drop = FALSE])
}

si789_run_s8 <- function(annotations = NULL, metadata = NULL, core_membership = NULL,
                         traits = NULL, output_dir = NULL) {
  if (!is.null(traits) && !is.null(annotations)) stop("Supply either pooled subtype traits or annotation records, not both.")
  prepared <- if (!is.null(traits)) {
    d <- if (all(names(si789_schemas()$s8_plot) %in% names(traits))) traits else si789_s8_source_plot_table(traits)
    list(plot_data = d)
  } else si789_prepare_s8(annotations, metadata, core_membership)
  fitted <- si789_fit_s8(prepared$plot_data)
  plots <- si789_plot_s8(prepared$plot_data, fitted, richness_display = "log10p1")
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    tables <- list(S8_subtype_metrics = prepared$plot_data,
      S8_Kruskal_Wallis = do.call(rbind, lapply(fitted$results, `[[`, "global")),
      S8_Dunn_BH = do.call(rbind, lapply(names(fitted$results), function(metric)
        data.frame(metric = metric, fitted$results[[metric]]$contrasts))),
      S8_letters = do.call(rbind, lapply(names(fitted$results), function(metric)
        data.frame(metric = metric, fitted$results[[metric]]$letters))),
      S8_excluded_undefined_metrics = fitted$excluded_missing)
    for (name in names(tables)) utils::write.csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")),
      row.names = FALSE, na = "", fileEncoding = "UTF-8")
    si789_require("patchwork")
    combined <- patchwork::wrap_plots(plots, ncol = 3L)
    si_save_plots(combined, output_dir, "S8", width = 14, height = 8)
  }
  list(prepared = prepared, fitted = fitted, plots = plots)
}

si789_fit_s9 <- function(fractions, metadata) {
  si789_check(fractions, c(sample_id = "character", panel = "character", numerator = "numeric",
                           denominator = "numeric", fraction = "numeric"), "S9 fractions", "fraction")
  if (!setequal(fractions$panel, letters[1:4]) || anyDuplicated(si789_key(fractions, c("sample_id", "panel"))))
    stop("S9 needs four panels with unique sample-panel rows.")
  if (any(!is.finite(fractions$denominator)) || any(!is.finite(fractions$numerator)) ||
      any(fractions$denominator < 0 | fractions$numerator < 0 | fractions$numerator > fractions$denominator) ||
      any(fractions$denominator != round(fractions$denominator) | fractions$numerator != round(fractions$numerator)))
    stop("Invalid S9 counts.")
  positive <- fractions$denominator > 0
  if (any(!is.na(fractions$fraction[!positive])) || anyNA(fractions$fraction[positive]) ||
      any(abs(fractions$fraction[positive] - fractions$numerator[positive] / fractions$denominator[positive]) > 1e-12))
    stop("S9 fractions disagree with their counts or undefined zero denominators.")
  m <- si789_metadata(metadata); d <- fractions; i <- match(d$sample_id, m$sample_id)
  if (anyNA(i)) stop("S9 contains samples outside metadata.")
  for (key in setdiff(names(si789_schemas()$metadata), "sample_id")) d[[key]] <- m[[key]][i]
  d$value <- d$fraction
  result <- lapply(letters[1:4], function(panel) {
    q <- d[d$panel == panel & positive, , drop = FALSE]
    r <- si789_lmm(q, "setting", si789_settings(), 6L)
    r$diagnostics$panel <- panel
    r$diagnostics$zero_denominator_samples <- sum(d$panel == panel & !positive)
    r
  })
  names(result) <- letters[1:4]
  list(results = result, model_input = d,
       excluded_zero_denominator = d[!positive, c("sample_id", "panel", "denominator"), drop = FALSE])
}

si789_half_violin <- function(d, groups, colours, y_label, cld = NULL, show_mean = FALSE) {
  si789_require("ggplot2")
  d <- d[!is.na(d$value), , drop = FALSE]
  d$x <- match(d$group, groups)
  if (!nrow(d) || anyNA(d$x)) stop("No usable points or unknown plot group.")
  polygons <- lapply(seq_along(groups), function(i) {
    v <- d$value[d$x == i]
    if (length(v) < 2L || diff(range(v)) == 0) return(NULL)
    den <- stats::density(v, n = 256, from = min(v), to = max(v))
    data.frame(group = groups[i], x = c(i + .42 * den$y / max(den$y), rep(i, length(den$x))),
               value = c(den$x, rev(den$x)))
  })
  polygons <- do.call(rbind, polygons)
  p <- ggplot2::ggplot()
  if (!is.null(polygons)) p <- p + ggplot2::geom_polygon(data = polygons,
    ggplot2::aes(x, value, group = group, fill = group), alpha = .65, colour = NA) +
    ggplot2::scale_fill_manual(values = colours)
  p <- p + ggplot2::geom_point(data = d, ggplot2::aes(x = x - .17, y = value, colour = group),
             position = ggplot2::position_jitter(width = .12, height = 0, seed = 197), size = .7, alpha = .55) +
    ggplot2::geom_boxplot(data = d, ggplot2::aes(x, value, group = group), width = .15,
                         fill = "white", outlier.shape = NA, linewidth = .4) +
    ggplot2::scale_x_continuous(breaks = seq_along(groups), labels = groups) +
    ggplot2::scale_colour_manual(values = colours) +
    ggplot2::labs(x = NULL, y = y_label) + ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), legend.position = "none",
                   axis.text.x = ggplot2::element_text(angle = 25, hjust = 1))
  span <- diff(range(d$value)); if (!is.finite(span) || span == 0) span <- 1
  top <- max(d$value) + .12 * span
  if (show_mean) {
    means <- stats::aggregate(d$value, list(group = d$group), mean); names(means)[2] <- "value"
    means$x <- match(means$group, groups); means$label <- sprintf("%.3f", means$value)
    p <- p + ggplot2::geom_text(data = means, ggplot2::aes(x = x, label = label), y = top, size = 3)
  }
  if (!is.null(cld)) {
    if (!all(groups %in% cld$group) || anyDuplicated(cld$group)) stop("Incomplete or duplicate model letters.")
    cld$x <- match(cld$group, groups)
    p <- p + ggplot2::geom_text(data = cld, ggplot2::aes(x = x, label = Letters),
                               y = top + if (show_mean) .1 * span else 0, size = 3)
  }
  p
}

si789_plot_s8 <- function(plot_data, fitted = NULL, richness_display = c("log10p1", "original")) {
  si789_check(plot_data, si789_schemas()$s8_plot, "S8 plot data", "value")
  si789_validate_s8_values(plot_data)
  if (anyDuplicated(si789_key(plot_data, c("arg_subtype", "metric")))) stop("S8 plot must have one point per subtype per metric.")
  richness_display <- match.arg(richness_display)
  colours <- stats::setNames(c("#E06C75", "#4E79A7", "#BDBDBD"), si789_groups())
  plots <- lapply(si789_metrics(), function(metric) {
    q <- plot_data[plot_data$metric == metric, , drop = FALSE]; q$group <- q$core_group
    label <- gsub("_", " ", metric)
    if (metric == "host_species_richness" && richness_display == "log10p1") {
      q$value <- log10(q$value + 1); label <- "Host species richness, log10(x+1)"
    }
    cld <- if (is.null(fitted)) NULL else fitted$results[[metric]]$letters
    si789_half_violin(q, si789_groups(), colours, label, cld)
  })
  stats::setNames(plots, si789_metrics())
}

si789_plot_s9 <- function(fractions, metadata, fitted = NULL) {
  m <- si789_metadata(metadata)
  si789_check(fractions, c(sample_id = "character", panel = "character", fraction = "numeric"), "S9 plot data", "fraction")
  fractions$group <- m$setting[match(fractions$sample_id, m$sample_id)]
  if (anyNA(fractions$group) || any(fractions$fraction < 0 | fractions$fraction > 1, na.rm = TRUE)) stop("Invalid S9 plot data.")
  labels <- c(a = "Clinical pathogen ARG-occurrence fraction", b = "Opportunistic/environment ARG-occurrence fraction",
              c = "Clinical pathogen unique-contig fraction", d = "Opportunistic/environment unique-contig fraction")
  colours <- stats::setNames(c("#E06C75", "#4E79A7", "#6EC5C1", "#F3A64A"), si789_settings())
  plots <- lapply(names(labels), function(panel) {
    q <- fractions[fractions$panel == panel, , drop = FALSE]; q$value <- q$fraction
    cld <- if (is.null(fitted)) NULL else fitted$results[[panel]]$letters
    si789_half_violin(q, si789_settings(), colours, labels[[panel]], cld, show_mean = TRUE) + ggplot2::labs(tag = panel)
  })
  stats::setNames(plots, names(labels))
}

si111213_need <- function(data, columns) {
  absent <- setdiff(columns, names(data))
  if (length(absent)) stop("Missing required columns: ", paste(absent, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

si111213_packages <- function(packages) {
  absent <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(absent)) stop("Required packages are not installed: ", paste(absent, collapse = ", "), call. = FALSE)
  invisible(TRUE)
}

si111213_mean <- function(x) {
  if (!length(x) || all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
}

si111213_keys <- function(data, columns, unique = TRUE) {
  si111213_need(data, columns)
  if (anyNA(data[columns]) || any(vapply(data[columns], function(x) any(!nzchar(as.character(x))), logical(1))))
    stop("Keys must be complete and nonempty: ", paste(columns, collapse = ", "), call. = FALSE)
  if (unique && anyDuplicated(data[columns])) stop("Duplicate rows at input unit: ", paste(columns, collapse = "/"), call. = FALSE)
  invisible(TRUE)
}

si111213_numeric <- function(data, columns, nonnegative = FALSE) {
  si111213_need(data, columns)
  for (column in columns) {
    x <- data[[column]]
    if (!is.numeric(x) || any(is.infinite(x)) || (nonnegative && any(x < 0, na.rm = TRUE)))
      stop("Invalid numeric values in ", column, call. = FALSE)
  }
  invisible(TRUE)
}

si111213_read <- function(path, sheet = NULL) {
  if (!file.exists(path)) stop("Input file does not exist: ", path, call. = FALSE)
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") return(utils::read.csv(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (ext %in% c("tsv", "txt")) return(utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE))
  if (ext == "rds") {
    out <- readRDS(path)
    if (!is.data.frame(out)) stop("RDS input must be a data frame.", call. = FALSE)
    return(out)
  }
  if (ext %in% c("xlsx", "xls")) {
    si111213_packages("readxl")
    if (is.null(sheet)) stop("Supply the input workbook sheet explicitly.", call. = FALSE)
    return(as.data.frame(readxl::read_excel(path, sheet = sheet)))
  }
  stop("Supported inputs: CSV, TSV, data-frame RDS or explicit Excel sheet.", call. = FALSE)
}

si111213_s11_activity <- function(activity,
    expected = NULL) {
  activity <- as.data.frame(activity)
  si111213_keys(activity, c("city", "hospital", "quarter"))
  si111213_keys(activity, c("hospital", "quarter"))
  si111213_numeric(activity, c("patient_days", "discharges"), nonnegative = TRUE)
  hospital_city <- unique(activity[c("hospital", "city")])
  if (anyDuplicated(hospital_city$hospital)) stop("One hospital maps to multiple cities.", call. = FALSE)
  observed <- c(hospital_quarters = nrow(activity), hospitals = length(unique(activity$hospital)),
                cities = length(unique(activity$city)))
  if (!is.null(expected) && any(observed[names(expected)] != expected))
    stop("Observed cohort dimensions differ from expected.", call. = FALSE)
  activity$patient_days_10000 <- activity$patient_days / 10000
  activity$length_of_stay <- ifelse(!is.na(activity$discharges) & activity$discharges > 0,
                                  activity$patient_days / activity$discharges, NA_real_)
  rows <- lapply(unique(as.character(activity$city)), function(city) {
    do.call(rbind, lapply(c("patient_days_10000", "length_of_stay"), function(metric) {
      value <- activity[[metric]][activity$city == city]
      n <- sum(!is.na(value))
      data.frame(city = city, metric = metric, mean = si111213_mean(value),
        SEM = if (n > 1L) stats::sd(value, na.rm = TRUE) / sqrt(n) else NA_real_, n = n,
        n_hospitals = length(unique(activity$hospital[activity$city == city])), stringsAsFactors = FALSE)
    }))
  })
  list(observations = activity, summary = do.call(rbind, rows), cohort = observed,
       weighting = "Equal hospital-quarter observations; SEM uses observed hospital-quarter count")
}

si111213_ddd_intensity <- function(ddd, activity) {
  si111213_keys(ddd, c("city", "hospital", "quarter", "route", "category"), unique = FALSE)
  si111213_keys(activity, c("city", "hospital", "quarter"))
  si111213_numeric(ddd, "ddd", nonnegative = TRUE)
  si111213_numeric(activity, "patient_days", nonnegative = TRUE)
  if (any(!ddd$route %in% c("oral", "injection"))) stop("Raw DDD routes must be oral or injection.", call. = FALSE)
  keys <- c("city", "hospital", "quarter", "route", "category")
  summed <- stats::aggregate(ddd["ddd"], ddd[keys], function(x) if (anyNA(x)) NA_real_ else sum(x))
  out <- merge(summed, activity[c("city", "hospital", "quarter", "patient_days")],
               by = c("city", "hospital", "quarter"), all.x = TRUE, sort = FALSE)
  if (nrow(out) != nrow(summed)) stop("Activity join multiplied DDD rows.", call. = FALSE)
  out$intensity <- ifelse(!is.na(out$patient_days) & out$patient_days > 0,
                          out$ddd / out$patient_days * 100, NA_real_)
  route_keys <- c("city", "hospital", "quarter", "category")
  oral <- out[out$route == "oral", c(route_keys, "intensity"), drop = FALSE]
  injection <- out[out$route == "injection", c(route_keys, "intensity"), drop = FALSE]
  combined <- merge(oral, injection, by = route_keys, all = TRUE, suffixes = c("_oral", "_injection"), sort = FALSE)
  combined$intensity <- combined$intensity_oral + combined$intensity_injection
  combined$route <- "combined"
  rbind(out[c(keys, "intensity")], combined[c(keys, "intensity")])
}

si111213_s11_amu <- function(intensities) {
  d <- as.data.frame(intensities)
  si111213_keys(d, c("city", "hospital", "quarter", "route", "category"))
  si111213_numeric(d, "intensity", nonnegative = TRUE)
  if (any(!d$route %in% c("oral", "injection", "combined"))) stop("Unknown administration route.", call. = FALSE)
  if (anyDuplicated(unique(d[c("hospital", "city")])$hospital)) stop("Hospital city mapping is inconsistent.", call. = FALSE)
  groups <- c("city", "hospital", "route", "category")
  hospital <- stats::aggregate(d["intensity"], d[groups], si111213_mean)
  counts <- stats::aggregate(list(valid_quarters = !is.na(d$intensity)), d[groups], sum)
  hospital <- merge(hospital, counts, by = groups, sort = FALSE)
  groups_city <- c("city", "route", "category")
  city <- stats::aggregate(hospital["intensity"], hospital[groups_city], si111213_mean)
  n <- stats::aggregate(list(n_hospitals = !is.na(hospital$intensity)), hospital[groups_city], sum)
  city <- merge(city, n, by = groups_city, sort = FALSE)
  total <- stats::aggregate(city["intensity"], city[c("city", "route")],
    function(x) if (anyNA(x)) NA_real_ else sum(x))
  list(hospital_means = hospital, city_contributions = city, city_totals = total,
       weighting = "Arithmetic mean of valid quarters in hospital; equal hospital means in city; explicit zeros retained")
}

si111213_plot_s11 <- function(activity_summary, amu_summary) {
  si111213_packages(c("ggplot2", "patchwork"))
  act <- activity_summary$summary
  amu <- amu_summary$city_contributions
  totals <- amu_summary$city_totals
  city_order <- unique(as.character(act$city))
  act$city <- factor(act$city, levels = city_order)
  amu$city <- factor(amu$city, levels = city_order)
  totals$city <- factor(totals$city, levels = city_order)
  if (anyNA(amu$city)) stop("AMU cities do not match the activity cohort.", call. = FALSE)
  theme <- ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  activity_plots <- lapply(c("patient_days_10000", "length_of_stay"), function(metric) {
    d <- act[act$metric == metric, , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(city, mean)) +
      ggplot2::geom_col(fill = "#6F95B8", width = .7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - SEM, ymax = mean + SEM), width = .18, na.rm = TRUE) +
      ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", mean), y = mean + ifelse(is.na(SEM), 0, SEM)), vjust = -.4, size = 3) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .17))) +
      ggplot2::labs(x = NULL, y = if (metric == "patient_days_10000") "Inpatient days (10,000 patient-days)" else "Length of stay (days)") + theme
  })
  amu_plots <- lapply(c("oral", "injection", "combined"), function(route) {
    d <- amu[amu$route == route, , drop = FALSE]
    tt <- totals[totals$route == route, , drop = FALSE]
    ggplot2::ggplot(d, ggplot2::aes(city, intensity, fill = category)) +
      ggplot2::geom_col(width = .7) +
      ggplot2::geom_text(data = tt, ggplot2::aes(city, intensity, label = sprintf("%.2f", intensity)), inherit.aes = FALSE, vjust = -.4, size = 3) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, .15))) +
      ggplot2::labs(x = NULL, y = "DDDs/100 patient-days", fill = NULL,
                    title = c(oral = "Oral", injection = "Injection", combined = "Oral + injection")[[route]]) + theme
  })
  panels <- c(activity_plots, amu_plots)
  names(panels) <- letters[1:5]
  list(panels = panels, figure = patchwork::wrap_plots(panels, ncol = 2, guides = "collect") +
         patchwork::plot_annotation(tag_levels = "a"))
}

si111213_config <- function() {
  if (!exists("fig6_network_config", mode = "function")) stop("Load the bundled Fig6 network helpers first.", call. = FALSE)
  config <- fig6_network_config()
  config$classes[["Physicochemical parameter"]] <- c("CODCr", "NH3_N", "TSS", "pH")
  config$labels <- c(config$labels, CODCr = "CODCr", NH3_N = "NH3-N", TSS = "TSS", pH = "pH",
                     flow = "Flow", wastewater_temperature = "Wastewater temperature")
  config$s12_predictors <- c("baidu.immigration", "Phenolic_compounds", "Aminoglycoside_class", "Macrolide_class", "pH",
    "baidu.emigration", "NSAIDs", "beta_Lactam_antibiotics", "Tetracycline_class", "CODCr")
  config$s13_network_predictors <- c("Sulfonamide_class_and_synergists", "Peptide_class", "Antihistamine_drugs", "PM10", "temp_mean",
    "Quinolone_class", "Other_types_of_antibiotics", "non_NSAIDs", "PM2.5", "prec_sum")
  config$s13_wwtp_predictors <- c("Sulfonamide_class_and_synergists", "Peptide_class", "Antihistamine_drugs", "PM10", "temp_mean", "NH3_N",
    "Quinolone_class", "Other_types_of_antibiotics", "non_NSAIDs", "PM2.5", "prec_sum", "TSS")

  config
}

si111213_quality_daily <- function(quality) {
  d <- as.data.frame(quality)
  si111213_keys(d, c("city", "date", "Index"), unique = FALSE)
  si111213_numeric(d, "value")
  d$date <- as.Date(d$date)
  if (anyNA(d$date)) stop("Unparseable quality dates.", call. = FALSE)
  allowed <- c("CODCr", "NH3_N", "TSS", "pH", "flow", "wastewater_temperature")
  if (any(!d$Index %in% allowed)) stop("Unknown quality Index; harmonize to the documented six canonical names.", call. = FALSE)
  d$value[d$Index == "flow" & !is.na(d$value) & d$value == 0] <- NA_real_
  stats::aggregate(d["value"], d[c("city", "date", "Index")], si111213_mean)
}

si111213_quality_lags <- function(samples, quality, config = si111213_config()) {
  si111213_keys(samples, c("city", "sample_date"), unique = FALSE)
  out <- as.data.frame(samples)
  out$sample_date <- as.Date(out$sample_date)
  if (anyNA(out$sample_date)) stop("Unparseable sample dates.", call. = FALSE)
  q <- si111213_quality_daily(quality)
  for (index in config$classes[["Physicochemical parameter"]]) {
    d <- q[q$Index == index, , drop = FALSE]
    for (lag in config$lags) {
      days <- unname(config$lag_days[lag])
      out[[paste0(index, "_", lag)]] <- vapply(seq_len(nrow(out)), function(i) {
        delta <- as.numeric(out$sample_date[i] - d$date)
        keep <- d$city == out$city[i] & delta >= 0 & delta <= if (days == 0L) 90L else days
        if (!any(keep)) return(NA_real_)
        if (days == 0L) d$value[which(keep)[which.min(delta[keep])]] else si111213_mean(d$value[keep])
      }, numeric(1))
    }
  }
  attr(out, "quality_city_date") <- q
  out
}

si111213_prepare_wwtp <- function(samples, daily, quality, config = si111213_config()) {
  si111213_keys(samples, "sample_id")
  si111213_need(samples, c("city", "site_nm", "setting", "sample_date", "season", config$outcome))
  si111213_numeric(samples, config$outcome, nonnegative = TRUE)
  si111213_need(daily, c("city", "date", "Index", "value"))
  daily_predictors <- unlist(config$classes[names(config$classes) != "Physicochemical parameter"], use.names = FALSE)
  absent <- setdiff(daily_predictors, unique(as.character(daily$Index)))
  if (length(absent)) stop("Missing entire daily predictor series: ", paste(absent, collapse = ", "), call. = FALSE)
  d <- samples[!is.na(samples$setting) & samples$setting == "WWTP", , drop = FALSE]
  if (!nrow(d)) stop("No WWTP samples supplied.", call. = FALSE)
  d$site_nm <- interaction(d$city, d$site_nm, drop = TRUE, sep = "::")
  d <- fig6_network_lags(d, daily[daily$Index %in% daily_predictors, , drop = FALSE], config)
  d <- si111213_quality_lags(d, quality, config)
  expected <- as.vector(outer(unlist(config$classes, use.names = FALSE), config$lags, paste, sep = "_"))
  si111213_need(d, expected)
  excluded <- sum(is.na(d[[config$outcome]]))
  d <- d[!is.na(d[[config$outcome]]), , drop = FALSE]
  if (!nrow(d)) stop("No WWTP samples have observed ARG outcomes.", call. = FALSE)
  prepared <- fig6_network_prepare_lmm(d, config)
  attr(prepared, "si111213_missing_outcomes_excluded") <- excluded
  attr(prepared, "si111213_scope") <- "WWTP"
  prepared
}

si111213_prepare_wwtp_rf <- function(prepared, config = si111213_config(), impute = TRUE) {
  if (!isTRUE(attr(prepared, "fig6_minmax")) || isTRUE(attr(prepared, "fig6_rf_imputed")))
    stop("RF requires the 0-1-scaled, pre-kNN input table.", call. = FALSE)
  predictors <- unlist(config$classes, use.names = FALSE)
  columns <- paste0(predictors, "_lag0")
  si111213_need(prepared, c(columns, "city", "sample_date", config$outcome))
  keep <- !is.na(prepared[[config$outcome]])
  out <- data.frame(city = as.factor(prepared$city[keep]),
                    time = as.factor(format(as.Date(prepared$sample_date[keep]), "%Y-%m")))
  out[[config$outcome]] <- prepared[[config$outcome]][keep]
  out[predictors] <- prepared[keep, columns, drop = FALSE]
  if (anyNA(out[c("city", "time")])) stop("RF context keys must be complete.", call. = FALSE)
  before <- vapply(out, function(x) sum(is.na(x)), integer(1))
  if (any(vapply(out[predictors], function(x) all(is.na(x)), logical(1))))
    stop("An entire lag0 predictor is unavailable; it cannot be imputed or silently dropped.", call. = FALSE)
  if (impute && anyNA(out)) {
    si111213_packages("VIM")
    out <- VIM::kNN(out, k = 5L, numFun = stats::median, imp_var = FALSE)
  }
  attr(out, "fig6_rf_imputed") <- impute && any(before > 0)
  attr(out, "si111213_rf_ready") <- !anyNA(out)
  attr(out, "fig6_rf_preprocessing") <- data.frame(Predictor = names(out), Missing_before = before,
    Missing_after = vapply(out, function(x) sum(is.na(x)), integer(1)))
  out
}

si111213_fit_wwtp_lmm <- function(prepared, outcome_label, config = si111213_config()) {
  if (isTRUE(attr(prepared, "fig6_rf_imputed"))) stop("Use pre-kNN data for LMMs.", call. = FALSE)
  expected <- as.vector(outer(unlist(config$classes, use.names = FALSE), config$lags, paste, sep = "_"))
  si111213_need(prepared, expected)
  fig6_network_fit_lmm(prepared, outcome_label, config)
}

si111213_fit_wwtp_rf <- function(prepared_rf, outcome_label, config = si111213_config()) {
  si111213_packages(c("rfPermute", "randomForest"))
  if (!isTRUE(attr(prepared_rf, "si111213_rf_ready")) || anyNA(prepared_rf)) stop("Complete the explicit WWTP RF preprocessing first.", call. = FALSE)
  if (config$ntree != 1000L || config$num.rep != 1000L) stop("S12 requires 1000 trees and 1000 permutations.", call. = FALSE)
  set.seed(config$seed)
  fit <- fig6_rf_permute(stats::reformulate(setdiff(names(prepared_rf), config$outcome), config$outcome),
    data = prepared_rf, importance = TRUE, ntree = 1000L, num.rep = 1000L, num.cores = config$num.cores)
  importance <- as.data.frame(randomForest::importance(fit, scale = TRUE), check.names = FALSE)
  importance$variable <- rownames(importance)
  table <- si111213_rf_table(importance, outcome_label, config)
  list(model = fit, importance = importance, plot_data = table,
    class_importance = si111213_rf_classes(table),
    metadata = data.frame(Outcome = outcome_label, N = nrow(prepared_rf), Trees = 1000L, Permutations = 1000L,
      Seed = config$seed, Scale_importance = TRUE, Significance = "nominal permutation P"),
    preprocessing = attr(prepared_rf, "fig6_rf_preprocessing"))
}

si111213_rf_table <- function(importance, outcome_label, config = si111213_config()) {
  d <- as.data.frame(importance, check.names = FALSE)
  si111213_need(d, "%IncMSE")
  key_col <- intersect(c("Predictor", "Index", "variable"), names(d))
  keys <- if (length(key_col)) as.character(d[[key_col[1]]]) else rownames(d)
  if (anyNA(keys) || any(!nzchar(keys)) || anyDuplicated(keys)) stop("RF predictor identifiers must be unique.", call. = FALSE)
  si111213_numeric(d, "%IncMSE")
  score <- d[["%IncMSE"]]
  if (anyNA(score) || sum(score) <= 0) stop("Signed RF importance sum must be positive for percent interpretation.", call. = FALSE)
  p_column <- intersect(c("%IncMSE.pval", "P", "p"), names(d))
  if (!length(p_column)) stop("Nominal permutation P values are required.", call. = FALSE)
  p <- d[[p_column[1]]]
  if (!is.numeric(p) || any(!is.na(p) & (!is.finite(p) | p < 0 | p > 1))) stop("Invalid RF permutation P values.", call. = FALSE)
  data.frame(Predictor = keys, Label = fig6_network_label(keys, config), Class = fig6_network_class(keys, config),
    Raw_importance = score, Importance = score / sum(score) * 100, P = p,
    Stars = fig6_network_stars(p), Outcome = outcome_label, stringsAsFactors = FALSE)
}

si111213_rf_classes <- function(table) {
  si111213_need(table, c("Class", "Importance"))
  out <- stats::aggregate(list(Percent = pmax(table$Importance, 0)), table["Class"], sum)
  total <- sum(out$Percent)
  if (!is.finite(total) || total <= 0) stop("Class percentages require a finite positive importance total.", call. = FALSE)
  out$Percent <- out$Percent / total * 100
  out
}

si111213_lmm_panels <- function(results, scope = c("wwtp", "network"), config = si111213_config()) {
  scope <- match.arg(scope)
  si111213_need(results, c("outcome", "class", "group", "lag", "coef", "ci_low", "ci_high", "p", "fdr"))
  si111213_keys(results, c("outcome", "group", "lag"))
  selected <- if (scope == "wwtp") config$s12_predictors else config$panel_predictors
  additional <- if (scope == "wwtp") config$s13_wwtp_predictors else config$s13_network_predictors
  universe <- unlist(config$classes[if (scope == "network") names(config$classes) != "Physicochemical parameter" else rep(TRUE, length(config$classes))], use.names = FALSE)
  if (length(intersect(selected, additional)) || !setequal(c(selected, additional), universe))
    stop("Panel selection does not form an exact complementary partition.", call. = FALSE)
  if (any(!results$group %in% universe)) stop("Unknown predictors in LMM results.", call. = FALSE)
  table <- fig6_network_lmm_table(results, config = config)
  first <- table[table$Predictor %in% selected, , drop = FALSE]
  second <- table[table$Predictor %in% additional, , drop = FALSE]
  first$Label <- factor(first$Label, levels = fig6_network_label(selected, config))
  second$Label <- factor(second$Label, levels = fig6_network_label(additional, config))
  list(selected = first, additional = second, all = table)
}

si111213_plot_lmm <- function(table, title = NULL, ncol = 5L) {
  si111213_packages("ggplot2")
  si111213_need(table, c("Predictor", "Label", "Lag", "Estimate", "CI_low", "CI_high", "Q"))
  si111213_keys(table, c("Predictor", "Lag"))
  if (!nrow(table)) stop("No additional/displayed estimates available for this panel.", call. = FALSE)
  if (anyNA(table[c("Estimate", "CI_low", "CI_high", "Q")])) stop("Plot requires complete model coefficients, intervals and Q values.", call. = FALSE)
  table$Direction <- ifelse(table$Q < .05 & table$Estimate > 0, "Positive",
    ifelse(table$Q < .05 & table$Estimate < 0, "Negative", "Not significant"))
  table$Stars <- ifelse(table$Q < .05, "*", "")
  table$star_y <- table$CI_high + .08 * pmax(table$CI_high - table$CI_low, .01)
  ggplot2::ggplot(table, ggplot2::aes(Lag, Estimate, colour = Direction)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = CI_low, ymax = CI_high), width = .17) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(y = star_y, label = Stars), colour = "black", size = 4) +
    ggplot2::facet_wrap(~Label, ncol = ncol, scales = "free_y", drop = TRUE) +
    ggplot2::scale_colour_manual(values = c(Positive = "#ca7579", Negative = "#4E79A7", "Not significant" = "grey55")) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.08, .18))) +
    ggplot2::labs(title = title, x = NULL, y = "Coefficient (0-1 scaled)", colour = NULL) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(), strip.background = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "bottom")
}

si111213_plot_rf <- function(table, title = NULL) {
  si111213_packages("ggplot2")
  si111213_need(table, c("Class", "Label", "Importance", "Stars"))
  if (!nrow(table)) stop("RF plot table is empty.", call. = FALSE)
  d <- table[order(table$Class, -table$Importance), , drop = FALSE]
  d$id <- seq_len(nrow(d))
  height <- log10(pmax(d$Importance, 0) + 1)
  span <- max(height)
  if (!is.finite(span) || span == 0) stop("Cannot plot empty importance variation.", call. = FALSE)
  baseline <- .82
  d$ymin <- baseline
  d$ymax <- baseline + height
  d$label <- paste0(d$Label, "\n", sprintf("%.1f%%", d$Importance), d$Stars)
  d$angle <- 90 - 360 * (d$id - .5) / nrow(d)
  d$hjust <- ifelse(d$angle < -90, 1, 0)
  d$angle <- ifelse(d$angle < -90, d$angle + 180, d$angle)
  classes <- si111213_rf_classes(d)
  classes$xmin <- vapply(classes$Class, function(x) min(d$id[d$Class == x]) - .5, numeric(1))
  classes$xmax <- vapply(classes$Class, function(x) max(d$id[d$Class == x]) + .5, numeric(1))
  ring <- max(d$ymax) + 1.1 * span
  classes$label <- paste0(classes$Class, "\n", sprintf("%.1f%%", classes$Percent))
  palette <- c("Meteorological factors" = "#AFCB86", "Air pollution" = "#8F9EAA", "Digital prescriptions" = "#B494C8",
    "Population migration" = "#6F95B8", "Physicochemical parameter" = "#C2A665", "Context factors" = "#C97868")
  ggplot2::ggplot() +
    ggplot2::geom_rect(data = d, ggplot2::aes(xmin = id - .5, xmax = id + .5, ymin = ymin, ymax = ymax, fill = Class), colour = "white", linewidth = .25) +
    ggplot2::geom_hline(yintercept = baseline, linetype = "dashed", colour = "grey30", linewidth = .4) +
    ggplot2::geom_text(data = d, ggplot2::aes(x = id, y = ymax + .06 * span, label = label, angle = angle, hjust = hjust), size = 2.4) +
    ggplot2::geom_rect(data = classes, ggplot2::aes(xmin = xmin, xmax = xmax, fill = Class), ymin = ring, ymax = ring + .18 * span, colour = "white") +
    ggplot2::geom_text(data = classes, ggplot2::aes(x = (xmin + xmax)/2, label = label), y = ring + .42 * span, size = 2.8) +
    ggplot2::scale_fill_manual(values = palette) +
    ggplot2::scale_x_continuous(limits = c(.5, nrow(d) + .5), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(limits = c(0, ring + .8 * span)) +
    ggplot2::coord_polar(theta = "x", clip = "off") +
    ggplot2::labs(title = title, caption = "Bar height: log10(max(relative importance, 0) + 1). Predictor labels retain signed percentages; class shares sum to 100%.") +
    ggplot2::theme_void(base_size = 10) + ggplot2::theme(legend.position = "none")
}

si111213_plot_s12 <- function(rf_total, rf_tier1, lmm_total, lmm_tier1, config = si111213_config()) {
  si111213_packages("patchwork")
  panels <- list(a = si111213_plot_rf(rf_total, "Total ARGs: WWTP"), b = si111213_plot_rf(rf_tier1, "Tier I ARGs: WWTP"),
    c = si111213_plot_lmm(si111213_lmm_panels(lmm_total, "wwtp", config)$selected, "Total ARGs: WWTP"),
    d = si111213_plot_lmm(si111213_lmm_panels(lmm_tier1, "wwtp", config)$selected, "Tier I ARGs: WWTP"))
  list(panels = panels, figure = patchwork::wrap_plots(panels, ncol = 2, guides = "collect") + patchwork::plot_annotation(tag_levels = "a"))
}

si111213_plot_s13 <- function(network_total, network_tier1, wwtp_total, wwtp_tier1, config = si111213_config()) {
  si111213_packages("patchwork")
  data <- list(network_total, network_tier1, wwtp_total, wwtp_tier1)
  scopes <- c("network", "network", "wwtp", "wwtp")
  titles <- c("Total ARGs: all wastewater", "Tier I ARGs: all wastewater", "Total ARGs: WWTP", "Tier I ARGs: WWTP")
  panels <- lapply(seq_along(data), function(i) si111213_plot_lmm(si111213_lmm_panels(data[[i]], scopes[i], config)$additional, titles[i], ncol = 6L))
  names(panels) <- letters[1:4]
  list(panels = panels, figure = patchwork::wrap_plots(panels, ncol = 1, guides = "collect") + patchwork::plot_annotation(tag_levels = "a"))
}

si111213_save <- function(plot, path, width, height, dpi = 300) {
  si111213_packages("ggplot2")
  if (!dir.exists(dirname(path))) dir.create(dirname(path), recursive = TRUE)
  ggplot2::ggsave(filename = path, plot = plot, width = width, height = height, units = "in", dpi = dpi, bg = "white")
  invisible(path)
}

si_packages <- function() c("dplyr", "tidyr", "ggplot2", "patchwork", "vegan", "permute",
                           "lme4", "lmerTest", "multcompView", "randomForest",
                           "rfPermute", "VIM", "nnls", "readxl")

si_check_packages <- function() {
  p <- si_packages()
  data.frame(Package = p, Installed = vapply(p, requireNamespace, logical(1), quietly = TRUE))
}

si_install_packages <- function(repos = "https://cloud.r-project.org") {
  status <- si_check_packages()
  missing <- status$Package[!status$Installed]
  if (length(missing)) utils::install.packages(missing, repos = repos)
  invisible(si_check_packages())
}

si_metadata <- function(metadata) {
  required <- c("sample_id", "city", "setting", "sample_date", "physical_site")
  if (!all(required %in% names(metadata))) stop("metadata requires: ", paste(required, collapse = ", "))
  m <- as.data.frame(metadata, stringsAsFactors = FALSE)
  m[required] <- lapply(m[required], as.character)
  if (anyNA(m[required]) || anyDuplicated(m$sample_id)) stop("Metadata identifiers must be complete and sample IDs unique.")
  if ("site_nm" %in% names(m)) {
    m$site_nm <- as.character(m$site_nm)
    if (anyNA(m$site_nm) || any(m$site_nm != m$physical_site))
      stop("site_nm must match physical_site for every sample.")
  } else m$site_nm <- m$physical_site
  m$sample_date <- as.Date(m$sample_date)
  if (anyNA(m$sample_date)) stop("sample_date must be YYYY-MM-DD or Date.")
  if (any(!m$setting %in% c("Hospital", "WWTP", "Community", "Wet market"))) stop("Use canonical setting names.")
  site_map <- unique(m[c("physical_site", "city", "setting")])
  if (anyDuplicated(site_map$physical_site)) stop("Physical-site IDs must be globally unique and map to one city/setting.")
  m$month <- format(m$sample_date, "%Y-%m")
  list(lower = m, upper = data.frame(Sample = m$sample_id, City = m$city,
    Setting = m$setting, Date = m$sample_date, PhysicalSite = m$physical_site,
    Month = m$month, stringsAsFactors = FALSE))
}

si_save_plots <- function(plots, directory, prefix, width = 12, height = 8) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Install ggplot2 first.")
  if (!dir.exists(directory)) dir.create(directory, recursive = TRUE)
  if (inherits(plots, "ggplot") || inherits(plots, "patchwork")) plots <- list(figure = plots)
  if (is.null(names(plots))) names(plots) <- seq_along(plots)
  paths <- character()
  for (key in names(plots)) {
    if (!inherits(plots[[key]], "ggplot") && !inherits(plots[[key]], "patchwork")) next
    path <- file.path(directory, paste0(prefix, "_", key, ".pdf"))
    ggplot2::ggsave(path, plots[[key]], width = width, height = height, bg = "white")
    paths <- c(paths, path)
  }
  invisible(paths)
}
