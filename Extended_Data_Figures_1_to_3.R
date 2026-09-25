
ed_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop('Install required packages: ', paste(missing, collapse = ', '))
}

ed_columns <- function(x, required, label) {
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(label, ' is missing columns: ', paste(missing, collapse = ', '))
  if (anyDuplicated(names(x))) stop(label, ' has duplicate column names.')
}

ed_read_table <- function(path, sheet = 1) {
  if (!file.exists(path)) stop('Input not found: ', path)
  extension <- tolower(tools::file_ext(path))
  if (extension %in% c('xlsx', 'xls')) {
    ed_require('readxl')
    return(as.data.frame(readxl::read_excel(path, sheet = sheet, col_types = 'text', .name_repair = 'minimal')))
  }
  if (!extension %in% c('csv', 'tsv', 'txt')) stop('Use a CSV, TSV or Excel input: ', path)
  utils::read.table(path, header = TRUE, sep = if (extension == 'csv') ',' else '\t',
    quote = '"', comment.char = '', check.names = FALSE, stringsAsFactors = FALSE,
    na.strings = c('', 'NA'), colClasses = 'character', fileEncoding = 'UTF-8-BOM')
}

ed_write_table <- function(x, path) {
  utils::write.table(x, path, sep = '\t', quote = FALSE, row.names = FALSE, na = 'NA',
                     fileEncoding = 'UTF-8')
}

ed_numeric <- function(x, label, allow_na = FALSE) {
  original <- x
  if (!is.numeric(x)) x <- suppressWarnings(as.numeric(as.character(x)))
  bad_text <- is.na(x) & !is.na(original) & trimws(as.character(original)) != ''
  if (any(bad_text)) stop(label, ' contains nonnumeric values.')
  if (any(!is.finite(x) & !is.na(x)) || (!allow_na && anyNA(x)))
    stop(label, ' contains missing or nonfinite values.')
  x
}

ed_family_names <- function(x) {
  x <- toupper(trimws(as.character(x)))
  sub('^BLA', '', x)
}

ed_settings <- function() c('Hospital', 'WWTP', 'Community', 'Wet market')

ed_setting_colours <- function() {
  c(Hospital = '#ED6A78', WWTP = '#4F7FB0', Community = '#69C6C3', 'Wet market' = '#EACA52')
}

ed_prepare_cohort <- function(abundance, metadata, families,
                              blank_policy = c('stop', 'nondetection'),
                              orientation = c('sample_rows', 'subtype_rows')) {
  blank_policy <- match.arg(blank_policy)
  orientation <- match.arg(orientation)
  required <- c('sample_id', 'setting', 'city', 'month', 'physical_site')
  ed_columns(metadata, required, 'Metadata')
  metadata <- as.data.frame(metadata[required], stringsAsFactors = FALSE)
  for (key in required) {
    metadata[[key]] <- trimws(as.character(metadata[[key]]))
    if (anyNA(metadata[[key]]) || any(metadata[[key]] == '')) stop('Missing metadata: ', key)
  }
  if (!nrow(metadata) || anyDuplicated(metadata$sample_id)) stop('Metadata sample IDs must be unique and nonempty.')
  metadata$setting[metadata$setting == 'Wet Market'] <- 'Wet market'
  if (any(!metadata$setting %in% ed_settings())) stop('Unknown wastewater setting.')
  site_map <- unique(metadata[c('physical_site', 'setting', 'city')])
  if (anyDuplicated(site_map$physical_site)) stop('A physical site has conflicting city or setting assignments.')
  id_column <- if (orientation == 'sample_rows') 'sample_id' else 'subtype'
  ed_columns(abundance, id_column, 'Abundance matrix')
  ids <- trimws(as.character(abundance[[id_column]]))
  if (!length(ids) || anyNA(ids) || any(ids == '') || anyDuplicated(ids)) stop('Abundance row identifiers must be unique.')
  values <- abundance[setdiff(names(abundance), id_column)]
  if (!ncol(values)) stop('Abundance matrix has no numeric columns.')
  values[] <- lapply(seq_along(values), function(i) ed_numeric(values[[i]], names(values)[i], allow_na = TRUE))
  x <- as.matrix(values)
  rownames(x) <- ids
  if (orientation == 'subtype_rows') x <- t(x)
  if (any(!is.finite(x) & !is.na(x)) || any(x < 0, na.rm = TRUE)) stop('Abundance must be finite and nonnegative.')
  if (anyNA(x)) {
    if (blank_policy == 'stop') stop('Missing abundance cells: explicitly set blank_policy = nondetection only when confirmed.')
    x[is.na(x)] <- 0
  }
  if (!setequal(rownames(x), metadata$sample_id)) stop('Abundance sample IDs and metadata sample IDs differ.')
  x <- x[match(metadata$sample_id, rownames(x)), , drop = FALSE]
  ed_columns(families, c('subtype', 'family'), 'Subtype-family map')
  families <- as.data.frame(families[c('subtype', 'family')], stringsAsFactors = FALSE)
  families$subtype <- trimws(as.character(families$subtype))
  families$family <- ed_family_names(families$family)
  if (anyNA(families) || any(families$subtype == '') || any(families$family == '') ||
      anyDuplicated(families$subtype)) stop('Subtype-family assignments must be complete and unique.')
  if (!all(families$subtype %in% colnames(x))) stop('Mapped family subtypes are absent from the abundance matrix.')
  list(abundance = x, metadata = metadata, families = families)
}

ed_selection <- function(x, mode = c('selected', 'screen'),
                         families = c('OXA', 'KPC', 'VIM', 'NDM', 'IMP'),
                         top_n = 20, columns = NULL) {
  mode <- match.arg(mode)
  if (!is.null(columns)) {
    if (is.null(names(columns)) || any(!unname(columns) %in% names(x))) stop('Invalid selection column map.')
    for (key in names(columns)) x[[key]] <- x[[columns[[key]]]]
  }
  aliases <- c(family = 'Family', subtype = 'Subtype', plot_order = 'Rank', display_label = 'Display_label')
  for (key in names(aliases)) if (!key %in% names(x) && aliases[[key]] %in% names(x)) x[[key]] <- x[[aliases[[key]]]]
  ed_columns(x, c('family', 'subtype'), 'Subtype selection')
  x$family <- ed_family_names(x$family)
  x$subtype <- trimws(as.character(x$subtype))
  x <- x[x$family %in% families, , drop = FALSE]
  if (!nrow(x) || anyNA(x$subtype) || any(x$subtype == '') || anyDuplicated(x$subtype))
    stop('Selected subtype identities must be nonempty and unique.')
  if (mode == 'selected') {
    ed_columns(x, 'plot_order', 'Original selected-subtype table')
    x$plot_order <- ed_numeric(x$plot_order, 'plot_order')
    if (any(x$plot_order <= 0 | x$plot_order != round(x$plot_order)) ||
        anyDuplicated(paste(x$family, x$plot_order, sep = ':'))) stop('Invalid original plotting order.')
    x <- x[order(match(x$family, families), x$plot_order), , drop = FALSE]
  } else {
    ed_columns(x, c('setting_r2', 'city_r2', 'setting_q'), 'Initial screen')
    for (key in c('setting_r2', 'city_r2', 'setting_q')) x[[key]] <- ed_numeric(x[[key]], key, allow_na = TRUE)
    if (any(x$setting_q < 0 | x$setting_q > 1, na.rm = TRUE)) stop('Initial adjusted P values must be in [0,1].')
    if (length(top_n) != 1 || !is.finite(top_n) || top_n < 1 || top_n != round(top_n)) stop('Invalid top_n.')
    eligible <- is.finite(x$setting_r2) & is.finite(x$city_r2) & is.finite(x$setting_q) & x$setting_q < .05
    x <- x[eligible, , drop = FALSE]
    if (!nrow(x)) stop('No eligible subtypes in the declared initial screen.')
    x <- do.call(rbind, lapply(families, function(fam) {
      d <- x[x$family == fam, , drop = FALSE]
      if (!nrow(d)) return(NULL)
      d <- d[order(-d$setting_r2, d$setting_q, d$subtype), , drop = FALSE]
      d <- head(d, top_n)
      d$plot_order <- seq_len(nrow(d))
      d
    }))
  }
  if (!'display_label' %in% names(x)) x$display_label <- x$subtype
  x$display_label <- as.character(x$display_label)
  missing_labels <- is.na(x$display_label) | trimws(x$display_label) == ''
  x$display_label[missing_labels] <- x$subtype[missing_labels]
  rownames(x) <- NULL
  x[c('family', 'subtype', 'plot_order', 'display_label')]
}

ed1_summary <- function(data, min_n = 10) {
  ed_columns(data, c('sample_id', 'country', 'region', 'source', 'abundance'), 'Global benchmark')
  d <- as.data.frame(data, stringsAsFactors = FALSE)
  for (key in c('sample_id', 'country', 'region', 'source')) {
    d[[key]] <- trimws(as.character(d[[key]]))
    if (anyNA(d[[key]]) || any(d[[key]] == '')) stop('Missing benchmark field: ', key)
  }
  d$source <- tolower(d$source)
  if (any(!d$source %in% c('study', 'reference'))) stop('Benchmark source must be study or reference.')
  if (any(d$source == 'reference' & tolower(d$country) %in% c('this study', 'present study', 'our study')))
    stop('Reference records cannot use a reserved study-group label.')
  if (anyDuplicated(d[c('source', 'sample_id')])) stop('Duplicate benchmark sample IDs within source.')
  d$abundance <- ed_numeric(d$abundance, 'Benchmark abundance')
  if (any(d$abundance < 0)) stop('Benchmark abundance must be nonnegative.')
  if (length(min_n) != 1 || !is.finite(min_n) || min_n < 1 || min_n != round(min_n)) stop('Invalid minimum display count.')
  d$country[d$country %in% c('China', 'China Mainland', 'China mainland', 'Mainland China')] <- 'Mainland China'
  d$country[d$country %in% c('Hong Kong', 'Hong Kong SAR', 'China Hong Kong SAR')] <- 'China Hong Kong SAR'
  d$country[d$country %in% c('Taiwan', 'China Taiwan')] <- 'China Taiwan'
  d$country[d$source == 'study'] <- 'This study'
  d$region[d$source == 'study'] <- 'This study'
  region_map <- unique(d[c('country', 'region')])
  if (anyDuplicated(region_map$country)) stop('A country/region has inconsistent geographic-region assignments.')
  groups <- split(d, factor(d$country, levels = unique(d$country)))
  all <- do.call(rbind, lapply(groups, function(g) {
    n <- nrow(g)
    data.frame(country = g$country[1], region = g$region[1], n = n,
      mean = mean(g$abundance), sd = if (n > 1) stats::sd(g$abundance) else NA_real_,
      sem = if (n > 1) stats::sd(g$abundance) / sqrt(n) else NA_real_)
  }))
  all <- all[order(-all$mean, all$country), , drop = FALSE]
  rownames(all) <- NULL
  list(all = all, display = all[all$n >= min_n, , drop = FALSE])
}

ed1_run <- function(data, output_dir, min_n = 10) {
  ed_require('ggplot2')
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  result <- ed1_summary(data, min_n)
  ed_write_table(result$all, file.path(output_dir, 'ED1_all_country_summaries.tsv'))
  ed_write_table(result$display, file.path(output_dir, 'ED1_display_summaries.tsv'))
  d <- result$display
  if (!nrow(d)) stop('No countries/regions meet the ED1 display threshold.')
  d$label <- paste0(d$country, ' (', d$n, ')')
  d$label <- factor(d$label, levels = d$label)
  palette <- c('This study' = '#ED6A78', Asia = '#EACA52', Europe = '#4F7FB0',
    Africa = '#A6BE78', Oceania = '#9988B9', 'South America' = '#D79569', 'North America' = '#69C6C3')
  extra <- setdiff(unique(d$region), names(palette))
  if (length(extra)) palette <- c(palette, setNames(grDevices::hcl.colors(length(extra), 'Dark 3'), extra))
  p <- ggplot2::ggplot(d, ggplot2::aes(label, mean, fill = region)) +
    ggplot2::geom_col(width = .78) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = mean - sem, ymax = mean + sem),
                          width = .18, linewidth = .4, na.rm = TRUE) +
    ggplot2::scale_fill_manual(values = palette, name = 'Geographic region') +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.02, .08))) +
    ggplot2::labs(x = NULL, y = 'Mean ARG abundance (copies per cell)') +
    ggplot2::theme_classic(base_size = 10) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 60, hjust = 1),
                   legend.position = 'right', panel.border = ggplot2::element_rect(fill = NA))
  width <- max(9, .36 * nrow(d) + 3)
  ggplot2::ggsave(file.path(output_dir, 'Extended_Data_Fig_1.pdf'), p, width = width, height = 5.5, limitsize = FALSE)
  ggplot2::ggsave(file.path(output_dir, 'Extended_Data_Fig_1.png'), p, width = width, height = 5.5, dpi = 300, limitsize = FALSE)
  invisible(list(summary = result, plot = p))
}

ed2_settings <- function() c("Hospital", "WWTP", "Community", "Wet market")

ed2_prepare <- function(cohort, families) {
  a <- cohort$abundance
  m <- cohort$metadata
  annotation <- cohort$families
  required <- c("sample_id", "setting", "city", "month", "physical_site")
  if (!is.matrix(a) || !is.numeric(a) || any(!is.finite(a)) || any(a < 0))
    stop("ED2 requires a finite nonnegative numeric abundance matrix")
  if (!all(required %in% names(m)) || nrow(m) != nrow(a) || anyNA(m[required]))
    stop("ED2 metadata are missing or not aligned with abundance")
  if (anyDuplicated(m$sample_id) || (!is.null(rownames(a)) &&
      !identical(rownames(a), as.character(m$sample_id))))
    stop("ED2 sample identifiers must be unique and aligned")
  if (!setequal(as.character(m$setting), ed2_settings()))
    stop("ED2 requires all four setting names")
  if (!all(c("subtype", "family") %in% names(annotation)) ||
      anyDuplicated(annotation$subtype) || anyDuplicated(colnames(a)))
    stop("ED2 subtype annotations must be unique")
  if (is.null(colnames(a)) || !all(annotation$subtype %in% colnames(a)))
    stop("ED2 requires every mapped subtype to occur in the abundance matrix")
  if (!length(families) || anyNA(families) || anyDuplicated(families))
    stop("ED2 families must be a nonempty unique vector")
  if (!all(families %in% annotation$family))
    stop("ED2 requested families must each have mapped subtypes")
  result <- lapply(families, function(family) {
    select <- match(annotation$subtype[annotation$family == family], colnames(a))
    total <- rowSums(a[, select, drop = FALSE])
    richness <- rowSums(a[, select, drop = FALSE] > 0)
    positive <- total[total > 0]
    pseudocount <- if (length(positive)) min(positive) / 2 else NA_real_
    data.frame(m[required], family = family, total = total, richness = richness,
               pseudocount = pseudocount, log_abundance = log10(total + pseudocount),
               stringsAsFactors = FALSE, row.names = NULL)
  })
  do.call(rbind, result)
}

ed2_design <- function(metadata) {
  d <- metadata
  d$setting <- factor(d$setting, levels = ed2_settings())
  for (v in c("city", "month", "physical_site")) d[[v]] <- factor(d[[v]])
  fail <- function(status, message) list(ok = FALSE, status = status, message = message, data = d)
  if (anyNA(d[c("setting", "city", "month", "physical_site")]) ||
      any(table(d$setting) == 0)) return(fail("invalid_design", "All four settings are required"))
  site_n <- table(d$physical_site)
  if (length(site_n) < 2 || !any(site_n > 1))
    return(fail("invalid_design", "At least two physical sites and repeated site observations are required"))
  x <- tryCatch(stats::model.matrix(~ setting + city + month, d), error = identity)
  if (inherits(x, "error")) return(fail("invalid_design", conditionMessage(x)))
  rank <- qr(x)$rank
  if (rank < ncol(x)) return(fail("rank_deficient", "The full fixed-effect design is rank deficient"))
  if (nrow(x) <= ncol(x)) return(fail("invalid_design", "The full design has no residual degrees of freedom"))
  list(ok = TRUE, status = "ok", message = "", data = d, rank = rank, columns = ncol(x))
}

ed2_empty_pairs <- function() {
  data.frame(family = character(), outcome = character(), group1 = character(), group2 = character(),
             estimate = double(), std_error = double(), z = double(), p_value = double(),
             p_adjusted = double(), ratio = double(), model_scale = character(),
             model_status = character(), usable_for_inference = logical())
}

ed2_letters <- function(pairwise, settings = ed2_settings(), alpha = 0.05) {
  if (nrow(pairwise) != choose(length(settings), 2) || any(!is.finite(pairwise$p_adjusted)))
    stop("ED2 letters require all finite pairwise adjusted P values")
  p <- matrix(1, length(settings), length(settings), dimnames = list(settings, settings))
  for (i in seq_len(nrow(pairwise))) {
    p[pairwise$group1[i], pairwise$group2[i]] <- pairwise$p_adjusted[i]
    p[pairwise$group2[i], pairwise$group1[i]] <- pairwise$p_adjusted[i]
  }
  letters <- multcompView::multcompLetters(p, threshold = alpha)$Letters[settings]
  for (i in seq_len(nrow(pairwise))) {
    shared <- length(intersect(strsplit(letters[pairwise$group1[i]], "")[[1]],
                               strsplit(letters[pairwise$group2[i]], "")[[1]])) > 0
    if (shared != (pairwise$p_adjusted[i] >= alpha))
      stop("ED2 compact-letter display failed its pairwise consistency check")
  }
  letters
}

ed2_contrasts <- function(fit, d, family, outcome, status) {
  settings <- ed2_settings()
  reference <- d[rep(1, length(settings)), , drop = FALSE]
  reference$setting <- factor(settings, levels = settings)
  fitted_x <- lme4::getME(fit, "X")
  x <- stats::model.matrix(~ setting + city + month, reference,
                          contrasts.arg = attr(fitted_x, "contrasts"))
  beta <- lme4::fixef(fit)
  x <- x[, names(beta), drop = FALSE]
  covariance <- as.matrix(stats::vcov(fit))[names(beta), names(beta), drop = FALSE]
  indices <- utils::combn(seq_along(settings), 2)
  contrast <- t(vapply(seq_len(ncol(indices)), function(i)
    x[indices[1, i], ] - x[indices[2, i], ], numeric(length(beta))))
  estimate <- as.vector(contrast %*% beta)
  se <- sqrt(rowSums((contrast %*% covariance) * contrast))
  if (any(!is.finite(se)) || any(se <= 0)) stop("Invalid Wald contrast standard errors")
  z <- estimate / se
  p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
  data.frame(family = family, outcome = outcome, group1 = settings[indices[1, ]],
             group2 = settings[indices[2, ]], estimate = estimate, std_error = se, z = z,
             p_value = p, p_adjusted = stats::p.adjust(p, method = "BH"),
             ratio = if (outcome == "richness") exp(estimate) else NA_real_,
             model_scale = if (outcome == "richness") "log expected richness" else "log10 abundance + c",
             model_status = status, usable_for_inference = status %in% c("ok", "singular"))
}

ed2_fit <- function(d, outcome) {
  family <- unique(d$family)
  design <- ed2_design(d)
  diagnostic <- data.frame(family = family, outcome = outcome, status = "not_fitted",
    n_samples = nrow(d), n_sites = length(unique(d$physical_site)),
    fixed_rank = if (design$ok) design$rank else NA_integer_,
    fixed_columns = if (design$ok) design$columns else NA_integer_,
    converged = FALSE, singular = NA, REML = outcome == "abundance", theta = NA_real_,
    sigma = NA_real_, log_likelihood = NA_real_, AIC = NA_real_,
    formula = "value ~ setting + city + month + (1 | physical_site)",
    message = "", warnings = "", letters_shown = FALSE, stringsAsFactors = FALSE)
  result <- list(model = NULL, diagnostic = diagnostic, pairwise = ed2_empty_pairs(),
                 letters = stats::setNames(rep("", 4), ed2_settings()))
  finish <- function(status, message = "") {
    result$diagnostic$status <- status
    result$diagnostic$message <- message
    result
  }
  if (all(d$total == 0)) return(finish("all_zero", "No positive family totals; no model or letters estimated"))
  if (!design$ok) return(finish(design$status, design$message))
  d <- design$data
  d$value <- if (outcome == "abundance") d$log_abundance else d$richness
  if (length(unique(d$value)) < 2) return(finish("constant_response", "The response is constant"))
  warnings <- character()
  fit <- tryCatch(withCallingHandlers({
    formula <- value ~ setting + city + month + (1 | physical_site)
    if (outcome == "abundance") {
      lme4::lmer(formula, data = d, REML = TRUE,
                 control = lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
    } else {
      lme4::glmer.nb(formula, data = d,
                    control = lme4::glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 200000)))
    }
  }, warning = function(w) {
    warnings <<- c(warnings, conditionMessage(w)); invokeRestart("muffleWarning")
  }, message = function(m) {
    warnings <<- c(warnings, conditionMessage(m)); invokeRestart("muffleMessage")
  }), error = identity)
  result$diagnostic$warnings <- paste(unique(warnings), collapse = " | ")
  if (inherits(fit, "error")) return(finish("fit_failed", conditionMessage(fit)))
  result$model <- fit
  singular <- lme4::isSingular(fit, tol = 1e-4)
  convergence <- unlist(fit@optinfo$conv$lme4$messages)
  convergence <- convergence[!grepl("boundary.*singular", convergence, ignore.case = TRUE)]
  optimizer <- unlist(fit@optinfo$conv$opt)
  bad_warning <- any(grepl("failed to converge|iteration limit|nearly unidentifiable|not positive definite",
                          warnings, ignore.case = TRUE))
  converged <- !length(convergence) && !any(optimizer != 0) && !bad_warning
  result$diagnostic$converged <- converged
  result$diagnostic$singular <- singular
  result$diagnostic$theta <- if (outcome == "richness") lme4::getME(fit, "glmer.nb.theta") else NA_real_
  result$diagnostic$sigma <- if (outcome == "abundance") stats::sigma(fit) else NA_real_
  result$diagnostic$log_likelihood <- as.numeric(stats::logLik(fit))
  result$diagnostic$AIC <- stats::AIC(fit)
  if (!converged) return(finish("nonconverged", paste(convergence, collapse = " | ")))
  status <- if (singular) "singular" else "ok"
  pairs <- tryCatch(ed2_contrasts(fit, d, family, outcome, status), error = identity)
  if (inherits(pairs, "error")) return(finish("contrast_failed", conditionMessage(pairs)))
  result$pairwise <- pairs
  letters <- tryCatch(ed2_letters(pairs), error = identity)
  if (inherits(letters, "error")) return(finish("letters_failed", conditionMessage(letters)))
  result$letters <- letters
  result$diagnostic$letters_shown <- TRUE
  finish(status, if (singular) "Singular random-effect fit; converged fixed-effect contrasts retained" else "")
}

ed2_panel <- function(d, outcome, letters, status = "ok") {
  settings <- ed2_settings()
  colors <- c(Hospital = "#ED6A78", WWTP = "#4F7FB0", Community = "#69C6C3", "Wet market" = "#EACA52")
  d$x <- match(d$setting, settings)
  d$value <- if (outcome == "abundance") d$log_abundance else d$richness
  y_label <- if (outcome == "abundance") "log10(family abundance + c)" else "Detected subtypes"
  base <- ggplot2::ggplot() + ggplot2::theme_classic(base_size = 10) +
    ggplot2::scale_x_continuous(breaks = 1:4, labels = settings, limits = c(.5, 4.6)) +
    ggplot2::labs(title = unique(d$family), x = NULL, y = y_label) +
    ggplot2::theme(legend.position = "none", axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                   plot.title = ggplot2::element_text(face = "bold", hjust = .5),
                   plot.margin = ggplot2::margin(8, 9, 6, 8))
  if (!any(is.finite(d$value))) return(base +
    ggplot2::annotate("text", x = 2.5, y = 0, label = "No positive totals\nTransformation undefined", size = 3) +
    ggplot2::scale_y_continuous(breaks = NULL))
  spread <- diff(range(d$value))
  if (spread == 0) spread <- max(1, abs(d$value[1]) * .2)
  ymax <- max(d$value)
  polygons <- lapply(seq_along(settings), function(i) {
    values <- d$value[d$x == i]
    if (length(unique(values)) < 2) return(NULL)
    density <- stats::density(values, from = min(values), to = max(values), n = 128)
    data.frame(x = c(i + .06, i + .06 + .35 * density$y / max(density$y), i + .06),
               y = c(density$x[1], density$x, tail(density$x, 1)), setting = settings[i])
  })
  polygon <- do.call(rbind, polygons)
  if (is.null(polygon)) polygon <- data.frame(x = numeric(), y = numeric(), setting = character())
  means <- data.frame(x = 1:4, setting = settings,
    label = paste0("Mean\n", formatC(vapply(settings, function(s) mean(d$value[d$setting == s]), numeric(1)),
                                     digits = 2, format = "f")), y = ymax + .12 * spread)
  letter_data <- data.frame(x = 1:4, y = ymax + .24 * spread, label = unname(letters[settings]))
  base + ggplot2::scale_fill_manual(values = colors, drop = FALSE) +
    ggplot2::geom_polygon(data = polygon, ggplot2::aes(x = x, y = y, group = setting, fill = setting),
                           alpha = .65, linewidth = .2, color = "grey40") +
    ggplot2::geom_point(data = d, ggplot2::aes(x = x - .2, y = value, fill = setting),
                        position = ggplot2::position_jitter(width = .11, height = 0, seed = 31),
                        shape = 21, stroke = .15, size = 1.1, alpha = .5) +
    ggplot2::geom_boxplot(data = d, ggplot2::aes(x = x, y = value, group = setting),
                          width = .10, outlier.shape = NA, fill = "white", linewidth = .4) +
    ggplot2::geom_text(data = means, ggplot2::aes(x = x, y = y, label = label), size = 2.5) +
    ggplot2::geom_text(data = letter_data, ggplot2::aes(x = x, y = y, label = label), size = 3.5) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(.04, .07))) +
    ggplot2::labs(subtitle = if (status == "ok") NULL else paste("Model:", status)) +
    ggplot2::coord_cartesian(clip = "off")
}

ed2_run <- function(cohort, output_dir, families = c("OXA", "KPC", "VIM", "NDM", "IMP")) {
  packages <- c("lme4", "ggplot2", "multcompView", "patchwork")
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("ED2 requires packages: ", paste(missing, collapse = ", "))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  totals <- ed2_prepare(cohort, families)
  fits <- list()
  diagnostics <- pairwise <- panels <- list()
  for (outcome in c("abundance", "richness")) for (family in families) {
    d <- totals[totals$family == family, , drop = FALSE]
    result <- ed2_fit(d, outcome)
    key <- paste(family, outcome, sep = "_")
    fits[key] <- list(result$model)
    diagnostics[[key]] <- result$diagnostic
    pairwise[[key]] <- result$pairwise
    panels[[key]] <- ed2_panel(d, outcome, result$letters, result$diagnostic$status)
  }
  diagnostics <- do.call(rbind, diagnostics)
  pairwise <- do.call(rbind, pairwise)
  rownames(diagnostics) <- rownames(pairwise) <- NULL
  pseudocounts <- do.call(rbind, lapply(families, function(family) {
    d <- totals[totals$family == family, ]
    data.frame(family = family, n_subtypes = sum(cohort$families$family == family &
      cohort$families$subtype %in% colnames(cohort$abundance)),
      n_samples = nrow(d), n_positive_samples = sum(d$total > 0),
      minimum_positive_family_total = if (any(d$total > 0)) min(d$total[d$total > 0]) else NA_real_,
      pseudocount = d$pseudocount[1])
  }))
  utils::write.csv(totals, file.path(output_dir, "ED2_family_totals.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(pseudocounts, file.path(output_dir, "ED2_pseudocounts.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(pairwise, file.path(output_dir, "ED2_pairwise.csv"), row.names = FALSE, na = "NA")
  utils::write.csv(diagnostics, file.path(output_dir, "ED2_model_diagnostics.csv"), row.names = FALSE, na = "NA")
  saveRDS(fits, file.path(output_dir, "ED2_model_fits.rds"))
  figure <- patchwork::wrap_plots(panels, ncol = length(families), nrow = 2, byrow = TRUE)
  ggplot2::ggsave(file.path(output_dir, "ED2.pdf"), figure, width = 3.2 * length(families), height = 8.5,
                  device = "pdf", limitsize = FALSE)
  ggplot2::ggsave(file.path(output_dir, "ED2.png"), figure, width = 3.2 * length(families), height = 8.5,
                  dpi = 300, bg = "white", limitsize = FALSE)
  invisible(list(totals = totals, pseudocounts = pseudocounts, diagnostics = diagnostics,
                  pairwise = pairwise, fits = fits, figure = figure))
}

ed3_setting_colours <- c(Hospital = "#ED6A78", WWTP = "#4F7FB0",
                        Community = "#69C6C3", `Wet market` = "#EACA52")

ed3_r2_components <- function(fixed_variance, site_variance, residual_variance) {
  components <- c(fixed_variance, site_variance, residual_variance)
  if (any(!is.finite(components)) || any(components < 0) || sum(components) <= 0) return(NA_real_)
  fixed_variance / sum(components)
}

ed3_model_c <- function(abundance, family_subtypes) {
  values <- abundance[, family_subtypes, drop = FALSE]
  positive <- values[is.finite(values) & values > 0]
  if (!length(positive)) return(NA_real_)
  min(positive) / 2
}

ed3_site_data <- function(cohort, selection, family) {
  selected <- selection[selection$family == family, , drop = FALSE]
  selected <- selected[order(selected$plot_order), , drop = FALSE]
  if (!nrow(selected)) stop("No selected subtypes for family ", family)
  meta <- cohort$metadata
  site_ids <- unique(as.character(meta$physical_site))
  groups <- lapply(site_ids, function(id) which(as.character(meta$physical_site) == id))
  consistent <- vapply(groups, function(i) all(vapply(c("setting", "city"),
    function(field) length(unique(as.character(meta[[field]][i]))) == 1L, logical(1))), logical(1))
  if (!all(consistent)) stop("Each physical site must have one setting and one city")
  means <- do.call(rbind, lapply(groups, function(i)
    colMeans(cohort$abundance[i, as.character(selected$subtype), drop = FALSE])))
  dimnames(means) <- list(site_ids, as.character(selected$subtype))
  site_meta <- meta[vapply(groups, `[`, integer(1), 1L), c("physical_site", "setting", "city"), drop = FALSE]
  site_meta$physical_site <- as.character(site_meta$physical_site)
  rownames(site_meta) <- NULL
  positive <- means[is.finite(means) & means > 0]
  plot_c <- if (length(positive)) min(positive) / 2 else NA_real_
  transformed <- if (is.finite(plot_c)) log10(means + plot_c) else means * NA_real_
  z <- transformed
  for (j in seq_len(ncol(transformed))) {
    profile <- transformed[, j]
    deviation <- stats::sd(profile)
    if (all(is.finite(profile)) && (length(profile) == 1L || deviation == 0)) {
      z[, j] <- 0
    } else if (is.finite(deviation) && deviation > 0) {
      z[, j] <- (profile - mean(profile)) / deviation
    }
  }
  list(means = means, log_abundance = transformed, z = z, plot_c = plot_c,
       metadata = site_meta, selection = selected)
}

ed3_branch_colour <- function(settings, palette = ed3_setting_colours,
                              min_descendants = 4L, min_purity = 0.70) {
  frequencies <- table(as.character(settings))
  if (!length(frequencies) || length(settings) < min_descendants) return("#B6B6B6")
  dominant <- names(frequencies)[which.max(frequencies)]
  if (max(frequencies) / length(settings) >= min_purity && dominant %in% names(palette)) palette[[dominant]] else "#B6B6B6"
}

ed3_cluster <- function(z, settings, palette = ed3_setting_colours) {
  if (!is.matrix(z) || nrow(z) < 1L || any(!is.finite(z))) stop("Clustering requires finite site profiles")
  if (length(settings) != nrow(z)) stop("Settings must align to site profiles")
  site_ids <- rownames(z)
  if (is.null(site_ids)) site_ids <- as.character(seq_len(nrow(z)))
  empty <- data.frame(x = numeric(), y = numeric(), xend = numeric(), yend = numeric(),
                      colour = character(), descendants = integer(), purity = numeric(),
                      dominant_setting = character(), parent = integer(), child = integer())
  if (nrow(z) == 1L) return(list(hclust = NULL, order = site_ids, segments = empty))
  tree <- stats::hclust(stats::dist(z, method = "euclidean"), method = "ward.D2")
  positions <- match(seq_len(nrow(z)), tree$order)
  nodes <- vector("list", nrow(tree$merge))
  segments <- list()
  child_node <- function(index) {
    if (index < 0L) list(x = positions[-index], y = 0, leaves = -index) else nodes[[index]]
  }
  for (i in seq_len(nrow(tree$merge))) {
    left <- child_node(tree$merge[i, 1L])
    right <- child_node(tree$merge[i, 2L])
    node <- list(x = mean(c(left$x, right$x)), y = tree$height[i], leaves = c(left$leaves, right$leaves))
    nodes[[i]] <- node
    for (k in seq_len(2L)) {
      child <- if (k == 1L) left else right
      descendant_settings <- as.character(settings[child$leaves])
      frequencies <- table(descendant_settings)
      annotation <- data.frame(colour = ed3_branch_colour(descendant_settings, palette),
        descendants = length(child$leaves), purity = max(frequencies) / length(child$leaves),
        dominant_setting = names(frequencies)[which.max(frequencies)], parent = i, child = tree$merge[i, k])
      vertical <- data.frame(x = child$x, y = child$y, xend = child$x, yend = node$y)
      horizontal <- data.frame(x = child$x, y = node$y, xend = node$x, yend = node$y)
      segments[[length(segments) + 1L]] <- cbind(rbind(vertical, horizontal), annotation[c(1L, 1L), ])
    }
  }
  list(hclust = tree, order = site_ids[tree$order], segments = do.call(rbind, segments))
}

ed3_marginal_r2 <- function(model) {
  variances <- as.data.frame(lme4::VarCorr(model))
  site_variance <- sum(variances$vcov[variances$grp == "physical_site"])
  residual_variance <- stats::sigma(model)^2
  fixed_variance <- stats::var(stats::predict(model, re.form = NA))
  ed3_r2_components(fixed_variance, site_variance, residual_variance)
}

ed3_fit_subtype <- function(abundance, metadata, pseudocount) {
  empty_stats <- data.frame(n_samples = 0L, n_sites = 0L, model_c = pseudocount,
    r2_full = NA_real_, r2_without_setting = NA_real_, r2_without_city = NA_real_,
    setting_r2 = NA_real_, city_r2 = NA_real_, setting_lrt = NA_real_, city_lrt = NA_real_,
    setting_df = NA_real_, city_df = NA_real_, setting_p = NA_real_, city_p = NA_real_, status = "not_fitted")
  names_model <- c("full", "without_setting", "without_city")
  fits <- stats::setNames(vector("list", length(names_model)), names_model)
  status <- data.frame(model = names_model, status = "not_fitted", usable = FALSE,
                       singular = NA, message = "", stringsAsFactors = FALSE)
  fail_all <- function(reason, message) {
    status$status <- reason
    status$message <- message
    empty_stats$status <- reason
    list(stats = empty_stats, status = status, models = fits, rows = integer())
  }
  if (!is.finite(pseudocount) || pseudocount <= 0) return(fail_all("invalid_pseudocount", "Family has no positive sample abundance"))
  data <- metadata[, c("sample_id", "setting", "city", "month", "physical_site"), drop = FALSE]
  data$response <- log10(abundance + pseudocount)
  included <- which(stats::complete.cases(data) & is.finite(data$response))
  data <- data[included, , drop = FALSE]
  for (column in c("setting", "city", "month", "physical_site")) data[[column]] <- droplevels(factor(data[[column]]))
  empty_stats$n_samples <- nrow(data)
  empty_stats$n_sites <- nlevels(data$physical_site)
  if (nrow(data) < 2L || stats::var(data$response) == 0) return(fail_all("constant_response", "Insufficient response variation"))
  formulas <- list(full = response ~ setting + city + month + (1 | physical_site),
                   without_setting = response ~ city + month + (1 | physical_site),
                   without_city = response ~ setting + month + (1 | physical_site))
  fixed_formulas <- list(full = response ~ setting + city + month,
                         without_setting = response ~ city + month,
                         without_city = response ~ setting + month)
  for (name in names_model) {
    index <- match(name, status$model)
    messages <- character()
    design <- tryCatch(stats::model.matrix(fixed_formulas[[name]], data), error = identity)
    if (inherits(design, "error")) {
      status$status[index] <- "invalid_design"
      status$message[index] <- conditionMessage(design)
      next
    }
    if (qr(design)$rank < ncol(design)) {
      status$status[index] <- "rank_deficient"
      status$message[index] <- "Fixed-effect design is not full rank"
      next
    }
    fit <- tryCatch(withCallingHandlers(
      lme4::lmer(formulas[[name]], data = data, REML = FALSE,
                 control = lme4::lmerControl(optimizer = "bobyqa", check.rankX = "stop.deficient",
                   optCtrl = list(maxfun = 200000))),
      warning = function(w) { messages <<- c(messages, conditionMessage(w)); invokeRestart("muffleWarning") },
      message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") }), error = identity)
    if (inherits(fit, "error")) {
      status$status[index] <- "fit_error"
      status$message[index] <- conditionMessage(fit)
      next
    }
    fits[[name]] <- fit
    convergence <- fit@optinfo$conv
    convergence_messages <- unlist(convergence$lme4$messages, use.names = FALSE)
    convergence_messages <- convergence_messages[!grepl("boundary.*singular", convergence_messages, ignore.case = TRUE)]
    failed <- (length(convergence$opt) && any(convergence$opt != 0)) || length(convergence_messages) > 0L
    status$singular[index] <- lme4::isSingular(fit)
    status$usable[index] <- !failed
    status$status[index] <- if (failed) "convergence_failure" else if (status$singular[index]) "ok_singular" else "ok"
    status$message[index] <- paste(unique(c(messages, convergence_messages)), collapse = " | ")
  }
  usable <- stats::setNames(status$usable, status$model)
  for (name in names_model) if (usable[[name]]) empty_stats[[paste0("r2_", name)]] <- ed3_marginal_r2(fits[[name]])
  for (effect in c("setting", "city")) {
    reduced <- paste0("without_", effect)
    if (!usable[["full"]] || !usable[[reduced]]) next
    ll_full <- stats::logLik(fits$full)
    ll_reduced <- stats::logLik(fits[[reduced]])
    degrees <- attr(ll_full, "df") - attr(ll_reduced, "df")
    statistic <- 2 * (as.numeric(ll_full) - as.numeric(ll_reduced))
    if (degrees <= 0 || !is.finite(statistic) || statistic < -1e-6) next
    empty_stats[[paste0(effect, "_lrt")]] <- max(statistic, 0)
    empty_stats[[paste0(effect, "_df")]] <- degrees
    empty_stats[[paste0(effect, "_p")]] <- stats::pchisq(max(statistic, 0), df = degrees, lower.tail = FALSE)
    empty_stats[[paste0(effect, "_r2")]] <- max(empty_stats$r2_full - empty_stats[[paste0("r2_", reduced)]], 0)
  }
  empty_stats$status <- if (all(usable)) "ok" else paste(unique(status$status[!status$usable]), collapse = ";")
  list(stats = empty_stats, status = status, models = fits, rows = included)
}

ed3_family_plot <- function(family, site_data, clustering, statistics,
                            palette = ed3_setting_colours) {
  selected <- site_data$selection
  n_sites <- nrow(site_data$means)
  n_subtypes <- nrow(selected)
  index <- match(clustering$order, rownames(site_data$means))
  site_meta <- site_data$metadata[index, , drop = FALSE]
  label <- if ("display_label" %in% names(selected)) as.character(selected$display_label) else as.character(selected$subtype)
  label[is.na(label) | !nzchar(label)] <- as.character(selected$subtype)[is.na(label) | !nzchar(label)]
  long <- expand.grid(x = seq_len(n_sites), subtype_index = seq_len(n_subtypes))
  long$y <- n_subtypes + 1L - long$subtype_index
  long$abundance <- as.vector(site_data$means[index, , drop = FALSE])
  long$z <- as.vector(site_data$z[index, , drop = FALSE])
  long$colour_z <- pmax(pmin(long$z, 3), -3)
  x_scale <- function() ggplot2::scale_x_continuous(limits = c(0.5, n_sites + 0.5), expand = c(0, 0))
  y_scale <- function(labels = TRUE) ggplot2::scale_y_continuous(
    breaks = seq_len(n_subtypes), labels = if (labels) rev(label) else NULL,
    limits = c(0.5, n_subtypes + 0.5), expand = c(0, 0))
  base_theme <- ggplot2::theme_minimal(base_size = 9) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(2, 4, 2, 4), legend.key.height = grid::unit(3, "mm"))
  bubble <- ggplot2::ggplot(long, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point(ggplot2::aes(size = abundance, colour = colour_z)) +
    ggplot2::scale_size_area(max_size = 4.3, name = "Site mean abundance\n(copies per cell)") +
    ggplot2::scale_colour_gradient2(low = "#376CAA", mid = "#F5F5F3", high = "#C94842",
      midpoint = 0, limits = c(-3, 3), breaks = c(-3, 0, 3), name = "Within-subtype Z") +
    x_scale() + y_scale() + ggplot2::labs(x = "Physical sites", y = NULL) + base_theme +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
  dendrogram <- ggplot2::ggplot(clustering$segments) +
    ggplot2::geom_segment(ggplot2::aes(x = x, y = y, xend = xend, yend = yend, colour = colour), linewidth = 0.35) +
    ggplot2::scale_colour_identity() + x_scale() + ggplot2::theme_void() +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 4, 0, 4))
  strip <- function(field, colours, title) {
    d <- data.frame(x = seq_len(n_sites), value = as.character(site_meta[[field]]))
    ggplot2::ggplot(d, ggplot2::aes(x = x, y = 1, fill = value)) +
      ggplot2::geom_tile(width = 1, height = 1) + ggplot2::scale_fill_manual(values = colours, name = title) +
      x_scale() + ggplot2::scale_y_continuous(breaks = 1, labels = title, expand = c(0, 0)) +
      base_theme + ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.text.y = ggplot2::element_text(size = 8))
  }
  cities <- sort(unique(as.character(site_meta$city)))
  city_colours <- stats::setNames(grDevices::hcl.colors(length(cities), palette = "Dark 3"), cities)
  setting_strip <- strip("setting", palette, "Setting")
  city_strip <- strip("city", city_colours, "City")
  statistics <- statistics[match(selected$subtype, statistics$subtype), , drop = FALSE]
  r2 <- data.frame(y = rep(rev(seq_len(n_subtypes)), each = 2),
    effect = rep(c("Setting R2", "City R2"), n_subtypes),
    r2 = as.vector(t(as.matrix(statistics[, c("setting_r2", "city_r2")]))))
  r2$effect <- factor(r2$effect, levels = c("Setting R2", "City R2"))
  r2$offset_y <- r2$y + ifelse(r2$effect == "Setting R2", 0.17, -0.17)
  bars <- ggplot2::ggplot(r2, ggplot2::aes(x = r2, y = offset_y, colour = effect)) +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = r2, yend = offset_y), linewidth = 2.2, na.rm = TRUE) +
    ggplot2::scale_colour_manual(values = c(`Setting R2` = "#D84E4B", `City R2` = "#3B73B9"), name = NULL) +
    ggplot2::scale_x_continuous(limits = c(0, max(0.05, r2$r2, na.rm = TRUE)), expand = ggplot2::expansion(mult = c(0, 0.08))) +
    y_scale(FALSE) + ggplot2::labs(x = expression(Delta * R[marginal]^2), y = NULL) + base_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 7))
  header <- ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0.5, label = family, hjust = 0, size = 5, fontface = "bold") +
    ggplot2::xlim(0, 1) + ggplot2::ylim(0, 1) + ggplot2::theme_void()
  pieces <- list(header, patchwork::plot_spacer(), dendrogram, patchwork::plot_spacer(),
                 setting_strip, patchwork::plot_spacer(), city_strip, patchwork::plot_spacer(), bubble, bars)
  patchwork::wrap_plots(pieces, ncol = 2, widths = c(8, 1.6),
    heights = c(0.4, 1.1, 0.22, 0.22, max(1.6, n_subtypes * 0.30)), guides = "collect") &
    ggplot2::theme(legend.position = "right", legend.text = ggplot2::element_text(size = 7),
                   legend.title = ggplot2::element_text(size = 8))
}

ed3_run <- function(cohort, selection, output_dir,
                    families = c("OXA", "KPC", "VIM", "NDM", "IMP")) {
  for (package in c("lme4", "ggplot2", "patchwork")) {
    if (!requireNamespace(package, quietly = TRUE)) stop("ED3 requires package: ", package)
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  write_csv <- function(data, filename) utils::write.csv(data, file.path(output_dir, filename), row.names = FALSE, na = "NA")
  write_matrix <- function(data, filename) write_csv(data.frame(physical_site = rownames(data), data, check.names = FALSE), filename)
  write_csv(selection, "ED3_display_selection.csv")
  results <- list()
  all_statistics <- list()
  all_model_status <- list()
  all_models <- list()
  family_status <- list()
  plots <- list()
  for (family in families) {
    selected <- selection[selection$family == family, , drop = FALSE]
    if (!nrow(selected)) {
      family_status[[family]] <- data.frame(family = family, status = "no_selected_subtypes", message = "")
      next
    }
    family_subtypes <- as.character(cohort$families$subtype[cohort$families$family == family])
    if (!length(family_subtypes) || !all(family_subtypes %in% colnames(cohort$abundance))) stop("Missing abundance columns for family ", family)
    model_c <- ed3_model_c(cohort$abundance, family_subtypes)
    fits <- lapply(as.character(selected$subtype), function(subtype)
      ed3_fit_subtype(cohort$abundance[, subtype], cohort$metadata, model_c))
    names(fits) <- as.character(selected$subtype)
    statistics <- do.call(rbind, lapply(names(fits), function(subtype)
      cbind(data.frame(family = family, subtype = subtype), fits[[subtype]]$stats)))
    statistics$setting_q <- stats::p.adjust(statistics$setting_p, method = "BH")
    statistics$city_q <- stats::p.adjust(statistics$city_p, method = "BH")
    model_status <- do.call(rbind, lapply(names(fits), function(subtype)
      cbind(data.frame(family = family, subtype = subtype), fits[[subtype]]$status)))
    all_statistics[[family]] <- statistics
    all_model_status[[family]] <- model_status
    all_models[[family]] <- fits
    write_csv(statistics, paste0("ED3_", family, "_statistics.csv"))
    site_data <- ed3_site_data(cohort, selection, family)
    write_matrix(site_data$means, paste0("ED3_", family, "_site_means.csv"))
    write_matrix(site_data$log_abundance, paste0("ED3_", family, "_site_log10.csv"))
    write_matrix(site_data$z, paste0("ED3_", family, "_site_z_unclipped.csv"))
    write_csv(data.frame(family = family, model_c = model_c, plot_c = site_data$plot_c,
      model_scope = "all sample abundances and all family subtypes",
      plot_scope = "site arithmetic means across displayed family subtypes"), paste0("ED3_", family, "_pseudocounts.csv"))
    results[[family]] <- list(sites = site_data, statistics = statistics, model_status = model_status, models = fits)
    if (any(!is.finite(site_data$z))) {
      family_status[[family]] <- data.frame(family = family, status = "plot_unavailable", message = "No positive site mean; log transform is undefined")
      next
    }
    clustering <- ed3_cluster(site_data$z, site_data$metadata$setting)
    ordered_meta <- site_data$metadata[match(clustering$order, site_data$metadata$physical_site), , drop = FALSE]
    write_csv(cbind(data.frame(plot_column = seq_along(clustering$order)), ordered_meta), paste0("ED3_", family, "_site_order.csv"))
    write_csv(clustering$segments, paste0("ED3_", family, "_dendrogram_edges.csv"))
    saveRDS(clustering$hclust, file.path(output_dir, paste0("ED3_", family, "_clustering.rds")))
    plot <- ed3_family_plot(family, site_data, clustering, statistics)
    width <- max(10, min(22, 6 + nrow(site_data$means) * 0.075))
    height <- max(4.5, 2.5 + nrow(selected) * 0.30)
    ggplot2::ggsave(file.path(output_dir, paste0("ED3_", family, ".pdf")), plot = plot,
      width = width, height = height, device = grDevices::pdf, useDingbats = FALSE, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(file.path(output_dir, paste0("ED3_", family, ".png")), plot = plot,
      width = width, height = height, dpi = 300, bg = "white", limitsize = FALSE)
    results[[family]]$clustering <- clustering
    results[[family]]$plot <- plot
    plots[[family]] <- plot
    family_status[[family]] <- data.frame(family = family, status = "ok", message = "")
  }
  statistics <- if (length(all_statistics)) do.call(rbind, all_statistics) else data.frame()
  model_status <- if (length(all_model_status)) do.call(rbind, all_model_status) else data.frame()
  statuses <- if (length(family_status)) do.call(rbind, family_status) else data.frame()
  write_csv(statistics, "ED3_statistics.csv")
  write_csv(model_status, "ED3_model_status.csv")
  write_csv(statuses, "ED3_family_status.csv")
  saveRDS(all_models, file.path(output_dir, "ED3_models.rds"))
  combined <- NULL
  if (length(plots)) {
    combined <- patchwork::wrap_plots(plots, ncol = 1)
    height <- sum(vapply(plots, function(x) 4.5, numeric(1))) +
      sum(vapply(results[names(plots)], function(x) max(0, nrow(x$sites$selection) * 0.30 - 2), numeric(1)))
    ggplot2::ggsave(file.path(output_dir, "ED3_combined.pdf"), combined, width = 16, height = height,
      device = grDevices::pdf, useDingbats = FALSE, bg = "white", limitsize = FALSE)
    ggplot2::ggsave(file.path(output_dir, "ED3_combined.png"), combined, width = 16, height = height,
      dpi = 200, bg = "white", limitsize = FALSE)
  }
  invisible(list(selection = selection, families = results, statistics = statistics,
                 model_status = model_status, family_status = statuses, combined_plot = combined))
}

extended_data_config <- function(input_dir = 'input', output_dir = 'extended_data_output') {
  list(
    benchmark_file = file.path(input_dir, 'global_benchmark.tsv'),
    abundance_file = file.path(input_dir, 'subtype_abundance.tsv'),
    metadata_file = file.path(input_dir, 'sample_metadata.tsv'),
    family_map_file = file.path(input_dir, 'subtype_family_map.tsv'),
    selection_file = file.path(input_dir, 'selected_subtypes.csv'),
    output_dir = output_dir,
    abundance_orientation = 'sample_rows',
    blank_policy = 'stop',
    metadata_columns = NULL,
    selection_mode = 'selected',
    selection_columns = NULL,
    families = c('OXA', 'KPC', 'VIM', 'NDM', 'IMP'),
    initial_screen_top_n = 20,
    minimum_country_samples = 10
  )
}

ed_input_templates <- function(input_dir = 'input') {
  dir.create(input_dir, recursive = TRUE, showWarnings = FALSE)
  headers <- list(
    global_benchmark.tsv = 'sample_id\tcountry\tregion\tsource\tabundance',
    subtype_abundance.tsv = 'sample_id\tsubtype_1',
    sample_metadata.tsv = 'sample_id\tsetting\tcity\tmonth\tphysical_site',
    subtype_family_map.tsv = 'subtype\tfamily',
    selected_subtypes.csv = 'family,subtype,plot_order,display_label'
  )
  paths <- file.path(input_dir, names(headers))
  if (any(file.exists(paths))) stop('Template generation will not overwrite existing input files.')
  for (i in seq_along(paths)) writeLines(headers[[i]], paths[i], useBytes = TRUE)
  invisible(paths)
}

run_extended_data <- function(config = extended_data_config(), figures = c(1, 2, 3)) {
  if (!length(figures) || anyNA(figures) || any(!figures %in% 1:3)) stop('figures must be selected from 1, 2 and 3.')
  figures <- unique(figures)
  required <- character()
  if (1 %in% figures) required <- c(required, config$benchmark_file)
  if (any(c(2, 3) %in% figures)) required <- c(required, config$abundance_file, config$metadata_file, config$family_map_file)
  if (3 %in% figures) required <- c(required, config$selection_file)
  if (any(!file.exists(required))) stop('Missing input files:\n', paste(required[!file.exists(required)], collapse = '\n'))
  if (!length(config$families) || any(!config$families %in% c('OXA', 'KPC', 'VIM', 'NDM', 'IMP')) || anyDuplicated(config$families))
    stop('Invalid family configuration.')
  cohort <- NULL
  if (any(c(2, 3) %in% figures)) {
    metadata <- ed_read_table(config$metadata_file)
    if (!is.null(config$metadata_columns)) {
      if (is.null(names(config$metadata_columns)) || any(!unname(config$metadata_columns) %in% names(metadata)))
        stop('Invalid metadata column map.')
      for (key in names(config$metadata_columns)) metadata[[key]] <- metadata[[config$metadata_columns[[key]]]]
    }
    cohort <- ed_prepare_cohort(ed_read_table(config$abundance_file), metadata,
      ed_read_table(config$family_map_file), blank_policy = config$blank_policy,
      orientation = config$abundance_orientation)
    absent <- setdiff(config$families, cohort$families$family)
    if (length(absent)) stop('Missing family assignments: ', paste(absent, collapse = ', '))
  }
  selection <- NULL
  if (3 %in% figures) {
    selection <- ed_selection(ed_read_table(config$selection_file), mode = config$selection_mode,
      families = config$families, top_n = config$initial_screen_top_n, columns = config$selection_columns)
    absent <- setdiff(config$families, selection$family)
    if (length(absent)) stop('No selected subtypes for: ', paste(absent, collapse = ', '))
    j <- match(selection$subtype, cohort$families$subtype)
    if (anyNA(j) || any(cohort$families$family[j] != selection$family)) stop('Selection and family map disagree.')
  }
  dir.create(config$output_dir, recursive = TRUE, showWarnings = FALSE)
  result <- list()
  if (1 %in% figures) {
    message('Extended Data Fig. 1')
    result$ED1 <- ed1_run(ed_read_table(config$benchmark_file), file.path(config$output_dir, 'ED1'), config$minimum_country_samples)
  }
  if (2 %in% figures) {
    message('Extended Data Fig. 2')
    result$ED2 <- ed2_run(cohort, file.path(config$output_dir, 'ED2'), families = config$families)
  }
  if (3 %in% figures) {
    message('Extended Data Fig. 3')
    result$ED3 <- ed3_run(cohort, selection, file.path(config$output_dir, 'ED3'), families = config$families)
  }
  writeLines(capture.output(utils::sessionInfo()), file.path(config$output_dir, 'sessionInfo.txt'))
  invisible(result)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 1L) {
    message('Usage: Rscript Extended_Data_Figures_1_to_3.R INPUT_DIR [OUTPUT_DIR] [1,2,3]')
  } else {
    cfg <- extended_data_config(args[1], if (length(args) >= 2L) args[2] else 'extended_data_output')
    requested <- if (length(args) >= 3L) as.integer(strsplit(args[3], ',', fixed = TRUE)[[1]]) else 1:3
    run_extended_data(cfg, requested)
  }
}
