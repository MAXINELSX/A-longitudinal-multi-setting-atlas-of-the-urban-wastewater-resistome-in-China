#!/usr/bin/env Rscript
# Complete self-contained Figure 4c-e plotting script.
# Source this file and call run_fig4_cde(DATA_DIRECTORY, OUTPUT_DIRECTORY).
# All figure data are read from external TSVs; no genomic analyses are rerun.
# Dependencies: ggplot2, patchwork, svglite.
# Individual panel scripts and input documentation accompany this file.

.fig4c <- local({
# Fig. 4c: cluster composition and wastewater representation.
# Plot styling recovered from plot_final17_mixed_clusters_25MAG_onefile.R.
# Reads the final plotting TSV directly without count adjustments.
# The separately supplied 233cluster.R is preserved in original_reference.
# Required package: ggplot2. Install dependencies separately.
#
# source("Fig4c.R")
# run_fig4c("local_inputs/Fig4c_cluster_data.tsv", "figures")
# Rscript Fig4c.R local_inputs/Fig4c_cluster_data.tsv figures
#
# Required TSV columns: Cluster, Clinical, Wastewater.
# Input row order and original audit columns are retained.

prepare_fig4c_data <- function(input_file) {
  if (missing(input_file) || length(input_file) != 1L ||
      is.na(input_file) || !file.exists(input_file)) {
    stop("Supply an existing Fig4c input TSV.", call. = FALSE)
  }
  df <- utils::read.delim(input_file, sep = "\t", header = TRUE,
                          stringsAsFactors = FALSE, check.names = FALSE)
  absent <- setdiff(c("Cluster", "Clinical", "Wastewater"), names(df))
  if (length(absent)) {
    stop("Fig4c input is missing: ", paste(absent, collapse = ", "), call. = FALSE)
  }
  df$Cluster <- as.character(df$Cluster)
  if (anyNA(df$Cluster) || any(!nzchar(df$Cluster)) || anyDuplicated(df$Cluster)) {
    stop("Cluster labels must be non-empty and unique.", call. = FALSE)
  }
  for (column in c("Clinical", "Wastewater")) {
    value <- suppressWarnings(as.numeric(df[[column]]))
    if (any(!is.finite(value)) || any(value < 1) || any(value != floor(value))) {
      stop(column, " must contain positive integer counts.", call. = FALSE)
    }
    df[[column]] <- value
  }
  totals <- df$Clinical + df$Wastewater
  percentages <- 100 * df$Wastewater / totals
  if ("Total" %in% names(df) && !isTRUE(all.equal(as.numeric(df$Total), totals))) {
    stop("Supplied Total does not equal Clinical + Wastewater.", call. = FALSE)
  }
  if ("WW_percent" %in% names(df) &&
      !isTRUE(all.equal(as.numeric(df$WW_percent), percentages))) {
    stop("Supplied WW_percent does not match the source counts.", call. = FALSE)
  }
  df$Total <- totals
  df$WW_percent <- percentages
  df$y <- rev(seq_len(nrow(df)))
  stopifnot(nrow(df) == 17L, sum(df$Clinical) == 221,
            sum(df$Wastewater) == 25, sum(df$Total) == 246)
  df
}

plot_fig4c <- function(df) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Install the required R package ggplot2 before running Fig4c.", call. = FALSE)
  }
  FONT <- "Arial"
  COL_CLINICAL <- "#B78387"
  COL_WASTEWATER <- "#6D82AC"
  if (.Platform$OS.type == "windows") {
    grDevices::windowsFonts(Arial = grDevices::windowsFont("Arial"))
  }
  n_clusters <- nrow(df)
  max_total <- max(df$Total)
  max_percent <- max(df$WW_percent)
  
  # Keep the left and right visual scales identical to the reference figure:
  # left values are genome/MAG counts; right values are wastewater percentages.
  left_limit <- -ceiling((max_total + 4) / 5) * 5
  right_limit <- max(60, ceiling((max_percent + 12) / 5) * 5)
  
  negative_breaks <- pretty(c(left_limit, 0), n = 4)
  negative_breaks <- negative_breaks[
    negative_breaks < 0 & negative_breaks >= left_limit
  ]
  positive_breaks <- c(25, 50, 75, 100)
  positive_breaks <- positive_breaks[positive_breaks <= right_limit]
  
  x_breaks <- c(negative_breaks, 0, positive_breaks)
  x_labels <- c(
    as.character(abs(negative_breaks)),
    "0",
    paste0(positive_breaks, "%")
  )
  
  bar_half_height <- 0.32
  
  # Position total-count labels inside larger bars and just outside very small bars.
  df$total_label_x <- ifelse(
    df$Total >= 6,
    -df$Total + 1.2,
    -df$Total - 1.3
  )
  df$total_label_hjust <- ifelse(df$Total >= 6, 0, 1)
  
  # Add extra room to the first low-percentage labels so they do not touch the dot.
  df$percent_label_x <- df$WW_percent + ifelse(df$WW_percent < 5, 4.5, 5.2)
  
  figure_cluster <- ggplot2::ggplot(df) +
    # Wastewater is the far-left blue portion.
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = -Total,
        xmax = -Clinical,
        ymin = y - bar_half_height,
        ymax = y + bar_half_height,
        fill = "Wastewater"
      ),
      linewidth = 0
    ) +
    # Clinical is the rose portion adjacent to zero.
    ggplot2::geom_rect(
      ggplot2::aes(
        xmin = -Clinical,
        xmax = 0,
        ymin = y - bar_half_height,
        ymax = y + bar_half_height,
        fill = "Clinical"
      ),
      linewidth = 0
    ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.85,
      color = "black"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = WW_percent,
        y = y,
        yend = y
      ),
      color = COL_WASTEWATER,
      linewidth = 1.0,
      lineend = "round"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(x = WW_percent, y = y),
      color = COL_WASTEWATER,
      size = 4.0
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = total_label_x,
        y = y,
        label = Total,
        hjust = total_label_hjust
      ),
      family = FONT,
      fontface = "bold",
      size = 4.7,
      color = "black"
    ) +
    ggplot2::geom_text(
      ggplot2::aes(
        x = percent_label_x,
        y = y,
        label = sprintf("%.1f%%", WW_percent)
      ),
      hjust = 0,
      family = FONT,
      fontface = "bold",
      size = 4.7,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = left_limit - 17,
      y = n_clusters + 1.38,
      label = "C",
      hjust = 0,
      vjust = 1,
      family = FONT,
      fontface = "bold",
      size = 19,
      color = "black"
    ) +
    ggplot2::annotate(
      "text",
      x = mean(c(left_limit, 0)),
      y = -0.42,
      label = "Number of genomes or MAGs",
      family = FONT,
      fontface = "bold",
      size = 4.8,
      vjust = 1
    ) +
    ggplot2::annotate(
      "text",
      x = mean(c(0, right_limit)),
      y = -0.42,
      label = "Wastewater members within each ANI cluster",
      family = FONT,
      fontface = "bold",
      size = 4.8,
      vjust = 1
    ) +
    ggplot2::scale_fill_manual(
      name = NULL,
      values = c(
        "Clinical" = COL_CLINICAL,
        "Wastewater" = COL_WASTEWATER
      ),
      breaks = c("Clinical", "Wastewater"),
      labels = c(
        paste0("Clinical (", sum(df$Clinical), ")"),
        paste0("Wastewater (", sum(df$Wastewater), ")")
      )
    ) +
    ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::scale_y_continuous(
      breaks = df$y,
      labels = df$Cluster,
      expand = ggplot2::expansion(mult = c(0, 0))
    ) +
    ggplot2::guides(
      fill = ggplot2::guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(linewidth = 0)
      )
    ) +
    ggplot2::coord_cartesian(
      xlim = c(left_limit, right_limit),
      ylim = c(0.5, n_clusters + 0.5),
      clip = "off"
    ) +
    ggplot2::theme_classic(base_family = FONT, base_size = 12) +
    ggplot2::theme(
      text = ggplot2::element_text(family = FONT, color = "black"),
      legend.position = "top",
      legend.justification = "center",
      legend.direction = "horizontal",
      legend.text = ggplot2::element_text(
        family = FONT,
        size = 15,
        color = "black",
        margin = ggplot2::margin(r = 24)
      ),
      legend.key.width = grid::unit(0.72, "cm"),
      legend.key.height = grid::unit(0.72, "cm"),
      legend.spacing.x = grid::unit(0.3, "cm"),
      axis.title = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(
        family = FONT,
        face = "bold",
        size = 11.3,
        color = "black",
        margin = ggplot2::margin(r = 7)
      ),
      axis.text.x = ggplot2::element_text(
        family = FONT,
        face = "bold",
        size = 10.8,
        color = "black"
      ),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line.y = ggplot2::element_blank(),
      axis.line.x = ggplot2::element_line(linewidth = 0.65, color = "black"),
      panel.grid.major.x = ggplot2::element_line(
        color = "#D9D9D9",
        linewidth = 0.5,
        linetype = "dashed"
      ),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(t = 4, r = 38, b = 82, l = 78)
    )
  

  figure_cluster
}

run_fig4c <- function(input_file, output_dir) {
  if (missing(output_dir) || length(output_dir) != 1L ||
      is.na(output_dir) || !nzchar(output_dir)) {
    stop("Supply an output directory for Fig4c.", call. = FALSE)
  }
  df <- prepare_fig4c_data(input_file)
  figure_cluster <- plot_fig4c(df)
  if (!capabilities("cairo")) {
    stop("This R installation needs Cairo support for PDF output.", call. = FALSE)
  }
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  paths <- c(pdf = file.path(output_dir, "Fig4c.pdf"),
             png = file.path(output_dir, "Fig4c.png"),
             data = file.path(output_dir, "Fig4c_plot_data.tsv"))
  utils::write.table(df[, setdiff(names(df), "y"), drop = FALSE],
                     file = paths[["data"]], sep = "\t", quote = FALSE, row.names = FALSE)
  ggplot2::ggsave(
    filename = paths[["pdf"]], plot = figure_cluster, device = grDevices::cairo_pdf,
    width = 11.0, height = 10.2, units = "in", family = "Arial", bg = "white"
  )
  ggplot2::ggsave(
    filename = paths[["png"]], plot = figure_cluster, width = 11.0, height = 10.2,
    units = "in", dpi = 600, bg = "white"
  )
  message("Fig4c exported to: ", normalizePath(output_dir, winslash = "/"))
  invisible(list(plot = figure_cluster, data = df, files = paths))
}
  list(run = run_fig4c)
})

.fig4d <- local({
#!/usr/bin/env Rscript
# Figure 4d: representative clinical-wastewater pairs.
# Recovered from the author's 26 August 2026 drawing script and supplied table.
# Plotting only: no ANI/SNP/read analysis or pair-selection analysis is rerun here.
# Blank specimen metadata stay blank; Deep sputum is displayed as Sputum.
# All supplied counts and pair identities are preserved.
# Usage: Rscript Fig4d.R <input.tsv> <output_dir> [bottom|right]
# In R: source("Fig4d.R"); run_fig4d("pairs.tsv", "figures")


required_columns <- c(
  "Pair_ID", "Clinical_display_ID", "Clinical_specimen", "Site",
  "Clinical_city", "Wastewater_city", "same_city", "Wastewater_setting",
  "ANI_min", "fixed_SNP_Mbp", "callable_fraction", "Direct_remap",
  "inStrain_partial", "Clinical_ARG_count", "Wastewater_ARG_count",
  "Shared_ARG_count", "High_priority_ARG"
)

numeric_columns <- c(
  "ANI_min", "fixed_SNP_Mbp", "callable_fraction", "Direct_remap",
  "inStrain_partial", "Clinical_ARG_count", "Wastewater_ARG_count",
  "Shared_ARG_count"
)

setting_order <- c("Community", "Hospital", "Wet market")

read_pairs <- function(path) {
  if (!file.exists(path)) stop("Input file does not exist: ", path)
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c(""),
    fileEncoding = "UTF-8"
  )
}

validate_pairs <- function(df) {
  missing_columns <- setdiff(required_columns, names(df))
  if (length(missing_columns) > 0L) {
    stop("Missing required columns: ", paste(missing_columns, collapse = ", "))
  }
  if (nrow(df) != 25L) stop("Expected 25 rows; found ", nrow(df))
  if (anyNA(df$Pair_ID) || any(trimws(df$Pair_ID) == "")) {
    stop("Pair_ID contains missing or blank values")
  }
  if (anyDuplicated(df$Pair_ID)) stop("Pair_ID values must be unique")

  bad_settings <- setdiff(unique(df$Wastewater_setting), setting_order)
  if (length(bad_settings) > 0L) {
    stop("Unexpected wastewater settings: ", paste(bad_settings, collapse = ", "))
  }
  observed_counts <- as.integer(table(factor(df$Wastewater_setting, levels = setting_order)))
  if (!identical(observed_counts, c(14L, 8L, 3L))) {
    stop(
      "Expected setting counts Community/Hospital/Wet market = 14/8/3; found ",
      paste(observed_counts, collapse = "/")
    )
  }

  for (column in numeric_columns) {
    converted <- suppressWarnings(as.numeric(df[[column]]))
    invalid <- is.na(converted) & !is.na(df[[column]])
    if (any(invalid)) stop("Non-numeric values in column ", column)
  }
  if (any(!(df$same_city %in% c("Yes", "No")))) {
    stop("same_city must contain only Yes or No")
  }
  required_numeric <- setdiff(numeric_columns, c("Direct_remap", "inStrain_partial"))
  for (column in required_numeric) {
    if (any(!is.finite(as.numeric(df[[column]])))) stop("Missing or non-finite values in ", column)
  }
  if (any(as.numeric(df$fixed_SNP_Mbp) <= 0)) stop("SNP density must be positive for the log-scale display")
  if (any(as.numeric(df$fixed_SNP_Mbp) < 7.5 | as.numeric(df$fixed_SNP_Mbp) > 1000)) {
    stop("SNP density is outside the original display range 7.5 to 1000; revise the axis explicitly")
  }
  if (any(as.numeric(df$callable_fraction) < 0.5 | as.numeric(df$callable_fraction) > 0.9)) {
    stop("Callable fraction is outside the original symbol-size range 0.5 to 0.9")
  }
  for (column in c("Clinical_ARG_count", "Wastewater_ARG_count", "Shared_ARG_count")) {
    values <- as.numeric(df[[column]])
    if (any(values < 0 | values != floor(values))) stop("ARG counts must be non-negative integers: ", column)
  }
  if (any(as.numeric(df$Shared_ARG_count) > pmin(as.numeric(df$Clinical_ARG_count),
                                               as.numeric(df$Wastewater_ARG_count)))) {
    stop("Shared ARG count exceeds one or both genome ARG counts")
  }
  for (column in c("Direct_remap", "inStrain_partial")) {
    if (any(!is.na(df[[column]]) & !(as.numeric(df[[column]]) %in% c(0, 1)))) {
      stop("Read-support fields must contain 0, 1 or missing values: ", column)
    }
  }
  invisible(TRUE)
}

extract_subtypes <- function(section) {
  section <- sub("^[^:]+:\\s*", "", section)
  entries <- trimws(unlist(strsplit(section, ",", fixed = TRUE)))
  entries <- entries[nzchar(entries)]
  if (length(entries) == 0L) return("")
  entries <- sub(".*--", "", entries)
  paste(entries, collapse = ", ")
}

parse_high_priority <- function(x) {
  if (length(x) == 0L || is.na(x) || !nzchar(trimws(x))) {
    return(list(shared = "", clinical_only = ""))
  }
  sections <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))
  shared_sections <- sections[grepl("^shared\\s*:", sections, ignore.case = TRUE)]
  clinical_sections <- sections[grepl("^clinical-only\\s*:", sections, ignore.case = TRUE)]
  shared <- if (length(shared_sections)) extract_subtypes(shared_sections[[1L]]) else ""
  clinical_only <- if (length(clinical_sections)) extract_subtypes(clinical_sections[[1L]]) else ""
  list(shared = shared, clinical_only = clinical_only)
}

normalize_pairs <- function(df) {
  validate_pairs(df)
  out <- df
  for (column in numeric_columns) out[[column]] <- as.numeric(out[[column]])

  specimen <- trimws(out$Clinical_specimen)
  specimen[is.na(out$Clinical_specimen) | specimen == ""] <- NA_character_
  specimen[specimen %in% c("Deep sputum", "Sputum")] <- "Sputum"
  allowed_specimens <- c("Oropharyngeal swab", "Sputum")
  unexpected <- setdiff(unique(stats::na.omit(specimen)), allowed_specimens)
  if (length(unexpected) > 0L) {
    stop("Unexpected clinical specimen values: ", paste(unexpected, collapse = ", "))
  }
  out$Clinical_specimen <- specimen
  out$Clinical_specimen_display <- ifelse(is.na(specimen), "", specimen)

  parsed <- lapply(out$High_priority_ARG, parse_high_priority)
  out$HP_shared <- vapply(parsed, `[[`, character(1), "shared")
  out$HP_clinical_only <- vapply(parsed, `[[`, character(1), "clinical_only")
  out$Wastewater_setting <- factor(out$Wastewater_setting, levels = setting_order)

  out <- out[order(out$Wastewater_setting, out$fixed_SNP_Mbp, out$Pair_ID), , drop = FALSE]
  rownames(out) <- NULL
  out$row_index <- seq_len(nrow(out))
  out
}

fig4d_palette <- list(
  Community = "#86BFC0",
  Hospital = "#C56D72",
  `Wet market` = "#E5BE4D"
)

fig4d_backgrounds <- list(
  Community = "#EEF7F7",
  Hospital = "#FBF1F1",
  `Wet market` = "#FFF9E8"
)

fig4d_shapes <- c(Community = 21L, Hospital = 22L, `Wet market` = 24L)

fig4d_typography <- list(
  section_header = 17.5,
  group_label = 15.0,
  column_header = 14.25,
  dense_header = 13.5,
  panel_header = 13.75,
  hp_key = 12.25,
  row_text = 13.375,
  row_meta = 13.125,
  count_value = 11.5,
  hp_annotation = 12.25,
  axis_tick = 14.0,
  axis_title = 15.25,
  legend_title = 15.0,
  legend_setting = 14.375,
  legend_value = 13.5,
  legend_specimen = 13.25
)

fig4d_export_geometry <- list(width_in = 17.5, height_in = 13, dpi = 300)

fig4d_layout <- list(
  data_top = 0.8625,
  data_bottom = 0.1700,
  axis_y = 0.145,
  axis_tick_y = 0.124,
  axis_title_y = 0.090,
  hp_line_offset = 0.24,
  read_sep = 0.522,
  table_right = 0.815,
  legend_left = 0.835
)

fig4d_header_x <- list(
  clinical = 0.113,
  specimen = 0.169,
  site = 0.215,
  remap = 0.445,
  instrain = 0.494,
  clinical_count = 0.550,
  wastewater_count = 0.615,
  shared_count = 0.680
)

fig4d_arg_x <- list(
  clinical_count = 0.550,
  wastewater_count = 0.615,
  shared_count = 0.680,
  hp_symbol = 0.711,
  hp_text = 0.721,
  hp_header = 0.765,
  key_shared_symbol = 0.713,
  key_shared_text = 0.721,
  key_clinical_symbol = 0.762,
  key_clinical_text = 0.770
)

fig4d_text <- function(label, x, y, size = 9, face = "plain", colour = "#111111",
                       just = "centre", rot = 0, lineheight = 0.9) {
  grid::grid.text(
    label,
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    just = just,
    rot = rot,
    gp = grid::gpar(
      fontfamily = "Arial",
      fontsize = size,
      fontface = face,
      col = colour,
      lineheight = lineheight
    )
  )
}

fig4d_point <- function(x, y, pch, fill, size_mm, colour = "#333333", lwd = 0.8) {
  grid::grid.points(
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    pch = pch,
    size = grid::unit(size_mm, "mm"),
    gp = grid::gpar(fill = fill, col = colour, lwd = lwd)
  )
}

fig4d_count_circle <- function(x, y, value, fill, border = "#596062", text_colour = "#202426",
                               radius_mm = 2.55) {
  grid::grid.circle(
    x = grid::unit(x, "npc"),
    y = grid::unit(y, "npc"),
    r = grid::unit(radius_mm, "mm"),
    gp = grid::gpar(fill = fill, col = border, lwd = 0.75)
  )
  fig4d_text(
    as.character(value), x, y,
    size = fig4d_typography$count_value,
    colour = text_colour
  )
}

draw_fig4d <- function(df, legend_position = c("bottom", "right")) {
  legend_position <- match.arg(legend_position)
  if (!inherits(df$Wastewater_setting, "factor")) df <- normalize_pairs(df)
  if (nrow(df) != 25L) stop("draw_fig4d requires exactly 25 normalized rows")

  table_left <- 0.025
  table_right <- fig4d_layout$table_right
  legend_left <- fig4d_layout$legend_left
  data_top <- fig4d_layout$data_top
  data_bottom <- fig4d_layout$data_bottom
  row_height <- (data_top - data_bottom) / nrow(df)
  row_y <- data_top - (seq_len(nrow(df)) - 0.5) * row_height

  x <- list(
    setting = 0.057,
    clinical = 0.116,
    specimen = 0.168,
    site = 0.209,
    source_sep = 0.226,
    same_city = 0.248,
    ani = 0.283,
    snp_left = 0.307,
    snp_right = 0.418,
    core_sep = 0.426,
    remap = 0.448,
    instrain = 0.486,
    read_sep = fig4d_layout$read_sep,
    clinical_count = fig4d_arg_x$clinical_count,
    wastewater_count = fig4d_arg_x$wastewater_count,
    shared_count = fig4d_arg_x$shared_count,
    hp_symbol = fig4d_arg_x$hp_symbol,
    hp_text = fig4d_arg_x$hp_text
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(x = 0, just = "left",
    width = if (legend_position == "bottom") 1.2 else 1,
    xscale = c(0, 1), yscale = c(0, 1)))
  grid::grid.rect(gp = grid::gpar(fill = "white", col = NA))

  for (setting in setting_order) {
    idx <- which(as.character(df$Wastewater_setting) == setting)
    if (length(idx) == 0L) next
    top <- row_y[min(idx)] + row_height / 2
    bottom <- row_y[max(idx)] - row_height / 2
    grid::grid.rect(
      x = grid::unit((table_left + table_right) / 2, "npc"),
      y = grid::unit((top + bottom) / 2, "npc"),
      width = grid::unit(table_right - table_left, "npc"),
      height = grid::unit(top - bottom, "npc"),
      gp = grid::gpar(fill = fig4d_backgrounds[[setting]], col = NA)
    )
    count <- length(idx)
    pair_word <- if (count == 1L) "pair" else "pairs"
    fig4d_text(
      paste0(setting, "\n", count, " ", pair_word),
      x$setting,
      (top + bottom) / 2,
      size = fig4d_typography$group_label,
      face = "bold",
      lineheight = 0.92
    )
  }

  section_box <- function(label, left, right) {
    grid::grid.rect(
      x = grid::unit((left + right) / 2, "npc"),
      y = grid::unit(0.965, "npc"),
      width = grid::unit(right - left, "npc"),
      height = grid::unit(0.055, "npc"),
      gp = grid::gpar(fill = "white", col = "#202020", lwd = 0.9)
    )
    fig4d_text(
      label, (left + right) / 2, 0.965,
      size = fig4d_typography$section_header,
      face = "bold"
    )
  }

  section_box("Source context", table_left, x$source_sep - 0.003)
  section_box("Core-genome\nrelatedness", x$source_sep + 0.002, x$core_sep - 0.003)
  section_box("Read-level\nsupport", x$core_sep + 0.002, x$read_sep - 0.003)
  section_box("ARG\nconcordance", x$read_sep + 0.002, table_right)

  header_y <- 0.916
  fig4d_text("Setting", x$setting, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("Clinical\nID", fig4d_header_x$clinical, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("Specimen", fig4d_header_x$specimen, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("Site", fig4d_header_x$site, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("Same\ncity", x$same_city, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("ANI\n(min)", x$ani, header_y, size = fig4d_typography$column_header, face = "bold")
  fig4d_text("Fixed core-genome\nSNPs per Mbp\n(log10 scale)",
             (x$snp_left + x$snp_right) / 2, header_y,
             size = fig4d_typography$dense_header, face = "bold")
  fig4d_text("Direct\nremap", fig4d_header_x$remap, header_y, size = fig4d_typography$panel_header, face = "bold")
  fig4d_text("inStrain\n(partial)", fig4d_header_x$instrain, header_y, size = fig4d_typography$panel_header, face = "bold")
  fig4d_text("Clinical\nARG count", fig4d_header_x$clinical_count, header_y, size = fig4d_typography$panel_header, face = "bold")
  fig4d_text("Wastewater\nARG count", fig4d_header_x$wastewater_count, header_y, size = fig4d_typography$panel_header, face = "bold")
  fig4d_text("Shared\nARG count", fig4d_header_x$shared_count, header_y, size = fig4d_typography$panel_header, face = "bold")
  fig4d_text("High-priority\nARG", fig4d_arg_x$hp_header, 0.916,
             size = fig4d_typography$panel_header, face = "bold")

  fig4d_point(fig4d_arg_x$key_shared_symbol, 0.891, 18, "#4F8888", 2.2,
              colour = "#4F8888")
  fig4d_text("Shared", fig4d_arg_x$key_shared_text, 0.891,
             size = fig4d_typography$hp_key, colour = "#4F8888", just = "left")
  fig4d_point(fig4d_arg_x$key_clinical_symbol, 0.891, 18, "#A6524C", 2.2,
              colour = "#A6524C")
  fig4d_text("Clinical only", fig4d_arg_x$key_clinical_text, 0.891,
             size = fig4d_typography$hp_key, colour = "#A6524C", just = "left")

  for (separator in c(x$source_sep, x$core_sep, x$read_sep)) {
    grid::grid.lines(
      x = grid::unit(c(separator, separator), "npc"),
      y = grid::unit(c(data_bottom, 0.937), "npc"),
      gp = grid::gpar(col = "#202020", lwd = 0.8, lty = "dashed")
    )
  }

  snp_min <- 7.5
  snp_max <- 1000
  snp_x <- function(value) {
    x$snp_left + (log10(value) - log10(snp_min)) /
      (log10(snp_max) - log10(snp_min)) * (x$snp_right - x$snp_left)
  }
  callable_size <- function(value) 2.6 + (value - 0.5) / 0.4 * 3.5

  for (i in seq_len(nrow(df))) {
    setting <- as.character(df$Wastewater_setting[i])
    colour <- fig4d_palette[[setting]]
    shape <- unname(fig4d_shapes[[setting]])
    y <- row_y[i]

    fig4d_text(df$Clinical_display_ID[i], x$clinical, y, size = fig4d_typography$row_text)
    fig4d_text(df$Site[i], x$site, y, size = fig4d_typography$row_text)

    if (!is.na(df$Clinical_specimen[i])) {
      if (df$Clinical_specimen[i] == "Oropharyngeal swab") {
        fig4d_point(x$specimen, y, 21, "white", 3.2, colour = "#111111", lwd = 0.8)
      } else if (df$Clinical_specimen[i] == "Sputum") {
        fig4d_point(x$specimen, y, 19, "#111111", 3.2, colour = "#111111", lwd = 0.8)
      }
    }

    same_colour <- if (df$same_city[i] == "Yes") "#C55359" else "#899093"
    fig4d_text(df$same_city[i], x$same_city, y, size = fig4d_typography$row_meta,
               face = "bold", colour = same_colour)
    fig4d_text(sprintf("%.1f", df$ANI_min[i]), x$ani, y, size = fig4d_typography$row_meta)

    target_x <- snp_x(df$fixed_SNP_Mbp[i])
    grid::grid.lines(
      x = grid::unit(c(x$snp_left, target_x), "npc"),
      y = grid::unit(c(y, y), "npc"),
      gp = grid::gpar(col = "#CBD1D1", lwd = 0.7)
    )
    fig4d_point(target_x, y, shape, colour, callable_size(df$callable_fraction[i]),
                colour = "#485052", lwd = 0.75)

    if (!is.na(df$Direct_remap[i]) && df$Direct_remap[i] > 0) {
      fig4d_point(x$remap, y, shape, colour, 3.2, colour = "#485052", lwd = 0.75)
    }
    if (!is.na(df$inStrain_partial[i]) && df$inStrain_partial[i] > 0) {
      fig4d_point(x$instrain, y, shape, colour, 3.2, colour = "#485052", lwd = 0.75)
    }

    fig4d_count_circle(x$clinical_count, y, df$Clinical_ARG_count[i], "#FFFFFF")
    fig4d_count_circle(x$wastewater_count, y, df$Wastewater_ARG_count[i], "#EFF5F5")
    fig4d_count_circle(x$shared_count, y, df$Shared_ARG_count[i], colour,
                       border = "#596062", text_colour = "#203234")

    has_shared <- nzchar(df$HP_shared[i])
    has_clinical <- nzchar(df$HP_clinical_only[i])
    if (has_shared) {
      shared_y <- if (has_clinical) y + row_height * fig4d_layout$hp_line_offset else y
      fig4d_point(x$hp_symbol, shared_y, 18, "#4F8888", 2.3, colour = "#4F8888")
      fig4d_text(paste(strwrap(df$HP_shared[i], width = 20), collapse = "\n"), x$hp_text, shared_y,
                 size = fig4d_typography$hp_annotation, face = "italic",
                 colour = "#4F8888", just = "left")
    }
    if (has_clinical) {
      clinical_y <- if (has_shared) y - row_height * fig4d_layout$hp_line_offset else y
      fig4d_point(x$hp_symbol, clinical_y, 18, "#A6524C", 2.3, colour = "#A6524C")
      fig4d_text(paste(strwrap(df$HP_clinical_only[i], width = 20), collapse = "\n"), x$hp_text, clinical_y,
                 size = fig4d_typography$hp_annotation,
                 face = "italic", colour = "#A6524C", just = "left")
    }
  }

  axis_y <- fig4d_layout$axis_y
  grid::grid.lines(
    x = grid::unit(c(x$snp_left, x$snp_right), "npc"),
    y = grid::unit(c(axis_y, axis_y), "npc"),
    gp = grid::gpar(col = "#202020", lwd = 0.8)
  )
  for (tick in c(10, 100, 1000)) {
    tick_x <- snp_x(tick)
    grid::grid.lines(
      x = grid::unit(c(tick_x, tick_x), "npc"),
      y = grid::unit(c(axis_y, axis_y - 0.006), "npc"),
      gp = grid::gpar(col = "#202020", lwd = 0.8)
    )
    fig4d_text(
      format(tick, big.mark = ",", scientific = FALSE), tick_x,
      fig4d_layout$axis_tick_y,
      size = fig4d_typography$axis_tick
    )
  }
  fig4d_text("Fixed core-genome SNPs/Mbp", (x$snp_left + x$snp_right) / 2,
             fig4d_layout$axis_title_y,
             size = fig4d_typography$axis_title, face = "bold")

  if (legend_position == "right") {
  fig4d_text("Wastewater setting", legend_left, 0.826, size = fig4d_typography$legend_title,
             face = "bold", just = "left")
  legend_y <- c(0.787, 0.752, 0.717)
  for (j in seq_along(setting_order)) {
    setting <- setting_order[j]
    fig4d_point(legend_left + 0.009, legend_y[j], unname(fig4d_shapes[[setting]]),
                fig4d_palette[[setting]], 4.2, colour = "#485052")
    fig4d_text(setting, legend_left + 0.028, legend_y[j],
               size = fig4d_typography$legend_setting, just = "left")
  }

  fig4d_text("Callable fraction", legend_left, 0.645, size = fig4d_typography$legend_title,
             face = "bold", just = "left")
  callable_values <- c(0.5, 0.7, 0.9)
  callable_y <- c(0.605, 0.557, 0.501)
  for (j in seq_along(callable_values)) {
    fig4d_point(legend_left + 0.012, callable_y[j], 21, "#111111",
                callable_size(callable_values[j]), colour = "#111111")
    fig4d_text(sprintf("%.1f", callable_values[j]), legend_left + 0.036,
               callable_y[j], size = fig4d_typography$legend_value, just = "left")
  }

  fig4d_text("Clinical specimen", legend_left, 0.430, size = fig4d_typography$legend_title,
             face = "bold", just = "left")
  fig4d_point(legend_left + 0.008, 0.389, 21, "white", 3.8, colour = "#111111")
  fig4d_text("Oropharyngeal swab", legend_left + 0.027, 0.389,
             size = fig4d_typography$legend_specimen, just = "left")
  fig4d_point(legend_left + 0.008, 0.350, 19, "#111111", 3.8, colour = "#111111")
  fig4d_text("Sputum", legend_left + 0.027, 0.350,
             size = fig4d_typography$legend_specimen, just = "left")

  }
  grid::popViewport()
  if (legend_position == "bottom") {
    # Bottom legends reproduce the assembly arrangement; right retains the recovered original.
    fig4d_text("Wastewater setting", 0.030, 0.055, size = 14.0, face = "bold", just = "left")
    for (j in seq_along(setting_order)) {
      setting <- setting_order[j]
      px <- c(0.225, 0.415, 0.595)[j]
      fig4d_point(px, 0.055, unname(fig4d_shapes[[setting]]), fig4d_palette[[setting]],
                  4.2, colour = "#485052")
      fig4d_text(setting, px + 0.018, 0.055, size = 13.5, just = "left")
    }
    fig4d_text("Callable fraction", 0.030, 0.020, size = 14.0, face = "bold", just = "left")
    for (j in seq_along(c(0.5, 0.7, 0.9))) {
      value <- c(0.5, 0.7, 0.9)[j]
      px <- c(0.225, 0.320, 0.415)[j]
      fig4d_point(px, 0.020, 21, "#111111", callable_size(value), colour = "#111111")
      fig4d_text(sprintf("%.1f", value), px + 0.018, 0.020, size = 13.0, just = "left")
    }
    fig4d_text("Clinical specimen", 0.505, 0.020, size = 14.0, face = "bold", just = "left")
    fig4d_point(0.690, 0.020, 21, "white", 3.8, colour = "#111111")
    fig4d_text("Oropharyngeal swab", 0.704, 0.020, size = 12.5, just = "left")
    fig4d_point(0.903, 0.020, 19, "#111111", 3.8, colour = "#111111")
    fig4d_text("Sputum", 0.918, 0.020, size = 12.5, just = "left")
  }
  invisible(TRUE)
}

export_fig4d <- function(df, output_stem, legend_position = c("bottom", "right")) {
  legend_position <- match.arg(legend_position)
  if (!requireNamespace("svglite", quietly = TRUE)) stop("Install the R package 'svglite' before running")
  width_in <- if (legend_position == "bottom") 14.5 else fig4d_export_geometry$width_in
  output_dir <- dirname(output_stem)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  png_path <- paste0(output_stem, ".png")
  pdf_path <- paste0(output_stem, ".pdf")
  svg_path <- paste0(output_stem, ".svg")

  grDevices::png(
    png_path,
    width = width_in,
    height = fig4d_export_geometry$height_in,
    units = "in",
    res = fig4d_export_geometry$dpi,
    bg = "white",
    type = "cairo"
  )
  tryCatch(draw_fig4d(df, legend_position), finally = grDevices::dev.off())

  grDevices::cairo_pdf(
    pdf_path,
    width = width_in,
    height = fig4d_export_geometry$height_in,
    bg = "white"
  )
  tryCatch(draw_fig4d(df, legend_position), finally = grDevices::dev.off())

  if (!requireNamespace("svglite", quietly = TRUE)) stop("R package 'svglite' is required")
  svglite::svglite(
    svg_path,
    width = width_in,
    height = fig4d_export_geometry$height_in,
    bg = "white"
  )
  tryCatch(draw_fig4d(df, legend_position), finally = grDevices::dev.off())

  c(png_path, pdf_path, svg_path)
}

run_fig4d <- function(input_file, output_dir, legend_position = c("bottom", "right")) {
  legend_position <- match.arg(legend_position)
  raw <- read_pairs(input_file)
  normalized <- normalize_pairs(raw)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  paths <- export_fig4d(normalized, file.path(output_dir, "Fig4d"), legend_position)
  utils::write.table(normalized, file.path(output_dir, "Fig4d_displayed_pairs.tsv"),
                     sep = "\t", quote = FALSE, row.names = FALSE, na = "", fileEncoding = "UTF-8")
  message("Figure 4d written: ", paste(paths, collapse = ", "))
  invisible(list(data = normalized, files = paths))
}
  list(run = run_fig4d)
})

.fig4e <- local({
# Fig. 4e: cross-source carbapenemase evidence and contig context.
# Recovered plotting source: plot_Figure_5_carpen.R (5 August 2026).
# Required packages: ggplot2, patchwork. No packages are installed automatically.
#
# These inputs are the recovered Figure_5 summary tables, NOT raw contig calls.
# The figure displays KPC-2 and OXA-48; the input also retains the original NDM-5.
# The context denominator is the sum of three displayed resolved classes: 11,986.
# The older five-class table included 6 unresolved and 2 conflicting records
# (11,994 total). The recovered three-class table excludes those 8 records.
# The taxonomic denominator is independently 6,920 host-assigned contigs.
# A recovered upstream contig table independently confirms these counts.
# Taxonomy includes one host-assigned mixed-context record omitted from the
# context panel; the two denominators intentionally describe different subsets.
# The recovered label "Aeromonas spp." groups A. media, A. veronii and
# A. hydrophila only (519 + 431 + 146 = 1,096), not all Aeromonas species.
# This portable plotting entry point uses the original summary tables; it does
# not rerun AMRFinder, taxonomic annotation, or the upstream contig merges.

prepare_fig4e <- function(data_dir) {
  read_input <- function(filename, required) {
    path <- file.path(data_dir, filename)
    if (!file.exists(path)) stop("Missing Fig. 4e input: ", path, call. = FALSE)
    x <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE,
                          quote = "", comment.char = "", fileEncoding = "UTF-8")
    if (!all(required %in% names(x)))
      stop("Missing columns in ", filename, ": ",
           paste(setdiff(required, names(x)), collapse = ", "), call. = FALSE)
    x
  }
  evidence <- read_input("e_carbapenemase_overlap.tsv",
                         c("Evidence", "Allele", "Supported", "Support_score"))
  context <- read_input("e_mobile_genetic_context.tsv", c("Context", "Contigs", "Percent"))
  taxa <- read_input("e_selected_taxa.tsv", c("Taxon_group", "Contigs", "Percent_of_classified"))

  evidence_order <- c("Clinical dataset", "Wastewater dataset", "Pair-level shared", "Candidate KP MAG host")
  display_alleles <- c("KPC-2", "OXA-48")
  displayed <- evidence[evidence$Allele %in% display_alleles, , drop = FALSE]
  expected_keys <- as.vector(outer(evidence_order, display_alleles, paste, sep = "|"))
  actual_keys <- paste(displayed$Evidence, displayed$Allele, sep = "|")
  if (nrow(displayed) != 8L || anyDuplicated(actual_keys) ||
      !setequal(actual_keys, expected_keys) ||
      anyNA(displayed$Supported) || !all(displayed$Supported %in% c("Yes", "No")))
    stop("The two displayed alleles must each have four unique Yes/No evidence records.", call. = FALSE)
  if (anyNA(displayed$Support_score) ||
      any(displayed$Support_score != as.integer(displayed$Supported == "Yes")))
    stop("Support_score disagrees with Supported.", call. = FALSE)

  context_names <- c("Chromosome associated", "Plasmid associated", "Phage associated")
  taxon_names <- c("Klebsiella pneumoniae", "Pseudomonas aeruginosa", "Enterobacter hormaechei",
                   "Klebsiella quasipneumoniae", "Aeromonas spp.", "Acinetobacter baumannii",
                   "Other identified taxa")
  check_counts <- function(x, name_col, expected_names) {
    if (nrow(x) != length(expected_names) || anyDuplicated(x[[name_col]]) ||
        !setequal(x[[name_col]], expected_names))
      stop("Unexpected or duplicate categories in ", name_col, call. = FALSE)
    if (!is.numeric(x$Contigs) || anyNA(x$Contigs) || any(!is.finite(x$Contigs)) ||
        any(x$Contigs < 0) || any(x$Contigs != floor(x$Contigs)) || sum(x$Contigs) == 0)
      stop("Contig counts must be nonnegative integers with a positive total.", call. = FALSE)
  }
  check_counts(context, "Context", context_names)
  check_counts(taxa, "Taxon_group", taxon_names)
  context <- context[match(context_names, context$Context), , drop = FALSE]
  taxa <- taxa[match(taxon_names, taxa$Taxon_group), , drop = FALSE]
  context$Percent_recalculated <- 100 * context$Contigs / sum(context$Contigs)
  taxa$Percent_recalculated <- 100 * taxa$Contigs / sum(taxa$Contigs)
  context$Context_plot <- factor(c("Chromosome", "Plasmid", "Phage"),
                                 levels = c("Phage", "Plasmid", "Chromosome"))
  taxa$Taxon_plot <- factor(c("K. pneumoniae", "P. aeruginosa", "E. hormaechei",
                              "K. quasipneumoniae", "Aeromonas spp.", "A. baumannii",
                              "Other identified taxa"),
                            levels = rev(c("K. pneumoniae", "P. aeruginosa", "E. hormaechei",
                                           "K. quasipneumoniae", "Aeromonas spp.", "A. baumannii",
                                           "Other identified taxa")))
  displayed$Evidence <- factor(displayed$Evidence, levels = evidence_order)
  displayed$Allele <- factor(displayed$Allele, levels = rev(display_alleles))
  displayed$Supported <- factor(displayed$Supported, levels = c("No", "Yes"))
  list(evidence_all = evidence, evidence = displayed, context = context, taxa = taxa,
       denominators = c(context = sum(context$Contigs), taxonomy = sum(taxa$Contigs)))
}

run_fig4e <- function(data_dir, output_dir) {
  dependencies <- c("ggplot2", "patchwork")
  missing <- dependencies[!vapply(dependencies, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  d <- prepare_fig4e(data_dir)
  comma <- function(x) format(x, big.mark = ",", scientific = FALSE, trim = TRUE)
  labels <- function(n, p) paste0(comma(n), " (", sprintf("%.1f", p), "%)")
  main_col <- "#66879B"
  text_col <- "#2B2B2B"
  theme_common <- ggplot2::theme_minimal(base_size = 9, base_family = "Arial") +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.major.x = ggplot2::element_line(colour = "#E4E8EB", linewidth = 0.32),
                   axis.title = ggplot2::element_text(size = 8.5, colour = text_col),
                   axis.text = ggplot2::element_text(size = 8.5, colour = text_col),
                   plot.title = ggplot2::element_text(size = 10.5, face = "bold", hjust = 0),
                   plot.subtitle = ggplot2::element_text(size = 8, colour = "#606060"),
                   plot.margin = ggplot2::margin(5, 8, 5, 5))
  p_evidence <- ggplot2::ggplot(d$evidence, ggplot2::aes(Evidence, Allele, fill = Supported)) +
    ggplot2::geom_point(shape = 21, size = 6.7, stroke = 0.9, colour = main_col) +
    ggplot2::scale_fill_manual(values = c(No = "white", Yes = main_col), drop = FALSE, name = NULL) +
    ggplot2::scale_x_discrete(labels = c("Clinical", "Wastewater", "Shared", "MAG host"),
                              expand = ggplot2::expansion(add = c(0.3, 0.3))) +
    ggplot2::scale_y_discrete(labels = c("OXA-48" = "blaOXA-48", "KPC-2" = "blaKPC-2")) +
    ggplot2::labs(title = "Cross-source carbapenemase evidence", x = NULL, y = NULL) +
    theme_common +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_line(colour = "#E4E8EB", linewidth = 0.32),
                   axis.text.y = ggplot2::element_text(face = "italic"),
                   legend.position = "top", legend.justification = "left",
                   legend.margin = ggplot2::margin(0, 0, 0, 0)) +
    ggplot2::guides(fill = ggplot2::guide_legend(override.aes = list(size = 3.2)))

  d$context$label <- labels(d$context$Contigs, d$context$Percent_recalculated)
  d$context$inside <- d$context$Contigs >= 5000
  d$context$label_colour <- ifelse(d$context$inside, "white", text_col)
  d$context$hjust <- ifelse(d$context$inside, 1.04, -0.08)
  d$context$fill <- ifelse(d$context$Context_plot == "Phage", "#EEF3F5", main_col)
  p_context <- ggplot2::ggplot(d$context, ggplot2::aes(Contigs, Context_plot)) +
    ggplot2::geom_col(ggplot2::aes(fill = fill), width = 0.55) +
    ggplot2::geom_text(ggplot2::aes(label = label, hjust = hjust, colour = label_colour),
                       size = 2.8, family = "Arial", fontface = "bold") +
    ggplot2::scale_fill_identity() + ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(labels = comma, breaks = c(0, 2000, 4000, 6000),
                               limits = c(0, 6300), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = "Predicted genetic context",
                   subtitle = paste(comma(d$denominators[["context"]]), "contigs"),
                   x = "Carbapenemase-carrying contigs", y = NULL) + theme_common

  d$taxa$label <- labels(d$taxa$Contigs, d$taxa$Percent_recalculated)
  d$taxa$inside <- d$taxa$Contigs >= 2500
  d$taxa$label_colour <- ifelse(d$taxa$inside, "white", text_col)
  d$taxa$hjust <- ifelse(d$taxa$inside, 1.03, -0.06)
  d$taxa$fill <- ifelse(d$taxa$Taxon_group == "Klebsiella pneumoniae", "#3F657C",
                        ifelse(d$taxa$Taxon_group == "Other identified taxa", "#8AA0AE", "#B8C7D0"))
  p_taxonomy <- ggplot2::ggplot(d$taxa, ggplot2::aes(Contigs, Taxon_plot)) +
    ggplot2::geom_col(ggplot2::aes(fill = fill), width = 0.56) +
    ggplot2::geom_text(ggplot2::aes(label = label, hjust = hjust, colour = label_colour),
                       size = 2.8, family = "Arial", fontface = "bold") +
    ggplot2::scale_fill_identity() + ggplot2::scale_colour_identity() +
    ggplot2::scale_x_continuous(labels = comma, breaks = c(0, 1000, 2000, 3000),
                               limits = c(0, 3300), expand = ggplot2::expansion(mult = c(0, 0))) +
    ggplot2::labs(title = "Taxonomic assignments",
                   subtitle = paste(comma(d$denominators[["taxonomy"]]), "species-assigned contigs"),
                   x = "Carbapenemase-carrying contigs", y = NULL) + theme_common

  plot <- patchwork::wrap_plots(p_evidence, p_context, p_taxonomy, ncol = 1,
                                heights = c(0.95, 1.15, 1.8))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  files <- file.path(output_dir, c("Fig4e.pdf", "Fig4e.png"))
  ggplot2::ggsave(files[[1L]], plot, width = 5.1, height = 8.2, units = "in", device = grDevices::cairo_pdf)
  ggplot2::ggsave(files[[2L]], plot, width = 5.1, height = 8.2, units = "in", dpi = 300, bg = "white")
  utils::write.table(d$context, file.path(output_dir, "Fig4e_context_recalculated.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(d$taxa, file.path(output_dir, "Fig4e_taxa_recalculated.tsv"),
                     sep = "\t", row.names = FALSE, quote = FALSE)
  invisible(list(plot = plot, data = d, files = files))
}
  list(run = run_fig4e)
})

run_fig4_cde <- function(data_dir, output_dir, panels = c("c", "d", "e"),
                        legend_position = c("bottom", "right")) {
  legend_position <- match.arg(legend_position)
  if (!dir.exists(data_dir)) stop("Input directory does not exist: ", data_dir)
  if (!length(panels) || any(!panels %in% c("c", "d", "e"))) stop("Panels must be c, d and/or e")
  panels <- unique(panels)
  results <- list()
  if ("c" %in% panels) results$c <- .fig4c$run(file.path(data_dir, "Fig4c_cluster_data.tsv"),
                                               file.path(output_dir, "c"))
  if ("d" %in% panels) results$d <- .fig4d$run(file.path(data_dir, "Fig4d_pairs.tsv"),
                                               file.path(output_dir, "d"), legend_position)
  if ("e" %in% panels) results$e <- .fig4e$run(data_dir, file.path(output_dir, "e"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  writeLines(capture.output(utils::sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  invisible(results)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args) < 2L || length(args) > 3L) {
    stop("Usage: Rscript Fig4_cde.R DATA_DIRECTORY OUTPUT_DIRECTORY [c,d,e]")
  }
  panels <- if (length(args) == 3L) strsplit(args[3], ",", fixed = TRUE)[[1]] else c("c", "d", "e")
  run_fig4_cde(args[1], args[2], panels)
}
