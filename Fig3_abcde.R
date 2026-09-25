
.fig3ab_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

.fig3a_pool_counts <- function(counts, sample_totals, metadata, group_order = NULL) {
  if (!is.matrix(counts) || !is.numeric(counts) || !nrow(counts) || !ncol(counts) ||
      any(!is.finite(counts)) || any(counts < 0 | counts != round(counts)))
    stop('Species input must contain nonnegative integer counts, not relative values.')
  if (is.null(rownames(counts)) || is.null(colnames(counts)) ||
      anyNA(rownames(counts)) || anyNA(colnames(counts)) ||
      any(!nzchar(rownames(counts))) || any(!nzchar(colnames(counts))) ||
      anyDuplicated(rownames(counts)) || anyDuplicated(colnames(counts)))
    stop('Species and sample count-matrix identifiers must be complete and unique.')
  if (!all(c('sample','group') %in% names(metadata))) stop('Metadata requires sample and group.')
  m <- as.data.frame(metadata[c('sample','group')], stringsAsFactors = FALSE)
  m[] <- lapply(m, as.character)
  if (!nrow(m) || anyNA(m) || any(!nzchar(m$sample)) || any(!nzchar(m$group)) || anyDuplicated(m$sample))
    stop('Metadata sample IDs must be unique and group assignments complete.')
  if (any(!m$sample %in% colnames(counts))) stop('Missing sample count columns.')
  if (!all(c('sample','total_species_contigs') %in% names(sample_totals)))
    stop('The denominator table requires sample and total_species_contigs.')
  totals <- as.data.frame(sample_totals[c('sample','total_species_contigs')], stringsAsFactors = FALSE)
  totals$sample <- as.character(totals$sample)
  if (anyNA(totals$sample) || any(!nzchar(totals$sample)) || anyDuplicated(totals$sample))
    stop('Sample denominator identifiers must be complete and unique.')
  if (!is.numeric(totals$total_species_contigs) || any(!is.finite(totals$total_species_contigs)) ||
      any(totals$total_species_contigs < 0 | totals$total_species_contigs != round(totals$total_species_contigs)))
    stop('Sample denominators must be nonnegative integer counts.')
  idx <- match(m$sample, totals$sample)
  if (anyNA(idx)) stop('Missing complete species-assigned contig denominator for a sample.')
  den <- totals$total_species_contigs[idx]
  x <- counts[,m$sample,drop=FALSE]
  if (any(colSums(x) > den)) stop('Displayed species counts exceed the complete sample denominator.')
  if (is.null(group_order)) group_order <- unique(m$group)
  if (anyNA(group_order) || anyDuplicated(group_order) || !setequal(group_order, m$group))
    stop('Group order must contain each observed group once.')
  pooled <- matrix(0, nrow(x), length(group_order), dimnames=list(rownames(x),group_order))
  group_den <- numeric(length(group_order))
  for (i in seq_along(group_order)) {
    rows <- which(m$group == group_order[i])
    pooled[,i] <- rowSums(x[,rows,drop=FALSE])
    group_den[i] <- sum(den[rows])
  }
  proportions <- sweep(pooled, 2, ifelse(group_den > 0, group_den, NA_real_), '/')
  list(counts=pooled, proportions=proportions,
       denominators=data.frame(group=group_order, total_species_contigs=group_den,
                               stringsAsFactors=FALSE))
}

run_fig3a <- function(data_dir, out_dir, render = TRUE) {
  .fig3ab_require(c("data.table", "dplyr", "tibble", "ComplexHeatmap", "circlize", "RColorBrewer"))

  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(RColorBrewer)
  })

  root_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  species_total_fp <- file.path(root_dir, "Fig3a_species_total_counts.tsv")
  mat_counts_fp    <- file.path(root_dir, "Fig3a_species_by_sample_counts.tsv")
  sample_totals_fp <- file.path(root_dir, "Fig3a_sample_species_contig_totals.tsv")
  arg_comp_fp      <- file.path(root_dir, "Fig3a_species_by_ARG_type_relative.tsv")
  meta_fp          <- file.path(root_dir, "Fig3a_sample_metadata.tsv")

  cat("Checking input files...\n")
  print(data.frame(
    file = c("Fig3a_species_total_counts.tsv",
             "Fig3a_species_by_sample_counts.tsv",
             "Fig3a_sample_species_contig_totals.tsv",
             "Fig3a_species_by_ARG_type_relative.tsv",
             "Fig3a_sample_metadata.tsv"),
    exists = c(file.exists(species_total_fp),
               file.exists(mat_counts_fp),
               file.exists(sample_totals_fp),
               file.exists(arg_comp_fp),
               file.exists(meta_fp))
  ))
  stopifnot(file.exists(species_total_fp),
            file.exists(mat_counts_fp),
            file.exists(sample_totals_fp),
            file.exists(arg_comp_fp),
            file.exists(meta_fp))

  species_total <- fread(species_total_fp)
  mat_counts    <- fread(mat_counts_fp)
  sample_totals <- fread(sample_totals_fp, data.table = FALSE)
  arg_comp      <- fread(arg_comp_fp)
  meta_df       <- fread(meta_fp)

  mat_counts_df <- as.data.frame(mat_counts)
  rownames(mat_counts_df) <- mat_counts_df[[1]]
  mat_counts_df[[1]] <- NULL
  mat_counts_df <- as.matrix(mat_counts_df)

  arg_comp_df <- as.data.frame(arg_comp)
  rownames(arg_comp_df) <- arg_comp_df[[1]]
  arg_comp_df[[1]] <- NULL
  arg_comp_df <- as.matrix(arg_comp_df)

  meta_df <- as.data.frame(meta_df)
  meta_df$sample <- as.character(meta_df$sample)

  site_col <- intersect(c("Setting", "setting", "Sample_Type", "sample_type", "Type", "type", "Site", "site"), colnames(meta_df))[1]
  city_col <- intersect(c("City", "city", "City_Code", "city_code"), colnames(meta_df))[1]

  if (is.na(site_col)) stop("Metadata lacks a sample-type column.")
  if (is.na(city_col)) stop("Metadata lacks a city column.")

  cat("Detected site column:", site_col, "\n")
  cat("Detected city column:", city_col, "\n")

  city_palette <- c(
    "BJ" = "#4E79A7",  "CQ" = "#E06C75",  "GY" = "#6CBF84",
    "GZ" = "#6EC5C1",  "HF" = "#F3A64A",  "HK" = "#E9C95B",
    "HR" = "#B38BCB",  "NJ" = "#C08F5A",  "SZ" = "#8FD0C7",
    "XA" = "#F6A6B2",  "XM" = "#9AD97C",  "ZB" = "#D8B4D8"
  )

  site_palette <- c(
    "Community"  = "#6EC5C1",
    "Hospital"   = "#E06C75",
    "Wet market" = "#E9C95B",
    "WWTP"       = "#4E79A7"
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

  meta_df[[site_col]] <- trimws(as.character(meta_df[[site_col]]))
  meta_df[[city_col]] <- trimws(as.character(meta_df[[city_col]]))
  meta_df$sample      <- trimws(as.character(meta_df$sample))

  meta_df[[site_col]][meta_df[[site_col]] %in% c("\u519c\u96c6\u8d38\u5e02\u573a", "WetMarket", "Wet_market")] <- "Wet market"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("\u793e\u533a")] <- "Community"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("\u533b\u9662")] <- "Hospital"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("\u6c61\u6c34\u5904\u7406\u5382", "Sewage treatment plant")] <- "WWTP"

  city_map <- c(
    "BJ" = "BJ", "Beijing" = "BJ",
    "CQ" = "CQ", "Chongqing" = "CQ",
    "GY" = "GY", "Guiyang" = "GY",
    "GZ" = "GZ", "Guangzhou" = "GZ",
    "HF" = "HF", "Hefei" = "HF",
    "HK" = "HK", "Haikou" = "HK",
    "HR" = "HR", "Harbin" = "HR", "Haerbin" = "HR",
    "NJ" = "NJ", "Nanjing" = "NJ",
    "SZ" = "SZ", "Shenzhen" = "SZ",
    "XA" = "XA", "Xian" = "XA", "Xi'an" = "XA",
    "XM" = "XM", "Xiamen" = "XM",
    "ZB" = "ZB", "Zibo" = "ZB"
  )
  city_map <- c(city_map, stats::setNames(
    c("BJ", "CQ", "GY", "GZ", "HF", "HK", "HR", "NJ", "SZ", "XA", "XM", "ZB"),
    c("\u5317\u4eac", "\u91cd\u5e86", "\u8d35\u9633", "\u5e7f\u5dde", "\u5408\u80a5", "\u6d77\u53e3",
      "\u54c8\u5c14\u6ee8", "\u5357\u4eac", "\u6df1\u5733", "\u897f\u5b89", "\u53a6\u95e8", "\u6dc4\u535a")))
  meta_df[[city_col]] <- ifelse(meta_df[[city_col]] %in% names(city_map),
                                unname(city_map[meta_df[[city_col]]]),
                                meta_df[[city_col]])

  site_order <- c("Hospital", "WWTP", "Community", "Wet market")
  city_order <- c("BJ", "CQ", "GY", "GZ", "HF", "HK", "HR", "NJ", "SZ", "XA", "XM", "ZB")

  species_order <- species_total %>%
    arrange(desc(total_count)) %>%
    pull(species)
  species_order <- intersect(species_order, rownames(mat_counts_df))
  mat_counts_df <- mat_counts_df[species_order, , drop = FALSE]
  arg_comp_df   <- arg_comp_df[species_order, , drop = FALSE]

  std_key <- function(x) {
    x <- tolower(trimws(x))
    x <- gsub("[ /-]+", "_", x)
    x
  }
  arg_type_map <- c(
    "aminoglycoside" = "Aminoglycoside",
    "antibacterial_fatty_acid" = "Antibacterial fatty acid",
    "bacitracin" = "Bacitracin",
    "bicyclomycin" = "Bicyclomycin",
    "bleomycin" = "Bleomycin",
    "chloramphenicol" = "Chloramphenicol",
    "defensin" = "Defensin",
    "edeine" = "Edeine",
    "factumycin" = "Factumycin",
    "florfenicol" = "Florfenicol",
    "fosfomycin" = "Fosfomycin",
    "fusidic_acid" = "Fusidic acid",
    "mls" = "MLS",
    "macrolide_lincosamide_streptogramin" = "MLS",
    "multidrug" = "Multidrug",
    "mupirocin" = "Mupirocin",
    "novobiocin" = "Novobiocin",
    "other_peptide_antibiotics" = "Other peptide antibiotics",
    "other_peptide_antibiotic" = "Other peptide antibiotics",
    "pleuromutilin_tiamulin" = "Pleuromutilin/Tiamulin",
    "polymyxin" = "Polymyxin",
    "puromycin" = "Puromycin",
    "quinolone" = "Quinolone",
    "rifamycin" = "Rifamycin",
    "streptothricin" = "Streptothricin",
    "sulfonamide" = "Sulfonamide",
    "tetracenomycin_c" = "Tetracenomycin C",
    "tetracycline" = "Tetracycline",
    "trimethoprim" = "Trimethoprim",
    "tunicamycin" = "Tunicamycin",
    "vancomycin" = "Vancomycin",
    "beta_lactam" = "beta_lactam",
    "betalactam" = "beta_lactam",
    "others" = "Others",
    "other" = "Others"
  )
  arg_keys <- std_key(colnames(arg_comp_df))
  arg_names_clean <- ifelse(arg_keys %in% names(arg_type_map),
                            arg_type_map[arg_keys],
                            colnames(arg_comp_df))
  arg_names_clean <- as.character(arg_names_clean)

  arg_comp_df2 <- do.call(cbind, lapply(unique(arg_names_clean), function(tp) {
    rowSums(arg_comp_df[, arg_names_clean == tp, drop = FALSE], na.rm = TRUE)
  }))
  arg_comp_df2 <- as.matrix(arg_comp_df2)
  rownames(arg_comp_df2) <- rownames(arg_comp_df)
  colnames(arg_comp_df2) <- unique(arg_names_clean)

  arg_cols_now <- colnames(arg_comp_df2)
  unknown_arg_cols <- setdiff(arg_cols_now, names(my_cols_cns))
  if (length(unknown_arg_cols) > 0) {
    other_sum <- rowSums(arg_comp_df2[, unknown_arg_cols, drop = FALSE], na.rm = TRUE)
    arg_comp_df2 <- arg_comp_df2[, setdiff(colnames(arg_comp_df2), unknown_arg_cols), drop = FALSE]
    if ("Others" %in% colnames(arg_comp_df2)) {
      arg_comp_df2[, "Others"] <- arg_comp_df2[, "Others"] + other_sum
    } else {
      arg_comp_df2 <- cbind(arg_comp_df2, Others = other_sum)
    }
  }
  arg_type_order <- colSums(arg_comp_df2, na.rm = TRUE) %>% sort(decreasing = TRUE) %>% names()
  arg_comp_df <- arg_comp_df2[, arg_type_order, drop = FALSE]
  arg_type_cols <- my_cols_cns[colnames(arg_comp_df)]

  sample_order_df <- meta_df %>%
    filter(sample %in% colnames(mat_counts_df)) %>%
    filter(!is.na(.data[[site_col]]), .data[[site_col]] != "",
           !is.na(.data[[city_col]]), .data[[city_col]] != "") %>%
    filter(.data[[site_col]] %in% site_order,
           .data[[city_col]] %in% city_order) %>%
    mutate(city_factor = factor(.data[[city_col]], levels = city_order),
           site_factor = factor(.data[[site_col]], levels = site_order)) %>%
    arrange(site_factor, city_factor, sample)

  sample_order <- sample_order_df$sample
  sample_order <- intersect(sample_order, colnames(mat_counts_df))
  mat_counts_df <- mat_counts_df[, sample_order, drop = FALSE]

  meta_plot <- sample_order_df %>%
    filter(sample %in% sample_order) %>%
    distinct(sample, .keep_all = TRUE)
  meta_plot <- meta_plot[match(sample_order, meta_plot$sample), , drop = FALSE]

  cat("Samples kept after filtering:", nrow(meta_plot), "\n")

  meta_plot$group <- paste(meta_plot[[city_col]], meta_plot[[site_col]], sep = " | ")
  group_order_df <- meta_plot %>%
    distinct(group, .keep_all = TRUE) %>%
    mutate(city_factor = factor(.data[[city_col]], levels = city_order),
           site_factor = factor(.data[[site_col]], levels = site_order)) %>%
    arrange(site_factor, city_factor)
  group_order <- group_order_df$group

  pooled <- .fig3a_pool_counts(mat_counts_df, sample_totals, meta_plot, group_order)

  meta_group <- group_order_df %>%
    dplyr::select(group, all_of(city_col), all_of(site_col))
  meta_group[[site_col]] <- factor(meta_group[[site_col]], levels = site_order)
  meta_group[[city_col]] <- factor(meta_group[[city_col]], levels = city_order)

  mat_rel_df <- pooled$proportions

  left_counts <- species_total %>%
    filter(species %in% species_order) %>%
    distinct(species, total_count) %>%
    as.data.frame()
  left_counts <- left_counts[match(species_order, left_counts$species), , drop = FALSE]
  left_bar <- left_counts$total_count
  names(left_bar) <- left_counts$species
  left_bar_rev <- -left_bar
  max_count <- max(left_bar, na.rm = TRUE)

  count_breaks <- c(0, 30000, 60000)
  count_breaks <- count_breaks[count_breaks <= max_count]

  highlight_species <- c("Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
                         "Enterococcus faecium", "Acinetobacter baumannii", "Staphylococcus aureus",
                         "Klebsiella quasipneumoniae", "Enterobacter cloacae")

  row_label_cols <- ifelse(species_order %in% highlight_species, "#D62728", "black")

  left_anno <- rowAnnotation(
    "Number of counts" = anno_barplot(
      left_bar_rev, baseline = 0,
      gp = gpar(fill = "#B7D7E8", col = "#9EC5DA"),
      border = FALSE, bar_width = 0.82,
      width = unit(5, "cm"),
      add_numbers = FALSE,
      axis_param = list(side = "top", at = -count_breaks, labels = count_breaks,
                        labels_rot = 0, gp = gpar(fontsize = 7))
    ),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 10, font = 1)
  )

  right_anno <- rowAnnotation(
    `ARG types` = anno_barplot(
      arg_comp_df, beside = FALSE,
      gp = gpar(fill = arg_type_cols[colnames(arg_comp_df)], col = NA),
      border = FALSE, bar_width = 1, width = unit(2.5, "cm"),
      axis = FALSE
    ),
    show_annotation_name = F
  )

  top_ha <- HeatmapAnnotation(
    `Sample type` = meta_group[[site_col]],
    City = meta_group[[city_col]],
    col = list(`Sample type` = site_palette, City = city_palette),
    annotation_legend_param = list(title_gp = gpar(fontsize = 10, font=1)),
    annotation_name_gp = gpar(fontsize = 10, font = 1),
    simple_anno_size = unit(4.2, "mm"), gap = unit(1.2, "mm"),
    show_annotation_name = TRUE
  )

  bubble_col_fun <- colorRamp2(c(0, 0.01, 0.02, 0.05),
                               c("#F7F7F7", "#F6D6B8", "#E98B4A", "#B30000"))

  cell_fun_bubble <- function(j, i, x, y, width, height, fill) {
    v <- mat_rel_df[i, j]
    if (is.na(v) || v == 0) return()
    max_radius <- min(unit.c(width, height)) * 1.0
    r <- sqrt(v / 0.10) * max_radius
    grid.circle(x = x, y = y, r = r,
                gp = gpar(fill = bubble_col_fun(v), col = NA, font = 1))
  }

  column_split_site <- factor(meta_group[[site_col]], levels = site_order)

  size_vals <- c(0.01, 0.02, 0.05)

  ref_width <- unit(1, "cm")
  ref_height <- unit(1, "cm")
  max_radius_ref <- min(ref_width, ref_height) * 0.5
  radii <- sqrt(size_vals / 0.10) * max_radius_ref

  bubble_col_fun <- colorRamp2(c(0, 0.01, 0.02, 0.05),
                               c("#F7F7F7", "#F6D6B8", "#E98B4A", "#B30000"))

  point_colors <- bubble_col_fun(size_vals)

  combined_legend <- Legend(
    title = "Species-assigned\ncontig fraction",
    type = "points",
    at = size_vals,
    labels = c("1%", "2%", "5%"),
    legend_gp = gpar(col = point_colors, fill=NULL),
    background = "transparent",
    size = unit(radii, "cm"),
    title_gp = gpar(fontsize = 10, font = 1),
    labels_gp = gpar(fontsize = 8, font = 1)
  )
  lgd_arg <- Legend(title = "ARG types", at = names(arg_type_cols),
                    labels = gsub("beta_lactam","Beta_lactam",gsub("Other peptide antibiotics","OPA",names(arg_type_cols))),
                    grid_height = unit(4, "mm"),tick_length = unit(0.4, "mm"),
                    legend_gp = gpar(fill = arg_type_cols),
                    ncol = 1, title_gp = gpar(fontsize = 10, font = 1)
  )
  all_legends <- list(combined_legend, lgd_arg)

  ht_bubble <- Heatmap(
    mat_rel_df,
    name = "Species-assigned\ncontig fraction",
    col = bubble_col_fun,
    rect_gp = gpar(type = "none"),
    cell_fun = cell_fun_bubble,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    heatmap_legend_param= list(title_gp = gpar(font=1),labels_gp = gpar(font=1)),
    column_split = column_split_site,
    column_gap = unit(1.3, "mm"),
    row_names_side = "left",
    show_row_names = FALSE,

    show_column_names = FALSE,
    column_labels = sub(" \\| .*", "", colnames(mat_rel_df)),
    column_names_rot = 0,
    column_names_gp = gpar(fontsize = 8),
    column_names_centered = FALSE,
    top_annotation = top_ha,
    left_annotation = left_anno,
    row_title = NULL,
    column_title = NULL,
    border = TRUE,
    use_raster = FALSE,
    show_heatmap_legend = FALSE
  ) + right_anno

  draw_bubble <- function(ht_obj, pdf_fp, png_fp, width = 11, height = 5.2) {

    add_row_names_to_left_bar <- function() {
      decorate_annotation("Number of counts", {
        n_rows <- length(species_order)

        x_pos_fixed <- unit(0.02, "npc")
        for (i in seq_along(species_order)) {
          species_name <- species_order[i]

          y_pos <- unit(n_rows - i + 1, "native")
          text_col <- ifelse(species_name %in% highlight_species, "#D62728", "black")
          grid.text(species_name, x = x_pos_fixed, y = y_pos,
                    just = "left",
                    gp = gpar(fontsize = 10, fontface = "italic", col = text_col))
        }
      })
    }

    pdf(pdf_fp, width = width, height = height)
    draw(ht_obj,
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         merge_legends = FALSE,
         annotation_legend_list = all_legends)
    add_row_names_to_left_bar()
    dev.off()

    png(png_fp, width = width, height = height, units = "in", res = 300)
    draw(ht_obj,
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         merge_legends = FALSE,
         annotation_legend_list = all_legends)
    add_row_names_to_left_bar()
    dev.off()
  }

  if (render) draw_bubble(ht_bubble,
              file.path(out_dir, "species_host_bubbleplot.pdf"),
              file.path(out_dir, "species_host_bubbleplot.png"))

  cat("Figure 3a completed.\n")
  cat("Output dir:", out_dir, "\n")
  cat("Generated files:\n")
  print(list.files(out_dir, full.names = TRUE))

  data.table::fwrite(meta_group, file.path(out_dir, "city_setting_groups.tsv"), sep = "\t")
  data.table::fwrite(pooled$denominators, file.path(out_dir, "city_setting_denominators.tsv"), sep = "\t")
  data.table::fwrite(data.frame(species=rownames(pooled$counts), pooled$counts, check.names=FALSE),
                    file.path(out_dir, "city_setting_species_counts.tsv"), sep = "\t")
  data.table::fwrite(data.frame(species = rownames(mat_rel_df), mat_rel_df, check.names = FALSE),
                    file.path(out_dir, "city_setting_relative_abundance.tsv"), sep = "\t")
  invisible(list(group_matrix = mat_rel_df, group_metadata = meta_group,
                 total_counts = left_bar, arg_composition = arg_comp_df,
                 group_counts = pooled$counts, group_denominators = pooled$denominators,
                 sample_count = nrow(meta_plot), heatmap = ht_bubble))

}

.fig3_annotation_count_columns <- function(data, aliases) {
  for (canonical in names(aliases)) {
    legacy <- unname(aliases[[canonical]])
    if (!canonical %in% names(data)) {
      if (!legacy %in% names(data)) stop("Missing ARG annotation-count field: ", canonical, call. = FALSE)
      names(data)[names(data) == legacy] <- canonical
    } else if (legacy %in% names(data)) {
      equal <- isTRUE(all.equal(suppressWarnings(as.numeric(data[[canonical]])),
                               suppressWarnings(as.numeric(data[[legacy]])),
                               tolerance = 0, check.attributes = FALSE))
      if (!equal) stop("Conflicting ARG annotation-count fields: ", canonical, " and ", legacy, call. = FALSE)
      data[[legacy]] <- NULL
    }
  }
  data
}

FIG3B_REFERENCE_TYPE_ORDER <- c("Multidrug", "MLS", "beta_lactam", "Bacitracin",
  "Aminoglycoside", "Tetracycline", "Polymyxin", "Others", "Mupirocin",
  "Chloramphenicol", "Trimethoprim")

.fig3b_palettes <- function() {
  list(type = c(
    "Aminoglycoside" = "#4E79A7", "Antibacterial fatty acid" = "#7FB8B5",
    "Bacitracin" = "#D6B97A", "Bicyclomycin" = "#8BCB88", "Bleomycin" = "#E39B76",
    "Chloramphenicol" = "#B5658D", "Defensin" = "#B7D989", "Edeine" = "#8FB6E8",
    "Factumycin" = "#D17C98", "Florfenicol" = "#AFC7E8", "Fosfomycin" = "#B7E0E5",
    "Fusidic acid" = "#E6BAC4", "MLS" = "#74C2B3", "Multidrug" = "#7E6BB5",
    "Mupirocin" = "#5F9EA0", "Novobiocin" = "#72B77E",
    "Other peptide antibiotics" = "#9FD58A", "Pleuromutilin/Tiamulin" = "#E4CF62",
    "Polymyxin" = "#75C0C1", "Puromycin" = "#9DBDE0", "Quinolone" = "#9A84D6",
    "Rifamycin" = "#A9C98A", "Streptothricin" = "#6DBDE3", "Sulfonamide" = "#C8A24B",
    "Tetracenomycin C" = "#D8BB92", "Tetracycline" = "#E8D98F",
    "Trimethoprim" = "#E39A9A", "Tunicamycin" = "#9AA3A8", "Vancomycin" = "#AFC9A0",
    "beta_lactam" = "#7DA34D", "Others" = "#D0D0D0"),
    context = c("Chromosome" = "#7DA34D", "Plasmid" = "#9A84D6", "Phage" = "#E39B76"),
    setting = c("Hospital" = "#E06C75", "WWTP" = "#4E79A7",
                "Community" = "#6EC5C1", "Wet market" = "#E9C95B"))
}

.fig3b_normalize <- function(values, aliases, field, allow_unknown = FALSE) {
  values <- trimws(as.character(values))
  if (anyNA(values) || any(!nzchar(values))) stop(field, " must not contain missing or blank values.", call. = FALSE)
  keys <- tolower(gsub("[[:space:]]+", " ", values))
  mapped <- unname(aliases[keys])
  unknown <- is.na(mapped)
  if (any(unknown) && !allow_unknown) {
    stop("Unsupported ", field, " value(s): ", paste(unique(values[unknown]), collapse = ", "), call. = FALSE)
  }
  mapped[unknown] <- values[unknown]
  mapped
}

prepare_fig3b_orf <- function(input, top_n = 10, type_order = NULL) {
  if (length(top_n) != 1L || !is.numeric(top_n) || !is.finite(top_n) || top_n < 1 || top_n != floor(top_n)) {
    stop("top_n must be a positive integer.", call. = FALSE)
  }
  if (is.character(input) && length(input) == 1L) {
    if (!file.exists(input)) stop("Input TSV does not exist: ", input, call. = FALSE)
    input <- read.delim(input, check.names = FALSE, stringsAsFactors = FALSE,
                        quote = "", comment.char = "", fileEncoding = "UTF-8-BOM")
  }
  required <- c("arg_type", "Context", "Setting", "n_ARG_ORFs")
  if (!is.data.frame(input) || !all(required %in% names(input))) {
    stop("Input must contain arg_type, Context, Setting, n_ARG_ORFs. Contig-count input is not accepted.", call. = FALSE)
  }
  if (anyDuplicated(names(input))) stop("Input column names must be unique.", call. = FALSE)
  input <- input[, required, drop = FALSE]
  counts <- input$n_ARG_ORFs
  if (!is.numeric(counts) || anyNA(counts) || any(!is.finite(counts)) || any(counts < 0) || any(counts != floor(counts))) {
    stop("n_ARG_ORFs must contain finite, non-negative numeric integers.", call. = FALSE)
  }
  total <- sum(counts)
  if (!is.finite(total) || total > 2^53 - 1) stop("ORF total exceeds the exact integer range supported by the plot.", call. = FALSE)
  if (total <= 0) stop("At least one positive ORF count is required; total is zero.", call. = FALSE)
  palettes <- .fig3b_palettes()
  type_aliases <- setNames(names(palettes$type), tolower(names(palettes$type)))
  type_aliases <- c(type_aliases, "macrolide-lincosamide-streptogramin" = "MLS",
                    "macrolide lincosamide streptogramin" = "MLS",
                    "beta-lactam" = "beta_lactam", "beta lactam" = "beta_lactam",
                    "beta-lactams" = "beta_lactam", "beta_lactams" = "beta_lactam")
  context_aliases <- c("chromosome" = "Chromosome", "chromosomal" = "Chromosome", "chr" = "Chromosome",
                       "plasmid" = "Plasmid", "plasmids" = "Plasmid",
                       "phage" = "Phage", "phages" = "Phage", "virus" = "Phage", "viral" = "Phage")
  setting_aliases <- c("hospital" = "Hospital", "wwtp" = "WWTP", "community" = "Community",
                       "wet market" = "Wet market", "wet_market" = "Wet market", "wetmarket" = "Wet market")
  input$arg_type <- .fig3b_normalize(input$arg_type, type_aliases, "arg_type", TRUE)
  input$Context <- .fig3b_normalize(input$Context, context_aliases, "Context")
  input$Setting <- .fig3b_normalize(input$Setting, setting_aliases, "Setting")
  input <- input[input$n_ARG_ORFs > 0, , drop = FALSE]
  type_totals <- aggregate(n_ARG_ORFs ~ arg_type, input, sum)
  type_totals <- type_totals[order(-type_totals$n_ARG_ORFs, type_totals$arg_type), , drop = FALSE]
  top_types <- head(type_totals$arg_type[type_totals$arg_type != "Others"], top_n)
  input$type_plot <- ifelse(input$arg_type %in% top_types, input$arg_type, "Others")
  paths <- aggregate(n_ARG_ORFs ~ type_plot + Context + Setting, input, sum)
  displayed_types <- c(top_types, if ("Others" %in% paths$type_plot) "Others")
  if (!is.null(type_order)) {
    if (!is.character(type_order) || !length(type_order)) stop("type_order must be a non-empty character vector or NULL.", call. = FALSE)
    type_order <- .fig3b_normalize(type_order, type_aliases, "type_order", TRUE)
    if (anyDuplicated(type_order)) stop("type_order contains duplicate canonical type names.", call. = FALSE)
    displayed_types <- c(intersect(type_order, displayed_types), setdiff(displayed_types, type_order))
  }
  type_order <- displayed_types
  context_order <- names(palettes$context)[names(palettes$context) %in% paths$Context]
  setting_order <- names(palettes$setting)[names(palettes$setting) %in% paths$Setting]
  paths <- paths[order(match(paths$type_plot, type_order), match(paths$Context, context_order),
                      match(paths$Setting, setting_order)), , drop = FALSE]
  rownames(paths) <- NULL
  node_names <- list(type_order, context_order, setting_order)
  fields <- c("type_plot", "Context", "Setting")
  prefixes <- c("type", "context", "setting")
  nodes <- do.call(rbind, lapply(seq_along(node_names), function(i) {
    nm <- node_names[[i]]
    totals <- vapply(nm, function(n) sum(paths$n_ARG_ORFs[paths[[fields[i]]] == n]), numeric(1))
    colors <- unname(palettes[[i]][nm]); colors[is.na(colors)] <- "#BDBDBD"
    data.frame(key = paste(prefixes[i], nm, sep = "::"), name = nm, layer = i - 1L,
               order = seq_along(nm) - 1L, total = unname(totals), pct = unname(totals) / total * 100,
               color = colors, stringsAsFactors = FALSE)
  }))
  rownames(nodes) <- NULL
  nodes$id <- seq_len(nrow(nodes)) - 1L
  nodes$group <- paste0(prefixes[nodes$layer + 1L], "__", nodes$id)
  label_names <- ifelse(nodes$name == "Wet market" & nodes$layer == 2L, "Wet Market", nodes$name)
  nodes$label <- sprintf("%s (%.1f%%)", label_names, nodes$pct)
  nodes <- nodes[, c("id", "key", "name", "layer", "order", "total", "pct", "color", "group", "label")]
  edge1 <- aggregate(n_ARG_ORFs ~ type_plot + Context, paths, sum)
  edge2 <- aggregate(n_ARG_ORFs ~ Context + Setting, paths, sum)
  links <- rbind(
    data.frame(source_key = paste0("type::", edge1$type_plot), target_key = paste0("context::", edge1$Context),
               value = edge1$n_ARG_ORFs, step = "ARG type -> Context", stringsAsFactors = FALSE),
    data.frame(source_key = paste0("context::", edge2$Context), target_key = paste0("setting::", edge2$Setting),
               value = edge2$n_ARG_ORFs, step = "Context -> Setting", stringsAsFactors = FALSE))
  links$source <- match(links$source_key, nodes$key) - 1L
  links$target <- match(links$target_key, nodes$key) - 1L
  links$group <- nodes$group[links$source + 1L]
  links$color <- nodes$color[links$source + 1L]
  links <- links[order(links$source, links$target), , drop = FALSE]
  rownames(links) <- NULL
  links <- links[, c("source", "target", "value", "group", "step", "source_key", "target_key", "color")]
  percentages <- expand.grid(Context = names(palettes$context), Setting = setting_order, stringsAsFactors = FALSE)
  percentages <- merge(percentages, edge2, by = c("Context", "Setting"), all.x = TRUE, sort = FALSE)
  percentages$n_ARG_ORFs[is.na(percentages$n_ARG_ORFs)] <- 0
  setting_totals <- setNames(nodes$total[nodes$layer == 2L], nodes$name[nodes$layer == 2L])
  percentages$total_setting <- unname(setting_totals[percentages$Setting])
  percentages$pct <- percentages$n_ARG_ORFs / percentages$total_setting * 100
  percentages <- percentages[order(match(percentages$Setting, setting_order),
                                    match(percentages$Context, names(palettes$context))), , drop = FALSE]
  rownames(percentages) <- NULL
  stopifnot(sum(paths$n_ARG_ORFs) == total,
            all(vapply(split(links$value, links$step), sum, numeric(1)) == total))
  list(nodes = nodes, links = links, path_counts = paths, setting_percentages = percentages,
       total = total, top_types = top_types)
}

.fig3b_xml <- function(text) {
  text <- gsub("&", "&amp;", as.character(text), fixed = TRUE)
  text <- gsub("<", "&lt;", text, fixed = TRUE)
  text <- gsub(">", "&gt;", text, fixed = TRUE)
  text <- gsub('"', "&quot;", text, fixed = TRUE)
  gsub("'", "&apos;", text, fixed = TRUE)
}

.fig3b_geometry <- function(prepared, width, height) {
  for (value in list(width, height)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value <= 0) {
      stop("width and height must be positive numbers.", call. = FALSE)
    }
  }
  if (width < 800 || height < 550) stop("Use width >= 800 and height >= 550 for the labels.", call. = FALSE)
  nodes <- prepared$nodes; links <- prepared$links
  margin <- c(top = 8, right = 78, bottom = 26, left = 32)
  inner_width <- width - margin["left"] - margin["right"]
  inner_height <- height - margin["top"] - margin["bottom"]
  gaps <- c(20, 44, 40)
  column_offsets <- c(0, 28, 28)
  column_counts <- as.numeric(table(factor(nodes$layer, levels = 0:2)))
  flow_height <- unname(min(inner_height - (column_counts - 1) * gaps - column_offsets))
  if (flow_height <= 0) stop("Increase height to accommodate the selected number of ARG types.", call. = FALSE)
  nodes$dx <- 36
  nodes$dy <- nodes$total / prepared$total * flow_height
  nodes$x <- NA_real_; nodes$y <- NA_real_
  x_positions <- unname(margin["left"] + c(0, (inner_width - 36) * 0.51, inner_width - 36))
  for (layer in 0:2) {
    idx <- which(nodes$layer == layer)
    idx <- idx[order(nodes$order[idx])]
    gap <- gaps[layer + 1L]
    y <- unname(margin["top"] + column_offsets[layer + 1L])
    for (i in idx) {
      nodes$x[i] <- x_positions[layer + 1L]; nodes$y[i] <- y
      y <- y + nodes$dy[i] + gap
    }
  }
  links$dy <- links$value / prepared$total * flow_height
  links$sy <- 0; links$ty <- 0
  for (id in nodes$id) {
    outgoing <- which(links$source == id)
    outgoing <- outgoing[order(nodes$y[links$target[outgoing] + 1L])]
    incoming <- which(links$target == id)
    incoming <- incoming[order(nodes$y[links$source[incoming] + 1L])]
    if (length(outgoing)) links$sy[outgoing] <- c(0, head(cumsum(links$dy[outgoing]), -1L))
    if (length(incoming)) links$ty[incoming] <- c(0, head(cumsum(links$dy[incoming]), -1L))
  }
  links$x0 <- nodes$x[links$source + 1L] + nodes$dx[links$source + 1L]
  links$x1 <- nodes$x[links$target + 1L]
  links$y0 <- nodes$y[links$source + 1L] + links$sy + links$dy / 2
  links$y1 <- nodes$y[links$target + 1L] + links$ty + links$dy / 2
  marks <- links[nodes$layer[links$target + 1L] == 2L, , drop = FALSE]
  marks$total_setting <- nodes$total[marks$target + 1L]
  marks$pct <- marks$value / marks$total_setting * 100
  marks$actual <- marks$y1
  marks <- marks[order(marks$actual, marks$source), , drop = FALSE]
  marks$y <- marks$actual
  for (target_id in unique(marks$target)) {
    idx <- which(marks$target == target_id)
    target <- nodes[target_id + 1L, ]
    inset <- min(2, target$dy / 2)
    lower <- target$y + inset
    upper <- target$y + target$dy - inset
    gap <- if (length(idx) > 1L) min(22, (upper - lower) / (length(idx) - 1L)) else 0
    y <- pmax(lower, pmin(upper, marks$actual[idx]))
    if (length(idx) > 1L) {
      for (i in 2:length(idx)) y[i] <- max(y[i], y[i - 1L] + gap)
      y[length(idx)] <- min(y[length(idx)], upper)
      for (i in seq.int(length(idx) - 1L, 1L)) y[i] <- min(y[i], y[i + 1L] - gap)
      y[1L] <- max(y[1L], lower)
      for (i in 2:length(idx)) y[i] <- max(y[i], y[i - 1L] + gap)
    }
    marks$y[idx] <- y
  }
  marks$text <- ifelse(marks$pct > 0 & marks$pct < 0.1, "<0.1%", sprintf("%.1f%%", marks$pct))
  list(nodes = nodes, links = links, labels = marks, flow_height = flow_height,
       width = width, height = height, node_width = 36, font_size = 20, pct_font_size = 17,
       margin = margin)
}

write_fig3b_svg <- function(prepared, file, width = 1000, height = 650, background = "transparent") {
  if (!is.list(prepared) || !all(c("nodes", "links", "total") %in% names(prepared))) {
    stop("prepared must be the result of prepare_fig3b_orf().", call. = FALSE)
  }
  if (length(file) != 1L || !is.character(file) || !nzchar(file)) stop("file must be an SVG output path.", call. = FALSE)
  if (file.exists(file)) stop("SVG output already exists; refusing to overwrite: ", file, call. = FALSE)
  if (!dir.exists(dirname(file))) stop("SVG parent directory does not exist: ", dirname(file), call. = FALSE)
  if (length(background) != 1L || !is.character(background) || is.na(background)) stop("background must be a CSS color.", call. = FALSE)
  geometry <- .fig3b_geometry(prepared, width, height)
  nodes <- geometry$nodes; links <- geometry$links; marks <- geometry$labels
  f <- function(x) formatC(x, format = "g", digits = 17, decimal.mark = ".")
  escape <- .fig3b_xml
  svg <- c('<?xml version="1.0" encoding="UTF-8"?>',
    sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s" role="img" data-total-ARG-ORFs="%s">', f(width), f(height), f(width), f(height), f(prepared$total)),
    '<title>ARG ORF type, genomic context, and sampling setting</title>',
    '<desc>Bar labels use all ARG ORFs as denominator. Blue percentages use ARG ORFs within the receiving setting. All positive flows are retained.</desc>',
    sprintf('<rect width="%s" height="%s" fill="%s"/>', f(width), f(height), escape(background)),
    '<g id="flows" fill="none">')
  for (i in order(-links$dy, links$source, links$target)) {
    d <- links[i, ]; midpoint <- (d$x0 + d$x1) / 2
    path <- sprintf("M%s,%s C%s,%s %s,%s %s,%s", f(d$x0), f(d$y0), f(midpoint), f(d$y0), f(midpoint), f(d$y1), f(d$x1), f(d$y1))
    title <- paste0(nodes$name[d$source + 1L], " -> ", nodes$name[d$target + 1L], ": ", format(d$value, scientific = FALSE, trim = TRUE), " ARG ORFs")
    svg <- c(svg, sprintf('<path class="fig3b-link" data-source="%d" data-target="%d" data-value="%s" d="%s" stroke="%s" stroke-width="%s" stroke-opacity="0.42"><title>%s</title></path>',
                         d$source, d$target, f(d$value), path, escape(d$color), f(d$dy), escape(title)))
  }
  svg <- c(svg, '</g>', '<g id="bars">')
  for (i in seq_len(nrow(nodes))) {
    d <- nodes[i, ]
    svg <- c(svg, sprintf('<rect class="fig3b-node" data-id="%d" data-layer="%d" data-order="%d" data-value="%s" x="%s" y="%s" width="36" height="%s" fill="%s" fill-opacity="0.9" stroke="#333333" stroke-width="0.7"><title>%s</title></rect>',
      d$id, d$layer, d$order, f(d$total), f(d$x), f(d$y), f(d$dy), escape(d$color), escape(d$label)))
  }
  svg <- c(svg, '</g>', '<g id="node-labels" font-family="Arial" font-size="20" fill="#000000">')
  for (i in seq_len(nrow(nodes))) {
    d <- nodes[i, ]; x <- d$x + if (d$layer == 0L) 44 else 18
    y <- d$y + if (d$layer == 0L) d$dy / 2 else d$dy + 22
    svg <- c(svg, sprintf('<text class="fig3b-node-label" x="%s" y="%s" dy="0.35em" text-anchor="%s">%s</text>',
      f(x), f(y), if (d$layer == 0L) "start" else "middle", escape(d$label)))
  }
  svg <- c(svg, '</g>', '<g id="setting-percentages" font-family="Arial" font-size="17" fill="#2864B4">')
  for (i in seq_len(nrow(marks))) {
    m <- marks[i, ]; x <- m$x1
    title <- sprintf("%s within %s: %.6f%%", nodes$name[m$source + 1L], nodes$name[m$target + 1L], m$pct)
    svg <- c(svg, sprintf('<text class="fig3b-setting-pct" data-source="%d" data-target="%d" data-denominator="%s" data-pct="%s" x="%s" y="%s" dy="0.35em" text-anchor="end">%s<title>%s</title></text>',
      m$source, m$target, f(m$total_setting), f(m$pct), f(x - 17), f(m$y), escape(m$text), escape(title)))
  }
  svg <- c(svg, '</g>', '</svg>')
  writeLines(enc2utf8(svg), file, useBytes = TRUE)
  invisible(geometry)
}

.fig3b_layout_js <- function() {
  'function(el, x, data) {
    var svg = d3.select(el).select("svg");
    var nodeSel = svg.selectAll(".node"), linkSel = svg.selectAll(".link");
    var nodes = [], links = [];
    nodeSel.each(function(d, i) {
      var m = data.nodes[i];
      d.fig3b = m; d.dx = m.dx; d.dy = m.dy;
      d.x = m.x - x.options.margin.left;
      d.y = m.y - x.options.margin.top;
      nodes.push(d);
    });
    linkSel.each(function(d) { d.dy = d.value / data.total * data.flow_height; links.push(d); });
    if (!nodes.length || !links.length) throw new Error("Fig3b: no rendered Sankey node/link data.");
    var graph = d3.select(nodeSel.node().parentNode);
    var margin = x.options.margin;
    var canvasWidth = data.width, canvasHeight = data.height;
    var innerWidth = canvasWidth - margin.left - margin.right;
    var innerHeight = canvasHeight - margin.top - margin.bottom;
    svg.attr("viewBox", "0 0 " + canvasWidth + " " + canvasHeight);
    if (typeof d3.sankey === "function") {
      d3.sankey().nodes(nodes).links(links).nodeWidth(36).relayout();
    } else {
      nodes.forEach(function(n) {
        var outgoing = links.filter(function(l) { return l.source === n; }).sort(function(a,b) { return a.target.y-b.target.y; });
        var incoming = links.filter(function(l) { return l.target === n; }).sort(function(a,b) { return a.source.y-b.source.y; });
        var sy = 0, ty = 0;
        outgoing.forEach(function(l) { l.sy = sy; sy += l.dy; });
        incoming.forEach(function(l) { l.ty = ty; ty += l.dy; });
      });
    }
    function path(d) {
      var x0 = d.source.x + d.source.dx, x1 = d.target.x;
      var x2 = (x0 + x1) / 2;
      var y0 = d.source.y + d.sy + d.dy / 2, y1 = d.target.y + d.ty + d.dy / 2;
      return "M" + x0 + "," + y0 + "C" + x2 + "," + y0 + " " + x2 + "," + y1 + " " + x1 + "," + y1;
    }
    svg.style("background", data.background);
    linkSel.attr("d", path).style("stroke-width", function(d) { return d.dy; }).style("stroke-opacity", 0.42)
      .on("mouseover", null).on("mouseout", null);
    nodeSel.attr("transform", function(d) { return "translate(" + d.x + "," + d.y + ")"; })
      .on(".drag", null).on("mouseover", null).on("mouseout", null).style("cursor", "default");
    nodeSel.select("rect").attr("width", 36).attr("height", function(d) { return d.dy; }).style("stroke", "#333").style("stroke-width", "0.7px")
      .style("cursor", "default");
    nodeSel.select("text").text(function(d) { return d.fig3b.label; })
      .attr("x", function(d) { return d.fig3b.layer === 0 ? 44 : 18; })
      .attr("y", function(d) { return d.fig3b.layer === 0 ? d.dy / 2 : d.dy + 22; })
      .attr("dy", "0.35em").attr("text-anchor", function(d) { return d.fig3b.layer === 0 ? "start" : "middle"; })
      .style("font-family", "Arial").style("font-size", "20px").style("fill", "#000");
    nodeSel.select("title").text(function(d) { return d.fig3b.name + "\\n" + d.fig3b.total.toLocaleString() + " ARG ORFs\\n" + d.fig3b.pct.toFixed(2) + "% of all ARG ORFs"; });
    linkSel.select("title").text(function(d) { return d.source.fig3b.name + " -> " + d.target.fig3b.name + "\\n" + d.value.toLocaleString() + " ARG ORFs"; });
    graph.selectAll(".fig3b-setting-pct").remove();
    var marks = links.filter(function(d) { return d.target.fig3b.layer === 2; }).map(function(d) {
      var saved = data.setting_labels.filter(function(m) {
        return m.source === d.source.fig3b.id && m.target === d.target.fig3b.id;
      })[0];
      return {link:d, actual:d.target.y + d.ty + d.dy / 2,
              pct:d.value / d.target.fig3b.total * 100, y:saved.y-margin.top};
    });
    marks.sort(function(a,b) { return a.actual-b.actual; });
    graph.selectAll(".fig3b-setting-pct").data(marks).enter().append("text")
      .attr("class", "fig3b-setting-pct").attr("x", function(m) { return m.link.target.x-17; })
      .attr("y", function(m) { return m.y; }).attr("dy", "0.35em").attr("text-anchor", "end")
      .text(function(m) { return m.pct > 0 && m.pct < 0.1 ? "<0.1%" : m.pct.toFixed(1)+"%"; })
      .style("font-family", "Arial").style("font-size", "17px").style("fill", "#2864B4")
      .append("title").text(function(m) { return m.link.source.fig3b.name+" within "+m.link.target.fig3b.name+": "+m.pct.toFixed(3)+"%"; });
    el.setAttribute("data-fig3b-layout", "fixed-three-column");
    el.setAttribute("data-fig3b-total", data.total);
  }'
}

plot_fig3b_orf <- function(input_file, out_dir, top_n = 10, width = 1000, height = 650,
                          export_static = TRUE, background = "transparent", type_order = NULL) {
  if (length(out_dir) != 1L || !is.character(out_dir) || !nzchar(out_dir)) stop("out_dir must be a directory path.", call. = FALSE)
  if (file.exists(out_dir) && !dir.exists(out_dir)) stop("out_dir already exists as a file.", call. = FALSE)
  if (dir.exists(out_dir) && length(list.files(out_dir, all.files = TRUE, no.. = TRUE))) {
    stop("Output directory is non-empty; choose a new directory to preserve existing results: ", out_dir, call. = FALSE)
  }
  p <- prepare_fig3b_orf(input_file, top_n, type_order = type_order)
  for (value in list(width, height)) {
    if (length(value) != 1L || !is.numeric(value) || !is.finite(value) || value <= 0) stop("width and height must be positive numbers.", call. = FALSE)
  }
  if (width < 800 || height < 550) stop("Use width >= 800 and height >= 550 to allow space for node and percentage labels.", call. = FALSE)
  if (length(export_static) != 1L || !is.logical(export_static) || is.na(export_static)) stop("export_static must be TRUE or FALSE.", call. = FALSE)
  if (length(background) != 1L || !is.character(background) || is.na(background)) stop("background must be a CSS color.", call. = FALSE)
  .fig3b_geometry(p, width, height)
  if (!dir.exists(out_dir) && !dir.create(out_dir, recursive = TRUE)) stop("Cannot create output directory: ", out_dir, call. = FALSE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = TRUE)
  svg <- file.path(out_dir, "Fig3b_ORF_original_style.svg")
  geometry <- write_fig3b_svg(p, svg, width = width, height = height, background = background)
  required_packages <- c("networkD3", "htmlwidgets", "jsonlite")
  missing <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
  widget <- NULL; html <- NULL
  if (length(missing)) {
    message("Optional HTML export skipped; unavailable packages: ", paste(missing, collapse = ", "), ". Native SVG and CSV tables are still generated. No packages were installed automatically.")
  } else {
  color_scale <- paste0("d3.scaleOrdinal().domain(", jsonlite::toJSON(p$nodes$group),
                        ").range(", jsonlite::toJSON(p$nodes$color), ")")
  widget <- networkD3::sankeyNetwork(Links = p$links, Nodes = p$nodes,
    Source = "source", Target = "target", Value = "value", NodeID = "label",
    NodeGroup = "group", LinkGroup = "group", units = "ARG ORFs", colourScale = color_scale,
    fontFamily = "Arial", fontSize = 20, nodeWidth = 36, nodePadding = 20,
    sinksRight = TRUE, iterations = 0, width = width, height = height,
    margin = as.list(geometry$margin))
  metadata <- list(nodes = lapply(seq_len(nrow(geometry$nodes)), function(i) as.list(geometry$nodes[i, ])),
                   total = p$total, width = width, height = height, background = background,
                   flow_height = geometry$flow_height,
                   setting_labels = lapply(seq_len(nrow(geometry$labels)), function(i) as.list(geometry$labels[i, ])))
  widget <- htmlwidgets::onRender(widget, .fig3b_layout_js(), data = metadata)
  html <- file.path(out_dir, "Fig3b_ORF_original_style.html")
  htmlwidgets::saveWidget(widget, html, selfcontained = FALSE,
                          libdir = "Fig3b_ORF_original_style_files", background = background,
                          title = "Fig. 3b | ARG ORF type, genomic context and setting")
  }
  tables <- c(nodes = "Fig3b_nodes.csv", links = "Fig3b_links.csv", path_counts = "Fig3b_path_counts.csv",
              setting_percentages = "Fig3b_setting_percentages.csv")
  for (name in names(tables)) write.csv(p[[name]], file.path(out_dir, tables[[name]]), row.names = FALSE, fileEncoding = "UTF-8")
  write.csv(data.frame(total_ARG_ORFs = p$total, displayed_ARG_types = sum(p$nodes$layer == 0L),
                       top_n = top_n, node_percentage_denominator = "all ARG ORFs",
                       right_link_percentage_denominator = "ARG ORFs within the receiving setting"),
            file.path(out_dir, "Fig3b_summary.csv"), row.names = FALSE, fileEncoding = "UTF-8")
  static_files <- svg
  static_errors <- character()
  if (export_static) {
    if (!requireNamespace("rsvg", quietly = TRUE)) {
      static_errors <- "Optional PNG/PDF conversion skipped: package rsvg is not installed. Native vector SVG and all CSV tables were saved."
      message(static_errors)
    } else {
      for (ext in c("png", "pdf")) {
        target <- file.path(out_dir, paste0("Fig3b_ORF_original_style.", ext))
        err <- tryCatch({
          if (ext == "png") rsvg::rsvg_png(svg, file = target, width = width * 2, height = height * 2)
          else rsvg::rsvg_pdf(svg, file = target, width = width, height = height)
          NULL
        }, error = identity)
        if (inherits(err, "error")) {
          msg <- paste0("Optional ", toupper(ext), " conversion failed: ", conditionMessage(err), ". Native SVG and CSV tables remain available.")
          static_errors <- c(static_errors, msg); message(msg)
        } else if (file.exists(target)) static_files <- c(static_files, target)
      }
    }
  }
  message("Saved Fig3b native SVG and audit tables. Total ARG ORFs: ", format(p$total, scientific = FALSE), ".")
  invisible(list(widget = widget, html = html, svg = svg, out_dir = out_dir, data = p, geometry = geometry,
                 static_files = static_files, static_errors = static_errors))
}

run_fig3b <- function(data_dir, out_dir, render = TRUE,
                      input_file = file.path(data_dir, "02_Fig3b_type_context_setting_counts.tsv")) {
  if (!render) return(invisible(prepare_fig3b_orf(input_file, top_n = 10,
                                                type_order = FIG3B_REFERENCE_TYPE_ORDER)))
  plot_fig3b_orf(input_file, out_dir, top_n = 10, background = "transparent",
                 type_order = FIG3B_REFERENCE_TYPE_ORDER)
}

fig3c_settings <- c('Hospital','WWTP','Community','Wet market')

prepare_fig3c_data <- function(points, metadata = NULL) {
  if (is.null(metadata)) metadata <- points
  aliases <- c(sample='Sample',total_arg_annotations='total_ARG_ORFs',
    mobile_arg_annotations='mobile_ARG_ORFs',plasmid_arg_annotations='plasmid_ARG_ORFs',
    virus_arg_annotations='phage_ARG_ORFs',chromosome_arg_annotations='chromosome_ARG_ORFs')
  for (key in names(aliases)) {
    alias <- aliases[[key]]
    if (!key %in% names(points) && alias %in% names(points)) points[[key]] <- points[[alias]]
    else if (key %in% names(points) && alias %in% names(points) &&
             !isTRUE(all.equal(as.character(points[[key]]),as.character(points[[alias]])))) {
      stop('Conflicting ORF input fields: ',key,' and ',alias)
    }
  }
  required <- c(names(aliases),'mobile_fraction')
  if (!all(required %in% names(points))) stop('Mobile-fraction input requires ORF counts or explicitly named annotation counts; contig counts are not ORF counts.')
  points$sample <- trimws(as.character(points$sample))
  if (anyNA(points$sample)||any(!nzchar(points$sample))||anyDuplicated(points$sample)) stop('Missing or duplicate sample ID in fraction input.')
  date_value <- function(x) {
    x <- as.character(x)
    x <- ifelse(grepl('^[0-9]{8}$',x),paste0(substr(x,1,4),'-',substr(x,5,6),'-',substr(x,7,8)),x)
    suppressWarnings(as.Date(x))
  }
  setting_value <- function(x) {
    x <- trimws(as.character(x))
    x[x %in% c('Wet Market','Wet_market','wet market')] <- 'Wet market'
    x
  }
  m <- as.data.frame(metadata,stringsAsFactors=FALSE)
  if (!'Sample' %in% names(m) && 'sample' %in% names(m)) m$Sample <- m$sample
  if (!all(c('Sample','City','Sample_Date') %in% names(m))) stop('Canonical metadata lacks required fields.')
  if (!'Setting' %in% names(m)) m$Setting <- m$Sample_Type
  if (!'PhysicalSite' %in% names(m)) {
    if (!'Sample_Rename' %in% names(m)) stop('Canonical metadata needs PhysicalSite or Sample_Rename.')
    m$PhysicalSite <- sub('_[0-9]{8}$','',as.character(m$Sample_Rename))
  }
  m$Setting <- setting_value(m$Setting)
  m$Sample_Date <- date_value(m$Sample_Date)
  if (!'SamplingMonth' %in% names(m)) m$SamplingMonth <- format(m$Sample_Date,'%Y-%m')
  keep <- c('Sample','City','Sample_Date','PhysicalSite','Setting','SamplingMonth')
  m <- m[keep]
  for (key in setdiff(keep,'Sample_Date')) m[[key]] <- trimws(as.character(m[[key]]))
  if (anyNA(m)||any(vapply(m[setdiff(keep,'Sample_Date')],function(x)any(!nzchar(x)),logical(1)))||
      anyDuplicated(m$Sample)||any(!m$Setting %in% fig3c_settings)||
      any(m$SamplingMonth != format(m$Sample_Date,'%Y-%m'))) stop('Incomplete, duplicate or inconsistent canonical metadata.')
  site_map <- unique(m[c('PhysicalSite','City','Setting')])
  if (anyDuplicated(site_map$PhysicalSite)) stop('Physical-site metadata maps to multiple cities or settings.')
  count_fields <- setdiff(names(aliases),'sample')
  for (key in c(count_fields,'mobile_fraction')) points[[key]] <- suppressWarnings(as.numeric(points[[key]]))
  if (any(points$total_arg_annotations < 0,na.rm=TRUE)) stop('Negative ORF denominator.')
  defined <- is.finite(points$total_arg_annotations) & points$total_arg_annotations > 0
  excluded <- points[!defined,,drop=FALSE]
  excluded$Sample <- excluded$sample
  excluded$exclusion_reason <- ifelse(is.na(excluded$total_arg_annotations),'Missing classified ORF denominator','Zero or nonfinite classified ORF denominator')
  points <- points[defined,,drop=FALSE]
  for (key in count_fields) {
    if (any(!is.finite(points[[key]]))||any(points[[key]]<0)||any(abs(points[[key]]-round(points[[key]]))>1e-8)) stop('Invalid ARG ORF count: ',key)
  }
  if (any(!is.finite(points$mobile_fraction))||any(points$mobile_fraction<0|points$mobile_fraction>1)) stop('Missing or invalid mobile fraction.')
  if (any(points$mobile_arg_annotations != points$plasmid_arg_annotations+points$virus_arg_annotations)||
      any(points$total_arg_annotations != points$mobile_arg_annotations+points$chromosome_arg_annotations)) stop('Classified ORF compartments do not sum to denominator.')
  if (any(abs(points$mobile_fraction-points$mobile_arg_annotations/points$total_arg_annotations)>1e-12)) stop('Mobile-fraction ratio differs from ORF counts.')
  if (all(c('unclassified_ARG_ORFs','final_ARG_ORFs') %in% names(points))) {
    unclassified <- suppressWarnings(as.numeric(points$unclassified_ARG_ORFs))
    final <- suppressWarnings(as.numeric(points$final_ARG_ORFs))
    if (any(!is.finite(unclassified))||any(unclassified<0)||any(unclassified!=round(unclassified))||
        any(!is.finite(final))||any(final != points$total_arg_annotations+unclassified)) stop('Final ORF total must equal classified plus unclassified ORFs.')
  }
  j <- match(m$Sample,points$sample)
  missing <- m[is.na(j),,drop=FALSE]
  outside <- points[!points$sample %in% m$Sample,,drop=FALSE]
  pp <- points[j[!is.na(j)],,drop=FALSE]
  d <- m[!is.na(j),,drop=FALSE]
  for (key in setdiff(names(pp),c('Sample','sample'))) {
    d[[if(key %in% c(keep,'Sample_Type')) paste0('Input_',key) else key]] <- pp[[key]]
  }
  disagreed <- rep(FALSE,nrow(d))
  for (key in intersect(c(keep,'Sample_Type'),names(pp))) {
    if (key=='Sample') next
    expected <- if(key=='Sample_Type') as.character(d$Setting) else as.character(d[[key]])
    actual <- if(key %in% c('Setting','Sample_Type')) setting_value(pp[[key]]) else if(key=='Sample_Date') as.character(date_value(pp[[key]])) else trimws(as.character(pp[[key]]))
    disagreed <- disagreed | is.na(actual) | actual != expected
  }
  disagreements <- d[disagreed,c(keep,grep('^Input_',names(d),value=TRUE)),drop=FALSE]
  d$Setting <- factor(d$Setting,levels=fig3c_settings)
  for (key in c('City','SamplingMonth','PhysicalSite')) d[[key]] <- factor(d[[key]])
  list(data=d,missing_samples=missing,outside_samples=outside,excluded_samples=excluded,
       metadata_disagreements=disagreements,metadata=m)
}

fit_fig3c_lmm <- function(d) {
  required <- c('mobile_fraction','Setting','City','SamplingMonth','PhysicalSite')
  if(!all(required %in% names(d))) stop('Missing model variables: ',paste(setdiff(required,names(d)),collapse=', '))
  if(anyNA(d[required])) stop('Incomplete model data; do not silently omit rows.')
  if(!setequal(as.character(d$Setting),fig3c_settings)) stop('All four settings are required for six contrasts.')
  d$Setting <- factor(d$Setting,levels=fig3c_settings)
  for(key in c('City','SamplingMonth','PhysicalSite')) d[[key]] <- factor(d[[key]])
  if(any(!is.finite(d$mobile_fraction))||any(d$mobile_fraction<0|d$mobile_fraction>1)) stop('Mobile fractions must be finite and between zero and one.')
  warnings_seen <- character(); messages_seen <- character()
  model_formula <- mobile_fraction ~ Setting + City + SamplingMonth + (1|PhysicalSite)
  fit <- withCallingHandlers(lme4::lmer(model_formula,data=d,REML=TRUE,
    control=lme4::lmerControl(optimizer='bobyqa'),na.action=na.fail),
    warning=function(w){warnings_seen<<-c(warnings_seen,conditionMessage(w));invokeRestart('muffleWarning')},
    message=function(m){messages_seen<<-c(messages_seen,conditionMessage(m));invokeRestart('muffleMessage')})
  beta <- lme4::fixef(fit);vc <- as.matrix(vcov(fit))
  if(length(attr(lme4::getME(fit,'X'),'col.dropped'))) stop('Rank-deficient fixed-effects design; contrast inference needs review.')
  nd <- d[rep(1,4),,drop=FALSE];nd$Setting <- factor(fig3c_settings,levels=fig3c_settings)
  X <- model.matrix(~Setting+City+SamplingMonth,nd)[,names(beta),drop=FALSE]
  contrasts <- do.call(rbind,lapply(combn(1:4,2,simplify=FALSE),function(j) {
    L <- X[j[1],]-X[j[2],];estimate <- sum(L*beta);se <- sqrt(drop(L%*%vc%*%L));z <- estimate/se
    data.frame(group1=fig3c_settings[j[1]],group2=fig3c_settings[j[2]],estimate=estimate,
      SE=se,Wald_z=z,CI_low=estimate-qnorm(.975)*se,CI_high=estimate+qnorm(.975)*se,
      P=2*pnorm(-abs(z)),Scale='Original mobile-fraction difference',stringsAsFactors=FALSE)
  }))
  contrasts$q_BH <- p.adjust(contrasts$P,'BH')
  contrasts$Comparison <- paste(contrasts$group1,contrasts$group2,sep='-')
  pmat <- matrix(1,4,4,dimnames=list(fig3c_settings,fig3c_settings))
  for(i in seq_len(nrow(contrasts))) pmat[contrasts$group1[i],contrasts$group2[i]] <- pmat[contrasts$group2[i],contrasts$group1[i]] <- contrasts$q_BH[i]
  ord <- fig3c_settings[order(drop(X%*%beta),decreasing=TRUE)]
  cld <- multcompView::multcompLetters(pmat[ord,ord],threshold=.05)$Letters
  letters <- data.frame(Setting=names(cld),Letters=unname(cld),stringsAsFactors=FALSE)
  for(i in seq_len(nrow(contrasts))) {
    shared <- length(intersect(strsplit(cld[contrasts$group1[i]],'')[[1]],strsplit(cld[contrasts$group2[i]],'')[[1]]))>0
    stopifnot(shared==(contrasts$q_BH[i]>=.05))
  }
  conv <- paste(unlist(fit@optinfo$conv$lme4$messages),collapse='; ')
  diagnostics <- data.frame(Outcome='Mobile ARG fraction',Family='Gaussian LMM',
    Formula=paste(deparse(model_formula),collapse=' '),Transform='None',Link='identity',REML=TRUE,
    n_samples=nrow(d),n_sites=length(unique(d$PhysicalSite)),n_cities=length(unique(d$City)),n_months=length(unique(d$SamplingMonth)),
    Inference='Two-sided asymptotic Wald z; 95% unadjusted contrast CI; BH across exactly six setting pairs',
    Singular=lme4::isSingular(fit),Convergence_messages=conv,Warnings=paste(unique(warnings_seen),collapse='; '),
    Messages=paste(unique(messages_seen),collapse='; '),
    Fitted_below_zero=sum(fitted(fit)<0),Fitted_above_one=sum(fitted(fit)>1),stringsAsFactors=FALSE)
  list(fit=fit,contrasts=contrasts,letters=letters,diagnostics=diagnostics)
}

plot_fig3c_lmm <- function(d,descriptive,letters) {
  palette <- c(Hospital='#E15759',WWTP='#4E79A7',Community='#76B7B2','Wet market'='#E3BA22')
  d$x_id <- as.integer(d$Setting)
  polygons <- do.call(rbind,lapply(seq_along(fig3c_settings),function(i) {
    values <- d$mobile_fraction[d$Setting==fig3c_settings[i]]
    den <- density(values,n=512,from=min(values),to=max(values))
    data.frame(Setting=fig3c_settings[i],x=c(i+.45*den$y/max(den$y),rep(i,length(den$x))),y=c(den$x,rev(den$x)))
  }))
  label <- merge(descriptive,letters,by='Setting',sort=FALSE)
  label$x_id <- match(label$Setting,fig3c_settings)
  label$mean_label <- sprintf('%.2f',label$Mean_mobile_fraction)
  span <- max(diff(range(d$mobile_fraction)),.1);top <- max(d$mobile_fraction)+.30*span
  set.seed(123)
  ggplot2::ggplot()+
    ggplot2::geom_polygon(data=polygons,ggplot2::aes(x=x,y=y,fill=Setting,group=Setting),alpha=.82,colour=NA)+
    ggplot2::geom_jitter(data=d,ggplot2::aes(x=x_id-.2,y=mobile_fraction,colour=Setting),width=.2,height=0,alpha=.82,size=1,shape=16,stroke=0)+
    ggplot2::geom_boxplot(data=d,ggplot2::aes(x=x_id,y=mobile_fraction,group=Setting,colour=Setting),width=.2,outlier.shape=NA,linewidth=.6,fill='white')+
    ggplot2::stat_summary(data=d,ggplot2::aes(x=x_id,y=mobile_fraction,group=Setting,colour=Setting),fun=median,geom='crossbar',width=.06,linewidth=.6)+
    ggplot2::geom_text(data=label,ggplot2::aes(x=x_id,y=top,label=mean_label,colour=Setting),size=4,vjust=2)+
    ggplot2::geom_text(data=label,ggplot2::aes(x=x_id,y=top,label=Letters,colour=Setting),size=4,vjust=0)+
    ggplot2::scale_fill_manual(values=palette,drop=FALSE)+ggplot2::scale_colour_manual(values=palette,drop=FALSE)+
    ggplot2::scale_x_continuous(breaks=1:4,labels=fig3c_settings)+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(.02,0)))+
    ggplot2::coord_cartesian(ylim=c(max(0,min(d$mobile_fraction)-.03*span),top+.1*span),clip='on')+
    ggplot2::labs(x=NULL,y='Mobile ARG fraction')+ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),panel.border=ggplot2::element_rect(colour='black',fill=NA,linewidth=.65),
      axis.line=ggplot2::element_blank(),axis.title.y=ggplot2::element_text(size=10,face='plain'),
      axis.text.x=ggplot2::element_text(colour='black',margin=ggplot2::margin(t=3)),axis.text.y=ggplot2::element_text(colour='black'),
      axis.ticks=ggplot2::element_line(colour='black',linewidth=.42),axis.ticks.length=grid::unit(.13,'cm'),
      legend.position='none',plot.margin=ggplot2::margin(5,3,3,5))
}

run_fig3c_lmm <- function(data_file,metadata_file=NULL,out_dir='results/panel_c') {
  required <- c('lme4','multcompView','ggplot2')
  if(!is.null(metadata_file) && grepl('\\.xlsx?$',metadata_file,ignore.case=TRUE)) required <- c(required,'readxl')
  missing_packages <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing_packages)) stop('Missing packages: ',paste(missing_packages,collapse=', '))

  points <- utils::read.delim(data_file,check.names=FALSE,stringsAsFactors=FALSE,encoding='UTF-8')
  meta <- if(is.null(metadata_file)) NULL else if(grepl('\\.xlsx?$',metadata_file,ignore.case=TRUE)) as.data.frame(readxl::read_excel(metadata_file)) else utils::read.csv(metadata_file,check.names=FALSE,stringsAsFactors=FALSE,encoding='UTF-8')
  prepared <- prepare_fig3c_data(points,meta)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  write_table <- function(x,name) utils::write.csv(x,file.path(out_dir,name),row.names=FALSE,na='')
  write_table(prepared$missing_samples,'canonical_samples_without_fraction.csv')
  write_table(prepared$outside_samples,'fraction_samples_outside_canonical_cohort.csv')
  write_table(prepared$excluded_samples,'excluded_undefined_denominator_samples.csv')
  write_table(prepared$metadata_disagreements,'metadata_disagreements_resolved_to_canonical.csv')
  write_table(prepared$data,'sample_model_input.csv')
  d <- prepared$data
  coverage <- data.frame(Canonical_samples=nrow(prepared$metadata),Input_fraction_samples=nrow(points),
    Analyzed_samples=nrow(d),Canonical_samples_without_fraction=nrow(prepared$missing_samples),
    Undefined_denominator_samples=nrow(prepared$excluded_samples),
    Input_samples_outside_canonical=nrow(prepared$outside_samples),Metadata_disagreement_rows=nrow(prepared$metadata_disagreements),
    Policy='Exact Sample ID join; canonical setting/city/date/site; no outcome or denominator imputation; observed 0/1 retained')
  write_table(coverage,'cohort_coverage.csv')
  result <- fit_fig3c_lmm(d)
  descriptive <- do.call(rbind,lapply(fig3c_settings,function(setting) {
    q <- d[d$Setting==setting,,drop=FALSE]
    data.frame(Setting=setting,Samples=nrow(q),Physical_sites=length(unique(q$PhysicalSite)),
      Mean_mobile_fraction=mean(q$mobile_fraction),Median_mobile_fraction=median(q$mobile_fraction),
      SD_mobile_fraction=sd(q$mobile_fraction),Min=min(q$mobile_fraction),Max=max(q$mobile_fraction),
      Total_classified_ARG_ORFs=sum(q$total_arg_annotations),Mobile_ARG_ORFs=sum(q$mobile_arg_annotations),
      Pooled_mobile_fraction=sum(q$mobile_arg_annotations)/sum(q$total_arg_annotations))
  }))
  write_table(descriptive,'setting_descriptive_means.csv')
  write_table(result$contrasts,'six_setting_contrasts_BH.csv')
  write_table(result$letters,'setting_letters.csv')
  write_table(result$diagnostics,'model_diagnostics.csv')
  write_table(as.data.frame(lme4::VarCorr(result$fit)),'variance_components.csv')
  beta <- lme4::fixef(result$fit)
  write_table(data.frame(Term=names(beta),Coefficient=unname(beta),SE=sqrt(diag(as.matrix(vcov(result$fit))))),'fixed_effects.csv')
  saveRDS(result$fit,file.path(out_dir,'fitted_model.rds'))
  restored_fit <- readRDS(file.path(out_dir,'fitted_model.rds'))
  stopifnot(inherits(restored_fit,'lmerMod'),nobs(restored_fit)==nrow(d),
    isTRUE(all.equal(lme4::fixef(restored_fit),lme4::fixef(result$fit),tolerance=0)),
    isTRUE(all.equal(as.numeric(logLik(restored_fit)),as.numeric(logLik(result$fit)),tolerance=0)))
  residuals <- data.frame(Sample=d$Sample,Fitted=fitted(result$fit),Residual=residuals(result$fit),
    Pearson_residual=residuals(result$fit,type='pearson'))
  write_table(residuals,'sample_residuals.csv')
  grDevices::pdf(file.path(out_dir,'residual_diagnostics.pdf'),width=8,height=4)
  old_par <- par(mfrow=c(1,2))
  plot(residuals$Fitted,residuals$Pearson_residual,xlab='Fitted mobile fraction',ylab='Pearson residual');abline(h=0,lty=2)
  qqnorm(residuals$Pearson_residual);qqline(residuals$Pearson_residual)
  par(old_par);grDevices::dev.off()
  plot <- plot_fig3c_lmm(d,descriptive,result$letters)
  ggplot2::ggsave(file.path(out_dir,'Fig3c.pdf'),plot,width=3.8,height=2.6)
  ggplot2::ggsave(file.path(out_dir,'Fig3c.png'),plot,width=3.8,height=2.6,dpi=300)
  writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))
  result$data <- d;result$descriptive <- descriptive;result$coverage <- coverage;result$plot <- plot
  message('Fig. 3c complete: ',nrow(d),' observed canonical samples; ',nrow(prepared$missing_samples),' canonical samples lack fraction data.')
  if(isTRUE(result$diagnostics$Singular)||nzchar(result$diagnostics$Convergence_messages)||nzchar(result$diagnostics$Warnings)) warning('Review exported model diagnostics.')
  invisible(result)
}

.fig3de_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages(invisible(lapply(packages, function(package) {
    library(package, character.only = TRUE)
  })))
}

.fig3de_remap_letters <- function(raw, reference = "a") {

  split_symbols <- function(x) unlist(strsplit(x[!is.na(x) & nzchar(x)], "", fixed = TRUE), use.names = FALSE)
  symbols <- unique(c(split_symbols(reference), split_symbols(raw)))
  alphabet <- c(letters, LETTERS)
  if (length(symbols) > length(alphabet)) stop("Too many compact-letter symbols", call. = FALSE)
  mapping <- stats::setNames(alphabet[seq_along(symbols)], symbols)
  vapply(raw, function(value) {
    if (is.na(value)) return(NA_character_)
    if (!nzchar(value)) return("")
    paste0(sort(unname(mapping[strsplit(value, "", fixed = TRUE)[[1]]])), collapse = "")
  }, character(1), USE.NAMES = FALSE)
}

.fig3de_dunn_two_groups <- function(values, groups) {

  groups <- droplevels(factor(groups))
  if (length(values) != length(groups) || anyNA(groups) || any(!is.finite(values)) || nlevels(groups) != 2L) {
    stop("Two-setting Dunn calculation requires finite observations in exactly two settings", call. = FALSE)
  }
  n_total <- length(values)
  ties <- as.numeric(table(values))
  rank_variance <- n_total * (n_total + 1) / 12 - sum(ties^3 - ties) / (12 * (n_total - 1))
  if (!is.finite(rank_variance) || rank_variance <= 0) {
    stop("Dunn calculation is undefined because the pooled ranks have zero variance", call. = FALSE)
  }
  mean_ranks <- tapply(rank(values, ties.method = "average"), groups, mean)
  n_group <- as.numeric(table(groups))
  z <- unname((mean_ranks[[1]] - mean_ranks[[2]]) / sqrt(rank_variance * sum(1 / n_group)))
  p <- 2 * stats::pnorm(-abs(z))
  data.frame(Comparison = paste(levels(groups), collapse = " - "), Z = z,
             P.unadj = p, P.adj = stats::p.adjust(p, method = "BH"))
}

.fig3d_hospital_flags <- function(site, letters) {
  reference <- letters[as.character(site) == 'Hospital']
  reference <- reference[!is.na(reference) & nzchar(reference)]
  if (length(reference) != 1L) return(rep('Not tested',length(site)))
  symbols <- strsplit(reference,'',fixed=TRUE)[[1]]
  vapply(seq_along(site),function(i) {
    if (is.na(letters[i]) || !nzchar(letters[i])) return('Not tested')
    if (as.character(site[i])=='Hospital') return('Same')
    if (length(intersect(symbols,strsplit(letters[i],'',fixed=TRUE)[[1]]))) 'Same' else 'Different'
  },character(1))
}

run_fig3d <- function(data_dir, out_dir) {
  .fig3de_packages(c("data.table", "dplyr", "ggplot2", "stringr", "FSA", "multcompView", "readr", "scales"))
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .figure_results_dir <- normalizePath(out_dir, mustWork = TRUE)
  fig3d_dunn_results <- list()
  fig3d_global_results <- list()
  fig3d_test_status <- list()
  fig3d_site_levels <- c("Hospital", "WWTP", "Community", "Wet market")

  fig3d_make_dunn_letters <- function(df, group_col, value_col, species) {
    dunn_df <- df %>%
      dplyr::mutate(
        .group = factor(.data[[group_col]], levels = fig3d_site_levels),
        .value = as.numeric(.data[[value_col]])
      ) %>%
      dplyr::filter(!is.na(.group), !is.na(.value))

    groups_present <- levels(droplevels(dunn_df$.group))
    if (length(groups_present) < 2L) {
      fig3d_test_status[[length(fig3d_test_status) + 1L]] <<- tibble::tibble(species = species, status = "not_tested", detail = "Fewer than two observed settings", n_groups = length(groups_present))
      return(tibble::tibble(
        site = factor(groups_present, levels = fig3d_site_levels),
        Letters = ""
      ))
    }

    tryCatch({
      global <- stats::kruskal.test(dunn_df$.value, droplevels(dunn_df$.group))
      fig3d_global_results[[length(fig3d_global_results) + 1L]] <<- data.frame(
        species = species, statistic = unname(global$statistic),
        df = unname(global$parameter), P = global$p.value, n = nrow(dunn_df),
        n_groups = length(groups_present), stringsAsFactors = FALSE)
      if (length(groups_present) == 2L) {
        dunn_result <- .fig3de_dunn_two_groups(dunn_df$.value, dunn_df$.group)
        engine_detail <- "Two-setting compatibility calculation: two-sided Dunn with pooled-rank tie correction; BH over one pair"
      } else {
        dunn_fit <- NULL
        invisible(utils::capture.output(
          dunn_fit <- FSA::dunnTest(
            x = dunn_df$.value,
            g = dunn_df$.group,
            method = "bh"
          )
        ))
        dunn_result <- dunn_fit$res
        engine_detail <- "FSA::dunnTest(method='bh')"
      }
      if (any(!is.finite(dunn_result$P.adj)) || any(dunn_result$P.adj < 0 | dunn_result$P.adj > 1)) {
        stop("Dunn-BH returned invalid adjusted p values", call. = FALSE)
      }
      dunn_result <- dunn_result %>%
        dplyr::mutate(Comparison = stringr::str_replace_all(Comparison, " - ", "-"))
      fig3d_dunn_results[[length(fig3d_dunn_results) + 1L]] <<- dplyr::mutate(dunn_result, species = species, .before = 1)
      adjusted_p <- dunn_result$P.adj
      names(adjusted_p) <- dunn_result$Comparison
      letter_result <- multcompView::multcompLetters(adjusted_p, threshold = 0.05)
      fig3d_test_status[[length(fig3d_test_status) + 1L]] <<- tibble::tibble(species = species, status = "ok", detail = engine_detail, n_groups = length(groups_present))
      tibble::tibble(
        site = factor(names(letter_result$Letters), levels = fig3d_site_levels),
        Letters = unname(letter_result$Letters)
      )
    }, error = function(e) {
      fig3d_test_status[[length(fig3d_test_status) + 1L]] <<- tibble::tibble(species = species, status = "error", detail = conditionMessage(e), n_groups = length(groups_present))
      readr::write_csv(dplyr::bind_rows(fig3d_test_status), file.path(.figure_results_dir, "Fig3d_Dunn_test_status.csv"))
      stop("Panel d Dunn test failed for ", species, ": ", conditionMessage(e), call. = FALSE)
    })
  }

  fig3d_mag_burden <- data.table::fread(
    file.path(data_dir, "Fig3d_MAG_total_burden.tsv"),
    sep = "\t",
    data.table = FALSE,
    encoding = "UTF-8"
  ) %>%
    dplyr::mutate(site = factor(site, levels = fig3d_site_levels))

  fig3d_site_species <- data.table::fread(
    file.path(data_dir, "Fig3d_species_burden_summary.tsv"),
    sep = "\t",
    data.table = FALSE,
    encoding = "UTF-8"
  ) %>%
    dplyr::mutate(site = factor(site, levels = fig3d_site_levels))

  fig3d_top_species_n <- 15L
  fig3d_species_order <- fig3d_site_species %>%
    dplyr::group_by(species) %>%
    dplyr::summarise(total_n = sum(n_mag, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(total_n)) %>%
    dplyr::slice_head(n = fig3d_top_species_n)

  fig3d_species_levels <- fig3d_species_order$species
  fig3d_letters <- fig3d_mag_burden %>%
    dplyr::filter(species %in% fig3d_species_levels) %>%
    dplyr::mutate(
      species = as.character(species),
      site = factor(site, levels = fig3d_site_levels)
    ) %>%
    dplyr::group_by(species) %>%
    dplyr::group_modify(~ fig3d_make_dunn_letters(.x, "site", "total_burden", .y$species)) %>%
    dplyr::ungroup()

  fig3d_plot_df <- fig3d_site_species %>%
    dplyr::filter(species %in% fig3d_species_levels) %>%
    dplyr::left_join(fig3d_species_order, by = "species") %>%
    dplyr::mutate(species = as.character(species)) %>%
    dplyr::left_join(fig3d_letters, by = c("species", "site")) %>%
    dplyr::mutate(Letters_raw = dplyr::coalesce(Letters, "")) %>%
    dplyr::group_by(species) %>%
    dplyr::group_modify(~ {
      species_df <- .x
      hospital_raw <- species_df$Letters_raw[species_df$site == "Hospital"]
      hospital_raw <- hospital_raw[!is.na(hospital_raw) & hospital_raw != ""]
      hospital_raw <- if (length(hospital_raw) == 0L) "a" else hospital_raw[[1]]
      species_df %>%
        dplyr::mutate(
          Letters = .fig3de_remap_letters(Letters_raw, reference = hospital_raw),
          diff_flag = .fig3d_hospital_flags(site, Letters)
        )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(species = factor(species, levels = rev(fig3d_species_levels)))

  readr::write_csv(
    fig3d_plot_df %>%
      dplyr::select(
        species, site, prop_pathogen_mags, mean_burden,
        Letters_raw, Letters, diff_flag
      ),
    file.path(.figure_results_dir, "Fig3d_Dunn_BH_letters.csv")
  )

  fig3d_fill_cols <- c("#F4FBF9", "#E2F3EE", "#C8E8E0", "#AEE0D6", "#8FD0C7", "#74BEB4")
  fig3d_plot <- ggplot2::ggplot(fig3d_plot_df, ggplot2::aes(x = site, y = species)) +
    ggplot2::geom_point(
      ggplot2::aes(fill = mean_burden, colour = diff_flag),
      shape = 21,
      size = 8.3,
      stroke = 0.6,
      alpha = 0.97
    ) +
    ggplot2::geom_text(
      ggplot2::aes(label = Letters, colour = diff_flag),
      size = 4.8,
      fontface = "plain"
    ) +
    ggplot2::scale_fill_gradientn(
      colours = fig3d_fill_cols,
      name = "Mean ARG burden\nper MAG",
      limits = range(fig3d_plot_df$mean_burden, na.rm = TRUE),
      oob = scales::squish
    ) +
    ggplot2::scale_colour_manual(
      values = c("Same" = "grey25", "Different" = "#D73027", "Not tested" = "grey25"),
      guide = "none"
    ) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      axis.text.y = ggplot2::element_text(face = "italic", color = "#C9332B", size = 13),
      axis.text.x = ggplot2::element_text(
        angle = 25, hjust = 1, vjust = 1, color = "black", size = 13
      ),
      panel.grid.major = ggplot2::element_line(color = "#eeeeee", linewidth = 0.35),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.direction = "horizontal",
      legend.title = ggplot2::element_text(size = 11, colour = "black"),
      legend.text = ggplot2::element_text(size = 10, colour = "black"),
      plot.title = ggplot2::element_blank()
    )

  ggplot2::ggsave(
    file.path(.figure_results_dir, "Fig3d.pdf"),
    fig3d_plot,
    width = 7,
    height = 6
  )
  ggplot2::ggsave(
    file.path(.figure_results_dir, "Fig3d.png"),
    fig3d_plot,
    width = 7,
    height = 6,
    dpi = 600
  )
  contrasts <- dplyr::bind_rows(fig3d_dunn_results)
  global_tests <- dplyr::bind_rows(fig3d_global_results)
  status <- dplyr::bind_rows(fig3d_test_status)
  readr::write_csv(contrasts, file.path(.figure_results_dir, "Fig3d_Dunn_BH_contrasts.csv"))
  readr::write_csv(global_tests, file.path(.figure_results_dir, "Fig3d_Kruskal_Wallis_tests.csv"))
  readr::write_csv(status, file.path(.figure_results_dir, "Fig3d_Dunn_test_status.csv"))
  readr::write_csv(fig3d_plot_df, file.path(.figure_results_dir, "Fig3d_plot_data.csv"))
  invisible(list(plot = fig3d_plot, plot_data = fig3d_plot_df, contrasts = contrasts,
                 global_tests = global_tests, test_status = status))
}

.fig3e_union_counts <- function(bubble, species, top_types, presence = NULL,
                               registry = NULL, premerged = FALSE) {
  group_fields <- c('site','species')
  key <- function(x,fields) do.call(paste,c(x[fields],sep='\r'))
  normalize_type <- function(x) {
    x <- trimws(as.character(x))
    aliases <- stats::setNames(c('beta-lactam','beta-lactam','MLS'),
      c('beta_lactam','\u03b2-lactam','macrolide-lincosamide-streptogramin'))
    mapped <- unname(aliases[tolower(x)]);x[!is.na(mapped)] <- mapped[!is.na(mapped)]
    known <- c(top_types,'Others','Rifamycin')
    idx <- match(tolower(x),tolower(known));x[!is.na(idx)] <- known[idx[!is.na(idx)]]
    x
  }
  if (!all(c(group_fields,'n_mag_site') %in% names(species))) stop('Species table requires site, species and n_mag_site.')
  species <- as.data.frame(species[c(group_fields,'n_mag_site')],stringsAsFactors=FALSE)
  for (field in group_fields) species[[field]] <- trimws(as.character(species[[field]]))
  species$n_mag_site <- as.numeric(species$n_mag_site)
  if (anyNA(species)||anyDuplicated(key(species,group_fields))||any(species$n_mag_site<=0)||
      any(species$n_mag_site!=round(species$n_mag_site))) stop('Invalid or duplicate species MAG denominator.')
  top_types <- setdiff(normalize_type(top_types),c('Others','Rifamycin'))
  if (is.null(presence) || is.null(registry)) {
    if (!is.null(presence)||!is.null(registry)) stop('Both MAG presence and MAG registry inputs are required.')
    if (!isTRUE(premerged)) stop('Cannot recover an Others MAG union from marginal counts. Supply MAG presence plus the unchanged MAG registry, or an explicitly premerged union table.')
    required <- c(group_fields,'arg_type_plot','n_type_carrying_mags','pct_mags')
    if (!all(required %in% names(bubble))) stop('Premerged union table lacks carrier counts or percentages.')
    result <- as.data.frame(bubble[required],stringsAsFactors=FALSE)
    result$arg_type_plot <- normalize_type(result$arg_type_plot)
    if (any(!result$arg_type_plot %in% c(top_types,'Others'))||anyDuplicated(key(result,c(group_fields,'arg_type_plot')))) stop('Premerged union table must have one row per displayed category and no Rifamycin or other unmerged categories.')
    result$n_type_carrying_mags <- as.numeric(result$n_type_carrying_mags)
    result$pct_mags <- as.numeric(result$pct_mags)
  } else {
    registry <- as.data.frame(registry,stringsAsFactors=FALSE)
    presence <- as.data.frame(presence,stringsAsFactors=FALSE)
    if (!all(c(group_fields,'mag_uid') %in% names(registry))) stop('MAG registry requires mag_uid, site and species.')
    registry <- registry[c('mag_uid',group_fields)]
    for (field in names(registry)) registry[[field]] <- trimws(as.character(registry[[field]]))
    if (anyNA(registry)||any(!nzchar(registry$mag_uid))||anyDuplicated(registry$mag_uid)) stop('Missing or duplicate MAG registry IDs.')
    registry <- registry[key(registry,group_fields) %in% key(species,group_fields),,drop=FALSE]
    observed <- table(factor(key(registry,group_fields),levels=key(species,group_fields)))
    if (any(as.numeric(observed)!=species$n_mag_site)) stop('MAG registry does not match the supplied species denominator; preserve the original quality subset.')
    if (!'mag_uid' %in% names(presence) && all(c('Sample','bin_id') %in% names(presence))) {
      if (anyNA(presence[c('Sample','bin_id')])) stop('Missing Sample/bin_id in MAG presence input.')
      presence$mag_uid <- paste(presence$Sample,presence$bin_id,sep='__')
    }
    if (!all(c('mag_uid','arg_type') %in% names(presence))) stop('MAG presence requires mag_uid (or Sample and bin_id) plus arg_type.')
    presence$mag_uid <- trimws(as.character(presence$mag_uid))
    matched <- match(presence$mag_uid,registry$mag_uid)
    presence <- presence[!is.na(matched),,drop=FALSE]
    matched <- matched[!is.na(matched)]
    types <- normalize_type(presence$arg_type)
    if (anyNA(types)||any(!nzchar(types))) stop('Missing ARG type in MAG presence input.')
    grouped <- ifelse(types %in% top_types,types,'Others')
    records <- unique(data.frame(registry[matched,,drop=FALSE],arg_type_plot=grouped,row.names=NULL))
    if (!nrow(records)) stop('MAG presence contains no ARG records for the selected MAG registry.')
    result <- aggregate(list(n_type_carrying_mags=rep(1,nrow(records))),
                        records[c(group_fields,'arg_type_plot')],sum)
    denom <- species$n_mag_site[match(key(result,group_fields),key(species,group_fields))]
    result$pct_mags <- 100*result$n_type_carrying_mags/denom
  }
  denom <- species$n_mag_site[match(key(result,group_fields),key(species,group_fields))]
  n <- result$n_type_carrying_mags
  if (anyNA(denom)||any(!is.finite(n))||any(n<0|n>denom)||any(n!=round(n))||
      any(!is.finite(result$pct_mags))||any(abs(result$pct_mags-100*n/denom)>1e-7)) {
    stop('Union carrier count or percentage is inconsistent with the supplied MAG denominator.')
  }
  result
}

run_fig3e <- function(data_dir, out_dir, mag_presence_file = NULL,
                     mag_registry_file = NULL, bubble_premerged = FALSE, render = TRUE) {
  .fig3de_packages(c("tidyverse", "scales", "grid", "data.table"))
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  out_dir <- normalizePath(out_dir, mustWork = TRUE)
  bubble_file  <- file.path(data_dir, "ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_bubble_pct.tsv")
  burden_file  <- file.path(data_dir, "ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_site_species_type_mean_burden.tsv")
  species_file <- file.path(data_dir, "ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_species_used_by_site.tsv")
  types_file   <- file.path(data_dir, "ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_argtypes_used.tsv")

  out_prefix <- "ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_R_yellow_noRifamycin"

  SITE_LEVELS <- c("Hospital", "WWTP", "Community", "Wet market")

  FIG_WIDTH <- 7.2
  BASE_HEIGHT <- 2.0
  ROW_HEIGHT <- 0.20
  HEIGHT_FACTOR <- 0.75
  MIN_HEIGHT <- 6.2

  BUBBLE_SIZE_RANGE <- c(0.25, 4.2)

  type_cols <- c(
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
    "OPA" = "#9FD58A",
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
    "beta-lactam" = "#7DA34D",
    "Others" = "#D0D0D0"
  )

  yellow_cols <- c(
    "#FFFDF2",
    "#FBF4D2",
    "#F5E7A8",
    "#EFD87D",
    "#E9C95B",
    "#D6AF3E",
    "#B9902A"
  )

  clean_chr <- function(x) trimws(as.character(x))

  fix_beta <- function(x) {
    x <- clean_chr(x)
    ifelse(x %in% c("beta_lactam", "\u03b2-lactam", "beta-lactam"), "beta-lactam", x)
  }

  wrap_type_label <- function(x) {
    x <- fix_beta(x)
    case_when(
      x == "OPA" ~ "OPA",
      x == "MLS" ~ "MLS",
      x == "beta-lactam" ~ "beta-lactam",
      x == "Others" ~ "Others",
      x == "Pleuromutilin/Tiamulin" ~ "Pleuro/\ntiamulin",
      x == "Chloramphenicol" ~ "Chloramphen\nicol",
      nchar(x) > 12 ~ str_wrap(x, width = 11),
      TRUE ~ x
    )
  }

  read_type_order <- function(types_file) {
    type_df <- fread(
      types_file,
      sep = "\t",
      header = TRUE,
      data.table = FALSE,
      encoding = "UTF-8"
    )

    if ("arg_type_plot" %in% names(type_df)) {
      out <- type_df$arg_type_plot
    } else if ("arg_type" %in% names(type_df)) {
      out <- type_df$arg_type
    } else {
      out <- type_df[[1]]
    }

    fix_beta(out)
  }

  files_need <- c(bubble_file, burden_file, species_file, types_file)

  cat("\nWorking directory:\n")
  print(getwd())

  cat("\nInput file check:\n")
  print(data.frame(
    file = files_need,
    exists = file.exists(files_need),
    size = ifelse(file.exists(files_need), file.info(files_need)$size, NA)
  ))

  if (!all(file.exists(files_need))) {
    stop("Some input files are missing.")
  }

  bubble_df <- fread(
    bubble_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    encoding = "UTF-8"
  )

  burden_raw <- fread(
    burden_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    encoding = "UTF-8"
  )

  species_df <- fread(
    species_file,
    sep = "\t",
    header = TRUE,
    data.table = FALSE,
    encoding = "UTF-8"
  )

  type_order_input <- read_type_order(types_file)

  type_order_input <- fix_beta(type_order_input)

  bubble_df <- bubble_df %>%
    mutate(
      site = clean_chr(site),
      species = clean_chr(species),
      arg_type_plot = fix_beta(arg_type_plot),
      pct_mags = as.numeric(pct_mags),
      n_type_carrying_mags = if ("n_type_carrying_mags" %in% names(.)) {
        as.numeric(n_type_carrying_mags)
      } else {
        NA_real_
      }
    ) %>%
    filter(site %in% SITE_LEVELS, arg_type_plot %in% type_order_input)

  species_df <- species_df %>%
    mutate(
      site = clean_chr(site),
      species = clean_chr(species),
      n_mag_site = as.numeric(n_mag_site)
    ) %>%
    filter(site %in% SITE_LEVELS, n_mag_site > 0)

  burden_raw <- burden_raw %>%
    mutate(
      site = clean_chr(site),
      species = clean_chr(species),
      mean_type_burden = as.numeric(mean_type_burden)
    ) %>%
    filter(site %in% SITE_LEVELS)

  type_order_input_no_rif <- setdiff(type_order_input, "Rifamycin")
  top_types_input <- setdiff(type_order_input_no_rif, "Others")

  if (is.null(mag_presence_file) && file.exists(file.path(data_dir,'Fig3e_MAG_ARG_presence.tsv'))) mag_presence_file <- file.path(data_dir,'Fig3e_MAG_ARG_presence.tsv')
  if (is.null(mag_registry_file) && file.exists(file.path(data_dir,'Fig3e_MAG_registry.tsv'))) mag_registry_file <- file.path(data_dir,'Fig3e_MAG_registry.tsv')
  presence <- if (is.null(mag_presence_file)) NULL else data.table::fread(mag_presence_file,data.table=FALSE)
  registry <- if (is.null(mag_registry_file)) NULL else data.table::fread(mag_registry_file,data.table=FALSE)
  bubble_df <- .fig3e_union_counts(bubble_df,species_df,top_types_input,presence,registry,premerged=bubble_premerged)

  if ("arg_type_plot" %in% names(burden_raw)) {
    burden_df <- burden_raw %>%
      mutate(
        arg_type_plot = fix_beta(arg_type_plot),
        arg_type_plot = if_else(arg_type_plot == "Rifamycin", "Others", arg_type_plot),
        arg_type_plot = if_else(arg_type_plot %in% top_types_input, arg_type_plot, "Others")
      ) %>%
      group_by(site, species, arg_type_plot) %>%
      summarise(
        mean_type_burden = sum(mean_type_burden, na.rm = TRUE),
        .groups = "drop"
      )
  } else if ("arg_type" %in% names(burden_raw)) {
    burden_df <- burden_raw %>%
      mutate(
        arg_type = fix_beta(arg_type),
        arg_type_plot = if_else(arg_type %in% top_types_input, arg_type, "Others")
      ) %>%
      group_by(site, species, arg_type_plot) %>%
      summarise(
        mean_type_burden = sum(mean_type_burden, na.rm = TRUE),
        .groups = "drop"
      )
  } else {
    stop("burden_file must contain either arg_type or arg_type_plot.")
  }

  burden_df <- burden_df %>%
    filter(arg_type_plot %in% c(top_types_input, "Others"))

  if (all(is.na(bubble_df$n_type_carrying_mags))) {
    type_order <- bubble_df %>%
      group_by(arg_type_plot) %>%
      summarise(type_abd = mean(pct_mags, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(type_abd)) %>%
      pull(arg_type_plot)
  } else {
    type_order <- bubble_df %>%
      group_by(arg_type_plot) %>%
      summarise(type_abd = sum(n_type_carrying_mags, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(type_abd)) %>%
      pull(arg_type_plot)
  }

  type_order <- c(setdiff(type_order, "Others"), "Others")
  type_order <- type_order[type_order %in% c(top_types_input, "Others")]

  cat("\nARG type order by abundance:\n")
  print(type_order)

  species_df <- species_df %>%
    mutate(site = factor(site, levels = SITE_LEVELS)) %>%
    arrange(site, desc(n_mag_site), species) %>%
    group_by(site) %>%
    mutate(
      label = paste0(species, " (", n_mag_site, ")"),
      row_id = row_number(),
      n_rows_site = n(),
      y_pos = n_rows_site - row_id + 1
    ) %>%
    ungroup()

  label_df <- species_df %>%
    transmute(
      site = factor(site, levels = SITE_LEVELS),
      species,
      label,
      y_pos,
      n_mag_site
    )

  bubble_df <- bubble_df %>%
    left_join(label_df, by = c("site", "species")) %>%
    filter(!is.na(y_pos)) %>%
    mutate(
      site = factor(site, levels = SITE_LEVELS),
      arg_type_plot = factor(arg_type_plot, levels = type_order)
    )

  burden_df <- burden_df %>%
    left_join(label_df, by = c("site", "species")) %>%
    filter(!is.na(y_pos)) %>%
    mutate(
      site = factor(site, levels = SITE_LEVELS),
      arg_type_plot = factor(arg_type_plot, levels = type_order)
    )

  burden_base <- burden_df %>%
    select(site, species, arg_type_plot, mean_type_burden) %>%
    group_by(site, species, arg_type_plot) %>%
    summarise(
      mean_type_burden = sum(mean_type_burden, na.rm = TRUE),
      .groups = "drop"
    )

  burden_grid <- label_df %>%
    select(site, species, label, y_pos, n_mag_site) %>%
    distinct() %>%
    crossing(arg_type_plot = type_order)

  burden_df <- burden_grid %>%
    left_join(
      burden_base,
      by = c("site", "species", "arg_type_plot")
    ) %>%
    mutate(
      mean_type_burden = replace_na(mean_type_burden, 0),
      site = factor(site, levels = SITE_LEVELS),
      arg_type_plot = factor(arg_type_plot, levels = type_order)
    )

  bubble_plot_df <- bubble_df %>%
    left_join(
      burden_df %>%
        select(site, species, arg_type_plot, mean_type_burden),
      by = c("site", "species", "arg_type_plot")
    ) %>%
    mutate(
      mean_type_burden = replace_na(mean_type_burden, 0)
    ) %>%
    filter(pct_mags > 0)

  n_type <- length(type_order)

  label_region_width <- 4.7

  bubble_gap_left <- 0.7
  bubble_start <- label_region_width + bubble_gap_left
  bubble_x <- bubble_start + seq_len(n_type) - 1
  names(bubble_x) <- type_order

  mid_gap <- 1.0

  bar_start <- max(bubble_x) + mid_gap + 1
  bar_width <- n_type

  max_burden <- burden_df %>%
    group_by(site, species, y_pos) %>%
    summarise(total_burden = sum(mean_type_burden, na.rm = TRUE), .groups = "drop") %>%
    summarise(max_burden = max(total_burden, na.rm = TRUE)) %>%
    pull(max_burden)

  if (is.na(max_burden) || max_burden <= 0) max_burden <- 1

  burden_tick_raw <- pretty(c(0, max_burden), n = 4)
  burden_tick_raw <- burden_tick_raw[burden_tick_raw >= 0 & burden_tick_raw <= max_burden]
  if (!0 %in% burden_tick_raw) burden_tick_raw <- c(0, burden_tick_raw)

  scale_burden_x <- function(v) {
    bar_start + (v / max_burden) * bar_width
  }

  label_x <- label_region_width - 0.12

  bubble_plot_df <- bubble_plot_df %>%
    mutate(
      x_pos = as.numeric(bubble_x[as.character(arg_type_plot)])
    )

  burden_rect_df <- burden_df %>%
    arrange(site, y_pos, arg_type_plot) %>%
    group_by(site, species, label, y_pos) %>%
    mutate(
      xmin_raw = cumsum(mean_type_burden) - mean_type_burden,
      xmax_raw = cumsum(mean_type_burden),
      xmin = scale_burden_x(xmin_raw),
      xmax = scale_burden_x(xmax_raw),
      ymin = y_pos - 0.31,
      ymax = y_pos + 0.31
    ) %>%
    ungroup() %>%
    filter(mean_type_burden > 0) %>%
    distinct(site, species, label, y_pos, arg_type_plot, .keep_all = TRUE)

  cat("Rows in burden_rect_df:", nrow(burden_rect_df), "\n")
  cat("Rows in bubble_plot_df:", nrow(bubble_plot_df), "\n")

  separator_x <- max(bubble_x) + mid_gap / 2 + 0.5

  x_breaks <- c(
    bubble_x,
    scale_burden_x(burden_tick_raw)
  )

  x_labels <- c(
    map_chr(type_order, wrap_type_label),
    as.character(burden_tick_raw)
  )

  type_cols_used <- type_cols[type_order]
  type_cols_used[is.na(type_cols_used)] <- "#D0D0D0"

  p <- ggplot() +
    geom_text(
      data = label_df,
      aes(x = label_x, y = y_pos, label = label),
      inherit.aes = FALSE,
      hjust = 1,
      size = 2.0,
      fontface = "italic",
      color = "#C00000"
    ) +
    geom_rect(
      data = burden_rect_df,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax,
        fill = arg_type_plot
      ),
      color = NA
    ) +
    geom_point(
      data = bubble_plot_df,
      aes(
        x = x_pos,
        y = y_pos,
        size = mean_type_burden,
        color = pct_mags
      ),
      alpha = 0.92
    ) +
    geom_vline(
      xintercept = separator_x,
      linewidth = 0.4,
      color = "black"
    ) +
    facet_grid(
      site ~ .,
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    coord_cartesian(
      xlim = c(0.2, scale_burden_x(max_burden) + 0.6),
      clip = "on"
    ) +
    scale_x_continuous(
      breaks = x_breaks,
      labels = x_labels,
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      breaks = NULL,
      labels = NULL,
      expand = expansion(add = c(0.45, 0.45))
    ) +
    scale_size_continuous(
      range = BUBBLE_SIZE_RANGE,
      breaks = pretty(bubble_plot_df$mean_type_burden, n = 4),
      name = "Mean type burden\nper MAG"
    ) +
    scale_color_gradientn(
      colours = yellow_cols,
      values = rescale(c(0, 25, 50, 75, 100)),
      limits = c(0, 100),
      oob = squish,
      name = "% MAGs"
    ) +
    scale_fill_manual(
      values = type_cols_used,
      name = "Burden type",
      drop = FALSE
    ) +
    guides(
      color = guide_colorbar(
        barwidth = unit(20, "mm"),
        barheight = unit(3, "mm"),
        order = 1
      ),
      size = guide_legend(
        override.aes = list(alpha = 0.85, color = "#666666"),
        keywidth = unit(3.4, "mm"),
        keyheight = unit(3.4, "mm"),
        order = 2
      ),
      fill = guide_legend(
        nrow = 2,
        byrow = TRUE,
        keywidth = unit(3.8, "mm"),
        keyheight = unit(3.0, "mm"),
        order = 3
      )
    ) +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = 8.6) +
    theme(
      panel.grid.major = element_line(color = "#ededed", linewidth = 0.28),
      panel.grid.minor = element_blank(),

      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.55),

      strip.background = element_rect(fill = "white", color = "black", linewidth = 0.55),
      strip.text.y.left = element_text(angle = 90, size = 8.0, face = "plain"),
      strip.placement = "outside",

      axis.title = element_blank(),
      axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.7),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),

      legend.position = "bottom",
      legend.direction = "horizontal",
      legend.title = element_text(size = 7.5, face = "plain"),
      legend.text = element_text(size = 6.8),
      legend.key.height = unit(3.0, "mm"),
      legend.key.width = unit(4.2, "mm"),
      legend.margin = margin(0, 0, 0, 0),
      legend.box.spacing = unit(0.02, "lines"),

      panel.spacing.y = unit(0.10, "lines"),
      panel.spacing.x = unit(0.08, "lines"),

      plot.margin = margin(2, 2, 2, 2, unit = "mm")
    )

  n_rows_total <- nrow(species_df)
  fig_height <- (BASE_HEIGHT + ROW_HEIGHT * n_rows_total) * HEIGHT_FACTOR
  fig_height <- max(fig_height, MIN_HEIGHT)
  if (render) ggplot2::ggsave(file.path(out_dir, "Fig3e.png"), p, width = FIG_WIDTH, height = fig_height, dpi = 300, bg = "white")
  if (render) ggplot2::ggsave(file.path(out_dir, "Fig3e.pdf"), p, width = FIG_WIDTH, height = fig_height, device = grDevices::cairo_pdf, bg = "white")
  readr::write_csv(bubble_plot_df, file.path(out_dir, "Fig3e_bubble_plot_data.csv"))
  readr::write_csv(burden_df, file.path(out_dir, "Fig3e_type_burden_plot_data.csv"))
  readr::write_csv(label_df, file.path(out_dir, "Fig3e_species_labels.csv"))
  invisible(list(plot = p, bubble_data = bubble_plot_df, burden_data = burden_df, bar_data = burden_rect_df, labels = label_df))
}

fig3_usage <- function() {
  cat(paste0(
    'Figure 3a-e\n',
    'Usage: Rscript Fig3_abcde.R --data-dir DATA --output-dir RESULTS [--panels a,b,c,d,e]\n',
    '       [--metadata-file CANONICAL_CSV_OR_XLSX]\n',
    '       [--mag-presence-file TSV --mag-registry-file TSV] [--bubble-premerged true|false]\n',
    'Defaults: DATA=data, RESULTS=results; c uses embedded metadata unless DATA/Fig3_sample_metadata.csv exists.\n',
    'Panel e requires MAG IDs and the unchanged registry, or an explicitly premerged union table.\n',
    'Example: Rscript Fig3_abcde.R --data-dir ../local_inputs --output-dir ../local_results --panels c\n',
    'No data are downloaded or imputed. Input files are never overwritten.\n'))
}

fig3_parse_args <- function(args) {
  opts <- list(data_dir='data',output_dir='results',panels=letters[1:5],metadata_file=NULL,
    mag_presence_file=NULL,mag_registry_file=NULL,bubble_premerged=FALSE,help=FALSE)
  if(any(args %in% c('--help','-h'))) {opts$help <- TRUE;return(opts)}
  if(length(args)%%2L!=0L) stop('Every option requires one value. Use --help.')
  if(!length(args)) return(opts)
  mapping <- c('--data-dir'='data_dir','--output-dir'='output_dir','--panels'='panels','--metadata-file'='metadata_file',
    '--mag-presence-file'='mag_presence_file','--mag-registry-file'='mag_registry_file','--bubble-premerged'='bubble_premerged')
  seen <- character()
  for(i in seq.int(1L,length(args),by=2L)) {
    key <- args[[i]]; value <- args[[i+1L]]
    if(!key %in% names(mapping)) stop('Unknown option: ',key)
    if(key %in% seen) stop('Repeated option: ',key)
    if(!nzchar(value)||startsWith(value,'--')) stop('Missing option value for ',key)
    seen <- c(seen,key)
    opts[[unname(mapping[[key]])]] <- value
  }
  if(is.character(opts$panels)&&length(opts$panels)==1L) opts$panels <- trimws(strsplit(opts$panels,',',fixed=TRUE)[[1]])
  if(!length(opts$panels)||any(!opts$panels %in% letters[1:5])||anyDuplicated(opts$panels)) stop('Panels must be unique letters a-e.')
  if(is.character(opts$bubble_premerged)) {
    if(!tolower(opts$bubble_premerged) %in% c('true','false')) stop('--bubble-premerged requires true or false.')
    opts$bubble_premerged <- tolower(opts$bubble_premerged)=='true'
  }
  opts
}

fig3_input_names <- function() list(
  a=c('Fig3a_species_total_counts.tsv','Fig3a_species_by_sample_counts.tsv','Fig3a_sample_species_contig_totals.tsv',
      'Fig3a_species_by_ARG_type_relative.tsv','Fig3a_sample_metadata.tsv'),
  b='02_Fig3b_type_context_setting_counts.tsv',
  c='08_Fig3c_sample_fractions.tsv',
  d=c('Fig3d_MAG_total_burden.tsv','Fig3d_species_burden_summary.tsv'),
  e=paste0('ESKAPEE_ARGTYPE_bubble_burdenStacked_top6Others_',
           c('bubble_pct','site_species_type_mean_burden','species_used_by_site','argtypes_used'),'.tsv'))

fig3_main <- function(args=commandArgs(trailingOnly=TRUE)) {
  opts <- fig3_parse_args(args)
  if(opts$help) {fig3_usage();return(invisible(NULL))}
  data_dir <- normalizePath(opts$data_dir,winslash='/',mustWork=TRUE)
  out_dir <- normalizePath(opts$output_dir,winslash='/',mustWork=FALSE)
  if(tolower(data_dir)==tolower(out_dir)) stop('The output directory must differ from the input directory.')
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  metadata <- opts$metadata_file
  if(is.null(metadata) && file.exists(file.path(data_dir,'Fig3_sample_metadata.csv'))) metadata <- file.path(data_dir,'Fig3_sample_metadata.csv')
  presence_file <- opts$mag_presence_file
  registry_file <- opts$mag_registry_file
  if(is.null(presence_file) && file.exists(file.path(data_dir,'Fig3e_MAG_ARG_presence.tsv'))) presence_file <- file.path(data_dir,'Fig3e_MAG_ARG_presence.tsv')
  if(is.null(registry_file) && file.exists(file.path(data_dir,'Fig3e_MAG_registry.tsv'))) registry_file <- file.path(data_dir,'Fig3e_MAG_registry.tsv')
  input_names <- fig3_input_names()
  inventory <- list(); status <- list(); results <- list()
  for(panel in opts$panels) {
    cat('\nRunning Figure 3',panel,'...\n',sep='')
    files <- file.path(data_dir,input_names[[panel]])
    if(panel=='c') files <- c(files,metadata)
    if(panel=='e') files <- c(files,presence_file,registry_file)
    started <- Sys.time()
    error_text <- ''
    value <- tryCatch({
      missing_files <- files[!file.exists(files)]
      if(length(missing_files)) stop('Missing input file(s): ',paste(basename(missing_files),collapse=', '))
      inventory[[panel]] <- data.frame(Panel=panel,File=basename(files),MD5=unname(tools::md5sum(files)))
      panel_out <- file.path(out_dir,paste0('panel_',panel))
      switch(panel,
        a=run_fig3a(data_dir,panel_out),
        b=run_fig3b(data_dir,panel_out),
        c=run_fig3c_lmm(files[[1]],metadata,panel_out),
        d=run_fig3d(data_dir,panel_out),
        e=run_fig3e(data_dir,panel_out,mag_presence_file=presence_file,
                   mag_registry_file=registry_file,bubble_premerged=opts$bubble_premerged))
    },error=function(e) {error_text <<- conditionMessage(e);NULL})
    results[[panel]] <- value
    status[[panel]] <- data.frame(Panel=panel,Success=!nzchar(error_text),
      Elapsed_seconds=as.numeric(difftime(Sys.time(),started,units='secs')),Error=error_text)
    utils::write.csv(do.call(rbind,status),file.path(out_dir,'run_status.csv'),row.names=FALSE)
    if(length(inventory)) utils::write.csv(do.call(rbind,inventory),file.path(out_dir,'input_file_hashes.csv'),row.names=FALSE)
    if(nzchar(error_text)) message('Figure 3',panel,' failed: ',error_text)
  }
  writeLines(capture.output(sessionInfo()),file.path(out_dir,'sessionInfo.txt'))
  state <- do.call(rbind,status)
  if(any(!state$Success)) stop('Panel failure(s): ',paste(state$Panel[!state$Success],collapse=', '),'. See run_status.csv. Successful panels were retained.')
  cat('\nCompleted panels: ',paste(opts$panels,collapse=', '),'\n',sep='')
  invisible(results)
}

if(sys.nframe()==0L) fig3_main()
