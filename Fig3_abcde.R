# Figure 3a-e. 
# 2026-09-22: panels b/c use author-confirmed ARG annotation-record units.
# Fig. 3b has no axis titles; counts, proportions and Fig. 3c model are unchanged.
# -------- fig3_ab.R --------

.fig3ab_require <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

run_fig3a <- function(data_dir, out_dir, render = TRUE) {
  .fig3ab_require(c("data.table", "dplyr", "tibble", "ComplexHeatmap", "circlize", "RColorBrewer"))
  # ========================================================= ----

  
  #rm(list = ls())
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tibble)
    library(ComplexHeatmap)
    library(circlize)
    library(grid)
    library(RColorBrewer)
  })
  
  # =========================================================
  # 0. PATHS (modify if needed) ----
  root_dir <- normalizePath(data_dir, winslash = "/", mustWork = TRUE)
  out_dir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  species_total_fp <- file.path(root_dir, "Fig3a_species_total_counts.tsv")
  mat_rel_fp       <- file.path(root_dir, "Fig3a_species_by_sample_relative.tsv")
  arg_comp_fp      <- file.path(root_dir, "Fig3a_species_by_ARG_type_relative.tsv")
  meta_fp          <- file.path(root_dir, "Fig3a_sample_metadata.tsv")
  
  cat("Checking input files...\n")
  print(data.frame(
    file = c("Fig3a_species_total_counts.tsv",
             "Fig3a_species_by_sample_relative.tsv",
             "Fig3a_species_by_ARG_type_relative.tsv",
             "Fig3a_sample_metadata.tsv"),
    exists = c(file.exists(species_total_fp),
               file.exists(mat_rel_fp),
               file.exists(arg_comp_fp),
               file.exists(meta_fp))
  ))
  stopifnot(file.exists(species_total_fp),
            file.exists(mat_rel_fp),
            file.exists(arg_comp_fp),
            file.exists(meta_fp))
  
  # =========================================================
  # 1. READ DATA ----
  species_total <- fread(species_total_fp)
  mat_rel       <- fread(mat_rel_fp)
  arg_comp      <- fread(arg_comp_fp)
  meta_df       <- fread(meta_fp)
  
  # heatmap matrix
  mat_rel_df <- as.data.frame(mat_rel)
  rownames(mat_rel_df) <- mat_rel_df[[1]]
  mat_rel_df[[1]] <- NULL
  mat_rel_df <- as.matrix(mat_rel_df)
  
  # ARG type composition matrix
  arg_comp_df <- as.data.frame(arg_comp)
  rownames(arg_comp_df) <- arg_comp_df[[1]]
  arg_comp_df[[1]] <- NULL
  arg_comp_df <- as.matrix(arg_comp_df)
  
  meta_df <- as.data.frame(meta_df)
  meta_df$sample <- as.character(meta_df$sample)
  
  # =========================================================
  # 2. DETECT METADATA COLUMNS ----
  site_col <- intersect(colnames(meta_df), c("Site", "site", "Sample_Type", "sample_type", "Type", "type"))[1]
  city_col <- intersect(colnames(meta_df), c("City", "city", "City_Code", "city_code"))[1]
  
  if (is.na(site_col)) stop("❌ metadata 中没找到 sample type 列")
  if (is.na(city_col)) stop("❌ metadata 中没找到 city 列")
  
  cat("Detected site column:", site_col, "\n")
  cat("Detected city column:", city_col, "\n")
  
  # =========================================================
  # 3. PALETTES ----
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
  
  # =========================================================
  # 4. CLEAN METADATA & UNIFY NAMES ----
  meta_df[[site_col]] <- trimws(as.character(meta_df[[site_col]]))
  meta_df[[city_col]] <- trimws(as.character(meta_df[[city_col]]))
  meta_df$sample      <- trimws(as.character(meta_df$sample))
  
  # unify site names
  meta_df[[site_col]][meta_df[[site_col]] %in% c("农集贸市场", "WetMarket", "Wet_market")] <- "Wet market"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("社区")] <- "Community"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("医院")] <- "Hospital"
  meta_df[[site_col]][meta_df[[site_col]] %in% c("污水处理厂", "Sewage treatment plant")] <- "WWTP"
  
  # city map to abbreviations
  city_map <- c(
    "BJ" = "BJ", "Beijing" = "BJ", "北京" = "BJ",
    "CQ" = "CQ", "Chongqing" = "CQ", "重庆" = "CQ",
    "GY" = "GY", "Guiyang" = "GY", "贵阳" = "GY",
    "GZ" = "GZ", "Guangzhou" = "GZ", "广州" = "GZ",
    "HF" = "HF", "Hefei" = "HF", "合肥" = "HF",
    "HK" = "HK", "Haikou" = "HK", "海口" = "HK",
    "HR" = "HR", "Harbin" = "HR", "Haerbin" = "HR", "哈尔滨" = "HR",
    "NJ" = "NJ", "Nanjing" = "NJ", "南京" = "NJ",
    "SZ" = "SZ", "Shenzhen" = "SZ", "深圳" = "SZ",
    "XA" = "XA", "Xian" = "XA", "Xi'an" = "XA", "西安" = "XA",
    "XM" = "XM", "Xiamen" = "XM", "厦门" = "XM",
    "ZB" = "ZB", "Zibo" = "ZB", "淄博" = "ZB"
  )
  meta_df[[city_col]] <- ifelse(meta_df[[city_col]] %in% names(city_map),
                                unname(city_map[meta_df[[city_col]]]),
                                meta_df[[city_col]])
  
  site_order <- c("Hospital", "WWTP", "Community", "Wet market")
  city_order <- c("BJ", "CQ", "GY", "GZ", "HF", "HK", "HR", "NJ", "SZ", "XA", "XM", "ZB")
  
  # =========================================================
  # 5. ORDER SPECIES BY TOTAL COUNT ----
  species_order <- species_total %>%
    arrange(desc(total_count)) %>%
    pull(species)
  species_order <- intersect(species_order, rownames(mat_rel_df))
  mat_rel_df    <- mat_rel_df[species_order, , drop = FALSE]
  arg_comp_df   <- arg_comp_df[species_order, , drop = FALSE]
  
  # =========================================================
  # 6. CLEAN ARG TYPE NAMES AND MERGE DUPLICATES ----
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
  
  arg_comp_df2 <- sapply(unique(arg_names_clean), function(tp) {
    rowSums(arg_comp_df[, arg_names_clean == tp, drop = FALSE], na.rm = TRUE)
  })
  arg_comp_df2 <- as.matrix(arg_comp_df2)
  rownames(arg_comp_df2) <- rownames(arg_comp_df)
  
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
  
  # =========================================================
  # 7. SAMPLE ORDER AND AGGREGATION TO GROUPS (city | site) ----
  sample_order_df <- meta_df %>%
    filter(sample %in% colnames(mat_rel_df)) %>%
    filter(!is.na(.data[[site_col]]), .data[[site_col]] != "",
           !is.na(.data[[city_col]]), .data[[city_col]] != "") %>%
    filter(.data[[site_col]] %in% site_order,
           .data[[city_col]] %in% city_order) %>%
    mutate(city_factor = factor(.data[[city_col]], levels = city_order),
           site_factor = factor(.data[[site_col]], levels = site_order)) %>%
    arrange(site_factor, city_factor, sample)
  
  sample_order <- sample_order_df$sample
  sample_order <- intersect(sample_order, colnames(mat_rel_df))
  mat_rel_df <- mat_rel_df[, sample_order, drop = FALSE]
  
  meta_plot <- sample_order_df %>%
    filter(sample %in% sample_order) %>%
    distinct(sample, .keep_all = TRUE)
  meta_plot <- meta_plot[match(sample_order, meta_plot$sample), , drop = FALSE]
  
  cat("Samples kept after filtering:", nrow(meta_plot), "\n")
  
  # Create groups (city | site) and aggregate by mean
  meta_plot$group <- paste(meta_plot[[city_col]], meta_plot[[site_col]], sep = " | ")
  group_order_df <- meta_plot %>%
    distinct(group, .keep_all = TRUE) %>%
    mutate(city_factor = factor(.data[[city_col]], levels = city_order),
           site_factor = factor(.data[[site_col]], levels = site_order)) %>%
    arrange(site_factor, city_factor)
  group_order <- group_order_df$group
  
  mat_rel_group <- sapply(group_order, function(g) {
    cols <- meta_plot$sample[meta_plot$group == g]
    cols <- intersect(cols, colnames(mat_rel_df))
    if (length(cols) == 0) return(rep(NA_real_, nrow(mat_rel_df)))
    if (length(cols) == 1) return(as.numeric(mat_rel_df[, cols]))
    rowMeans(mat_rel_df[, cols, drop = FALSE], na.rm = TRUE)
  })
  mat_rel_group <- as.matrix(mat_rel_group)
  rownames(mat_rel_group) <- rownames(mat_rel_df)
  colnames(mat_rel_group) <- group_order
  
  meta_group <- group_order_df %>%
    dplyr::select(group, all_of(city_col), all_of(site_col))
  meta_group[[site_col]] <- factor(meta_group[[site_col]], levels = site_order)
  meta_group[[city_col]] <- factor(meta_group[[city_col]], levels = city_order)
  
  mat_rel_df <- mat_rel_group   # now matrix of groups
  
  # =========================================================
  # 8. LEFT COUNTS BAR (增加宽度，用于放置物种名) ###***
  left_counts <- species_total %>%
    filter(species %in% species_order) %>%
    distinct(species, total_count) %>%
    as.data.frame()
  left_counts <- left_counts[match(species_order, left_counts$species), , drop = FALSE]
  left_bar <- left_counts$total_count
  names(left_bar) <- left_counts$species
  left_bar_rev <- -left_bar
  max_count <- max(left_bar, na.rm = TRUE)
  
  # 自定义 breaks (可修改)
  count_breaks <- c(0, 30000, 60000)
  count_breaks <- count_breaks[count_breaks <= max_count]
  
  # 高亮物种（仅用于颜色，不再用于行名显示）
  highlight_species <- c("Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
                         "Enterococcus faecium", "Acinetobacter baumannii", "Staphylococcus aureus",
                         "Klebsiella quasipneumoniae", "Enterobacter cloacae")
  # 注意：row_label_cols 仍用于后续手动添加的行名颜色 ###***
  row_label_cols <- ifelse(species_order %in% highlight_species, "#D62728", "black")
  
  # 左侧条形图宽度增加到 5 cm，以便容纳物种名 ###***
  left_anno <- rowAnnotation(
    "Number of counts" = anno_barplot(
      left_bar_rev, baseline = 0,
      gp = gpar(fill = "#B7D7E8", col = "#9EC5DA"),
      border = FALSE, bar_width = 0.82,
      width = unit(5, "cm"),                # 宽度增加到 5 cm ###***
      add_numbers = FALSE,
      axis_param = list(side = "top", at = -count_breaks, labels = count_breaks,
                        labels_rot = 0, gp = gpar(fontsize = 7))
    ),
    annotation_name_side = "top",
    annotation_name_rot = 0,
    annotation_name_gp = gpar(fontsize = 10, font = 1)
  )
  
  # =========================================================
  # 9. RIGHT ARG TYPE ANNOTATION (移除原有轴线，事后添加内嵌刻度尺) ###***
  right_anno <- rowAnnotation(
    `ARG types` = anno_barplot(
      arg_comp_df, beside = FALSE,
      gp = gpar(fill = arg_type_cols[colnames(arg_comp_df)], col = NA),
      border = FALSE, bar_width = 1, width = unit(2.5, "cm"),
      axis = FALSE
    ),
    show_annotation_name = F
  )
  
  # =========================================================
  # 10. TOP ANNOTATION (Sample type and City) ----
  top_ha <- HeatmapAnnotation(
    `Sample type` = meta_group[[site_col]],
    City = meta_group[[city_col]],
    col = list(`Sample type` = site_palette, City = city_palette),
    annotation_legend_param = list(title_gp = gpar(fontsize = 10, font=1)),
    annotation_name_gp = gpar(fontsize = 10, font = 1),
    simple_anno_size = unit(4.2, "mm"), gap = unit(1.2, "mm"),
    show_annotation_name = TRUE
  )
  
  # =========================================================
  # 11. BUBBLE PLOT VIA CELL_FUN (replaces heatmap) ----
  bubble_col_fun <- colorRamp2(c(0, 0.01, 0.02, 0.05),
                               c("#F7F7F7", "#F6D6B8", "#E98B4A", "#B30000"))
  
  cell_fun_bubble <- function(j, i, x, y, width, height, fill) {
    v <- mat_rel_df[i, j]
    if (is.na(v) || v == 0) return()
    max_radius <- min(unit.c(width, height)) * 1.0 ###***0.5
    r <- sqrt(v / 0.10) * max_radius
    grid.circle(x = x, y = y, r = r,
                gp = gpar(fill = bubble_col_fun(v), col = NA, font = 1))
  }
  
  column_split_site <- factor(meta_group[[site_col]], levels = site_order)
  
  # 合并图例（离散的点，不同颜色和大小）
  size_vals <- c(0.01, 0.02, 0.05)   # 相对丰度值
  # 计算半径（与气泡图中 cell_fun 的比例一致，max_radius 系数为 0.5）
  ref_width <- unit(1, "cm")
  ref_height <- unit(1, "cm")
  max_radius_ref <- min(ref_width, ref_height) * 0.5
  radii <- sqrt(size_vals / 0.10) * max_radius_ref
  
  # 颜色映射函数（与气泡图保持一致）
  bubble_col_fun <- colorRamp2(c(0, 0.01, 0.02, 0.05),
                               c("#F7F7F7", "#F6D6B8", "#E98B4A", "#B30000"))
  # 为每个 size_vals 获取颜色
  point_colors <- bubble_col_fun(size_vals)
  
  combined_legend <- Legend(
    title = "Relative\nabundance",
    type = "points",
    at = size_vals,
    labels = c("1%", "2%", "5%"),
    legend_gp = gpar(col = point_colors, fill=NULL),  # fill 为内部色，col 为边框（黑色）
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
  # 关键修改：隐藏原有行名，因为我们会在左侧条形图中手动添加 ###***
  ht_bubble <- Heatmap(
    mat_rel_df,
    name = "Relative\nabundance",
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
    show_row_names = FALSE,                # 隐藏原有行名 ###***
    # row_names_gp 不再需要，因为手动添加 ###***
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
  
  # =========================================================
  # 12. DRAW AND SAVE (手动添加行名到左侧条形图内部，左对齐) ###***
  # Extra canvas height keeps all 28 original ARG-class legend entries visible.
  draw_bubble <- function(ht_obj, pdf_fp, png_fp, width = 11, height = 5.2) {
    # 辅助函数：将物种名添加到左侧条形图内部，固定左对齐 ###***
    add_row_names_to_left_bar <- function() {
      decorate_annotation("Number of counts", {
        n_rows <- length(species_order)
        # 固定 x 坐标：左侧偏移 0.02 npc（可根据需要调整，确保在条形内且不重叠）
        x_pos_fixed <- unit(0.02, "npc")
        for (i in seq_along(species_order)) {
          species_name <- species_order[i]
          # y 坐标：顶部行 (i=1) 对应 y = n_rows，底部行 (i=n_rows) 对应 y = 1
          y_pos <- unit(n_rows - i + 1, "native")
          text_col <- ifelse(species_name %in% highlight_species, "#D62728", "black")
          grid.text(species_name, x = x_pos_fixed, y = y_pos,
                    just = "left", 
                    gp = gpar(fontsize = 10, fontface = "italic", col = text_col))
        }
      })
    }
    
    # PDF
    pdf(pdf_fp, width = width, height = height)
    draw(ht_obj,
         heatmap_legend_side = "right",
         annotation_legend_side = "right",
         merge_legends = FALSE,
         annotation_legend_list = all_legends)
    add_row_names_to_left_bar()
    dev.off()
    
    # PNG
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
  
  cat("✅ Bubble plot generated!\n")
  cat("Output dir:", out_dir, "\n")
  cat("Generated files:\n")
  print(list.files(out_dir, full.names = TRUE))
  
  
  data.table::fwrite(meta_group, file.path(out_dir, "city_setting_groups.tsv"), sep = "\t")
  data.table::fwrite(data.frame(species = rownames(mat_rel_df), mat_rel_df, check.names = FALSE),
                    file.path(out_dir, "city_setting_relative_abundance.tsv"), sep = "\t")
  audit_warnings <- data.frame(
    issue = c("retained_source_scope", "retained_denominator", "retained_bubble_geometry"),
    detail = c(
      paste0("The supplied study total_count and ARG-class composition retain the original 2652-sample scope; the study matrix has 2607 available cohort samples. This run retains ", nrow(meta_plot), " samples. No counts were recalculated."),
      "Input relative abundances are preserved. Source normalization includes higher taxonomic assignments and is not species-only; displaying 20 species does not change that denominator.",
      "Author's bubble radius and legend scaling are preserved. Some bubbles can extend beyond the smaller cell dimension; scientific values are unchanged."
    )
  )
  data.table::fwrite(audit_warnings, file.path(out_dir, "audit_warnings.tsv"), sep = "\t")
  invisible(list(group_matrix = mat_rel_df, group_metadata = meta_group,
                 total_counts = left_bar, arg_composition = arg_comp_df,
                 sample_count = nrow(meta_plot), heatmap = ht_bubble,
                 audit_warnings = audit_warnings))
  
}

# Fig. 3b,c count ARG annotations, one final best-hit annotation per ORF.
# Different ARG-bearing ORFs on the same contig contribute separate records.
# Legacy "*contigs" input headers are retained as accepted aliases only.
# These aggregate tables cannot re-check ORF-level deduplication.
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

run_fig3b <- function(data_dir, out_dir, render = TRUE) {
  .fig3ab_require(c("dplyr", "readr", "ggplot2", "ggalluvial", "forcats", "stringr"))
  
  suppressPackageStartupMessages({
    library(dplyr)
    library(readr)
    library(ggplot2)
    library(ggalluvial)
    library(forcats)
    library(stringr)
  })
  
  infile <- file.path(data_dir, "Fig3b_FigS7_figure_input_type.tsv")
  outdir <- normalizePath(out_dir, winslash = "/", mustWork = FALSE)
  dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
  
  site_order <- c("Hospital", "WWTP", "Community", "Wet market")
  carrier_order <- c( "chromosome","plasmid", "virus")
  
  type_palette <- c(
    "Aminoglycoside" = "#4E79A7",
    "Antibacterial fatty acid" = "#2F5D62",
    "Bacitracin" = "#C7A76C",
    "beta_lactam" = "#556B2F",
    "Bicyclomycin" = "#7BAE7F",
    "Bleomycin" = "#D97A5C",
    "Chloramphenicol" = "#A23B72",
    "MLS" = "#6FB1A0",
    "Multidrug" = "#5B3F8C",
    "Tetracycline" = "#E8D98F",
    "Polymyxin" = "#75C0C1",
    "Mupirocin" = "#5F9EA0",
    "Others" = "#BDBDBD"
  )
  
  df <- .fig3_annotation_count_columns(
    readr::read_tsv(infile, show_col_types = FALSE),
    c(n_arg_annotations = "n_contigs")
  ) %>%
    mutate(
      sample = as.character(sample),
      type_plot = as.character(type_plot),
      predicted_class = as.character(predicted_class),
      Site = as.character(Site),
      n_arg_annotations = as.numeric(n_arg_annotations)
    ) %>%
    filter(
      !is.na(type_plot), type_plot != "",
      !is.na(predicted_class), predicted_class != "",
      !is.na(Site), Site != "",
      !is.na(n_arg_annotations), n_arg_annotations > 0
    )
  
  df <- df %>%
    mutate(
      predicted_class = str_trim(predicted_class),
      predicted_class = case_when(
        predicted_class %in% c("Plasmid", "plasmid") ~ "plasmid",
        predicted_class %in% c("Chromosome", "chromosome") ~ "chromosome",
        predicted_class %in% c("Virus", "virus") ~ "virus",
        TRUE ~ predicted_class
      ),
      Site = str_trim(Site),
      Site = case_when(
        Site %in% c("Hospital") ~ "Hospital",
        Site %in% c("WWTP") ~ "WWTP",
        Site %in% c("Community") ~ "Community",
        Site %in% c("Wet market", "Wet_market", "WetMarket") ~ "Wet market",
        TRUE ~ Site
      )
    ) %>%
    filter(
      predicted_class %in% carrier_order,
      Site %in% site_order
    )
  
  top9_types <- df %>%
    group_by(type_plot) %>%
    summarise(total_n = sum(n_arg_annotations, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_n)) %>%
    slice_head(n = 9) %>%
    pull(type_plot)
  
  df2 <- df %>%
    mutate(
      type_top = if_else(type_plot %in% top9_types, type_plot, "Others")
    )
  
  plot_dat <- df2 %>%
    group_by(type_top, predicted_class, Site) %>%
    summarise(n = sum(n_arg_annotations, na.rm = TRUE), .groups = "drop") %>%
    filter(n > 0)
  
  type_order <- plot_dat %>%
    group_by(type_top) %>%
    summarise(total_n = sum(n), .groups = "drop") %>%
    arrange(desc(total_n)) %>%
    pull(type_top)
  
  type_order <- c(setdiff(type_order, "Others"), intersect("Others", type_order))
  
  plot_dat <- plot_dat %>%
    mutate(
      type_top = factor(type_top, levels = type_order),
      predicted_class = factor(predicted_class, levels = carrier_order),
      Site = factor(Site, levels = site_order)
    )
  
  readr::write_tsv(plot_dat, file.path(outdir, "type_top9_others_carrier_site_counts.tsv"))
  
  p <- ggplot(
    plot_dat,
    aes(axis1 = type_top, axis2 = predicted_class, axis3 = Site, y = n)
  ) +
    geom_alluvium(
      aes(fill = type_top),
      width = 0.16,
      alpha = 0.7,
      knot.pos = 0.35,
      decreasing = FALSE
    ) +
    geom_stratum(
      width = 0.22,
      fill = "grey95",
      colour = "black",
      linewidth = 0.45
    ) +
    geom_text(
      stat = "stratum",
      aes(label = after_stat(stratum)),
      size = 3,
      family = "sans"
    ) +
    scale_x_discrete(
      limits = c("ARG family", "Carrier", "Site"),
      expand = c(0.06, 0.06)
    ) +
    scale_fill_manual(
      values = type_palette,
      breaks = type_order,
      drop = FALSE
    ) +
    labs(
      title = "Type -> carrier -> site",
      x = NULL,
      y = NULL,
      fill = "ARG family"
    ) +
    theme_bw(base_size = 13) +
    theme(
      panel.grid = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 0.4),
      axis.ticks = element_line(colour = "black", linewidth = 0.4),
      plot.title = element_text(face = "bold", hjust = 0, size = 15),
      axis.title.y = element_blank(),
      axis.text.x = element_text(face = "bold", size = 12),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.height = unit(0.55, "cm"),
      legend.key.width = unit(0.45, "cm"),
      legend.position = "right",
      plot.margin = margin(8, 12, 8, 8)
    )
  
  if (render) ggsave(
    filename = file.path(outdir, "Type_to_carrier_to_site_top9_others.pdf"),
    plot = p,
    width = 13.5,
    height = 8.5,
    units = "in"
  )
  
  if (render) ggsave(
    filename = file.path(outdir, "Type_to_carrier_to_site_top9_others.png"),
    plot = p,
    width = 13.5,
    height = 8.5,
    units = "in",
    dpi = 600
  )
  
  cat("Done.\n")
  cat("Output dir:", outdir, "\n")
  
  mobile_unit_fraction <- sum(df$n_arg_annotations[df$predicted_class %in% c("plasmid", "virus")]) / sum(df$n_arg_annotations)
  audit_summary <- data.frame(
    metric = c("retained_samples", "retained_rows", "ARG_annotations", "mobile_ARG_annotations", "mobile_ARG_annotation_fraction"),
    value = c(dplyr::n_distinct(df$sample), nrow(df), sum(df$n_arg_annotations),
              sum(df$n_arg_annotations[df$predicted_class %in% c("plasmid", "virus")]), mobile_unit_fraction)
  )
  readr::write_tsv(audit_summary, file.path(outdir, "count_unit_summary.tsv"))
  audit_warnings <- data.frame(
    issue = c("counting_unit", "input_scope", "disconnected_preparation_removed"),
    detail = c(
      "Flows sum ARG annotation records: one final best-hit annotation per ORF, with separate ORFs on the same contig counted separately. The legacy n_contigs input field is an alias for n_arg_annotations. The pooled mobile fraction is the plasmid- plus virus-associated ARG annotation count divided by all included ARG annotation counts, not a unique-contig fraction.",
      "The supplied plotting table and the original disconnected preparation table have different sample scopes. This function preserves the plotting input and its original filters; missing/excluded samples are not imputed.",
      "The disconnected preparation branch was omitted. Its saved table excludes all Wet market samples and does not feed this plotted result."
    )
  )
  readr::write_tsv(audit_warnings, file.path(outdir, "audit_warnings.tsv"))
  invisible(list(plot_data = plot_dat, type_palette = type_palette, plot = p,
                 mobile_unit_fraction = mobile_unit_fraction,
                 audit_summary = audit_summary, audit_warnings = audit_warnings))
  
}


# -------- fig3c_lmm.R --------

fig3c_settings <- c('Hospital','WWTP','Community','Wet market')

prepare_fig3c_data <- function(points, metadata) {
  # One row per sample; counts refer to final ARG annotations, not unique contigs.
  points <- .fig3_annotation_count_columns(points, c(
    total_arg_annotations='total_arg_contigs', mobile_arg_annotations='mobile_arg_contigs',
    plasmid_arg_annotations='plasmid_arg_contigs', virus_arg_annotations='virus_arg_contigs',
    chromosome_arg_annotations='chromosome_arg_contigs'))
  required_points <- c('sample','Sample_Type','Sample_Date','total_arg_annotations',
    'mobile_arg_annotations','plasmid_arg_annotations','virus_arg_annotations','chromosome_arg_annotations','mobile_fraction')
  if(!all(required_points %in% names(points))) stop('Mobile-fraction input lacks required fields.')
  points$sample <- trimws(as.character(points$sample))
  if(anyNA(points$sample)||any(!nzchar(points$sample))||anyDuplicated(points$sample)) stop('Missing or duplicate sample ID in fraction input.')
  if(!all(c('Sample','City','Sample_Date') %in% names(metadata))) stop('Canonical metadata lacks required fields.')
  m <- as.data.frame(metadata,stringsAsFactors=FALSE)
  if(!'Setting' %in% names(m)) m$Setting <- m$Sample_Type
  m$Setting <- trimws(as.character(m$Setting))
  m$Setting[m$Setting %in% c('Wet Market','Wet_market','wet market')] <- 'Wet market'
  if(!'PhysicalSite' %in% names(m)) {
    if(!'Sample_Rename' %in% names(m)) stop('Canonical metadata needs PhysicalSite or Sample_Rename.')
    m$PhysicalSite <- sub('_[0-9]{8}$','',as.character(m$Sample_Rename))
  }
  date_text <- as.character(m$Sample_Date)
  dates <- as.Date(ifelse(grepl('^[0-9]{8}$',date_text),
    paste0(substr(date_text,1,4),'-',substr(date_text,5,6),'-',substr(date_text,7,8)),date_text))
  m$Sample_Date <- dates
  if(!'SamplingMonth' %in% names(m)) m$SamplingMonth <- format(dates,'%Y-%m')
  keep <- c('Sample','City','Sample_Date','PhysicalSite','Setting','SamplingMonth')
  m <- m[keep]
  for(key in setdiff(keep,'Sample_Date')) m[[key]] <- trimws(as.character(m[[key]]))
  if(anyNA(m)||any(vapply(m[setdiff(keep,'Sample_Date')],function(x)any(!nzchar(x)),logical(1)))||
     anyDuplicated(m$Sample)||any(!m$Setting %in% fig3c_settings)||
     any(m$SamplingMonth != format(m$Sample_Date,'%Y-%m'))) stop('Incomplete, duplicate or inconsistent canonical metadata.')
  site_map <- unique(m[c('PhysicalSite','City','Setting')])
  if(anyDuplicated(site_map$PhysicalSite)) stop('Physical-site metadata maps to multiple cities or settings.')
  numeric_fields <- c('total_arg_annotations','mobile_arg_annotations','plasmid_arg_annotations','virus_arg_annotations','chromosome_arg_annotations','mobile_fraction')
  for(key in numeric_fields) points[[key]] <- suppressWarnings(as.numeric(points[[key]]))
  if(any(!is.finite(points$total_arg_annotations))||any(points$total_arg_annotations<=0)) stop('Missing or nonpositive denominator: do not impute mobile fraction.')
  if(any(!is.finite(points$mobile_fraction))||any(points$mobile_fraction<0|points$mobile_fraction>1)) stop('Missing or invalid mobile fraction.')
  for(key in setdiff(numeric_fields,'mobile_fraction')) {
    if(any(!is.finite(points[[key]]))||any(points[[key]]<0)||any(abs(points[[key]]-round(points[[key]]))>1e-8)) stop('Invalid ARG annotation count: ',key)
  }
  if(any(points$mobile_arg_annotations != points$plasmid_arg_annotations+points$virus_arg_annotations)||
     any(points$total_arg_annotations != points$mobile_arg_annotations+points$chromosome_arg_annotations)) stop('ARG annotation counts by compartment do not sum to denominator.')
  if(any(abs(points$mobile_fraction-points$mobile_arg_annotations/points$total_arg_annotations)>1e-12)) stop('Mobile-fraction ratio differs from input counts.')
  j <- match(m$Sample,points$sample)
  missing <- m[is.na(j),,drop=FALSE]
  outside <- points[!points$sample %in% m$Sample,,drop=FALSE]
  pp <- points[j[!is.na(j)],,drop=FALSE]
  # Preserve input covariates as audit columns; inference uses canonical metadata.
  d <- m[!is.na(j),,drop=FALSE]
  for(key in setdiff(names(pp),'sample')) {
    d[[if(key %in% c('City','Sample_Date','Sample_Type')) paste0('Input_',key) else key]] <- pp[[key]]
  }
  input_setting <- as.character(d$Input_Sample_Type)
  input_setting[input_setting=='Wet Market'] <- 'Wet market'
  disagreed <- input_setting != d$Setting | as.character(d$Input_Sample_Date) != format(d$Sample_Date,'%Y%m%d')
  disagreements <- d[disagreed,c('Sample','Input_Sample_Type','Setting','Input_Sample_Date','Sample_Date'),drop=FALSE]
  d$Setting <- factor(d$Setting,levels=fig3c_settings)
  for(key in c('City','SamplingMonth','PhysicalSite')) d[[key]] <- factor(d[[key]])
  list(data=d,missing_samples=missing,outside_samples=outside,metadata_disagreements=disagreements,metadata=m)
}

fit_fig3c_lmm <- function(d) {
  if(length(unique(d$Setting))!=4L) stop('All four settings are required for six contrasts.')
  required <- c('mobile_fraction','Setting','City','SamplingMonth','PhysicalSite')
  if(anyNA(d[required])) stop('Incomplete model data; do not silently omit rows.')
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
    den <- density(d$mobile_fraction[d$Setting==fig3c_settings[i]],n=512)
    data.frame(Setting=fig3c_settings[i],x=c(i+.45*den$y/max(den$y),rep(i,length(den$x))),y=c(den$x,rev(den$x)))
  }))
  label <- merge(descriptive,letters,by='Setting',sort=FALSE)
  label$x_id <- match(label$Setting,fig3c_settings)
  label$mean_label <- sprintf('%.2f',label$Mean_mobile_fraction)
  span <- diff(range(d$mobile_fraction));top <- max(d$mobile_fraction)+.16*span
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
    ggplot2::coord_cartesian(ylim=c(min(d$mobile_fraction),top+.1*span),clip='off')+
    ggplot2::labs(x=NULL,y='Mobile ARG fraction')+ggplot2::theme_bw(base_size=12)+
    ggplot2::theme(panel.grid=ggplot2::element_blank(),panel.border=ggplot2::element_rect(colour='black',fill=NA,linewidth=.65),
      axis.line=ggplot2::element_blank(),axis.title.y=ggplot2::element_text(size=10,face='plain'),
      axis.text.x=ggplot2::element_text(colour='black',margin=ggplot2::margin(t=3)),axis.text.y=ggplot2::element_text(colour='black'),
      axis.ticks=ggplot2::element_line(colour='black',linewidth=.42),axis.ticks.length=grid::unit(.13,'cm'),
      legend.position='none',plot.margin=ggplot2::margin(5,3,3,5))
}

run_fig3c_lmm <- function(data_file,metadata_file,out_dir) {
  required <- c('lme4','multcompView','ggplot2')
  if(grepl('\\.xlsx?$',metadata_file,ignore.case=TRUE)) required <- c(required,'readxl')
  missing_packages <- required[!vapply(required,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing_packages)) stop('Missing packages: ',paste(missing_packages,collapse=', '))
  # Mark UTF-8 input strings explicitly: native-encoding strings saved under the C
  # locale can fail to deserialize under UTF-8 when sample IDs contain Chinese text.
  points <- utils::read.delim(data_file,check.names=FALSE,stringsAsFactors=FALSE,encoding='UTF-8')
  meta <- if(grepl('\\.xlsx?$',metadata_file,ignore.case=TRUE)) as.data.frame(readxl::read_excel(metadata_file)) else utils::read.csv(metadata_file,check.names=FALSE,stringsAsFactors=FALSE,encoding='UTF-8')
  prepared <- prepare_fig3c_data(points,meta)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  write_table <- function(x,name) utils::write.csv(x,file.path(out_dir,name),row.names=FALSE,na='')
  write_table(prepared$missing_samples,'canonical_samples_without_fraction.csv')
  write_table(prepared$outside_samples,'fraction_samples_outside_canonical_cohort.csv')
  write_table(prepared$metadata_disagreements,'metadata_disagreements_resolved_to_canonical.csv')
  write_table(prepared$data,'sample_model_input.csv')
  d <- prepared$data
  coverage <- data.frame(Canonical_samples=nrow(prepared$metadata),Input_fraction_samples=nrow(points),
    Analyzed_samples=nrow(d),Canonical_samples_without_fraction=nrow(prepared$missing_samples),
    Input_samples_outside_canonical=nrow(prepared$outside_samples),Metadata_disagreement_rows=nrow(prepared$metadata_disagreements),
    Policy='Exact Sample ID join; canonical setting/city/date/site; no outcome or denominator imputation; observed 0/1 retained')
  write_table(coverage,'cohort_coverage.csv')
  result <- fit_fig3c_lmm(d)
  descriptive <- do.call(rbind,lapply(fig3c_settings,function(setting) {
    q <- d[d$Setting==setting,,drop=FALSE]
    data.frame(Setting=setting,Samples=nrow(q),Physical_sites=length(unique(q$PhysicalSite)),
      Mean_mobile_fraction=mean(q$mobile_fraction),Median_mobile_fraction=median(q$mobile_fraction),
      SD_mobile_fraction=sd(q$mobile_fraction),Min=min(q$mobile_fraction),Max=max(q$mobile_fraction),
      Total_ARG_annotations=sum(q$total_arg_annotations),Mobile_ARG_annotations=sum(q$mobile_arg_annotations),
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


# -------- fig3_de.R --------

.fig3de_packages <- function(packages) {
  missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "), call. = FALSE)
  suppressPackageStartupMessages(invisible(lapply(packages, function(package) {
    library(package, character.only = TRUE)
  })))
}

.fig3de_remap_letters <- function(raw, reference = "a") {
  # Rename individual symbols bijectively, preserving every shared-letter relation.
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
  # FSA/dunn.test versions may fail for a single pair under method="bh".
  # This is the same two-sided Dunn rank comparison, including tie correction.
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

run_fig3d <- function(data_dir, out_dir) {
  .fig3de_packages(c("data.table", "dplyr", "ggplot2", "stringr", "FSA", "multcompView", "readr", "scales"))
  data_dir <- normalizePath(data_dir, mustWork = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  .figure_results_dir <- normalizePath(out_dir, mustWork = TRUE)
  fig3d_dunn_results <- list()
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
        Letters = "a"
      ))
    }
  
    tryCatch({
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
          non_h_n_letters = dplyr::n_distinct(Letters[site != "Hospital" & Letters != ""]),
          non_h_has_diff = non_h_n_letters > 1L,
          diff_flag = dplyr::case_when(
            site == "Hospital" ~ "Same",
            Letters == "" ~ "Same",
            non_h_has_diff ~ "Different",
            Letters != "a" ~ "Different",
            TRUE ~ "Same"
          )
        )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(species = factor(species, levels = rev(fig3d_species_levels)))
  
  readr::write_csv(
    fig3d_plot_df %>%
      dplyr::select(
        species, site, prop_pathogen_mags, mean_burden,
        Letters_raw, Letters, non_h_n_letters, non_h_has_diff, diff_flag
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
      values = c("Same" = "grey25", "Different" = "#D73027"),
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
  status <- dplyr::bind_rows(fig3d_test_status)
  readr::write_csv(contrasts, file.path(.figure_results_dir, "Fig3d_Dunn_BH_contrasts.csv"))
  readr::write_csv(status, file.path(.figure_results_dir, "Fig3d_Dunn_test_status.csv"))
  readr::write_csv(fig3d_plot_df, file.path(.figure_results_dir, "Fig3d_plot_data.csv"))
  invisible(list(plot = fig3d_plot, plot_data = fig3d_plot_df, contrasts = contrasts, test_status = status))
}

run_fig3e <- function(data_dir, out_dir) {
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
    ifelse(x %in% c("beta_lactam", "β-lactam", "beta-lactam"), "beta-lactam", x)
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
  
  bubble_df <- bubble_df %>%
    mutate(
      arg_type_plot = if_else(arg_type_plot == "Rifamycin", "Others", arg_type_plot),
      arg_type_plot = if_else(arg_type_plot %in% top_types_input, arg_type_plot, "Others")
    ) %>%
    group_by(site, species, arg_type_plot) %>%
    summarise(
      pct_mags = max(pct_mags, na.rm = TRUE),
      n_type_carrying_mags = if (all(is.na(n_type_carrying_mags))) {
        NA_real_
      } else {
        sum(n_type_carrying_mags, na.rm = TRUE)
      },
      .groups = "drop"
    )
  
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
  ggplot2::ggsave(file.path(out_dir, "Fig3e.png"), p, width = FIG_WIDTH, height = fig_height, dpi = 300, bg = "white")
  ggplot2::ggsave(file.path(out_dir, "Fig3e.pdf"), p, width = FIG_WIDTH, height = fig_height, device = grDevices::cairo_pdf, bg = "white")
  readr::write_csv(bubble_plot_df, file.path(out_dir, "Fig3e_bubble_plot_data.csv"))
  readr::write_csv(burden_df, file.path(out_dir, "Fig3e_type_burden_plot_data.csv"))
  readr::write_csv(label_df, file.path(out_dir, "Fig3e_species_labels.csv"))
  invisible(list(plot = p, bubble_data = bubble_plot_df, burden_data = burden_df, bar_data = burden_rect_df, labels = label_df))
}


# -------- fig3_entrypoint.R --------
# Command-line entrypoint. This file is appended to the a-e modules by the bundler.
fig3_usage <- function() {
  cat(paste0(
    'Figure 3a-e (panel f is deliberately excluded)\n',
    'Usage: Rscript Fig3_abcde.R --data-dir DATA --output-dir RESULTS [--panels a,b,c,d,e]\n',
    '       [--metadata-file CANONICAL_CSV_OR_XLSX]\n',
    'Defaults: DATA=data, RESULTS=results; canonical metadata=DATA/Fig3_sample_metadata.csv\n',
    'Example: Rscript Fig3_abcde.R --data-dir ../local_inputs --output-dir ../local_results --panels c\n',
    'No data are downloaded or imputed. Input files are never overwritten.\n'))
}

fig3_parse_args <- function(args) {
  opts <- list(data_dir='data',output_dir='results',panels=letters[1:5],metadata_file=NULL,help=FALSE)
  if(any(args %in% c('--help','-h'))) {opts$help <- TRUE;return(opts)}
  if(length(args)%%2L!=0L) stop('Every option requires one value. Use --help.')
  if(!length(args)) return(opts)
  mapping <- c('--data-dir'='data_dir','--output-dir'='output_dir','--panels'='panels','--metadata-file'='metadata_file')
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
  if(!length(opts$panels)||any(!opts$panels %in% letters[1:5])||anyDuplicated(opts$panels)) stop('Panels must be unique letters a-e; panel f is not part of this code.')
  opts
}

fig3_input_names <- function() list(
  a=c('Fig3a_species_total_counts.tsv','Fig3a_species_by_sample_relative.tsv',
      'Fig3a_species_by_ARG_type_relative.tsv','Fig3a_sample_metadata.tsv'),
  b='Fig3b_FigS7_figure_input_type.tsv',
  c='Fig3c_mobile_fraction_violin_points.tsv',
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
  metadata <- if(is.null(opts$metadata_file)) file.path(data_dir,'Fig3_sample_metadata.csv') else opts$metadata_file
  input_names <- fig3_input_names()
  inventory <- list(); status <- list(); results <- list()
  for(panel in opts$panels) {
    cat('\nRunning Figure 3',panel,'...\n',sep='')
    files <- file.path(data_dir,input_names[[panel]])
    if(panel=='c') files <- c(files,metadata)
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
        e=run_fig3e(data_dir,panel_out))
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

# Sourcing this script defines functions only; Rscript executes the selected panels.
if(sys.nframe()==0L) fig3_main()


