if (.Platform$OS.type == "windows") invisible(suppressWarnings(Sys.setlocale("LC_ALL", "English_United States.utf8")))
fig5_extra_library <- Sys.getenv("FIG5_R_LIB", "")
if (nzchar(fig5_extra_library)) .libPaths(c(fig5_extra_library, .libPaths()))

prepare_fig5_data <- function(long, meta, expected_samples=NULL, expected_subtypes=NULL) {
  required <- c('Sample','Setting','City','PhysicalSite','Sample_Date','SamplingMonth')
  if(!all(required %in% names(meta))) stop('Metadata need: ',paste(required,collapse=', '))
  if(!all(c('Sample','Subtype','Abundance') %in% names(long))) stop('Abundance table needs Sample, Subtype, Abundance')
  meta <- as.data.frame(meta,stringsAsFactors=FALSE)
  meta$Sample <- as.character(meta$Sample)
  if(anyNA(meta[required]) || anyDuplicated(meta$Sample) || any(!nzchar(meta$Sample))) stop('Missing or duplicated sample metadata')
  for(v in c('Setting','City','PhysicalSite','SamplingMonth')) {
    meta[[v]] <- as.character(meta[[v]])
    if(any(!nzchar(meta[[v]]))) stop('Blank metadata field: ',v)
  }
  if(any(!meta$Setting %in% c('Hospital','WWTP','Community','Wet market'))) stop('Unknown setting')
  dates <- as.Date(meta$Sample_Date)
  if(anyNA(dates) || any(format(dates,'%Y-%m')!=meta$SamplingMonth)) stop('Date and sampling-month mismatch')
  meta$Sample_Date <- dates
  site_map <- unique(meta[c('PhysicalSite','Setting','City')])
  if(anyDuplicated(site_map$PhysicalSite)) stop('A physical site maps to more than one setting/city')
  long <- as.data.frame(long,stringsAsFactors=FALSE)
  long$Sample <- as.character(long$Sample); long$Subtype <- as.character(long$Subtype)
  if(anyNA(long[c('Sample','Subtype')]) || any(!nzchar(long$Subtype))) stop('Missing sample/subtype identifier')
  if(anyDuplicated(paste(long$Sample,long$Subtype,sep='\034'))) stop('Duplicate sample-subtype records')
  if(!setequal(long$Sample,meta$Sample)) stop('Abundance and metadata sample sets differ')
  if(!is.null(expected_samples) && nrow(meta)!=expected_samples) stop('Unexpected sample count: ',nrow(meta))
  subtypes <- unique(long$Subtype)
  if(!is.null(expected_subtypes) && length(subtypes)!=expected_subtypes) stop('Unexpected Tier I subtype count: ',length(subtypes))
  if(nrow(long)!=nrow(meta)*length(subtypes) || any(table(long$Sample)!=length(subtypes))) stop('Incomplete sample-subtype matrix')
  raw <- as.character(long$Abundance)
  blank <- is.na(raw) | trimws(raw)==''
  y <- suppressWarnings(as.numeric(raw))
  if(any(is.na(y)&!blank) || any(!is.finite(y[!blank])) || any(y<0,na.rm=TRUE)) stop('Invalid abundance: use non-negative numbers or confirmed blank nondetections')
  y[blank] <- 0
  long$Abundance <- y
  long$ARG_type <- sub('__.*$','',long$Subtype)
  if(any(long$ARG_type==long$Subtype)) stop('Subtype IDs must retain the ARG-type__subtype prefix')
  idx <- match(long$Sample,meta$Sample)
  if('Setting' %in% names(long) && any(as.character(long$Setting)!=meta$Setting[idx])) stop('Workbook and canonical setting disagree')
  for(v in setdiff(required,'Sample')) long[[v]] <- meta[[v]][idx]
  totals <- rowsum(y,group=long$Sample,reorder=FALSE)
  richness <- rowsum(as.integer(y>0),group=long$Sample,reorder=FALSE)
  blank_counts <- rowsum(as.integer(blank),group=long$Sample,reorder=FALSE)
  samples <- meta
  samples$Total_abundance <- as.numeric(totals[match(meta$Sample,rownames(totals)),1])
  samples$Richness <- as.integer(richness[match(meta$Sample,rownames(richness)),1])
  samples$Raw_blank_count <- as.integer(blank_counts[match(meta$Sample,rownames(blank_counts)),1])
  list(long=long,meta=meta,samples=samples,subtypes=subtypes,
       blank_cells=long[blank,c('Sample','Subtype','Setting')],
       zero_profiles=samples[samples$Total_abundance==0,,drop=FALSE])
}

read_fig5_data <- function(data_root, metadata_file, audit_dir=NULL) {
  if(!requireNamespace('readxl',quietly=TRUE)) stop('Install R package readxl')
  suppressWarnings(Sys.setlocale('LC_ALL','English_United States.utf8'))
  if(!file.exists(metadata_file)) stop('Canonical sample metadata not found: ',metadata_file)
  if(grepl('\\.xlsx$',metadata_file,ignore.case=TRUE)) {
    raw_meta <- readxl::read_excel(metadata_file)
    if(!all(c('Sample','City','Sample_Date','Sample_Rename') %in% names(raw_meta))) stop('Unexpected source metadata schema')
    meta <- data.frame(Sample=raw_meta$Sample,City=raw_meta$City,
      Sample_Date=as.Date(as.character(raw_meta$Sample_Date),'%Y%m%d'),
      PhysicalSite=sub('_[0-9]{8}$','',raw_meta$Sample_Rename),stringsAsFactors=FALSE)
    meta$SamplingMonth <- format(meta$Sample_Date,'%Y-%m')
  } else meta <- utils::read.csv(metadata_file,stringsAsFactors=FALSE,check.names=FALSE,fileEncoding='UTF-8')
  locations <- c(Hospital='\u533b\u9662',WWTP='\u6c61\u6c34\u5904\u7406\u5382',Community='\u793e\u533a','Wet market'='\u519c\u96c6\u8d38\u5e02\u573a')
  suffix <- Sys.getenv('FIG5_SUBTYPE_SUFFIX', '.xlsx')
  folder <- if(dir.exists(file.path(data_root,'filtered_subtype_xlsx'))) file.path(data_root,'filtered_subtype_xlsx') else data_root
  blocks <- list(); manifests <- list(); ref <- NULL; seen <- character()
  for(setting in names(locations)) {
    path <- file.path(folder,paste0('subtype_by_site_',locations[[setting]],suffix))
    if(!file.exists(path)) stop('Missing Tier I workbook: ',path)
    for(sheet in readxl::excel_sheets(path)) {
      z <- readxl::read_excel(path,sheet=sheet,.name_repair='minimal')
      if(ncol(z)<2L) stop('No sample columns in ',basename(path),' / ',sheet)
      ids <- names(z)[-1]; subtypes <- as.character(z[[1]])
      if(anyNA(subtypes) || anyDuplicated(subtypes) || anyDuplicated(ids) || any(ids %in% seen)) stop('Duplicate subtype/sample keys in input workbook')
      if(is.null(ref)) ref <- subtypes
      if(!setequal(ref,subtypes)) stop('Inconsistent Tier I subtype sets across sheets')
      z <- z[match(ref,subtypes),,drop=FALSE]
      b <- data.frame(Sample=rep(ids,each=length(ref)),Subtype=rep(ref,length(ids)),
        Abundance=as.vector(as.matrix(z[-1])),Setting=setting,stringsAsFactors=FALSE)
      blocks[[length(blocks)+1L]] <- b
      seen <- c(seen,ids)
      manifests[[length(manifests)+1L]] <- data.frame(Workbook=basename(path),Sheet=sheet,Setting=setting,Samples=length(ids),Subtypes=length(ref))
    }
  }
  long <- do.call(rbind,blocks)
  if(!'Setting' %in% names(meta)) {
    setting_map <- unique(long[c('Sample','Setting')])
    meta$Setting <- setting_map$Setting[match(meta$Sample,setting_map$Sample)]
  }
  result <- prepare_fig5_data(long,meta)
  result$manifest <- do.call(rbind,manifests)
  if(!is.null(audit_dir)) {
    dir.create(audit_dir,recursive=TRUE,showWarnings=FALSE)
    for(name in c('meta','samples','manifest','blank_cells','zero_profiles'))
      utils::write.csv(result[[name]],file.path(audit_dir,paste0(name,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
    writeLines('Blank subtype abundance cells and numeric zeros denote nondetection. Missing sample metadata are not imputed.',file.path(audit_dir,'blank_policy.txt'))
  }
  result
}

fig5_bc_settings <- c("Hospital", "WWTP", "Community", "Wet market")
fig5_bc_site_columns <- c("Site.Hospital", "Site.WWTP", "Site.Community", "Site.Wet.market")
fig5_bc_levels <- paste("Level", c("I", "II", "III", "IV"))
fig5_bc_sites_palette <- c(Community="#6EC5C1", Hospital="#E06C75", `Wet market`="#E9C95B", WWTP="#4E79A7")
fig5_bc_tier_palette <- c(`Level I`="#B2182B", `Level II`="#EF8A62", `Level III`="#BDBDBD", `Level IV`="#E0E0E0")
fig5_bc_type_palette <- c(Aminoglycoside="#4E79A7", Bleomycin="#E39B76", Chloramphenicol="#B5658D", Florfenicol="#AFC7E8", Fosfomycin="#B7E0E5", MLS="#74C2B3", Multidrug="#7E6BB5", Polymyxin="#75C0C1", Quinolone="#9A84D6", Rifamycin="#A9C98A", Streptothricin="#6DBDE3", Sulfonamide="#C8A24B", Tetracycline="#E8D98F", Trimethoprim="#E39A9A", Vancomycin="#AFC9A0", beta_lactam="#7DA34D", Bacitracin="#D6B97A", Others="#CFCFCF")

fig5_bc_type_name <- function(x) {
  x <- as.character(x)
  x[x=="macrolide-lincosamide-streptogramin"] <- "MLS"
  x[x=="glycopeptide"] <- "vancomycin"
  idx <- !x %in% c("MLS", "beta_lactam")
  x[idx] <- paste0(toupper(substr(x[idx], 1, 1)), substring(x[idx], 2))
  x
}

validate_fig5_bc_table <- function(x, scope=c("risk", "tier_i")) {
  scope <- match.arg(scope)
  need <- c("arg_type", "arg_subtype", "core_flag", "risk_level", "host_breadth", "mobility_ratio", fig5_bc_site_columns)
  miss <- setdiff(need, names(x))
  if(length(miss)) stop("Missing source fields: ", paste(miss, collapse=", "))
  if(anyNA(x$arg_subtype) || any(!nzchar(trimws(x$arg_subtype))) || anyDuplicated(x$arg_subtype)) stop("Missing or duplicate subtype IDs")
  if(anyNA(x$arg_type) || any(!nzchar(x$arg_type))) stop("Missing ARG type")
  if(anyNA(x$core_flag) || any(toupper(as.character(x$core_flag))!="TRUE")) stop("Input contains non-core or undefined core records")
  for(i in seq_along(fig5_bc_settings)) {
    value <- trimws(as.character(x[[fig5_bc_site_columns[i]]]))
    absent <- is.na(value) | value %in% c("", "NA", "0", "FALSE", "False", "false", "No", "NO", "no")
    if(any(!absent & value!=fig5_bc_settings[i])) stop("Unknown CORE membership token in ", fig5_bc_site_columns[i])
    x[[fig5_bc_settings[i]]] <- !absent
  }
  x$n_settings <- rowSums(x[fig5_bc_settings])
  if(any(x$n_settings==0)) stop("Core records must have at least one setting membership")
  if(anyNA(x$risk_level) || any(!x$risk_level %in% fig5_bc_levels)) stop("Unexpected tier label")
  if(scope=="tier_i" && any(x$risk_level!="Level I")) stop("Tier I input contains other tier labels")
  mob <- suppressWarnings(as.numeric(x$mobility_ratio)); breadth <- suppressWarnings(as.numeric(x$host_breadth))
  if(any(!is.finite(mob) | mob<0 | mob>1) || any(!is.finite(breadth) | breadth<1 | breadth!=floor(breadth))) stop("Invalid mobility fraction or host breadth")
  expected_tier <- ifelse(mob<0.1,"Level IV",ifelse(mob<0.5,"Level III",ifelse(mob>=0.7 & breadth>=3,"Level I","Level II")))
  if(any(expected_tier!=x$risk_level)) stop("Recorded tier disagrees with the supplied threshold definitions")
  x$arg_type_label <- fig5_bc_type_name(x$arg_type)
  x$region <- apply(x[fig5_bc_settings],1,function(z) paste(fig5_bc_settings[as.logical(z)],collapse=" + "))
  x$region_code <- apply(x[fig5_bc_settings],1,function(z) paste(c("H","W","C","M")[as.logical(z)],collapse="+"))
  x
}

read_fig5_bc_inputs <- function(data_dir) {
  read_input <- function(name) {
    fp <- file.path(data_dir,name)
    if(!file.exists(fp)) stop("Missing required Fig5 b/c input: ",fp)
    read.csv(fp, check.names=FALSE, stringsAsFactors=FALSE, fileEncoding="UTF-8")
  }
  risk <- validate_fig5_bc_table(read_input("subtype_risk_4level_core_only_final.csv"),"risk")
  tier_i <- validate_fig5_bc_table(read_input("subtype_level_I_core_only.csv"),"tier_i")
  core <- read_input("core_subtype_list.csv")
  panel <- read_input("panel_subtype_mobility_summary.csv")
  if(!"arg_subtype" %in% names(core) || !"subtype_raw" %in% names(panel)) stop("Invalid core/mobility provenance schema")
  if(anyDuplicated(core$arg_subtype) || anyDuplicated(panel$subtype_raw)) stop("Duplicate provenance subtype IDs")
  supported <- intersect(core$arg_subtype,panel$subtype_raw)
  if(any(!risk$arg_subtype %in% supported)) stop("Tier table includes a subtype outside contig-supported core")
  if(!setequal(tier_i$arg_subtype,risk$arg_subtype[risk$risk_level=="Level I"])) stop("Tier I file is not the exact Tier I subset")
  matched <- risk[match(tier_i$arg_subtype,risk$arg_subtype),]
  if(any(matched$region!=tier_i$region) || any(matched$arg_type!=tier_i$arg_type)) stop("Tier I identities or memberships disagree with full tier file")
  if(any(!fig5_bc_site_columns %in% names(core))) stop("Core provenance is missing setting flags")
  core_match <- core[match(risk$arg_subtype,core$arg_subtype),]
  for(i in seq_along(fig5_bc_settings)) {
    z <- core_match[[fig5_bc_site_columns[i]]]
    flag <- !is.na(z) & z==fig5_bc_settings[i]
    if(any(flag!=risk[[fig5_bc_settings[i]]])) stop("Core membership disagrees with core provenance")
  }
  counts <- c(nrow(core),length(setdiff(core$arg_subtype,panel$subtype_raw)),length(supported),length(setdiff(supported,risk$arg_subtype)),nrow(risk))
  long <- do.call(rbind,lapply(fig5_bc_settings,function(s) data.frame(arg_subtype=risk$arg_subtype[risk[[s]]],arg_type=risk$arg_type_label[risk[[s]]],risk_level=risk$risk_level[risk[[s]]],setting=s)))
  nt <- sort(table(long$arg_type),decreasing=TRUE)
  top <- names(nt)[seq_len(min(9L,length(nt)))]
  long$type_display <- ifelse(long$arg_type %in% top,long$arg_type,"Others")
  flow <- aggregate(rep(1L,nrow(long)),long[c("type_display","risk_level","setting")],sum)
  names(flow)[4] <- "n_subtypes"
  flow$type_display <- factor(flow$type_display,levels=c(top,"Others"))
  flow$risk_level <- factor(flow$risk_level,levels=fig5_bc_levels)
  flow$setting <- factor(flow$setting,levels=fig5_bc_settings)
  regions <- aggregate(rep(1L,nrow(tier_i)),tier_i[c("region","region_code","n_settings")],sum)
  names(regions)[4] <- "n_subtypes"
  regions <- regions[order(-regions$n_settings,-regions$n_subtypes,regions$region),]
  donuts <- aggregate(rep(1L,nrow(tier_i)),tier_i[c("region","region_code","arg_type_label")],sum)
  names(donuts)[4] <- "n_subtypes"
  donuts$fraction <- donuts$n_subtypes / ave(donuts$n_subtypes,donuts$region,FUN=sum)
  provenance <- data.frame(stage=c("Core subtypes","Core without contig support","Contig-supported core","Contig-supported core without host-set assignment","Prioritized subtypes"),n_subtypes=as.integer(counts))
  list(risk=risk,tier_i=tier_i,flow=flow,regions=regions,donuts=donuts,provenance=provenance)
}

run_fig5_bc <- function(data_dir, out_dir) {
  packages <- c("ggplot2","ggalluvial","VennDiagram","patchwork")
  absent <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(absent)) stop("Missing R packages: ",paste(absent,collapse=", "),". Install them before running this module.")
  x <- read_fig5_bc_inputs(data_dir)
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  write_out <- function(d,name) write.csv(d,file.path(out_dir,name),row.names=FALSE,na="",fileEncoding="UTF-8")
  write_out(x$provenance,"Fig5bc_provenance_counts.csv")
  write_out(x$flow,"Fig5b_alluvial_subtype_setting_occurrences.csv")
  write_out(as.data.frame(table(factor(x$risk$risk_level,levels=fig5_bc_levels))),"Fig5b_unique_tier_counts.csv")
  write_out(x$tier_i[c("arg_subtype","arg_type",fig5_bc_settings,"n_settings","region")],"Fig5c_core_memberships.csv")
  write_out(x$regions,"Fig5c_exact_regions.csv")
  write_out(x$donuts,"Fig5c_region_type_counts.csv")
  types <- union(x$risk$arg_type_label,"Others")
  palette <- fig5_bc_type_palette
  extra <- setdiff(types,names(palette))
  if(length(extra)) palette <- c(palette,setNames(grDevices::hcl.colors(length(extra),"Dark 3"),extra))
  p_b <- ggplot2::ggplot(x$flow,ggplot2::aes(axis1=type_display,axis2=risk_level,axis3=setting,y=n_subtypes)) +
    ggalluvial::geom_alluvium(ggplot2::aes(fill=risk_level),width=.16,alpha=.95,knot.pos=.35,color="white",linewidth=.15) +
    ggalluvial::geom_stratum(ggplot2::aes(fill=ggplot2::after_stat(stratum)),width=.16,color="white",linewidth=.4) +
    ggalluvial::stat_stratum(geom="text",ggplot2::aes(label=sub("Level","Tier",ggplot2::after_stat(stratum))),size=3.3,color="black") +
    ggplot2::scale_x_discrete(limits=c("ARG type","Tier","Core in setting"),expand=c(.08,.08)) +
    ggplot2::scale_fill_manual(values=c(fig5_bc_tier_palette,palette,fig5_bc_sites_palette),breaks=fig5_bc_levels,labels=sub("Level","Tier",fig5_bc_levels)) +
    ggplot2::labs(x=NULL,y="Subtype-setting occurrences",fill="Prioritization tier",caption="A subtype contributes once for each setting in which it is core.") +
    ggplot2::theme_minimal(base_size=12) + ggplot2::theme(panel.grid=ggplot2::element_blank(),axis.text.y=ggplot2::element_blank(),legend.position="right",plot.caption=ggplot2::element_text(hjust=0,size=9))
  sets <- lapply(fig5_bc_settings,function(s) x$tier_i$arg_subtype[x$tier_i[[s]]]); names(sets) <- fig5_bc_settings
  n_int <- function(i) length(Reduce(intersect,sets[i]))
  venn <- VennDiagram::draw.quad.venn(area1=length(sets[[1]]),area2=length(sets[[2]]),area3=length(sets[[3]]),area4=length(sets[[4]]),n12=n_int(c(1,2)),n13=n_int(c(1,3)),n14=n_int(c(1,4)),n23=n_int(c(2,3)),n24=n_int(c(2,4)),n34=n_int(c(3,4)),n123=n_int(c(1,2,3)),n124=n_int(c(1,2,4)),n134=n_int(c(1,3,4)),n234=n_int(c(2,3,4)),n1234=n_int(1:4),category=fig5_bc_settings,fill=fig5_bc_sites_palette[fig5_bc_settings],alpha=rep(.4,4),col=rep("grey65",4),lwd=rep(1.35,4),cex=1.7,fontface="plain",fontfamily="sans",label.col="black",cat.cex=1.35,cat.fontface="plain",cat.fontfamily="sans",cat.col=fig5_bc_sites_palette[fig5_bc_settings],cat.pos=c(180,110,70,0),cat.dist=rep(.05,4),margin=.03,scaled=FALSE,ind=FALSE)

  label_xy <- list(Hospital=c(.12,.35), WWTP=c(.89,.35), Community=c(.29,.94), `Wet market`=c(.66,.94))
  for(i in seq_along(venn)) {
    if(inherits(venn[[i]],"text") && length(venn[[i]]$label)==1L && venn[[i]]$label %in% names(label_xy)) {
      xy <- label_xy[[venn[[i]]$label]]
      venn[[i]]$x <- grid::unit(xy[1],"npc"); venn[[i]]$y <- grid::unit(xy[2],"npc")
    }
  }
  p_venn <- patchwork::wrap_elements(full=grid::grobTree(children=venn))
  dd <- x$donuts
  dd$region_code <- factor(dd$region_code,levels=x$regions$region_code)
  totals <- x$regions; totals$region_code <- factor(totals$region_code,levels=x$regions$region_code)
  p_donut <- ggplot2::ggplot(dd,ggplot2::aes(x=2,y=fraction,fill=arg_type_label)) +
    ggplot2::geom_col(width=.85,color="white",linewidth=.25) + ggplot2::coord_polar(theta="y") +
    ggplot2::xlim(.5,2.5) + ggplot2::facet_wrap(~region_code,ncol=3) +
    ggplot2::geom_text(data=totals,ggplot2::aes(x=.5,y=0,label=n_subtypes),inherit.aes=FALSE,size=3.7) +
    ggplot2::scale_fill_manual(values=palette,drop=FALSE,breaks=sort(unique(x$tier_i$arg_type_label))) +
    ggplot2::labs(fill="ARG type",caption="Exact core-membership regions. H: hospital; W: WWTP; C: community; M: wet market.\nDonut areas are not scaled by subtype count; centre labels give region totals.") +
    ggplot2::theme_void(base_size=11) + ggplot2::theme(strip.text=ggplot2::element_text(size=10),legend.position="bottom",plot.caption=ggplot2::element_text(hjust=0,size=8)) +
    ggplot2::guides(fill=ggplot2::guide_legend(ncol=3))
  p_c <- patchwork::wrap_plots(p_venn,p_donut,ncol=2,widths=c(1,1.3))
  save_plot <- function(plot,name,w,h) {
    ggplot2::ggsave(file.path(out_dir,paste0(name,".pdf")),plot,width=w,height=h,device=grDevices::cairo_pdf)
    ggplot2::ggsave(file.path(out_dir,paste0(name,".png")),plot,width=w,height=h,dpi=300,bg="white")
  }
  save_plot(p_b,"Fig5b_alluvial",10.5,6)
  save_plot(p_venn,"Fig5c_core_Venn",5.2,5.2)
  save_plot(p_donut,"Fig5c_type_donuts",7.5,8)
  save_plot(p_c,"Fig5c_Venn_and_type_donuts",13.5,8)
  invisible(x)
}

fig5d_require_packages <- function() {
  needed <- c("vegan", "ggplot2")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install the declared dependencies first: ",
                            paste(missing, collapse = ", "), call. = FALSE)
}

fig5d_prepare_data <- function(long, meta, zero_policy = c("error", "drop")) {
  zero_policy <- match.arg(zero_policy)
  required_long <- c("Sample", "Subtype", "Abundance")
  required_meta <- c("Sample", "Setting", "City", "PhysicalSite")
  if (!all(required_long %in% names(long))) {
    stop("long must contain Sample, Subtype, and Abundance.", call. = FALSE)
  }
  if (!all(required_meta %in% names(meta))) {
    stop("meta must contain Sample, Setting, City, and PhysicalSite.", call. = FALSE)
  }
  long <- as.data.frame(long[, required_long, drop = FALSE])
  meta <- as.data.frame(meta)
  for (nm in c("Sample", "Subtype")) long[[nm]] <- as.character(long[[nm]])
  for (nm in required_meta) meta[[nm]] <- as.character(meta[[nm]])
  if (!nrow(long) || !nrow(meta)) stop("Input data are empty.", call. = FALSE)
  if (anyNA(long[c("Sample", "Subtype")]) ||
      any(!nzchar(trimws(long$Sample))) || any(!nzchar(trimws(long$Subtype)))) {
    stop("Sample and Subtype cannot be blank or missing.", call. = FALSE)
  }
  if (anyNA(meta[required_meta]) ||
      any(vapply(meta[required_meta], function(x) any(!nzchar(trimws(x))), logical(1)))) {
    stop("Required sample metadata cannot be blank or missing.", call. = FALSE)
  }
  if (anyDuplicated(meta$Sample)) stop("Metadata Sample IDs must be unique.", call. = FALSE)
  if (!is.numeric(long$Abundance) || any(!is.finite(long$Abundance))) {
    stop("Abundance contains non-numeric or non-finite values; resolve them upstream.",
         call. = FALSE)
  }
  if (any(long$Abundance < 0)) stop("Abundance contains negative values.", call. = FALSE)
  if (anyDuplicated(long[c("Sample", "Subtype")])) {
    stop("Duplicate sample-subtype rows must be resolved explicitly upstream.", call. = FALSE)
  }
  unknown <- setdiff(long$Sample, meta$Sample)
  absent <- setdiff(meta$Sample, long$Sample)
  if (length(unknown)) stop("Abundances have samples absent from metadata: ",
                            paste(unknown, collapse = ", "), call. = FALSE)
  if (length(absent)) stop("Metadata samples have no abundance records: ",
                           paste(absent, collapse = ", "), call. = FALSE)
  setting_levels <- c("Community", "Hospital", "Wet market", "WWTP")
  if (any(!meta$Setting %in% setting_levels)) {
    stop("Setting must be Community, Hospital, Wet market, or WWTP.", call. = FALSE)
  }
  subtype_levels <- sort(unique(long$Subtype))
  if (nrow(long) != nrow(meta) * length(subtype_levels)) {
    stop("Incomplete sample-subtype matrix: every sample must have an explicit ",
         "record for every subtype. Resolve missing records upstream; ",
         "they are not silently replaced with zero.", call. = FALSE)
  }
  mat <- matrix(0, nrow(meta), length(subtype_levels),
                dimnames = list(meta$Sample, subtype_levels))
  mat[cbind(match(long$Sample, meta$Sample), match(long$Subtype, subtype_levels))] <-
    long$Abundance
  zero <- rowSums(mat) == 0
  excluded <- meta[zero, , drop = FALSE]
  if (any(zero) && zero_policy == "error") {
    stop("Zero-total profiles cannot enter Bray-Curtis PCoA: ",
         paste(meta$Sample[zero], collapse = ", "),
         ". Use zero_policy='drop' only after reviewing them; exclusions are exported.",
         call. = FALSE)
  }
  if (any(zero)) {
    warning("Explicitly dropping ", sum(zero), " zero-total sample(s); identities are retained.",
            call. = FALSE)
  }
  mat <- mat[!zero, , drop = FALSE]
  meta <- meta[!zero, , drop = FALSE]
  if (nrow(mat) < 3L) stop("At least three nonzero samples are required.", call. = FALSE)
  meta$Setting <- factor(meta$Setting, levels = setting_levels)
  rownames(meta) <- NULL
  list(matrix = mat, metadata = meta, excluded_zero_samples = excluded)
}

fig5d_ordination <- function(mat) {
  fig5d_require_packages()
  dissimilarity <- vegan::vegdist(mat, method = "bray")
  pcoa <- stats::cmdscale(dissimilarity, eig = TRUE, k = 2)
  if (ncol(pcoa$points) < 2L || any(!is.finite(pcoa$points))) {
    stop("This abundance matrix does not yield two finite PCoA axes.", call. = FALSE)
  }
  positive_sum <- sum(pcoa$eig[pcoa$eig > 0])
  points <- as.data.frame(pcoa$points)
  names(points) <- c("PCoA1", "PCoA2")
  points$Sample <- rownames(points)
  rownames(points) <- NULL
  list(distance = dissimilarity, points = points, eigenvalues = pcoa$eig,
       axis_percent = 100 * pcoa$eig[1:2] / positive_sum,
       negative_eigen_fraction = sum(abs(pcoa$eig[pcoa$eig < 0])) / sum(abs(pcoa$eig)))
}

run_fig5d <- function(long, meta, output_dir, ellipse_level = 0.68,
                      permutations = 999L, seed = 123L,
                      zero_policy = c("error", "drop"),
                      width = 7.2, height = 5.5, dpi = 600,
                      compute_dispersion = TRUE) {
  fig5d_require_packages()
  zero_policy <- match.arg(zero_policy)
  if (length(ellipse_level) != 1L || !is.finite(ellipse_level) ||
      ellipse_level <= 0 || ellipse_level >= 1) stop("ellipse_level must be between 0 and 1.")
  if (length(permutations) != 1L || !is.finite(permutations) ||
      permutations < 1L || permutations != as.integer(permutations)) {
    stop("permutations must be a positive integer.")
  }
  if (!is.logical(compute_dispersion) || length(compute_dispersion) != 1L ||
      is.na(compute_dispersion)) stop("compute_dispersion must be TRUE or FALSE.")
  prepared <- fig5d_prepare_data(long, meta, zero_policy)
  message("Fig5d: data validated; calculating Bray-Curtis PCoA.")
  ordination <- fig5d_ordination(prepared$matrix)
  metadata <- prepared$metadata
  if (nlevels(droplevels(metadata$Setting)) < 2L ||
      nrow(metadata) <= nlevels(droplevels(metadata$Setting))) {
    stop("PERMANOVA and PERMDISP require multiple settings and residual degrees of freedom.")
  }
  had_seed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = .GlobalEnv)
  on.exit({
    if (had_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
    else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
      rm(".Random.seed", envir = .GlobalEnv)
  }, add = TRUE)
  dissimilarity <- ordination$distance
  message("Fig5d: ordination complete; running unrestricted sample-level PERMANOVA.")
  set.seed(seed)
  permanova <- vegan::adonis2(dissimilarity ~ Setting, data = metadata,
                             permutations = permutations)
  statistics <- data.frame(
    Test = "PERMANOVA", Term = "Setting",
    Df = permanova$Df[1], R2 = permanova$R2[1],
    F = permanova$F[1], P = permanova$`Pr(>F)`[1],
    PermutationsRequested = as.integer(permutations), Seed = as.integer(seed),
    PermutationScheme = "unrestricted samples",
    Inference = "exploratory; repeated sampling not adjusted")
  plot_data <- cbind(ordination$points,
                     metadata[match(ordination$points$Sample, metadata$Sample),
                              setdiff(names(metadata), "Sample"), drop = FALSE])
  plotted <- fig5d_plot(plot_data, statistics, ordination$axis_percent, ellipse_level)
  p <- plotted$plot
  ellipse_ok <- plotted$ellipse_ok
  audit <- data.frame(
    n_samples_input = nrow(meta), n_samples_used = nrow(metadata),
    n_zero_samples_removed = nrow(prepared$excluded_zero_samples),
    n_physical_sites = length(unique(metadata$PhysicalSite)),
    n_features = ncol(prepared$matrix),
    abundance_scale = "original copies per cell; no transformation or normalization",
    distance = "Bray-Curtis", ordination = "stats::cmdscale", correction = "none",
    axis_denominator = "sum of positive eigenvalues",
    PCoA1_percent = ordination$axis_percent[1], PCoA2_percent = ordination$axis_percent[2],
    negative_eigen_fraction = ordination$negative_eigen_fraction,
    ellipse_type = "normal data ellipse; not mean confidence interval",
    ellipse_level = ellipse_level,
    plot_limits_include_ellipses = TRUE,
    settings_without_ellipse = paste(names(ellipse_ok)[!ellipse_ok], collapse = ";"),
    permdisp_requested = compute_dispersion, permdisp_completed = FALSE,
    permdisp_center = if (compute_dispersion) "spatial median" else "not run",
    permdisp_bias_adjust = if (compute_dispersion) FALSE else NA,
    inference = "exploratory sample-level tests; repeated sampling not adjusted",
    vegan_version = as.character(utils::packageVersion("vegan")),
    ggplot2_version = as.character(utils::packageVersion("ggplot2")))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  ggplot2::ggsave(file.path(output_dir, "Fig5d_PCoA.png"), p,
                  width = width, height = height, dpi = dpi, bg = "white")
  pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf
  ggplot2::ggsave(file.path(output_dir, "Fig5d_PCoA.pdf"), p,
                  width = width, height = height, device = pdf_device, bg = "white")
  utils::write.csv(statistics, file.path(output_dir, "Fig5d_PCoA_stats.csv"), row.names = FALSE)
  utils::write.csv(audit, file.path(output_dir, "Fig5d_PCoA_audit.csv"), row.names = FALSE)
  utils::write.csv(plot_data, file.path(output_dir, "Fig5d_PCoA_scores.csv"), row.names = FALSE)
  utils::write.csv(prepared$excluded_zero_samples,
                   file.path(output_dir, "Fig5d_PCoA_excluded_zero_samples.csv"), row.names = FALSE)
  utils::write.csv(data.frame(Axis = seq_along(ordination$eigenvalues),
                              Eigenvalue = ordination$eigenvalues),
                   file.path(output_dir, "Fig5d_PCoA_eigenvalues.csv"), row.names = FALSE)
  message("Fig5d: main plot, scores, audit, and PERMANOVA results saved.")
  dispersion <- permdisp <- NULL
  if (compute_dispersion) {
    message("Fig5d: running PERMDISP; main outputs are already available.")
    dispersion <- vegan::betadisper(dissimilarity, metadata$Setting,
                                    type = "median", bias.adjust = FALSE)
    set.seed(seed)
    permdisp <- vegan::permutest(dispersion, permutations = permutations)
    dispersion_statistics <- statistics[1, , drop = FALSE]
    dispersion_statistics$Test <- "PERMDISP"
    dispersion_statistics$Df <- permdisp$tab$Df[1]
    dispersion_statistics$R2 <- NA_real_
    dispersion_statistics$F <- permdisp$tab$F[1]
    dispersion_statistics$P <- permdisp$tab$`Pr(>F)`[1]
    statistics <- rbind(statistics, dispersion_statistics)
    audit$permdisp_completed <- TRUE
    utils::write.csv(statistics, file.path(output_dir, "Fig5d_PCoA_stats.csv"), row.names = FALSE)
    utils::write.csv(audit, file.path(output_dir, "Fig5d_PCoA_audit.csv"), row.names = FALSE)
  }
  list(plot = p, data = plot_data, stats = statistics, audit = audit,
       matrix = prepared$matrix, ordination = ordination,
       excluded_zero_samples = prepared$excluded_zero_samples,
       permanova = permanova, dispersion = dispersion, permdisp = permdisp)
}

fig5d_plot <- function(plot_data, statistics, axis_percent, ellipse_level = 0.68) {
  fig5d_require_packages()
  plot_data$Setting <- factor(plot_data$Setting,
                              levels = c("Community", "Hospital", "Wet market", "WWTP"))

  ellipse_ok <- vapply(split(plot_data, plot_data$Setting, drop = TRUE), function(x) {
    nrow(x) >= 4L && qr(stats::cov(x[c("PCoA1", "PCoA2")]))$rank == 2L
  }, logical(1))
  ellipse_data <- plot_data[as.character(plot_data$Setting) %in% names(ellipse_ok)[ellipse_ok], ]
  palette <- c("Community" = "#6EC5C1", "Hospital" = "#E06C75",
               "Wet market" = "#E9C95B", "WWTP" = "#4E79A7")
  ellipse_plot <- ggplot2::ggplot(plot_data,
                                  ggplot2::aes(PCoA1, PCoA2, colour = Setting)) +
    ggplot2::stat_ellipse(data = ellipse_data, ggplot2::aes(group = Setting),
                         type = "norm", level = ellipse_level, linewidth = 0.6,
                         alpha = 0.18, show.legend = FALSE)
  ellipse_coords <- if (nrow(ellipse_data)) ggplot2::ggplot_build(ellipse_plot)$data[[1]]
                    else data.frame(x = numeric(), y = numeric())
  format_p <- function(p) if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
  annotation <- paste0("PERMANOVA\nSetting: ", format_p(statistics$P[1]),
                        sprintf(", R2 = %.3f", statistics$R2[1]))

  x_range <- range(c(plot_data$PCoA1, ellipse_coords$x), na.rm = TRUE)
  y_range <- range(c(plot_data$PCoA2, ellipse_coords$y), na.rm = TRUE)
  x_pad <- max(diff(x_range) * 0.08, .Machine$double.eps)
  y_pad <- max(diff(y_range) * 0.12, .Machine$double.eps)
  p <- ellipse_plot +
    ggplot2::geom_point(size = 2.9, alpha = 0.92) +
    ggplot2::annotate("text", x = x_range[1] + x_pad,
                      y = y_range[2] + 0.2 * y_pad, label = annotation,
                      hjust = 0, vjust = 1, size = 4) +
    ggplot2::scale_colour_manual(values = palette, drop = FALSE, name = "Site type") +
    ggplot2::labs(x = sprintf("PCoA1 (%.1f%%)", axis_percent[1]),
                  y = sprintf("PCoA2 (%.1f%%)", axis_percent[2])) +
    ggplot2::coord_fixed(xlim = x_range + c(-x_pad, x_pad),
                         ylim = y_range + c(-y_pad, 0.35 * y_pad), clip = "off") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(axis.line = ggplot2::element_line(linewidth = 0.7, colour = "black"),
                   axis.ticks = ggplot2::element_line(linewidth = 0.6, colour = "black"),
                   axis.ticks.length = grid::unit(0.18, "cm"),
                   axis.title = ggplot2::element_text(size = 15),
                   axis.text = ggplot2::element_text(size = 14),
                   legend.position = "right", legend.title = ggplot2::element_text(size = 11),
                   legend.text = ggplot2::element_text(size = 10),
                   legend.key = ggplot2::element_blank(),
                   plot.margin = ggplot2::margin(5, 18, 5, 5)) +
    ggplot2::guides(colour = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  list(plot = p, ellipse_ok = ellipse_ok)
}

validate_nb_response <- function(y) {
  if (anyNA(y) || any(!is.finite(y)) || any(y < 0) || any(abs(y-round(y)) > 1e-8))
    stop("Negative-binomial responses must be nonnegative integer counts.")
  TRUE
}

fit_fig5_model <- function(samples, outcome, family, out_dir) {
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if(!family %in% c('negative_binomial','gaussian_original')) stop('Unsupported or unconfirmed model family.')
  if(!outcome %in% c('Richness','Total_abundance')) stop('Unknown outcome.')
  if(outcome=='Total_abundance' && family=='negative_binomial')
    stop('Tier I abundance is continuous copies/cell; negative-binomial fitting or rounding is not permitted.')
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

fig5_half_violin <- function(mapping=NULL, data=NULL, ...) {
  geom <- ggplot2::ggproto("GeomFig5HalfViolin", ggplot2::Geom,
    required_aes=c("x","y"),
    default_aes=ggplot2::aes(weight=1,colour=NA,fill="grey80",linewidth=.3,alpha=NA,linetype="solid"),
    setup_data=function(data,params) {
      if(is.null(data$width)) data$width <- if(is.null(params$width)) .6 else params$width
      data$ymin <- ave(data$y,data$group,FUN=min)
      data$ymax <- ave(data$y,data$group,FUN=max)
      data$xmin <- data$x;data$xmax <- data$x+data$width/2
      data
    },
    draw_group=function(data,panel_params,coord) {
      data$xmaxv <- data$x+data$violinwidth*(data$xmax-data$x)
      up <- data[order(data$y),];up$x <- up$xmaxv
      down <- data[order(-data$y),];down$x <- down$xmin
      shape <- rbind(up,down,up[1,])
      ggplot2::GeomPolygon$draw_panel(shape,panel_params,coord)
    },
    draw_key=ggplot2::draw_key_polygon)
  ggplot2::layer(data=data,mapping=mapping,stat="ydensity",geom=geom,
    position="identity",show.legend=FALSE,inherit.aes=TRUE,
    params=list(trim=TRUE,scale="area",na.rm=FALSE,...))
}

plot_fig5_violin <- function(samples, outcome, model_result, out_dir, basename) {
  settings <- c("Hospital","WWTP","Community","Wet market")
  palette <- c(Hospital="#E06C75",WWTP="#4E79A7",Community="#6EC5C1","Wet market"="#E9C95B")
  d <- samples
  d$Setting <- factor(d$Setting,levels=settings)
  d$x_id <- as.numeric(d$Setting);d$value <- d[[outcome]]
  means <- aggregate(value~Setting,d,mean)
  means$x_id <- match(means$Setting,settings)
  means$label <- sprintf(if(outcome=="Richness") "%.1f" else "%.2f",means$value)
  letters <- model_result$letters
  letters$x_id <- match(letters$Site,settings)
  top <- max(d$value); span <- max(diff(range(d$value)),top*.1,1e-6)
  mean_y <- top+.055*span; letter_y <- top+.13*span
  plot <- ggplot2::ggplot()+
    fig5_half_violin(ggplot2::aes(x=x_id,y=value,fill=Setting),data=d,width=.54,alpha=.82,colour=NA)+
    ggplot2::geom_point(data=d,ggplot2::aes(x=x_id-.15,y=value,colour=Setting),
      position=ggplot2::position_jitter(width=.04,height=0,seed=123),alpha=.15,size=.9,shape=16)+
    ggplot2::geom_boxplot(data=d,ggplot2::aes(x=x_id,y=value,group=Setting,colour=Setting),
      width=.082,outlier.shape=NA,linewidth=.6,fill="white")+
    ggplot2::geom_text(data=means,ggplot2::aes(x=x_id,label=label,colour=Setting),
      y=mean_y,size=4.2,vjust=0)+
    ggplot2::geom_text(data=letters,ggplot2::aes(x=x_id,label=Letters,colour=Site),
      y=letter_y,size=6.2,vjust=0)+
    ggplot2::scale_fill_manual(values=palette)+ggplot2::scale_colour_manual(values=palette)+
    ggplot2::scale_x_continuous(breaks=1:4,labels=settings,expand=ggplot2::expansion(mult=c(.04,.06)))+
    ggplot2::scale_y_continuous(expand=ggplot2::expansion(mult=c(0,0)))+
    ggplot2::coord_cartesian(ylim=c(0,top+.24*span),clip="off")+
    ggplot2::labs(x=NULL,y=if(outcome=="Richness") "Tier I ARG subtype richness" else "Tier I ARG abundance (copies/cell)")+
    ggplot2::theme_classic(base_size=12)+
    ggplot2::theme(panel.border=ggplot2::element_rect(colour="black",fill=NA,linewidth=.65),
      axis.line=ggplot2::element_blank(),axis.text.x=ggplot2::element_text(angle=25,hjust=1,size=12,colour="black"),
      axis.text.y=ggplot2::element_text(size=12,colour="black"),legend.position="none",
      plot.margin=ggplot2::margin(12,14,12,12))
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(file.path(out_dir,paste0(basename,".pdf")),plot,width=5.3,height=4.8,device=grDevices::cairo_pdf)
  ggplot2::ggsave(file.path(out_dir,paste0(basename,".png")),plot,width=5.3,height=4.8,dpi=400,bg="white")
  means$Letters <- letters$Letters[match(as.character(means$Setting),letters$Site)]
  names(means)[names(means)=="value"] <- "Arithmetic_mean"
  utils::write.csv(means,file.path(out_dir,paste0(basename,"_means_letters.csv")),row.names=FALSE)
  plot
}

fig5_model_dependencies <- function() {
  packages <- c("lme4","dplyr","multcompView","ggplot2")
  missing <- packages[!vapply(packages,requireNamespace,logical(1),quietly=TRUE)]
  if(length(missing)) stop("Install required R packages: ",paste(missing,collapse=", "))
}
run_fig5e <- function(samples, output_dir) {
  fig5_model_dependencies()
  model <- fit_fig5_model(samples,"Total_abundance","gaussian_original",output_dir)
  model$plot <- plot_fig5_violin(samples,"Total_abundance",model,output_dir,"Fig5e")
  invisible(model)
}
run_tier1_richness <- function(samples, output_dir) {
  fig5_model_dependencies()
  model <- fit_fig5_model(samples,"Richness","negative_binomial",output_dir)
  model$plot <- plot_fig5_violin(samples,"Richness",model,output_dir,"TierI_richness")
  invisible(model)
}

fig5_fg_months <- function(months) {
  dates <- as.Date(paste0(months, "-01"))
  if (!length(dates) || anyNA(dates)) stop("months must contain valid YYYY-MM values.")
  format(seq(min(dates), max(dates), by = "month"), "%Y-%m")
}

fig5_fg_palette <- function() {
  c("Aminoglycoside" = "#4E79A7", "Antibacterial fatty acid" = "#7FB8B5",
    "Bacitracin" = "#D6B97A", "Bicyclomycin" = "#8BCB88", "Bleomycin" = "#E39B76",
    "Chloramphenicol" = "#B5658D", "Defensin" = "#B7D989", "Edeine" = "#8FB6E8",
    "Factumycin" = "#D17C98", "Florfenicol" = "#AFC7E8", "Fosfomycin" = "#B7E0E5",
    "Fusidic acid" = "#E6BAC4", "MLS" = "#74C2B3", "Multidrug" = "#7E6BB5",
    "Mupirocin" = "#5F9EA0", "Novobiocin" = "#72B77E",
    "Other peptide antibiotics" = "#9FD58A", "Pleuromutilin/Tiamulin" = "#E4CF62",
    "Polymyxin" = "#75C0C1", "Puromycin" = "#9DBDE0", "Quinolone" = "#9A84D6",
    "Rifamycin" = "#A9C98A", "Streptothricin" = "#6DBDE3", "Sulfonamide" = "#C8A24B",
    "Tetracenomycin C" = "#D8BB92", "Tetracycline" = "#E8D98F",
    "Trimethoprim" = "#E39A9A", "Tunicamycin" = "#9AA3A8",
    "Vancomycin" = "#AFC9A0", "beta_lactam" = "#7DA34D", "Others" = "#CFCFCF")
}

fig5_fg_type_labels <- function(x) {
  palette <- fig5_fg_palette()
  dictionary <- setNames(names(palette), tolower(gsub(" ", "_", names(palette), fixed = TRUE)))
  dictionary["macrolide-lincosamide-streptogramin"] <- "MLS"
  dictionary["pleuromutilin_tiamulin"] <- "Pleuromutilin/Tiamulin"
  dictionary["pleuromutilin/tiamulin"] <- "Pleuromutilin/Tiamulin"
  dictionary["unknown"] <- "Others"
  dictionary["glycopeptide"] <- "Vancomycin"
  result <- unname(dictionary[tolower(gsub(" ", "_", x, fixed = TRUE))])
  if (anyNA(result)) stop("Unmapped ARG types: ", paste(unique(x[is.na(result)]), collapse=", "))
  result
}

prepare_fig5_fg <- function(long, meta) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  long_required <- c("Sample", "Subtype", "ARG_type", "Abundance")
  meta_required <- c("Sample", "Setting", "City", "SamplingMonth")
  if (!all(long_required %in% names(long))) stop("long requires: ", paste(long_required, collapse = ", "))
  if (!all(meta_required %in% names(meta))) stop("meta requires: ", paste(meta_required, collapse = ", "))
  long <- as.data.frame(long[, long_required], stringsAsFactors = FALSE)
  meta <- as.data.frame(meta[, meta_required], stringsAsFactors = FALSE)
  long[1:3] <- lapply(long[1:3], as.character)
  meta[] <- lapply(meta, as.character)
  if (anyDuplicated(meta$Sample)) stop("Duplicate Sample in metadata.")
  if (anyDuplicated(long[c("Sample", "Subtype")])) stop("Duplicate Sample and Subtype rows in long input.")
  if (anyNA(meta) || any(vapply(meta, function(x) any(!nzchar(x)), logical(1)))) stop("Metadata keys cannot be missing.")
  if (anyNA(long[1:3]) || any(vapply(long[1:3], function(x) any(!nzchar(x)), logical(1)))) stop("Abundance keys cannot be missing.")
  if (!is.numeric(long$Abundance) || any(!is.na(long$Abundance) & (!is.finite(long$Abundance) | long$Abundance < 0))) {
    stop("Abundance must be numeric, finite and nonnegative, or NA for missing measurements.")
  }
  if (any(!long$Sample %in% meta$Sample)) stop("Some abundance samples have no metadata.")
  classes_per_subtype <- dplyr::distinct(long, Subtype, ARG_type)
  if (anyDuplicated(classes_per_subtype$Subtype)) stop("A Subtype maps to multiple ARG_type values.")
  settings <- c("Hospital", "WWTP", "Community", "Wet market")
  if (any(!meta$Setting %in% settings)) stop("Metadata contains an unsupported Setting.")
  if (any(!grepl("^[0-9]{4}-(0[1-9]|1[0-2])$", meta$SamplingMonth))) stop("SamplingMonth must use YYYY-MM.")
  months <- fig5_fg_months(meta$SamplingMonth)
  meta <- meta[meta$SamplingMonth %in% months, , drop = FALSE]
  long <- long[long$Sample %in% meta$Sample, , drop = FALSE]
  if (!nrow(meta) || !nrow(long)) stop("No abundance samples in the metadata interval.")
  classes <- sort(unique(long$ARG_type))
  cities <- sort(unique(meta$City))
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)
  expected <- dplyr::count(dplyr::distinct(long, Subtype, ARG_type), ARG_type, name = "n_expected_subtypes")
  observed <- long |>
    dplyr::group_by(Sample, ARG_type) |>
    dplyr::summarise(n_subtypes = dplyr::n(), n_missing = sum(is.na(Abundance)),
                     sample_abundance = if (anyNA(Abundance)) NA_real_ else sum(Abundance), .groups = "drop") |>
    dplyr::left_join(expected, by = "ARG_type")
  observed$sample_abundance[observed$n_subtypes != observed$n_expected_subtypes] <- NA_real_
  sample_class <- expand.grid(Sample = meta$Sample, ARG_type = classes, stringsAsFactors = FALSE) |>
    dplyr::left_join(meta, by = "Sample") |>
    dplyr::left_join(observed, by = c("Sample", "ARG_type"))

  summarize_by <- function(keys) {
    sample_class |>
      dplyr::group_by(dplyr::across(dplyr::all_of(keys))) |>
      dplyr::summarise(mean_abundance = safe_mean(sample_abundance),
                       n_samples = dplyr::n(), n_measured = sum(!is.na(sample_abundance)), .groups = "drop")
  }
  complete_summary <- function(grid, observed, keys) {
    out <- dplyr::left_join(grid, observed, by = keys)
    out$n_samples[is.na(out$n_samples)] <- 0L
    out$n_measured[is.na(out$n_measured)] <- 0L
    out$status <- ifelse(out$n_samples == 0, "unsampled",
                         ifelse(is.na(out$mean_abundance), "missing_measurement",
                                ifelse(out$mean_abundance == 0, "observed_zero", "observed_positive")))
    out$ARG_label <- fig5_fg_type_labels(out$ARG_type)
    out
  }
  f_keys <- c("City", "Setting", "SamplingMonth", "ARG_type")
  g_keys <- c("Setting", "SamplingMonth", "ARG_type")
  f <- complete_summary(expand.grid(City = cities, Setting = settings, SamplingMonth = months,
                                    ARG_type = classes, stringsAsFactors = FALSE), summarize_by(f_keys), f_keys)
  city_keys <- c("City", "Setting", "ARG_type")
  f_city <- complete_summary(expand.grid(City = cities, Setting = settings,
                                         ARG_type = classes, stringsAsFactors = FALSE),
                             summarize_by(city_keys), city_keys)
  g <- complete_summary(expand.grid(Setting = settings, SamplingMonth = months,
                                    ARG_type = classes, stringsAsFactors = FALSE), summarize_by(g_keys), g_keys)
  g <- g |>
    dplyr::group_by(Setting, ARG_type) |>
    dplyr::mutate(row_mean = safe_mean(mean_abundance),
                   relative_deviation = ifelse(is.finite(row_mean) & row_mean > 0,
                                              (mean_abundance - row_mean) / row_mean, NA_real_)) |>
    dplyr::ungroup()
  deviations <- g$relative_deviation[is.finite(g$relative_deviation)]
  color_cap <- if (length(deviations)) max(abs(stats::quantile(deviations, c(0.05, 0.95), names = FALSE))) else 0.30
  if (!is.finite(color_cap) || color_cap <= 0) color_cap <- 0.30
  g$fill_value <- pmax(-color_cap, pmin(color_cap, g$relative_deviation))
  positive <- g$mean_abundance[is.finite(g$mean_abundance) & g$mean_abundance > 0]
  size_limits <- if (length(positive) >= 2) stats::quantile(positive, c(0.05, 0.95), names = FALSE) else c(1e-6, 1e-2)
  if (any(!is.finite(size_limits)) || size_limits[1] <= 0 || size_limits[1] >= size_limits[2]) size_limits <- c(1e-6, 1e-2)
  g$size_abundance <- ifelse(g$mean_abundance > 0,
                            pmax(size_limits[1], pmin(size_limits[2], g$mean_abundance)), NA_real_)
  size_breaks <- 10^seq(-6, 0)
  size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]
  if (length(size_breaks) < 3) {
    size_breaks <- 10^pretty(log10(size_limits), n = 4)
    size_breaks <- size_breaks[size_breaks >= size_limits[1] & size_breaks <= size_limits[2]]
  }
  if (!length(size_breaks)) size_breaks <- size_limits
  class_order <- f_city |>
    dplyr::group_by(ARG_label) |>
    dplyr::summarise(total = sum(mean_abundance, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(total))
  f_class_order <- c(setdiff(class_order$ARG_label, "Others"), intersect("Others", class_order$ARG_label))
  g_class_order <- g |>
    dplyr::group_by(ARG_type) |>
    dplyr::summarise(overall_mean = safe_mean(mean_abundance), .groups = "drop") |>
    dplyr::arrange(dplyr::desc(overall_mean))
  counts <- dplyr::distinct(f, City, Setting, SamplingMonth, n_samples)
  list(f_city = as.data.frame(f_city), f = as.data.frame(f), g = as.data.frame(g), sample_class = as.data.frame(sample_class),
       sampling_counts = as.data.frame(counts), months = months, settings = settings, cities = cities,
       f_class_order = f_class_order, g_class_order = g_class_order$ARG_type,
       size_color = list(size_limits = size_limits, size_breaks = size_breaks, color_cap = color_cap,
                         size_transform = "log10", color_statistic = "(monthly_mean - row_mean) / row_mean"))
}

plot_fig5_fg <- function(result) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  f <- result$f
  f_city <- result$f_city
  g <- result$g
  f$SamplingMonth <- factor(f$SamplingMonth, levels = result$months)
  g$SamplingMonth <- factor(g$SamplingMonth, levels = result$months)
  f$Setting <- factor(f$Setting, levels = result$settings)
  g$Setting <- factor(g$Setting, levels = result$settings)
  f$ARG_label <- factor(f$ARG_label, levels = result$f_class_order)
  f_city$ARG_label <- factor(f_city$ARG_label, levels = result$f_class_order)
  f_city$Setting <- factor(f_city$Setting, levels = result$settings)
  f_city$City <- factor(f_city$City, levels = result$cities)
  g$ARG_type <- factor(g$ARG_type, levels = result$g_class_order)
  clean_theme <- ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   axis.text.x = ggplot2::element_text(angle = 50, hjust = 1, vjust = 1),
                   legend.position = "bottom", legend.key = ggplot2::element_blank())
  palette <- fig5_fg_palette()
  city_zeros <- f_city |>
    dplyr::group_by(Setting, City) |>
    dplyr::summarise(all_zero = all(status == "observed_zero"), .groups = "drop")
  city_zeros <- city_zeros[city_zeros$all_zero, , drop = FALSE]
  f_main <- ggplot2::ggplot(f_city, ggplot2::aes(City, mean_abundance, fill = ARG_label)) +
    ggplot2::geom_col(width = 0.74, na.rm = TRUE) +
    ggplot2::geom_point(data = city_zeros, ggplot2::aes(City, 0), inherit.aes = FALSE,
                        shape = 4, size = 1.5, color = "grey35") +
    ggplot2::facet_wrap(~Setting, ncol = 2, drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.015, 0.04))) +
    ggplot2::scale_fill_manual(values = palette, limits = result$f_class_order, drop = FALSE, name = "ARG class") +
    ggplot2::labs(x = "City", y = "Mean Tier I abundance (copies/cell)",
                  caption = paste0("Arithmetic means across individual samples within each city and setting, November 2024-December 2025.\n",
                                   "Blank city-setting: no usable abundance data. Cross: sampled, all class means zero.")) +
    clean_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 0),
                                plot.caption = ggplot2::element_text(size = 8, hjust = 0)) +
    ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4))
  f_plots <- setNames(lapply(result$cities, function(city) {
    dat <- f[f$City == city, , drop = FALSE]
    zeros <- dat |>
      dplyr::group_by(Setting, SamplingMonth) |>
      dplyr::summarise(all_zero = all(status == "observed_zero"), .groups = "drop")
    zeros <- zeros[zeros$all_zero, , drop = FALSE]
    ggplot2::ggplot(dat, ggplot2::aes(SamplingMonth, mean_abundance, fill = ARG_label)) +
      ggplot2::geom_col(width = 0.74, na.rm = TRUE) +
      ggplot2::geom_point(data = zeros, ggplot2::aes(SamplingMonth, 0), inherit.aes = FALSE,
                          shape = 4, size = 1.5, color = "grey35") +
      ggplot2::facet_wrap(~Setting, ncol = 2, drop = FALSE) +
      ggplot2::scale_x_discrete(drop = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.015, 0.04))) +
      ggplot2::scale_fill_manual(values = palette, limits = result$f_class_order, drop = FALSE, name = "ARG class") +
      ggplot2::labs(title = city, x = NULL, y = "Mean Tier I abundance (copies/cell)",
                    caption = "Blank month: no usable abundance data. Cross: sampled, all class means zero.") +
      clean_theme + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 8, angle = 50, hjust = 1),
                                  plot.caption = ggplot2::element_text(size = 8)) +
      ggplot2::guides(fill = ggplot2::guide_legend(ncol = 4))
  }), result$cities)
  positive <- g[!is.na(g$size_abundance), , drop = FALSE]
  zeros <- g[g$status == "observed_zero", , drop = FALSE]
  g_plot <- ggplot2::ggplot(g, ggplot2::aes(SamplingMonth, ARG_type)) +
    ggplot2::geom_blank() +
    ggplot2::geom_point(data = positive, ggplot2::aes(size = size_abundance, fill = fill_value),
                        shape = 21, color = "grey68", stroke = 0.22, na.rm = TRUE) +
    ggplot2::geom_point(data = zeros, shape = 4, size = 1.3, color = "grey45") +
    ggplot2::facet_grid(rows = ggplot2::vars(Setting), scales = "free_y", space = "free_y", switch = "y", drop = FALSE) +
    ggplot2::scale_x_discrete(drop = FALSE) +
    ggplot2::scale_y_discrete(position = "right", labels = fig5_fg_type_labels, drop = FALSE) +
    ggplot2::scale_fill_gradient2(low = "#4E79A7", mid = "#F7F7F7", high = "#D7301F", midpoint = 0,
                                  limits = c(-1, 1) * result$size_color$color_cap, na.value = "white",
                                  name = "Relative deviation\nfrom row mean") +
    ggplot2::scale_size_continuous(trans = "log10", range = c(0.5, 4.6),
                                   limits = result$size_color$size_limits,
                                   breaks = result$size_color$size_breaks,
                                   labels = function(x) format(x, scientific = TRUE, digits = 1),
                                   name = "Mean abundance\n(copies/cell; log10 size)") +
    ggplot2::labs(x = NULL, y = NULL,
                  caption = paste0("Bubble size: abundance clipped at global positive 5th/95th percentiles, then log10 scaled.\n",
                                   "Cross: observed zero; blank: no usable abundance data. Color = (monthly mean - row mean)/row mean;\n",
                                   "symmetric 5th/95th percentile color cap. Row mean includes sampled zeros.")) +
    clean_theme +
    ggplot2::theme(axis.text.y = ggplot2::element_text(size = 8),
                   strip.placement = "outside", strip.text.y.left = ggplot2::element_text(angle = 90),
                   panel.spacing.y = grid::unit(0.18, "cm"),
                   plot.caption = ggplot2::element_text(size = 7.5, hjust = 0),
                   legend.title = ggplot2::element_text(size = 9), legend.text = ggplot2::element_text(size = 8)) +
    ggplot2::guides(fill = ggplot2::guide_colorbar(order = 1), size = ggplot2::guide_legend(order = 2))
  list(f = f_main, f_monthly = f_plots, g = g_plot)
}

run_fig5_fg <- function(long, meta, output_dir, export_monthly_si = FALSE) {
  result <- prepare_fig5_fg(long, meta)
  plots <- plot_fig5_fg(result)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tables <- list(Fig5f_city_setting_class = result$f_city,
                 SI_monthly_city_setting_class = result$f,
                 Fig5g_monthly_setting_class = result$g,
                 Fig5fg_sample_class = result$sample_class,
                 Fig5fg_sampling_counts = result$sampling_counts)
  for (name in names(tables)) utils::write.csv(tables[[name]], file.path(output_dir, paste0(name, ".csv")),
                                              row.names = FALSE, na = "", fileEncoding = "UTF-8")
  utils::write.csv(data.frame(size_lower = result$size_color$size_limits[1],
                              size_upper = result$size_color$size_limits[2],
                              color_cap = result$size_color$color_cap, size_transform = "log10",
                              color_statistic = result$size_color$color_statistic),
                   file.path(output_dir, "Fig5g_display_parameters.csv"), row.names = FALSE)
  saveRDS(result, file.path(output_dir, "Fig5fg_prepared.rds"))
  ggplot2::ggsave(file.path(output_dir, "Fig5f_TierI_city_composition.pdf"), plots$f, width = 11, height = 8, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "Fig5f_TierI_city_composition.png"), plots$f, width = 11, height = 8, dpi = 300, bg = "white")
  if (isTRUE(export_monthly_si)) {
    f_dir <- file.path(output_dir, "SI_monthly_by_city")
    dir.create(f_dir, recursive = TRUE, showWarnings = FALSE)
    for (city in names(plots$f_monthly)) {
      stem <- file.path(f_dir, paste0("SI_monthly_", gsub("[^[:alnum:]_-]", "_", city)))
      ggplot2::ggsave(paste0(stem, ".pdf"), plots$f_monthly[[city]], width = 9, height = 7, bg = "white")
      ggplot2::ggsave(paste0(stem, ".png"), plots$f_monthly[[city]], width = 9, height = 7, dpi = 300, bg = "white")
    }
  }
  ggplot2::ggsave(file.path(output_dir, "Fig5g_TierI_monthly_bubble.pdf"), plots$g, width = 9.4, height = 10.5, bg = "white")
  ggplot2::ggsave(file.path(output_dir, "Fig5g_TierI_monthly_bubble.png"), plots$g, width = 9.4, height = 10.5, dpi = 300, bg = "white")
  invisible(list(data = result, plots = plots))
}

run_fig5 <- function(data_root, metadata_file, annotation_dir,
                     output_dir = "Fig5_results", panels = letters[2:7],
                     include_richness = FALSE, prepared_data = NULL,
                     export_monthly_si = FALSE) {
  if (!length(panels) || any(!panels %in% letters[2:7]))
    stop("panels must contain one or more of b, c, d, e, f, g.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  tier <- read_fig5_bc_inputs(annotation_dir)
  needs_samples <- any(panels %in% c("d", "e", "f", "g")) || include_richness
  result <- list()
  if (needs_samples) {
    cohort <- if (is.null(prepared_data)) {
      read_fig5_data(data_root, metadata_file, file.path(output_dir, "input_audit"))
    } else {
      checked <- prepare_fig5_data(prepared_data$long, prepared_data$meta)
      checked
    }
    if (!setequal(cohort$subtypes, tier$tier_i$arg_subtype))
      stop("The abundance matrix is not the exact Tier I roster used in b/c.")
    classes <- sub("__.*$", "", tier$tier_i$arg_subtype)
    if (any(classes != tier$tier_i$arg_type))
      stop("ARG-type prefix and the Tier I roster annotation disagree.")
    result$cohort_counts <- data.frame(samples=nrow(cohort$meta),
      sites=length(unique(cohort$meta$PhysicalSite)),
      cities=length(unique(cohort$meta$City)),
      months=length(unique(cohort$meta$SamplingMonth)), subtypes=length(cohort$subtypes))
    utils::write.csv(result$cohort_counts, file.path(output_dir,"cohort_counts.csv"), row.names=FALSE)
  }
  if (any(panels %in% c("b", "c")))
    result$bc <- run_fig5_bc(annotation_dir, file.path(output_dir,"bc"))
  if ("d" %in% panels)
    result$d <- run_fig5d(cohort$long, cohort$meta, file.path(output_dir,"d_PCoA"),
      ellipse_level=.68, permutations=999L, seed=123L, zero_policy="error",
      compute_dispersion=TRUE)
  if ("e" %in% panels)
    result$e <- run_fig5e(cohort$samples, file.path(output_dir,"e_LMM"))
  if (any(panels %in% c("f", "g")))
    result$fg <- run_fig5_fg(cohort$long, cohort$meta, file.path(output_dir,"fg"),
      export_monthly_si=export_monthly_si)
  if (isTRUE(include_richness))
    result$richness <- run_tier1_richness(cohort$samples, file.path(output_dir,"S4a_richness_NB"))
  writeLines(capture.output(sessionInfo()), file.path(output_dir,"R_session_info.txt"))
  invisible(result)
}

fig5_command_line <- function() {
  args <- commandArgs(trailingOnly=TRUE)
  if (!length(args)) {
    message("Source this file and call run_fig5(), or use: Rscript Fig5_b_to_g.R DATA_ROOT METADATA ANNOTATION_DIR OUTPUT_DIR [PANELS] [--richness]")
    return(invisible(NULL))
  }
  if (length(args)<4L) stop("Four paths are required: DATA_ROOT METADATA ANNOTATION_DIR OUTPUT_DIR")
  panel_arg <- args[seq_along(args)>4L & !grepl("^--",args)]
  panels <- if (length(panel_arg)) strsplit(panel_arg[1], "", fixed=TRUE)[[1]] else letters[2:7]
  run_fig5(args[1],args[2],args[3],args[4],panels=panels,
    include_richness="--richness" %in% args, export_monthly_si="--monthly-si" %in% args)
}

if (sys.nframe() == 0L) fig5_command_line()
