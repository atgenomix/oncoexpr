#' @title prepare data
#' @description precessing the gtf file into expression profile
#' @param entirepath entire path for gtf files
#' @param groups_list the group name of samples
#' @return expression profile
#' @export
input_and_prep_file <-  function(entirepath , groups_list ){
  groups_path_variable <- paste0(groups_list, c("_gtf_file"))
  groups_data_variable <- paste0(groups_list, c("_gtf_data"))
  groups_expr_variable <- paste0(groups_list, c("_rna_expr"))
  print(groups_expr_variable)
  
  names(groups_path_variable) <- groups_list
  names(groups_data_variable) <- groups_list
  names(groups_expr_variable) <- groups_list
  
  print(groups_expr_variable)
  lapply(seq_len(length(groups_list)), function(i) assign(groups_data_variable[[groups_list[i]]], read_tsv(get(entirepath[i]), comment = "#", col_names = FALSE), envir=.GlobalEnv) )
  lapply(seq_len(length(groups_list)), function(i) assign(groups_expr_variable[[groups_list[i]]], gtf2expr(get(groups_data_variable[[groups_list[i]]])), envir=.GlobalEnv) )
  
  print(dim(groups_expr_variable[[groups_list[1]]]))
  print(dim(groups_expr_variable[[groups_list[2]]]))
  list_ <- list("path" = groups_path_variable, "gtf" = groups_data_variable, "expr" = groups_expr_variable)
  return(list_)
  
}

#' @title prepare data
#' @description precessing the gtf file into expression profile
#' @param entirepath entire path for gtf files
#' @param groups_list the group name of samples
#' @return expression profile
#' @export

input_and_prep_file_directly <-  function(entirepath, groups_list){
  groups_path_variable <- paste0(groups_list, c("_gtf_file"))
  groups_data_variable <- paste0(groups_list, c("_gtf_data"))
  groups_expr_variable <- paste0(groups_list, c("_rna_expr"))
  print(groups_expr_variable)
  names(groups_path_variable) <- groups_list
  names(groups_data_variable) <- groups_list
  names(groups_expr_variable) <- groups_list
  
  print(groups_expr_variable)
  assign(groups_path_variable[[groups_list[1]]], entirepath[1], envir=.GlobalEnv)
  assign(groups_path_variable[[groups_list[2]]], entirepath[2], envir=.GlobalEnv)
  
  assign(groups_data_variable[[groups_list[1]]], read_tsv(get(groups_path_variable[[groups_list[1]]]), comment = "#", col_names = FALSE), envir=.GlobalEnv)
  assign(groups_data_variable[[groups_list[2]]], read_tsv(get(groups_path_variable[[groups_list[2]]]), comment = "#", col_names = FALSE), envir=.GlobalEnv)
  
  
  assign(groups_expr_variable[[groups_list[1]]], gtf2expr(get(groups_data_variable[[groups_list[1]]])), envir=.GlobalEnv)
  assign(groups_expr_variable[[groups_list[2]]], gtf2expr(get(groups_data_variable[[groups_list[2]]])), envir=.GlobalEnv)
  print(dim(groups_expr_variable[[groups_list[1]]]))
  print(dim(groups_expr_variable[[groups_list[2]]]))
  list_ <- list("path" = groups_path_variable, "gtf" = groups_data_variable, "expr" = groups_expr_variable)
  return(list_)
  
}




#' @title gtf to expr profile
#' @description 輸入RNA-seq pipline生成的 .gtf file
#' @param gtf_data gtf_data的絕對路徑
#' @return 經過處理後的表現量表格，包括Cov, FPKM, TPM三種量化指標
#' @export
gtf2expr <- function(gtf_data) {
  colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  
  gtf_data <- gtf_data %>%
    filter(feature == "transcript") %>%
    mutate(
      gene_id = ifelse(grepl('gene_id', attribute), sub('.*gene_id "([^"]+)".*', '\\1', attribute), NA),
      Gene = ifelse(grepl('gene_name', attribute), sub('.*gene_name "([^"]+)".*', '\\1', attribute), NA),
      cov = ifelse(grepl('cov', attribute), as.numeric(sub('.*cov "([^"]+)".*', '\\1', attribute)), NA),
      FPKM = ifelse(grepl('FPKM', attribute), as.numeric(sub('.*FPKM "([^"]+)".*', '\\1', attribute)), NA),
      TPM = ifelse(grepl('TPM', attribute), as.numeric(sub('.*TPM "([^"]+)".*', '\\1', attribute)), NA)
    )
  
  gene_expression <- gtf_data %>%
    group_by(Gene) %>%
    summarise(
      max_cov = ifelse(all(is.na(cov)), NA, max(cov, na.rm = TRUE)),
      max_FPKM = ifelse(all(is.na(FPKM)), NA, max(FPKM, na.rm = TRUE)),
      max_TPM = ifelse(all(is.na(TPM)), NA, max(TPM, na.rm = TRUE))
    )
  
  return(gene_expression)
}

#' @title prepare for expression pattern
#' @description compare two samples' expression pattern
#' @param group1 first sample 
#' @param group2 sencond sample
#' @param groups_name sample type name e.g. tissue & ctc
#' @return a list with dataframe and specific gene list 
#' @export

Prep4ExprComp <- function(group1, group2, groups_name = c("tissue", "ctc")){
  group2_matched <- group2[match(group1$"Gene", group2$"Gene"),]
  data_ <- data.frame(group1$max_TPM, group2_matched$max_TPM)
  rownames(data_) <- group1$Gene
  print(colnames(data_))
  rm_noise_ <- data_[,1] > 1 | data_[,2] > 1
  
  only_group1_high_gene_profile <- data_[data_[,"group1.max_TPM"] < 100 & data_[,"group2_matched.max_TPM"] < 0.001,]
  only_group2_high_gene_profile <- data_[data_[,"group1.max_TPM"] < 0.001 & data_[,"group2_matched.max_TPM"] < 100,]
  
  special_gene_list <- c(rownames(only_group1_high_gene_profile), rownames(only_group2_high_gene_profile))
  
  data_v2 <- as.data.frame(data_[rm_noise_,])
  data_v2 <- data_v2 + 0.0001
  colnames(data_v2) <- groups_name
  list_ <- list("special_gene_list" = special_gene_list, "result" = data_v2)
  return(list_)
  
}

#' @title target gene expression
#' @description 針對目標基因回傳兩個樣本的基因表現量
#' @param selectMode select or exclude to gene list
#' @param geneList_ interesting genes
#' @param groups_list_ term of group
#' @param save_path_ save path
#' @param save_filename_  save name
#' @param expr_profile_ target gene expression profile
#' @return the expression of genes which you are interesting
#' @export

target_exprofile <- function(selectMode = TRUE, geneList_, groups_list_ , save_path_ = NULL, save_filename_ = "target_gene_expr.csv", expr_profile_){
  
  target_exprprfile_ <-  expr_profile_ %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  
  colnames(target_exprprfile_) <- groups_list_
  
  if(!is.null(save_path_)){
    #collect from instance to local 
    write.csv(target_exprprfile_, paste0(save_path_, save_filename_))
    print(paste("save at the path:", paste0(save_path_, save_filename_)))
    
  }else(
    
    print("do not save file")
    
  )
  
  return(target_exprprfile_ )
  
}

#' @title expression pattern analysis
#' @description expression pattern analysis for two samples by scatter plot 
#' @param data_v2 dataframe derived from Prep4ExprComp function
#' @param special_gene_list gene list derived from Prep4ExprComp function
#' @param save_path_ where you want to save file
#' @param save_filename_ plot filename
#' @param cutoff_tpm what target genes would be shown in red color
#' @param x_col x axis tile of scatter plot
#' @param y_col y axis tile of scatter plot
#' @return scatter plot for expression pattern analysis
#' @export

expr_pattern <- function(data_v2, special_gene_list, save_path_ = NULL, save_filename_ = "scatter_plot.png", cutoff_tpm = 7, x_col = "tissue", y_col = "ctc"){
  data_v2$distance_to_diagonal <- abs(log2(data_v2[,x_col]) - log2(data_v2[,y_col]))
  color_ <- ifelse(data_v2$distance_to_diagonal>cutoff_tpm, "outlier", "consistent")
  data_v2$color_ <- color_
  genes_beyond <- data_v2[data_v2$distance_to_diagonal>cutoff_tpm, ]
  genes_beyond_v2 <- genes_beyond[!is.element(rownames(genes_beyond), special_gene_list), ]
  genes_beyond_v2$label_ <- rownames(genes_beyond_v2)
  scatter_plot <- ggplot(data_v2, aes(x = .data[[x_col]], y = .data[[y_col]]))+
    geom_point(aes(color = color_), size=1.5, alpha=0.4) +  
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    geom_text_repel(data = genes_beyond_v2, aes(label = .data[["label_"]]), size = 1.5, force = 1, max.overlaps = Inf, segment.size=0.1) +
    labs(x = x_col, y =  y_col ) +
    scale_x_continuous(limits = c(NA, 200000),expand = c(0, 0), trans = "log10", 
                       labels = scales::number_format(accuracy = 0.0001)) +
    scale_y_continuous(limits = c(NA, 200000),expand = c(0, 0), trans = "log10", 
                       labels = scales::number_format(accuracy = 0.0001)) +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = c("blue", "red")) 
  save_parameter <- paste0(save_path_, save_filename_)
  ggsave(save_parameter, plot = scatter_plot , width = 28, height = 28, units="cm",dpi = 300)
  return(scatter_plot)
  
}

#' @title generate pseudo colData with random classification
#' @description generate random classification for colData
#' @param normcount normCount table to get colnames (sample id)
#' @param genecol to remove gene column
#' @return pseudo colData
#' @export

generate_colData_random <- function(normcount = normCount, genecol = "gene_id"){
    rownames(normcount) <- normcount[, genecol]  # 假設第一列是基因名稱
    normcount <- normcount[, -which(colnames(normcount) == genecol)]  # 去除基因名稱列，確保只有數值矩陣
    classification <- sample(c("T","NT"), ncol(normcount), replace = TRUE)
    colData <- data.frame("mainCode" =colnames(normcount) , "subCode" = classification)
    rownames(colData) <- colnames(normcount)
    return(colData)
}


#' @title generate colData for RNAseqShinyApp
#' @description generate classification for colData
#' @param normcount normCount table to get colnames (sample id)
#' @param grouplist group list
#' @param genecol to remove gene column
#' @return colData
#' @export

generate_colData <- function(normcount = normCount, grouplist  , genecol = "GeneSymbol" ){
    rownames(normcount) <- normcount[, genecol]  # 假設第一列是基因名稱
    normcount <- normcount[, -which(colnames(normcount) == genecol)]  # 去除基因名稱列，確保只有數值矩陣
    classification <-  grouplist
    colData <- data.frame("mainCode" = colnames(normcount) , "subCode" = classification)
    rownames(colData) <- colnames(normcount)
    return(colData)
}


#' @title time filter for filename selection
#' @description extract the latest file from filename list with time information
#' @param filename_list list of filenames
#' @return filename list with time information
#' @export

get_latest_file_group_df <- function(filename_list) {

  filename_split <- strsplit(filename_list, "_")
  
  extract_time <- function(parts) {
    time_str <- parts[grepl("^\\d{8}t\\d{6}z$", parts)]
    if (length(time_str) == 1) return(time_str) else return(NA)
  }
  
  time_strs <- sapply(filename_split, extract_time)
  time_posix <- as.POSIXct(time_strs, format = "%Y%m%dt%H%M%Sz", tz = "UTC")
  
  df <- data.frame(
    file = filename_list,
    time_str = time_strs,
    time = time_posix,
    stringsAsFactors = FALSE
  )
  
  latest_time <- max(df$time, na.rm = TRUE)
  print("latest_time")
  print(latest_time)
  df$is_latest <- !is.na(df$time) & (df$time == latest_time)
  
  #df <- df[order(df$time, decreasing = TRUE, na.last = TRUE), ]
  
  return(df)
}


#' @title limma analysis with MAE input
#' @description differential expression by limma
#' @param mae multiple assays experiment
#' @param assayName assay
#' @param subCodeCol which col for comparasion
#' @param coef coefficient for limma
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @return limma analysis table
#' @export

LimmaMAE <- function(mae,
                     assayName = "RNAseq",
                     # 假設分組資訊在 colData 中的 subCode 欄位
                     subCodeCol = "subCode",
                     coef = 2, # 要比較的係數 (group 的哪一層)
                     pval_cut = 0.05,
                     lfc_cut = 1,
                     useAdjP = FALSE) {
  rna_counts <- assays(mae[[assayName]])$max_TPM
  sample_info <- colData(mae[[assayName]])
  print(sample_info)
  group <- factor(sample_info[[subCodeCol]])
  print(group)

  dge <- DGEList(counts = rna_counts, group = group)
  dge <- calcNormFactors(dge)

  design <- model.matrix(~group)

  v <- voom(dge, design)

  fit <- lmFit(v, design)
  fit <- eBayes(fit)

  dt <- decideTests(fit, p.value = pval_cut, lfc = lfc_cut)
  sumDT <- summary(dt)

  topTab <- topTable(fit, coef = coef, number = Inf)

  return(list(
    fit            = fit,
    topTable       = topTab,
    decideTests    = dt,
    upDownSummary  = sumDT
  ))
}


#' @title se to longformat
#' @description se from MAE package transfered to long format dataframe
#' @param se_ se object 
#' @param assay_name selecting assay from se object
#' @param selected_row_ select gene
#' @param selected_col_ select sample
#' @param selected_colData_ filter colData column
#' @return pivot long format from se object
#' @export
#' 
se2longformat <- function(se_ = exp_brca, 
                          assay_name = "fpkm_unstrand", 
                          selected_row_ = c(1:5), 
                          selected_col_ = c("TCGA-5L-AAT0-01A-12R-A41B-07", "TCGA-A2-A04U-01A-11R-A115-07"), 
                          selected_colData_ = "tissue_type") {
  
  assay_mat <- assay(se_, assay_name)
  
  if (is.null(colnames(assay_mat))) {
    stop("Assay matrix does not have column names!")
  }
  
  new_colData <- colData(se_)
  rownames(new_colData) <- colnames(se_)
  
  se_new <- SummarizedExperiment(
    assays = list(data = assay_mat),
    colData = new_colData,
    rowRanges = rowRanges(se_)
  )
  
  mae_new <- MultiAssayExperiment(setNames(list(se_new), assay_name))
  colData(mae_new) <- colData(se_new)
  if (is.null(selected_row_)) {
    selected_row_ <- rownames(se_new)
  } else if (is.numeric(selected_row_)) {
    selected_row_ <- rownames(se_new)[selected_row_]
  }
  if (is.null(selected_col_)) {
    selected_col_ <- seq_len(ncol(se_new))
  }
  
  subset_mae <- mae_new[selected_row_, selected_col_, assay_name]
  
  se2longformat_ <- longFormat(subset_mae, colDataCols = selected_colData_)
  
  return(se2longformat_)
}



#' @title generate expression matrix with multiple samples
#' @description convert long to wide for SQL dataframe from spark database
#' @param data fileterd data from spark database 
#' @param mainCode_ mainCode selection
#' @param subCode_ subCode selection
#' @return a expression matrix with gene row and sample column
#' @export
#' 
convert_long_to_wide <- function(data, mainCode_, subCode_) {
  data <- data %>%
    mutate(sample_name = paste(!!sym(mainCode_), !!sym(subCode_), sep = "_"))
  
  rna_expr <- data %>%
    select(Gene, sample_name, value) %>%
    pivot_wider(names_from = sample_name, values_from = value, values_fill = NULL)
  
  return(rna_expr)
}

#' @title sample information table
#' @description generate sample information table for MAE colData 
#' @param data fileterd data from spark database 
#' @param mainCode_ mainCode selection
#' @param subCode_ subCode selection
#' @return a table about sample codes
#' @export
#' 
generate_sample_info <- function(data, mainCode_, subCode_) {
  sample_info <- data %>%
    distinct(!!sym(mainCode_), !!sym(subCode_)) %>%
    mutate(sample_id = paste0("S", sprintf("%03d", row_number())))
  
  sample_info_table <- sample_info %>%
    select(sample_id, !!sym(mainCode_), !!sym(subCode_))
  
  return(sample_info_table)
}
#' Extract the three key tables from a SummarizedExperiment
#'
#' Fast helper that returns the expression matrix, gene-level metadata and
#' sample-level metadata from a single `SummarizedExperiment` object.
#'
#' @param se A `SummarizedExperiment` (or `RangedSummarizedExperiment`).
#' @param assay_name `character(1)` Name of the assay to extract.  
#'   If `NULL` (default) the first assay is used.
#'
#' @return A named `list` with three elements  
#' * **`exp_data`** – expression matrix (`matrix` or `DelayedMatrix`)  
#' * **`gene_info`** – gene-level metadata (`data.frame`)  
#' * **`sample_info`** – sample-level metadata (`data.frame`)
#'
#' @examples
#' \dontrun{
#' se <- readRDS("query_TCGA_BRCA.rds")
#' tbls <- query_se_tables(se, assay_name = "tpm_unstrand")
#' str(tbls)
#' }
#' @export
query_se_tables <- function(se, assay_name = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (is.null(assay_name)) {
    assay_name <- names(assays(se))[1]
  } else if (!assay_name %in% names(assays(se))) {
    stop("`assay_name` not found in this SummarizedExperiment: ",
         paste(names(assays(se)), collapse = ", "))
  }

  expr_mat    <- assay(se, assay_name)
  gene_info   <- as.data.frame(rowData(se))
  sample_info <- as.data.frame(colData(se))

  list(
    exp_data    = expr_mat,
    gene_info   = gene_info,
    sample_info = sample_info
  )
}


#' Filter gene metadata by biotype (or any column)
#'
#' Convenience wrapper around **dplyr** for sub-setting a gene annotation
#' table and returning only the desired columns.
#'
#' @param gene_info A `data.frame` – typically `rowData()` converted to a
#'   data frame.
#' @param gene_biotype Character vector of biotypes to keep (default:
#'   `"protein_coding"`).  Use `NULL` to skip biotype filtering.
#' @param gene_id_col Column name that contains the stable gene ID.
#' @param gene_name_col Column name that contains the gene symbol.
#' @param gene_type_col Column name that stores the biotype information.
#'
#' @return A `data.frame` with three columns: `gene_id`, `gene_name`,
#'   `gene_type`.
#'
#' @examples
#' \dontrun{
#' se  <- readRDS("query_TCGA_BRCA.rds")
#' gi  <- query_se_tables(se)$gene_info
#' pcs <- geneinfo_filter(gi, gene_biotype = "protein_coding")
#' head(pcs)
#' }
#' @export
geneinfo_filter <- function(gene_info,
                            gene_biotype  = "protein_coding",
                            gene_id_col   = "gene_id",
                            gene_name_col = "gene_name",
                            gene_type_col = "gene_type") {

  requireNamespace("dplyr")
  requireNamespace("rlang")

  out <- gene_info |>
    dplyr::filter(
      if (!is.null(gene_biotype))
        .data[[gene_type_col]] %in% gene_biotype
      else
        TRUE
    ) |>
    dplyr::transmute(
      gene_id   = .data[[gene_id_col]],
      gene_name = .data[[gene_name_col]],
      gene_type = .data[[gene_type_col]]
    )

  out
}

#' Sub-set an expression matrix to a selected gene list
#'
#' Aligns a gene list produced by \code{\link{geneinfo_filter}} (or any
#' similar data frame) with an expression matrix returned by
#' \code{\link{query_se_tables}}, optionally replacing the rownames by gene
#' symbols and guaranteeing matrix dimensions.
#'
#' @param tables A list created by \code{\link{query_se_tables}}.
#' @param genes_df A `data.frame` containing at least the `gene_id` column.
#' @param id_col Name of the gene-ID column in `genes_df`.
#' @param symbol_col Name of the symbol column in `genes_df`.
#' @param use_symbol Logical; if `TRUE` (default) replace rownames by gene
#'   symbols (duplicates are made unique via `make.unique()`).
#' @param drop Passed to the matrix sub-setting (`drop = FALSE` keeps
#'   2-dimensional structure even if only one row is selected).
#'
#' @return A list with  
#' * **`expr_mat`** – sub-set expression matrix  
#' * **`gene_info`** – the supplied `genes_df` (for convenience)  
#' * **`sample_info`** – carried over from `tables`
#'
#' @examples
#' \dontrun{
#' se  <- readRDS("query_TCGA_BRCA.rds")
#' tbl <- query_se_tables(se, "tpm_unstrand")
#'
#' pcs <- geneinfo_filter(tbl$gene_info)      # protein-coding by default
#' res <- subset_expr(tbl, pcs)
#'
#' dim(res$expr_mat)          # genes × samples
#' head(rownames(res$expr_mat))
#' }
#' @export
subset_expr <- function(tables,
                        genes_df,
                        id_col      = "gene_id",
                        symbol_col  = "gene_name",
                        use_symbol  = TRUE,
                        drop        = FALSE) {

  idx <- match(genes_df[[id_col]], rownames(tables$exp_data))
  idx <- idx[!is.na(idx)]

  if (!length(idx))
    stop("None of the selected genes are present in `tables$exp_data`.")

  expr_mat <- tables$exp_data[idx, , drop = drop]

  if (use_symbol) {
    sym <- genes_df[[symbol_col]][match(rownames(expr_mat), genes_df[[id_col]])]
    rownames(expr_mat) <- make.unique(sym)
  }

  list(
    expr_mat    = expr_mat,
    gene_info   = genes_df,
    sample_info = tables$sample_info
  )
}
