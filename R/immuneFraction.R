
#' @title Immune Cell Fraction Prediction
#' @description predict the composition of immune cell types for your samples
#' @param group1_ sample 1
#' @param group2_ sample 2
#' @param groups_list sample types
#' @param method quantiseq, xCell, CIBERSORT
#' @param tumorStatus if your sample is tumor TRUE or FALSE
#' @param save_path what the path you want to save the files
#' @return plot for immune cell faction of your samples
#' @export
ImmuneFractionPlot <- function(group1_, group2_, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL){
  
  tryCatch({  
    print("ensuring the dataframe be correct")
    if (!all(c("Gene", "max_TPM") %in% colnames(group1_)) ||
        !all(c("Gene", "max_TPM") %in% colnames(group2_))) {
      stop("Input data frames must contain 'Gene' and 'max_TPM' columns")
    }
    # 2. 匹配基因並創建數據框
    print("match genes position")
    group2_matched <- group2_[match(group1_$"Gene", group2_$"Gene"), ]
    
    # 3. 創建並處理合併後的數據框
    print("creat merged dataframe")
    df <- data.frame(
      Gene = group1_$"Gene",
      group1 = group1_$"max_TPM",
      group2 = group2_matched$"max_TPM",
      stringsAsFactors = FALSE
    )
    # 4. 移除 NA 值的行
    print("remove na value")
    df <- df[complete.cases(df), ]
    
    print("set colnames")        
    colnames(df) <- c("Gene", groups_list)
    
    # 6. 準備用於 deconvolution 的矩陣
    print("transfer Gene column to be rownames")
    expression_matrix <- df %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    
    print(head(expression_matrix))
    # 7. 運行 deconvolution
    res_quantiseq <- immunedeconv::deconvolute(
      expression_matrix,
      method = method_,
      tumor = tumorStatus 
    )
    
    
  }, error = function(e) {
    message("Error in immuneFractionPlot: ", e$message)
    # 診斷信息
    message("\nDiagnostic information:")
    message("group1 structure:")
    print(str(group1_))
    message("\ngroup2 structure:")
    print(str(group2_))
    message("\ngroups_list:")
    print(groups_list)
    return(NULL)
  })
  my_colors <- c(
    # T細胞
    "T cell CD4+ (non-regulatory)" = "#D62728",     # 深紅色
    "T cell CD8+" = "#FF7F0E",                      # 橙色
    "T cell regulatory (Tregs)" = "#FF9896",        # 浅红色
    
    # Macrophage M1 和 M2
    "Macrophage M1" = "#2CA02C",                    # 綠色
    "Macrophage M2" = "#98DF8A",                    # 淺綠色
    # NK cell 和 B cell
    "NK cell" = "#1F77B4",                          # 藍色
    "B cell" = "#AEC7E8",                           # 淺藍色
    # Neutrophil 和 Monocyte
    "Neutrophil" = "#9467BD",                       # 深紫色
    "Monocyte" = "#8C564B",                         # 棕色
    # Myeloid dendritic cell
    "Myeloid dendritic cell" = "#BC80BD",           # 淺紫色
    # 未分類细胞
    "uncharacterized cell" = "#D0E1F9"              # 淺粉藍色
  )
  
  p <-  res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    #scale_fill_manual(values = my_colors)+
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  return(p)
  ggsave(past0(method_, "_", "immuneFractionPlot.jpg"), path = save_path ,plot = p, device = "png", width = 35, height = 25, units = "cm", dpi = 300)
  
  
}

#' @title Immune Cell Fraction Prediction
#' @description predict the composition of immune cell types for your samples
#' @param mae multiple assays experiments
#' @param groups_list sample types
#' @param method quantiseq, xCell, CIBERSORT
#' @param tumorStatus if your sample is tumor TRUE or FALSE
#' @param save_path what the path you want to save the files
#' @return plot for immune cell faction of your samples
#' @export
ImmuneFractionPlot_mae <- function(mae, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL){
  
  rna_counts  <- assays(mae[["RNAseq"]])$max_TPM
  # 6. 準備用於 deconvolution 的矩陣
  print("transfer Gene column to be rownames")
  
  expression_matrix <- rna_counts %>%
    #tibble::column_to_rownames("Gene") %>% #rownames 要是gene symbol
    as.matrix()
  
  print(head(expression_matrix))
  # 7. 運行 deconvolution
  res_quantiseq <- immunedeconv::deconvolute(
    expression_matrix,
    method = method_,
    tumor = tumorStatus 
  )
  
  
  
  my_colors <- c(
    # T細胞
    "T cell CD4+ (non-regulatory)" = "#D62728",     # 深紅色
    "T cell CD8+" = "#FF7F0E",                      # 橙色
    "T cell regulatory (Tregs)" = "#FF9896",        # 浅红色
    
    # Macrophage M1 和 M2
    "Macrophage M1" = "#2CA02C",                    # 綠色
    "Macrophage M2" = "#98DF8A",                    # 淺綠色
    # NK cell 和 B cell
    "NK cell" = "#1F77B4",                          # 藍色
    "B cell" = "#AEC7E8",                           # 淺藍色
    # Neutrophil 和 Monocyte
    "Neutrophil" = "#9467BD",                       # 深紫色
    "Monocyte" = "#8C564B",                         # 棕色
    # Myeloid dendritic cell
    "Myeloid dendritic cell" = "#BC80BD",           # 淺紫色
    # 未分類细胞
    "uncharacterized cell" = "#D0E1F9"              # 淺粉藍色
  )
  
  p <-  res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    #scale_fill_manual(values = my_colors)+
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  return(p)
  #ggsave(past0(method_, "_", "immuneFractionPlot.jpg"), path = save_path ,plot = p, device = "png", width = 35, height = 25, units = "cm", dpi = 300)
  
  
}
