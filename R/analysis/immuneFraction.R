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
ImmuneFractionPlot <- function(group1_, group2_, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL) {
  tryCatch({
    print("Ensuring the dataframe is correct")
    if (!all(c("Gene", "max_TPM") %in% colnames(group1_)) ||
        !all(c("Gene", "max_TPM") %in% colnames(group2_))) {
      stop("Input data frames must contain 'Gene' and 'max_TPM' columns")
    }
    
    print("Matching genes position")
    group2_matched <- group2_[match(group1_$"Gene", group2_$"Gene"), ]
    
    print("Creating merged dataframe")
    df <- data.frame(
      Gene = group1_$"Gene",
      group1 = group1_$"max_TPM",
      group2 = group2_matched$"max_TPM",
      stringsAsFactors = FALSE
    )
    
    print("Removing NA values")
    df <- df[complete.cases(df), ]
    
    print("Setting column names")
    colnames(df) <- c("Gene", groups_list)
    
    print("Transferring Gene column to rownames")
    expression_matrix <- df %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    
    print(head(expression_matrix))
    
    print("Running deconvolution")
    res_quantiseq <- immunedeconv::deconvolute(
      expression_matrix,
      method = method_,
      tumor = tumorStatus 
    )
    
  }, error = function(e) {
    message("Error in immuneFractionPlot: ", e$message)
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
    "T cell CD4+ (non-regulatory)" = "#D62728",     # Dark red
    "T cell CD8+" = "#FF7F0E",                      # Orange
    "T cell regulatory (Tregs)" = "#FF9896",        # Light red
    "Macrophage M1" = "#2CA02C",                    # Green
    "Macrophage M2" = "#98DF8A",                    # Light green
    "NK cell" = "#1F77B4",                          # Blue
    "B cell" = "#AEC7E8",                           # Light blue
    "Neutrophil" = "#9467BD",                       # Dark purple
    "Monocyte" = "#8C564B",                         # Brown
    "Myeloid dendritic cell" = "#BC80BD",           # Light purple
    "uncharacterized cell" = "#D0E1F9"              # Light pink blue
  )
  
  p <- res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  return(p)
  ggsave(paste0(method_, "_", "immuneFractionPlot.jpg"), path = save_path, plot = p, device = "png", width = 35, height = 25, units = "cm", dpi = 300)
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
ImmuneFractionPlot_mae <- function(mae, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL) {
  rna_counts <- assays(mae[["RNAseq"]])$max_TPM
  
  print("Transferring Gene column to rownames")
  expression_matrix <- rna_counts %>%
    as.matrix()
  
  print(head(expression_matrix))
  
  print("Running deconvolution")
  res_quantiseq <- immunedeconv::deconvolute(
    expression_matrix,
    method = method_,
    tumor = tumorStatus 
  )
  
  my_colors <- c(
    "T cell CD4+ (non-regulatory)" = "#D62728",     # Dark red
    "T cell CD8+" = "#FF7F0E",                      # Orange
    "T cell regulatory (Tregs)" = "#FF9896",        # Light red
    "Macrophage M1" = "#2CA02C",                    # Green
    "Macrophage M2" = "#98DF8A",                    # Light green
    "NK cell" = "#1F77B4",                          # Blue
    "B cell" = "#AEC7E8",                           # Light blue
    "Neutrophil" = "#9467BD",                       # Dark purple
    "Monocyte" = "#8C564B",                         # Brown
    "Myeloid dendritic cell" = "#BC80BD",           # Light purple
    "uncharacterized cell" = "#D0E1F9"              # Light pink blue
  )
  
  p <- res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  
  return(p)
}