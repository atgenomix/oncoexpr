#' @title TissueTyper
#' @description predict the original tissue what the sample from 
#' @param group1_ sample 1
#' @param group2_ sample 2
#' @param save_path what the path you want to save the files
#' @param gtexFilePath "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
#' @param groups_list sample types
#' @param sample_id samply id
#' @return heatmap plot for correlation coefficient between multiple tissues of GTEx and your samples
#' @export
TissueTyper <- function(selectMode = FALSE, geneList_ = NULL , group1_ , group2_ , save_path, gtexFilePath, groups_list = c("tissue", "ctc"), sample_id ="BC021-RNA"){
  print("讀取檔案")
  date_ <- Sys.Date() 
  print(date_)
  gtex_tissue_gct_data <- read.delim(gtexFilePath, skip=2)
  head(gtex_tissue_gct_data)
  #gtex_gene_expression <- gtex_tissue_gct_data %>%
  #                          dplyr::group_by(Description) %>%
  #                          dplyr::summarise(across(where(is.numeric), max))%>%
  #                          ungroup()                   
  print("將表現量型別轉成數值並取最大值來合併重複基因")
  
  gtex_gene_expression <- gtex_tissue_gct_data %>%
    group_by(Description) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
    ungroup()
  print("取最大值")                                    
  colnames(gtex_gene_expression)[2] <- "ENSG_id"
  colnames(gtex_gene_expression)[1] <- "Gene"
  print("修改colnames")
  group1_tpm <- group1_[,c("Gene","max_TPM")]
  group2_tpm <- group2_[,c("Gene","max_TPM")]
  print("取得樣本資訊")
  
  colnames(group1_tpm) <- c("Gene", paste0(groups_list[1],"_TPM"))
  colnames(group2_tpm) <- c("Gene", paste0(groups_list[2],"_TPM"))
  print("修改colnames")
  
  merged_df <- gtex_gene_expression %>%
    right_join(group1_tpm, by = "Gene") %>%
    right_join(group2_tpm, by = "Gene")
  print("合併GTEx和兩樣本表現量")
  merged_df <- as.data.frame(na.omit(merged_df))
  print(dim(merged_df))
  #20241114 add new function: select genes
  rownames(merged_df) <- merged_df$"Gene"
  merged_df <-  merged_df %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  print(dim(merged_df))
  cor_matrix<- cor.fk(merged_df[,-c(1:2)]) #"kendall 用cor.fk"
  print("計算kendall's tau")
  cor_data <- melt(cor_matrix)
  colnames(cor_data) <- c("Group1", "Group2", "Correlation")
  
  group1_cor_result <- cor_data[grepl(paste0(groups_list[1],"_TPM"), cor_data$'Group1'),]
  group2_cor_result <- cor_data[grepl(paste0(groups_list[2],"_TPM"), cor_data$'Group1'),]
  
  group1_cor_rank <- group1_cor_result[order(group1_cor_result$'Correlation', decreasing=TRUE),]
  group2_cor_rank <- group2_cor_result[order(group2_cor_result$'Correlation', decreasing=TRUE),]
  print("完成相關係數排序，準備輸出成檔案")
  write.csv(group1_cor_rank, paste0(save_path, paste(c(sample_id, date_ , groups_list[1], "cor_rank.csv"), collapse="_")))
  write.csv(group2_cor_rank, paste0(save_path, paste(c(sample_id, date_ , groups_list[2], "cor_rank.csv"), collapse="_")))
  print("準備繪製成heatmap")
  color_palette <- colorRampPalette( c("white","orange", "red"))(100)
  breaks <- seq( 0, 1, length.out = length(color_palette) + 1)
  heatmap_cor_plot <- pheatmap(
    cor_matrix,
    color = color_palette,
    breaks = breaks,
    na_col = "grey", main = "TissueTyper"
  )
  print("輸出heatmap")
  print(paste0(sample_id, "_", date_, "_GTEx_cor_analysis.png"))
  ggsave(   paste0(sample_id, "_", date_ ,"_GTEx_cor_analysis.png"), 
            plot = heatmap_cor_plot, 
            path = save_path, 
            width = 30, 
            height = 28, 
            units="cm", 
            dpi = 300
  )
  list_ <- list("plot"= heatmap_cor_plot, "rank1"= group1_cor_rank, "rank2"= group2_cor_rank)
  return(list_ )
  
  
  
}



#' @title multiple samples TissueTyper
#' @description predict the original tissue what the sample from 
#' @param selectMode select or excluding
#' @param geneList_ target gene set
#' @param sample_list_ multiple sample id  
#' @param save_path what the path you want to save the files
#' @param gtexFilePath "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
#' @param groups_list sample types
#' @return heatmap plot for correlation coefficient between multiple tissues of GTEx and your samples
#' @export
#' 
#' 

mTissueTyper <- function(selectMode = FALSE, geneList_ = NULL , sample_list_, save_path, gtexFilePath, groups_list = c("tissue", "ctc")){
  print(sample_list_)
  rna_expr <- paste0(unlist(sample_list_), c("_rna_expr"))
  samples_data <- mget(rna_expr, envir = .GlobalEnv)  # 防止變數不存在
  print(samples_data[[1]])
  print("讀取檔案")
  date_ <- Sys.Date() 
  print(date_)
  gtex_tissue_gct_data <- read.delim(gtexFilePath, skip=2)
  head(gtex_tissue_gct_data)
  #gtex_gene_expression <- gtex_tissue_gct_data %>%
  #                          dplyr::group_by(Description) %>%
  #                          dplyr::summarise(across(where(is.numeric), max))%>%
  #                          ungroup()                   
  print("將表現量型別轉成數值並取最大值來合併重複基因")
  
  gtex_gene_expression <- gtex_tissue_gct_data %>%
    group_by(Description) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
    ungroup()
  print("取最大值")                                    
  colnames(gtex_gene_expression)[2] <- "ENSG_id"
  colnames(gtex_gene_expression)[1] <- "Gene"
  print("修改colnames")
  library(purrr)
  processed_list <- lapply(seq_along(samples_data), function(i) {
    samples_data[[i]] %>%
      select(Gene, max_TPM) %>%
      rename(!!unlist(sample_list)[i] := max_TPM)
  })
  
  
  # 使用 reduce 合併處理後的資料框
  merged_df <- purrr::reduce(processed_list, function(df1, df2) {
    right_join(df1, df2, by = "Gene")
  }, .init = gtex_gene_expression)
  
  write.csv(merged_df, paste0(save_path, paste0(paste(patient_code, collapse="_"), "_", date_, "_GTEx_merged_df.csv")))
  
  print("合併GTEx和兩樣本表現量")
  merged_df <- as.data.frame(na.omit(merged_df))
  print(dim(merged_df))
  #20241114 add new function: select genes
  rownames(merged_df) <- merged_df$"Gene"
  merged_df <-  merged_df %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  print(dim(merged_df))
  cor_matrix <- cor.fk(merged_df[,-c(1:2)]) #"kendall 用cor.fk"
  print("計算kendall's tau")
  cor_data <- melt(cor_matrix)
  colnames(cor_data) <- c("Group1", "Group2", "Correlation")
  print("處理數據表格")
  cor_results <- paste0(unlist(sample_list_), "_cor_result")
  print(cor_results)
  cor_results_order <- paste0(unlist(sample_list_), "_cor_order")
  lapply(1:length(unlist(sample_list_)), function(s) assign(cor_results[s], cor_data[grepl(unlist(sample_list_)[s], cor_data$'Group1'),], envir = .GlobalEnv) )
  print(get(cor_results[1]))
  mget_cor_results <- mget(cor_results, envir = .GlobalEnv)
  
  lapply(seq_along(mget_cor_results), function(s) {
    assign(
      cor_results_order[s], 
      mget_cor_results[[s]][order(mget_cor_results[[s]]$Correlation, decreasing = TRUE), ], 
      envir = .GlobalEnv
    )
    write.csv(get(cor_results_order[s]), paste0(save_path, paste(c(unlist(sample_list_)[s], date_ , "cor_rank.csv"), collapse="_")))
  })
  print("完成相關係數排序，準備輸出成檔案")
  
  print("準備繪製成heatmap")
  color_palette <- colorRampPalette( c("white","orange", "red"))(100)
  breaks <- seq( 0, 1, length.out = length(color_palette) + 1)
  heatmap_cor_plot <- pheatmap(
    cor_matrix,
    color = color_palette,
    breaks = breaks,
    na_col = "grey", main = "TissueTyper"
  )
  print("輸出heatmap")
  print(paste0(paste(patient_code, collapse="_"), "_", date_, "_GTEx_cor_analysis.png"))
  ggsave( paste0(paste(patient_code, collapse="_"), "_", date_, "_GTEx_cor_analysis.png"), 
          plot = heatmap_cor_plot, 
          path = save_path, 
          width = 30, 
          height = 28, 
          units="cm", 
          dpi = 300
  )
  list_ <- list("plot"= heatmap_cor_plot, "merged_df" = merged_df, "cor_matrix" = cor_matrix)
  return(list_ )
  
  
  
}