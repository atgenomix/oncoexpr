#' @title Create a MultiAssayExperiment object from airway data using DESeq2 analysis.
#' @description This function loads the airway dataset, performs DESeq2 analysis,
#' and creates a MultiAssayExperiment (MAE) object. The DESeq2 results 
#' (p-value and log2FoldChange) are attached to the rowData of the RNAseq assay.
#' @return A MultiAssayExperiment object.
#' @export
create_mae_airway <- function() {
  # 載入必要套件
  library(airway)
  library(DESeq2)
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  
  data(airway)
  se <- airway
  dds <- DESeqDataSet(se, design = ~ dex)
  keep <- rowSums(counts(dds)) >= 10
  dds <- dds[keep, ]
  
  dds$dex <- relevel(dds$dex, ref = "untrt")
  dds <- DESeq(dds)
  res <- results(dds)
  res <- as.data.frame(res)
  
  # 取得 normalized counts 作為表現矩陣
  assay_data <- as.matrix(counts(dds, normalized = TRUE))
  # 取得 sample 資訊，並確保 rownames 與 assay_data 的 colnames 一致
  sample_info <- as.data.frame(colData(dds))
  rownames(sample_info) <- colnames(assay_data)
  # 若原本沒有 group 欄位，這邊以 dex 作為分組資訊
  sample_info$subCode <- sample_info$dex
  
  # 建立 SummarizedExperiment，並將 DESeq2 結果附加到 rowData 中
  se_expression_matrix <- SummarizedExperiment(
    assays = list(normCount = assay_data),
    colData = sample_info,
    rowData = S4Vectors::DataFrame(res[rownames(assay_data), ])
  )
  
  # 建立 MultiAssayExperiment 物件
  mae <- MultiAssayExperiment(
    experiments = list(RNAseq = se_expression_matrix),
    colData = sample_info
  )
  
  mae
}

#' Create a MultiAssayExperiment object from airway data using edgeR analysis.
#'
#' This function loads the airway dataset, performs edgeR differential expression analysis,
#' and creates a MultiAssayExperiment (MAE) object. The edgeR results (including p-value,
#' log2FoldChange, and FDR) are attached to the rowData of the RNAseq assay.
#'
#' @return A MultiAssayExperiment object.
#' @export
create_mae_airway_edgeR <- function() {
  # 載入必要套件
  library(airway)
  library(edgeR)
  library(SummarizedExperiment)
  library(MultiAssayExperiment)
  
  # 載入 airway 資料
  data(airway)
  se <- airway
  
  # 建立 DGEList 物件
  dge <- DGEList(counts = assay(se))
  
  # 過濾低表現基因 (使用 filterByExpr)
  keep <- filterByExpr(dge)
  dge <- dge[keep, , keep.lib.sizes = FALSE]
  
  # 計算 normalization factors
  dge <- calcNormFactors(dge)
  
  # 建立設計矩陣，使用 colData(se)$dex 作為分組依據
  design <- model.matrix(~ dex, data = as.data.frame(colData(se)))
  
  # 估計 dispersion
  dge <- estimateDisp(dge, design)
  
  # 使用 glmQLFit 與 glmQLFTest 進行差異分析
  fit <- glmQLFit(dge, design)
  qlf <- glmQLFTest(fit, coef = 2)
  
  # 擷取所有基因的 edgeR 結果
  edgeR_res <- topTags(qlf, n = Inf)$table
  edgeR_res <- as.data.frame(edgeR_res)
  
  # 取得 normalized counts (以 CPM 表示，不取 log 轉換)
  norm_counts <- as.matrix(cpm(dge, normalized.lib.sizes = TRUE))
  
  # 取得 sample 資訊，並確保 rownames 與 norm_counts 的 colnames 一致
  sample_info <- as.data.frame(colData(se))
  rownames(sample_info) <- colnames(norm_counts)
  sample_info$subCode <- sample_info$dex
  
  # 取出在 norm_counts 與 edgeR_res 中皆存在的基因
  common_genes <- intersect(rownames(norm_counts), rownames(edgeR_res))
  norm_counts <- norm_counts[common_genes, , drop = FALSE]
  edgeR_res_sub <- edgeR_res[common_genes, , drop = FALSE]
  
  # 建立 SummarizedExperiment 物件，並將 edgeR 結果附加到 rowData 中
  se_expression_matrix <- SummarizedExperiment(
    assays = list(normCount = norm_counts),
    colData = sample_info,
    rowData = S4Vectors::DataFrame(edgeR_res_sub)
  )
  
  # 建立 MultiAssayExperiment 物件
  mae <- MultiAssayExperiment(
    experiments = list(RNAseq = se_expression_matrix),
    colData = sample_info
  )
  
  mae
}
#' Generate a ComplexHeatmap heatmap based on a MAE object with separate RNAseq and DEG assays.
#'
#' This function extracts normalized counts from the RNAseq assay and differential expression
#' results (p-value and log2FoldChange) from the DEG assay of a MAE object. It subsets the data
#' based on a provided gene list (if given) and generates a heatmap with three tracks:
#' the z-score of expression, -log10(p-value), and log2FoldChange.
#'
#' It assumes that the RNAseq SummarizedExperiment has an assay named "normCount" and that
#' the DEG SummarizedExperiment contains the differential expression results in its rowData.
#'
#' @param mae A MultiAssayExperiment object containing at least two experiments:
#' one for RNAseq expression (with an assay "normCount") and one for DEG results.
#' @param geneList A character vector of gene identifiers to include. If NULL, all genes are used.
#'
#' @return A \code{ComplexHeatmap} object.
#' @export
#' 
make_heatmap_mae <- function(mae, geneList = NULL) {
  # 載入必要套件
  library(InteractiveComplexHeatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(GetoptLong)
  
  # 從 MAE 中取得 RNAseq 的 SummarizedExperiment，
  # DEG 補充資訊已整合在 rowData 中
  rna_se <- experiments(mae)[["RNAseq"]]
  
  # 取得 normalized counts (assay 名稱為 "normCount")
  m <- assay(rna_se, "normCount")
  print("assay normCount ok")
  # 根據 geneList 進行子集化，若 geneList 為 NULL 則使用所有基因
  if (!is.null(geneList)) {
    valid_genes <- intersect(geneList, rownames(m))
    if (length(valid_genes) == 0) {
      stop("None of the genes in geneList were found in the MAE RNAseq assay.")
    }
    m_sub <- m[valid_genes, , drop = FALSE]
    rd_sub <- rowData(rna_se)[valid_genes, , drop = FALSE]
  } else {
    m_sub <- m
    rd_sub <- rowData(rna_se)
  }

  # 從 rowData 中提取 pvalue 與 log2FoldChange 資料
  if (!("PValue" %in% colnames(rd_sub)) || !("logFC" %in% colnames(rd_sub))) {
    warning("The RNAseq assay rowData does not contain 'PValue' and 'logFC'. Only z-score heatmap is displayed.")
    # 僅顯示 z-score heatmap
    ht <- Heatmap(t(scale(t(m_sub))), name = "z-score",
                  top_annotation = HeatmapAnnotation(
                    group = as.data.frame(colData(rna_se))$group
                  ),
                  show_row_names = TRUE, show_column_names = FALSE,
                  column_title = paste0(nrow(m_sub), " selected genes"),
                  show_row_dend = FALSE)
  } else {
    pval <- as.numeric(unlist(rd_sub$PValue))
    pval[pval == 0] <- .Machine$double.eps
    pval[is.na(pval)] <- 1  # NA 值轉換為 1 (-log10(1)=0)
    log2fc <- as.numeric(unlist(rd_sub$logFC))
    
    # 取得 sample 資訊（假設 RNAseq 與 DEG 使用相同的 colData）
    sample_info <- as.data.frame(colData(rna_se))
    
    # 主要 heatmap：以 z-score 呈現 normalized counts
    ht1 <- Heatmap(t(scale(t(m_sub))), name = "z-score",
                   top_annotation = HeatmapAnnotation(
                     group = sample_info$"subCode"
                   ),
                   show_row_names = FALSE, show_column_names = FALSE,
                   column_title = paste0(nrow(m_sub), " selected genes"),
                   show_row_dend = FALSE)
    
    # 第二 heatmap：呈現 -log10(p-value)
    ht2 <- Heatmap(-log10(pval), name = "-log10(p-value)",
                   show_row_names = FALSE, show_column_names = FALSE,
                   width = unit(5, "mm"))
    
    # 第三 heatmap：呈現 log2FoldChange
    ht3 <- Heatmap(log2fc, name = "log2FoldChange",
                   show_row_names = FALSE, show_column_names = FALSE,
                   width = unit(5, "mm"),
                   col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")))
    
    ht <- ht1 + ht2 + ht3
  }
  
  ht <- draw(ht, merge_legend = TRUE)
  ht
}




#' Run a Shiny App to display the interactive heatmap.
#'
#' This function launches a Shiny app that displays the heatmap generated from a MAE object
#' using a specified gene list. The UI is minimal and only shows the heatmap output.
#'
#' @param mae A MultiAssayExperiment object.
#' @param geneList A character vector of gene identifiers to include in the heatmap. If NULL, all genes will be used.
#'
#' @export
run_heatmap_app <- function(mae, geneList = NULL) {
  # 載入必要套件
  library(shiny)
  library(InteractiveComplexHeatmap)
  
  ui <- fluidPage(
    originalHeatmapOutput("ht", height = 800, containment = TRUE)
  )
  
  server <- function(input, output, session) {
    ht <- make_heatmap_mae(mae, geneList)
    if (!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    } else {
      output$ht_heatmap <- renderPlot({
        grid::grid.newpage()
        grid::grid.text("No data available.")
      })
    }
  }
  
  shinyApp(ui, server)
}

#' Generate pheatmap from a MAE Object
#'
#' This function extracts RNA-seq normalized counts (assay "normCount") from the RNAseq experiment 
#' in a MAE object, optionally subsets the data based on a provided gene list, and retrieves sample 
#' annotation from the colData (defaulting to the "group" column). It then uses the pheatmap package 
#' to generate a heatmap that displays both gene names and sample names.
#'
#' @param mae A MultiAssayExperiment object that must contain an experiment named "RNAseq" with a "normCount" assay.
#' @param geneList Optional character vector specifying the genes to subset. If NULL, all genes are used.
#'
#' @return A pheatmap object.
#' @export

make_heatmap_mae_pheatmap <- function(mae, geneList = NULL) {
  # 若尚未安裝 pheatmap，請先安裝
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    install.packages("pheatmap")
  }
  library(pheatmap)
  
  # 從 MAE 中取得 RNAseq 的 SummarizedExperiment
  if (!"RNAseq" %in% names(experiments(mae))) {
    stop("RNAseq experiment not found in the provided MAE object.")
  }
  rna_se <- experiments(mae)[["RNAseq"]]
  
  # 取得 normalized counts (assay 名稱為 "normCount")
  if (!"normCount" %in% assayNames(rna_se)) {
    stop("normCount assay not found in RNAseq experiment.")
  }
  m <- assay(rna_se, "normCount")
  
  # 若有提供 geneList，則以 geneList 為依據子集化基因
  if (!is.null(geneList)) {
    valid_genes <- intersect(geneList, rownames(m))
    if (length(valid_genes) == 0) {
      stop("None of the genes in geneList were found in the MAE RNAseq assay.")
    }
    m <- m[valid_genes, , drop = FALSE]
  }
  
  # 從 colData 取得樣本的 annotation，預設使用 "group" 欄位，若無則嘗試 "subCode"
  sample_info <- as.data.frame(colData(rna_se))
  if ("group" %in% colnames(sample_info)) {
    annotation_col <- data.frame(Group = sample_info$group)
  } else if ("subCode" %in% colnames(sample_info)) {
    annotation_col <- data.frame(Group = sample_info$subCode)
  } else {
    annotation_col <- data.frame(Group = rep("Unknown", nrow(sample_info)))
  }
  rownames(annotation_col) <- colnames(m)
  
  # 繪製 heatmap，scale = "row" 可將每個基因做 z-score 正規化，並顯示基因與樣本名稱
  p <- pheatmap(m,
                scale = "row",
                annotation_col = annotation_col,
                show_rownames = TRUE,
                show_colnames = TRUE,
                main = paste0(nrow(m), " genes (pheatmap)"))
  
  return(p)
}
