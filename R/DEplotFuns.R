
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
                     coef = 2,            # 要比較的係數 (group 的哪一層)
                     pval_cut = 0.05,
                     lfc_cut  = 1,
                     useAdjP  = FALSE) {
  
  # 從 mae 中擷取資料：這裡示範從名為 "RNAseq" 的 experiment 拿資料
  # 且假設其中的 assays 包含一個名為 "max_TPM" 的矩陣
  rna_counts  <- assays(mae[[assayName]])$max_TPM
  sample_info <- colData(mae[[assayName]])
  print(sample_info)
  # 構建分組因子 (這裡以 subCode 作為分組資訊)
  group <- factor(sample_info[[subCodeCol]])
  print(group)
  # 建立 DGEList，並進行標準化
  dge <- DGEList(counts = rna_counts, group = group)
  dge <- calcNormFactors(dge)  # TMM 等方法
  
  # 建置設計矩陣 (含截距項)
  design <- model.matrix(~ group)
  
  # voom 轉換
  v <- voom(dge, design)
  
  # 線性模型擬合 + Bayes校正
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # 使用 decideTests() 判定哪些基因顯著差異
  # 在此設定 p.value 與 lfc 閾值，作為上/下調判定標準
  dt <- decideTests(fit, p.value = pval_cut, lfc = lfc_cut)
  sumDT <- summary(dt)  # 統計上調/下調/不顯著的數目
  
  # 用 topTable() 匯出完整結果 (logFC, P.Value, adj.P.Val 等)
  # coef = 2 代表第二個參數 (group 的對比)，請根據實際情況調整
  topTab <- topTable(fit, coef = coef, number = Inf)
  
  # 最後以 list 形式回傳，以便後續使用
  return(list(
    fit            = fit,      # 可給 ggvolcano_limma() 繪圖
    topTable       = topTab,   # 所有基因的差異分析統計
    decideTests    = dt,       # limma 對基因上調/下調/不顯著的判定
    upDownSummary  = sumDT     # 上下調總結
  ))
}

#' @title limma analysis with MAE input
#' @description differential expression by limma 
#' @param fit limma fit result
#' @param assayName assay 
#' @param subCodeCol which col for comparasion
#' @param coef coefficient for limma 
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @param ptAlpha transparent
#' @param labelSize gene symbol text size
#' @param pointSize point plot size
#' @param topN show a number of genes symbol with high log fold-change
#' @return limma analysis table 
#' @export


ggvolcano_limma <- function(fit,
                            coef      = 2,
                            lfc_cut   = 1,
                            pval_cut  = 0.05,
                            useAdjP   = FALSE,
                            title     = "Volcano Plot",
                            topN      = 0,
                            geneCol   = NULL,
                            pointSize = 2,       # 加入點大小參數
                            ptAlpha   = 0.6,     # 加入點透明度參數
                            labelSize = 3,       # 加入基因標籤字型大小參數
                            ...) {
  
  # 1) 從 limma 的 fit 物件取出指定係數的結果
  tt <- topTable(
    fit      = fit,
    coef     = coef,
    number   = Inf,
    ...
  )
  
  # 如果沒有 geneCol，則使用 rownames 作為 gene 名稱
  if (is.null(geneCol)) {
    tt$gene <- rownames(tt)
  } else {
    if (!geneCol %in% colnames(tt)) {
      stop("您指定的 geneCol 不存在於 topTable() 輸出的欄位中，請確認。")
    }
    tt$gene <- tt[[geneCol]]
  }
  
  # 2) 選擇 p-value 或 adj.P.Val
  pvalCol <- if (useAdjP) "adj.P.Val" else "P.Value"
  if (!pvalCol %in% colnames(tt)) {
    stop(paste0("無法在 topTable 結果中找到欄位: ", pvalCol))
  }
  
  # 3) 建立繪圖用的資料框
  plotData <- data.frame(
    gene  = tt$gene,
    logFC = tt$logFC,
    pval  = tt[[pvalCol]]
  )
  
  # 4) 根據閾值新增顏色欄位（上調：red，下調：blue，其餘 grey）
  plotData$color <- "grey"
  plotData$color[plotData$logFC >=  lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"
  
  # 5) 計算 -log10(p-value)
  plotData$negLogP <- -log10(plotData$pval)
  
  # 6) 繪製基本圖形，將點大小、透明度參數納入
  p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) +
    geom_point(size = pointSize, alpha = ptAlpha) +
    scale_color_identity() +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
    theme_classic() +
    labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = bquote(-log[10]~.(pvalCol))
    )
  
  # 7) 標註前 topN 個最顯著基因（若有指定 topN > 0）
  if (topN > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
    p <- p + 
      ggrepel::geom_text_repel(
        data = topGenesDF,
        aes(label = gene),
        size = labelSize,
        max.overlaps = Inf
      )
  }
  
  return(p)
}



#' @title ggvolcano_custom
#' @description ggvolcano_custom for RNAseq analysis with custom parameters
#' @param df DEG table
#' @param geneName gene symbol
#' @param pValCol which col for p-value
#' @param logFCCol which col for log fold-change
#' @param coef coefficient for limma 
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @param title plot title
#' @param ptAlpha transparent
#' @param labelSize gene symbol text size
#' @param geneCol gene symbol column
#' @param pointSize point plot size
#' @param topN show a number of genes symbol with high log fold-change
#' @return ggplot object
#' @export


ggvolcano_custom <- function (df, geneName, pValCol = "PValue", logFCCol = "logFC", 
    coef = 2, lfc_cut = 1, pval_cut = 0.05, useAdjP = FALSE, 
    title = "Volcano Plot", topN = 20, geneCol = NULL, pointSize = 2, 
    ptAlpha = 0.6, labelSize = 3) {
    plotData <- data.frame(gene = geneName, logFC = df[[logFCCol]], pval = df[[pValCol]])
    plotData$color <- "grey"
    plotData$color[plotData$logFC >= lfc_cut & plotData$pval <= 
        pval_cut] <- "red"
    plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= 
        pval_cut] <- "blue"
    plotData$negLogP <- -log10(plotData$pval)
    p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) + 
        geom_point(size = pointSize, alpha = ptAlpha) + scale_color_identity() + 
        geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", 
            color = "black") + geom_hline(yintercept = -log10(pval_cut), 
        linetype = "dashed", color = "black") + theme_classic() + 
        labs(title = title, x = expression(log[2] ~ "Fold Change"), 
            y = bquote(-log[10] ~ .(pValCol)))
    if (topN > 0) {
        topGenesDF <- head(plotData[order(plotData$pval), ], 
            n = topN)
        p <- p + ggrepel::geom_text_repel(data = topGenesDF, 
            aes(label = gene), size = labelSize, max.overlaps = Inf)
    }
    return(p)
}


View(volcanoData)

ggvolcano_custom_interactive <- function(df, geneName, pValCol = "PValue", logFCCol = "logFC", 
                               coef = 2, lfc_cut = 1, pval_cut = 0.05, useAdjP = FALSE, 
                               title = "Volcano Plot", topN = 20, geneCol = NULL, 
                               pointSize = 2, ptAlpha = 0.6, labelSize = 3) {
  
  # 整理資料
  plotData <- data.frame(gene = geneName, 
                         logFC = df[[logFCCol]], 
                         pval = df[[pValCol]])
  plotData$color <- "grey"
  plotData$color[plotData$logFC >= lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"
  plotData$negLogP <- -log10(plotData$pval)
  
  # 使用互動式的 geom_point_interactive 取代原本的 geom_point
  p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) +
    ggiraph::geom_point_interactive(aes(tooltip = gene, data_id = gene), 
                                     size = pointSize, alpha = ptAlpha) +
    scale_color_identity() +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
    theme_classic() +
    labs(title = title, 
         x = expression(log[2] ~ "Fold Change"), 
         y = bquote(-log[10] ~ .(pValCol)))
  
  # 標註 topN 基因
  if (topN > 0&& nrow(plotData) > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
    p <- p + ggrepel::geom_text_repel(data = topGenesDF, 
                                      aes(label = gene), 
                                      size = labelSize, 
                                      max.overlaps = Inf)
  }
  
  return(p)
}
