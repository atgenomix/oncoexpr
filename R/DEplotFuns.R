
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
  
  rna_counts  <- assays(mae[[assayName]])$max_TPM
  sample_info <- colData(mae[[assayName]])
  print(sample_info)
  group <- factor(sample_info[[subCodeCol]])
  print(group)

  dge <- DGEList(counts = rna_counts, group = group)
  dge <- calcNormFactors(dge) 
  
  design <- model.matrix(~ group)
  
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
                            pointSize = 2,      
                            ptAlpha   = 0.6,     
                            labelSize = 3,       
                            ...) {
  
  tt <- topTable(
    fit      = fit,
    coef     = coef,
    number   = Inf,
    ...
  )
  

  if (is.null(geneCol)) {
    tt$gene <- rownames(tt)
  } else {
    if (!geneCol %in% colnames(tt)) {
      stop("您指定的 geneCol 不存在於 topTable() 輸出的欄位中，請確認。")
    }
    tt$gene <- tt[[geneCol]]
  }
  
  pvalCol <- if (useAdjP) "adj.P.Val" else "P.Value"
  if (!pvalCol %in% colnames(tt)) {
    stop(paste0("無法在 topTable 結果中找到欄位: ", pvalCol))
  }

  plotData <- data.frame(
    gene  = tt$gene,
    logFC = tt$logFC,
    pval  = tt[[pvalCol]]
  )
  
  plotData$color <- "grey"
  plotData$color[plotData$logFC >=  lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"

  plotData$negLogP <- -log10(plotData$pval)

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


ggvolcano_custom_interactive <- function(df, geneName, pValCol = "PValue", logFCCol = "logFC", 
                               coef = 2, lfc_cut = 1, pval_cut = 0.05, useAdjP = FALSE, 
                               title = "Volcano Plot", topN = 20, geneCol = NULL, 
                               pointSize = 2, ptAlpha = 0.6, labelSize = 3) {
  
  plotData <- data.frame(gene = geneName, 
                         logFC = df[[logFCCol]], 
                         pval = df[[pValCol]])
  plotData$color <- "grey"
  plotData$color[plotData$logFC >= lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"
  plotData$negLogP <- -log10(plotData$pval)
  
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
  
  if (topN > 0&& nrow(plotData) > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
    p <- p + ggrepel::geom_text_repel(data = topGenesDF, 
                                      aes(label = gene), 
                                      size = labelSize, 
                                      max.overlaps = Inf)
  }
  
  return(p)
}
