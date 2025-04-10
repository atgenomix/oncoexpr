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
                            coef = 2,
                            lfc_cut = 1,
                            pval_cut = 0.05,
                            useAdjP = FALSE,
                            title = "Volcano Plot",
                            topN = 0,
                            geneCol = NULL,
                            pointSize = 2,
                            ptAlpha = 0.6,
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
  plotData$color[plotData$logFC >= lfc_cut & plotData$pval <= pval_cut] <- "red"
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
      x     = expression(log[2] ~ "Fold Change"),
      y     = bquote(-log[10] ~ .(pvalCol))
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


ggvolcano_custom <- function(
    df, geneName, pValCol = "PValue", logFCCol = "logFC",
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
    geom_point(size = pointSize, alpha = ptAlpha) +
    scale_color_identity() +
    geom_vline(
      xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed",
      color = "black"
    ) +
    geom_hline(
      yintercept = -log10(pval_cut),
      linetype = "dashed", color = "black"
    ) +
    theme_classic() +
    labs(
      title = title, x = expression(log[2] ~ "Fold Change"),
      y = bquote(-log[10] ~ .(pValCol))
    )
  if (topN > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ],
      n = topN
    )
    p <- p + ggrepel::geom_text_repel(
      data = topGenesDF,
      aes(label = gene), size = labelSize, max.overlaps = Inf
    )
  }
  return(p)
}

#' @title ggvolcano_custom_interactive
#' @description ggvolcano_custom with custom parameters and interactive plot for RNAseq analysis
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
#' @param highlight gene with highlight
#' @return ggplot object
#' @export

ggvolcano_custom_interactive <- function(df, geneName, pValCol = "PValue", logFCCol = "logFC",
                                         coef = 2, lfc_cut = 1, pval_cut = 0.05, useAdjP = FALSE,
                                         title = "Volcano Plot", topN = 20, geneCol = NULL,
                                         pointSize = 2, ptAlpha = 0.6, labelSize = 3,
                                         highlight = NULL) {
  plotData <- data.frame(
    gene = geneName,
    logFC = df[[logFCCol]],
    pval = df[[pValCol]]
  )

  plotData$negLogP <- -log10(plotData$pval)

  plotData$color <- "grey"
  plotData$color[plotData$logFC >= lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"

  if (!is.null(highlight)) {
    plotData$color[plotData$gene %in% highlight] <- "orange"
  }

  p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) +
    ggiraph::geom_point_interactive(aes(tooltip = gene, data_id = gene),
      size = pointSize, alpha = ptAlpha
    ) +
    scale_color_identity() +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
    theme_classic() +
    labs(
      title = title,
      x = expression(log[2] ~ "Fold Change"),
      y = bquote(-log[10] ~ .(pValCol))
    )

  if (topN > 0 && nrow(plotData) > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
    p <- p + ggrepel::geom_text_repel(
      data = topGenesDF,
      aes(label = gene),
      size = labelSize,
      max.overlaps = Inf
    )
  }

  return(p)
}

#' @title PCAplot
#' @description PCA plot for RNAseq expression
#' @param pcaResult pca table
#' @param colData gene symbol
#' @param pcX x axis PC dimension
#' @param pcY y axis PC dimension
#' @return PCA plot
#' @export
#'
createPCAPlot <- function(pcaResult, colData, pcX, pcY) {
  # Convert PCA results to data frame
  pcaData <- as.data.frame(pcaResult$x)
  pcaData$Sample <- rownames(pcaData)

  # Standardize sample names by replacing dots with dashes
  pcaData$Sample <- gsub("\\.", "-", pcaData$Sample)

  # Merge PCA results with colData based on sample IDs
  mergedData <- merge(pcaData, colData, by.x = "Sample", by.y = "mainCode", all.x = TRUE)
  missing <- mergedData[is.na(mergedData$subCode), "Sample"]
  if (length(missing) > 0) {
    message("Missing group annotation for samples: ", paste(missing, collapse = ", "))
  }
  # Create the PCA plot with ggplot2, mapping color to the group (subCode)
  ggplot(mergedData, aes_string(x = pcX, y = pcY, label = "Sample", color = "subCode")) +
    geom_point(size = 3) +
    geom_text(vjust = -0.5) +
    labs(
      title = paste("PCA Plot:", pcX, "vs", pcY),
      x = pcX, y = pcY, color = "Group"
    ) +
    theme_minimal()
}

#' @title PCAplot with arrow and cluster
#' @description PCA plot for RNAseq expression
#' @param pcaResult pca table
#' @param colData gene symbol
#' @param pcX x axis PC dimension
#' @param pcY y axis PC dimension
#' @return PCA plot
#' @export

createPCAPlot_withArrows <- function(pcaResult, colData, pcX, pcY) {
  pcaData <- as.data.frame(pcaResult$x)
  pcaData$Sample <- rownames(pcaData)
  mergedData <- merge(pcaData, colData, by.x = "Sample", by.y = "mainCode", all.x = TRUE)

  centroids <- mergedData %>%
    group_by(subCode) %>%
    summarise(
      cx = mean(.data[[pcX]]),
      cy = mean(.data[[pcY]]),
      .groups = "drop"
    )

  mergedData <- mergedData %>%
    left_join(centroids, by = "subCode") %>%
    mutate(
      distToCenter = sqrt((.data[[pcX]] - cx)^2 + (.data[[pcY]] - cy)^2),
      dx = .data[[pcX]] - cx,
      dy = .data[[pcY]] - cy,
      arrowRatio = 0.9,
      xend = cx + dx * arrowRatio,
      yend = cy + dy * arrowRatio
    )

  circleData <- mergedData %>%
    group_by(subCode) %>%
    summarise(
      radius = max(distToCenter),
      cx = first(cx), cy = first(cy),
      .groups = "drop"
    )

  var_percent <- round(100 * pcaResult$sdev^2 / sum(pcaResult$sdev^2), 1)
  names(var_percent) <- colnames(pcaResult$x)
  xlab <- paste0(pcX, " (", var_percent[pcX], "%)")
  ylab <- paste0(pcY, " (", var_percent[pcY], "%)")

  ggplot(mergedData, aes(x = !!sym(pcX), y = !!sym(pcY), color = subCode)) +
    # Group background circles
    geom_circle(data = circleData,
                aes(x0 = cx, y0 = cy, r = radius, fill = subCode),
                inherit.aes = FALSE,
                color = NA, alpha = 0.15) +

    # Sample → center arrows
    geom_segment(
      aes(x = cx, y = cy, xend = xend, yend = yend),
      linetype = "dashed",
      arrow = arrow(length = unit(0.3, "cm"), type = "closed"),
      size = 1.2
    ) +

    # Sample points with border
    geom_point(shape = 21, size = 4, stroke = 1, fill = "white") +

    # Labels (smart position)
    ggrepel::geom_text_repel(aes(label = Sample),
                             size = 3.5,
                             max.overlaps = 50,
                             box.padding = 0.3) +

    # Group center
    geom_point(data = centroids,
               aes(x = cx, y = cy),
               inherit.aes = FALSE,
               shape = 4, size = 6, stroke = 2, color = "black") +

    labs(title = "PCA: Group Structure with Arrows",
         x = xlab, y = ylab, color = "Group", fill = "Group") +
    theme_minimal(base_size = 14)
}
