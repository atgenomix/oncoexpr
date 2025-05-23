#' Compute Expression Data in Long Format
#'
#' This function converts a given expression matrix and associated sample metadata 
#' (colData) into a long-format data frame suitable for downstream plotting.
#' The function extracts gene expression values from the matrix and assigns each 
#' sample a group based on the matching of sample names with the metadata.
#'
#' @param exprMatrix A numeric matrix or data frame of expression data, where rows represent genes
#'   and columns represent samples. It is expected that one column is named "GeneSymbol" (or that row names
#'   contain gene symbols) which will be used to annotate the genes.
#' @param colData A data frame containing sample metadata. It must include a column named 
#'   \code{subCode} representing the group labels. This function will also add a \code{mainCode} column 
#'   to \code{colData} based on the column names of \code{exprMatrix} (excluding the "GeneSymbol" column).
#'
#' @return A data frame in long format with the following columns:
#' \describe{
#'   \item{GeneSymbol}{Gene names.}
#'   \item{sample}{Sample identifiers.}
#'   \item{expression}{Gene expression values.}
#'   \item{group}{Group label corresponding to each sample, taken from \code{colData$subCode}.}
#' }
#'
#' @examples
#' \dontrun{
#'   # Assuming normCount is an expression matrix and colData is a data frame with a 'subCode' column:
#'   long_expr <- transfExprFormat(normCount, colData)
#' }
#'
#' @export
transfExprFormat <- function(exprMatrix = normCount, colData = colData) {
  df <- as.data.frame(exprMatrix)
  colData$mainCode <- colnames(df)[-which(colnames(df) == "GeneSymbol")]
  print(colnames(df))
  long_df <- tidyr::pivot_longer(df, cols = -GeneSymbol, names_to = "sample", values_to = "expression")
  long_df$group <- colData$subCode[match(long_df$sample, colData$mainCode)]
  
  return(long_df)
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

# #' @title PCAplot
# #' @description PCA plot for RNAseq expression
# #' @param pcaResult pca table
# #' @param colData gene symbol
# #' @param pcX x axis PC dimension
# #' @param pcY y axis PC dimension
# #' @return PCA plot
# #' @export
# #'
# createPCAPlot <- function(pcaResult, colData, pcX, pcY) {
#   # Convert PCA results to data frame
#   pcaData <- as.data.frame(pcaResult$x)
#   pcaData$Sample <- rownames(pcaData)

#   # Standardize sample names by replacing dots with dashes
#   #pcaData$Sample <- gsub("\\.", "-", pcaData$Sample)

#   # Merge PCA results with colData based on sample IDs
#   mergedData <- merge(pcaData, colData, by.x = "Sample", by.y = "mainCode", all.x = TRUE)
#   missing <- mergedData[is.na(mergedData$subCode), "Sample"]
#   if (length(missing) > 0) {
#     message("Missing group annotation for samples: ", paste(missing, collapse = ", "))
#   }
#   # Create the PCA plot with ggplot2, mapping color to the group (subCode)
#   ggplot(mergedData, aes_string(x = pcX, y = pcY, label = "Sample", color = "subCode")) +
#     geom_point(size = 3) +
#     geom_text(vjust = -0.5) +
#     labs(
#       title = paste("PCA Plot:", pcX, "vs", pcY),
#       x = pcX, y = pcY, color = "Group"
#     ) +
#     theme_minimal()
# }
