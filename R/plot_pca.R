
#' @title PCAplot with distance and cluster labeling
#' @description PCA plot for RNAseq expression
#' @param pcaResult pca table
#' @param colData gene symbol
#' @param enableClustering set clustering or not
#' @param centers centers number (default = 2)
#' @return PCA plot
#' @export
createPCAPlot <- function(pcaResult, colData, enableClustering = FALSE, centers = 2) {
  pcaData <- as.data.frame(pcaResult$x)
  pcaData$Sample <- rownames(pcaData)
  
  explained <- summary(pcaResult)$importance["Proportion of Variance", ]
  xlab <- paste0("PC1 (", round(explained["PC1"] * 100, 1), "%)")
  ylab <- paste0("PC2 (", round(explained["PC2"] * 100, 1), "%)")

  x_range <- range(pcaData$PC1)
  y_range <- range(pcaData$PC2)

  x_padding <- 0.15 * diff(x_range)
  y_padding <- 0.15 * diff(y_range)
  x_lim <- c(x_range[1] - x_padding, x_range[2] + x_padding)
  y_lim <- c(y_range[1] - y_padding, y_range[2] + y_padding)

  if (enableClustering) {

    clusterData <- pcaData[, c("PC1", "PC2")]
    set.seed(123)
    km.res <- kmeans(clusterData, centers = centers, nstart = 25)
    
    plot <- fviz_cluster(km.res,
                        data = clusterData,
                        stand = FALSE,
                        ellipse.type = "norm", 
                        star.plot = TRUE,
                        repel = TRUE,
                        geom = "point",
                        show.clust.cent = TRUE,
                        labelsize = 4,
                        ggtheme = theme_minimal()
    )
    
    plot <- plot + coord_cartesian(xlim = x_lim, ylim = y_lim) + labs(title = "K‑Means Clustering",
          x = xlab,
          y = ylab,
          color = "Cluster") + 
          guides(shape = guide_none())
    
  } else {
    mergedData <- merge(pcaData, colData, by.x = "Sample", by.y = "mainCode", all.x = TRUE)
    
    plot <- ggplot(mergedData, aes(x = PC1, y = PC2, label = Sample, color = subCode)) +

      coord_cartesian(xlim = x_lim, ylim = y_lim) +
      geom_point(size = 3) +

      ggrepel::geom_text_repel(size = 4, max.overlaps = 200) +
      labs(title = "Cases vs. Controls",
          x = xlab,
          y = ylab,
          color = "Group") +
      theme_minimal()
  }
  
  return(plot)
}