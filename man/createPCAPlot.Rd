% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_pca.R
\name{createPCAPlot}
\alias{createPCAPlot}
\title{PCAplot with distance and cluster labeling}
\usage{
createPCAPlot(pcaResult, colData, enableClustering = FALSE, centers = 2)
}
\arguments{
\item{pcaResult}{pca table}

\item{colData}{gene symbol}

\item{enableClustering}{set clustering or not}

\item{centers}{centers number (default = 2)}
}
\value{
PCA plot
}
\description{
PCA plot for RNAseq expression
}
