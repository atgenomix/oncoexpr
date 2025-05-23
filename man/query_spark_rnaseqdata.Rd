% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spark_query.R
\name{query_spark_rnaseqdata}
\alias{query_spark_rnaseqdata}
\title{expression data}
\usage{
query_spark_rnaseqdata(
  sc,
  use_table,
  mainCode_ = "BC026",
  subCode_ = "PBMC",
  normMethod = "max_TPM"
)
}
\arguments{
\item{sc}{spark connection}

\item{use_table}{selecting table to show tables}

\item{mainCode_}{patient code or treatment code}

\item{subCode_}{tissue type or drug dose}

\item{normMethod}{expression normalization methods (max TPM, FPKM, read count)}
}
\value{
a list of related information and the expression data
}
\description{
get the certain samples rnaseq profile
}
