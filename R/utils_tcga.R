#' Extract the three key tables from a SummarizedExperiment
#'
#' Fast helper that returns the expression matrix, gene-level metadata and
#' sample-level metadata from a single `SummarizedExperiment` object.
#'
#' @param se A `SummarizedExperiment` (or `RangedSummarizedExperiment`).
#' @param assay_name `character(1)` Name of the assay to extract.  
#'   If `NULL` (default) the first assay is used.
#'
#' @return A named `list` with three elements  
#' * **`exp_data`** – expression matrix (`matrix` or `DelayedMatrix`)  
#' * **`gene_info`** – gene-level metadata (`data.frame`)  
#' * **`sample_info`** – sample-level metadata (`data.frame`)
#'
#' @examples
#' \dontrun{
#' se <- readRDS("query_TCGA_BRCA.rds")
#' tbls <- query_se_tables(se, assay_name = "tpm_unstrand")
#' str(tbls)
#' }
#' @export
query_se_tables <- function(se, assay_name = NULL) {
  stopifnot(inherits(se, "SummarizedExperiment"))

  if (is.null(assay_name)) {
    assay_name <- names(assays(se))[1]
  } else if (!assay_name %in% names(assays(se))) {
    stop("`assay_name` not found in this SummarizedExperiment: ",
         paste(names(assays(se)), collapse = ", "))
  }

  expr_mat    <- assay(se, assay_name)
  gene_info   <- as.data.frame(rowData(se))
  sample_info <- as.data.frame(colData(se))

  list(
    exp_data    = expr_mat,
    gene_info   = gene_info,
    sample_info = sample_info
  )
}


#' Filter gene metadata by biotype (or any column)
#'
#' Convenience wrapper around **dplyr** for sub-setting a gene annotation
#' table and returning only the desired columns.
#'
#' @param gene_info A `data.frame` – typically `rowData()` converted to a
#'   data frame.
#' @param gene_biotype Character vector of biotypes to keep (default:
#'   `"protein_coding"`).  Use `NULL` to skip biotype filtering.
#' @param gene_id_col Column name that contains the stable gene ID.
#' @param gene_name_col Column name that contains the gene symbol.
#' @param gene_type_col Column name that stores the biotype information.
#'
#' @return A `data.frame` with three columns: `gene_id`, `gene_name`,
#'   `gene_type`.
#'
#' @examples
#' \dontrun{
#' se  <- readRDS("query_TCGA_BRCA.rds")
#' gi  <- query_se_tables(se)$gene_info
#' pcs <- geneinfo_filter(gi, gene_biotype = "protein_coding")
#' head(pcs)
#' }
#' @export
geneinfo_filter <- function(gene_info,
                            gene_biotype  = "protein_coding",
                            gene_id_col   = "gene_id",
                            gene_name_col = "gene_name",
                            gene_type_col = "gene_type") {

  requireNamespace("dplyr")
  requireNamespace("rlang")

  out <- gene_info |>
    dplyr::filter(
      if (!is.null(gene_biotype))
        .data[[gene_type_col]] %in% gene_biotype
      else
        TRUE
    ) |>
    dplyr::transmute(
      gene_id   = .data[[gene_id_col]],
      gene_name = .data[[gene_name_col]],
      gene_type = .data[[gene_type_col]]
    )

  out
}

#' Sub-set an expression matrix to a selected gene list
#'
#' Aligns a gene list produced by \code{\link{geneinfo_filter}} (or any
#' similar data frame) with an expression matrix returned by
#' \code{\link{query_se_tables}}, optionally replacing the rownames by gene
#' symbols and guaranteeing matrix dimensions.
#'
#' @param tables A list created by \code{\link{query_se_tables}}.
#' @param genes_df A `data.frame` containing at least the `gene_id` column.
#' @param id_col Name of the gene-ID column in `genes_df`.
#' @param symbol_col Name of the symbol column in `genes_df`.
#' @param use_symbol Logical; if `TRUE` (default) replace rownames by gene
#'   symbols (duplicates are made unique via `make.unique()`).
#' @param drop Passed to the matrix sub-setting (`drop = FALSE` keeps
#'   2-dimensional structure even if only one row is selected).
#'
#' @return A list with  
#' * **`expr_mat`** – sub-set expression matrix  
#' * **`gene_info`** – the supplied `genes_df` (for convenience)  
#' * **`sample_info`** – carried over from `tables`
#'
#' @examples
#' \dontrun{
#' se  <- readRDS("query_TCGA_BRCA.rds")
#' tbl <- query_se_tables(se, "tpm_unstrand")
#'
#' pcs <- geneinfo_filter(tbl$gene_info)      # protein-coding by default
#' res <- subset_expr(tbl, pcs)
#'
#' dim(res$expr_mat)          # genes × samples
#' head(rownames(res$expr_mat))
#' }
#' @export
subset_expr <- function(tables,
                        genes_df,
                        id_col      = "gene_id",
                        symbol_col  = "gene_name",
                        use_symbol  = TRUE,
                        drop        = FALSE) {

  idx <- match(genes_df[[id_col]], rownames(tables$exp_data))
  idx <- idx[!is.na(idx)]

  if (!length(idx))
    stop("None of the selected genes are present in `tables$exp_data`.")

  expr_mat <- tables$exp_data[idx, , drop = drop]

  if (use_symbol) {
    sym <- genes_df[[symbol_col]][match(rownames(expr_mat), genes_df[[id_col]])]
    rownames(expr_mat) <- make.unique(sym)
  }

  list(
    expr_mat    = expr_mat,
    gene_info   = genes_df,
    sample_info = tables$sample_info
  )
}
