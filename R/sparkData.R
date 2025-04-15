
#' @title se to longformat
#' @description se from MAE package transfered to long format dataframe
#' @param se_ se object 
#' @param assay_name selecting assay from se object
#' @param selected_row_ select gene
#' @param selected_col_ select sample
#' @param selected_colData_ filter colData column
#' @return pivot long format from se object
#' @export
#' 
se2longformat <- function(se_ = exp_brca, 
                          assay_name = "fpkm_unstrand", 
                          selected_row_ = c(1:5), 
                          selected_col_ = c("TCGA-5L-AAT0-01A-12R-A41B-07", "TCGA-A2-A04U-01A-11R-A115-07"), 
                          selected_colData_ = "tissue_type") {
  
  assay_mat <- assay(se_, assay_name)
  
  if (is.null(colnames(assay_mat))) {
    stop("Assay matrix does not have column names!")
  }
  
  new_colData <- colData(se_)
  rownames(new_colData) <- colnames(se_)
  
  se_new <- SummarizedExperiment(
    assays = list(data = assay_mat),
    colData = new_colData,
    rowRanges = rowRanges(se_)
  )
  
  mae_new <- MultiAssayExperiment(setNames(list(se_new), assay_name))
  colData(mae_new) <- colData(se_new)
  if (is.null(selected_row_)) {
    selected_row_ <- rownames(se_new)
  } else if (is.numeric(selected_row_)) {
    selected_row_ <- rownames(se_new)[selected_row_]
  }
  if (is.null(selected_col_)) {
    selected_col_ <- seq_len(ncol(se_new))
  }
  
  subset_mae <- mae_new[selected_row_, selected_col_, assay_name]
  
  se2longformat_ <- longFormat(subset_mae, colDataCols = selected_colData_)
  
  return(se2longformat_)
}


#' @title view spark
#' @description show spark datasets and tables
#' @param sc spark connection
#' @param use_dataset selecting dataset to show tables
#' @return tables in the datasets
#' @export
#' 
show_spark_tables <- function(sc, use_dataset){
  show_tables_input <- paste0("SHOW TABLES IN spark_catalog.", use_dataset)
  show_tables <- DBI::dbGetQuery(sc, show_tables_input)
  print(show_tables)
  use_tables_input <- paste0("USE spark_catalog.", use_dataset)
  print(use_tables_input)
  DBI::dbExecute(sc, use_tables_input)
  return(show_tables)
  
}

#' @title expression data
#' @description get the certain samples rnaseq profile
#' @param sc spark connection
#' @param use_table selecting table to show tables
#' @param mainCode_ patient code or treatment code
#' @param subCode_ tissue type or drug dose
#' @param normMethod expression normalization methods (max TPM, FPKM, read count)
#' @return a list of related information and the expression data
#' @export
#' 
query_spark_rnaseqdata <- function(sc, use_table, mainCode_ = "BC026", subCode_ = "PBMC", normMethod = "max_TPM"){  
  df <- tbl(sc, use_table)
  df_filtering <- df %>%
    filter(
      main_code %in% mainCode_,
      sub_code  %in% subCode_,
      normalization %in% normMethod
    )
  
  df_filtering_local <- df_filtering %>% collect() %>% tibble::as_tibble()
  list_return <- list ("mainCode" = mainCode_, "subCode" = subCode_, "normMethod" = normMethod,  "data" = df_filtering_local )
  return(list_return)
}


#' @title generate expression matrix with multiple samples
#' @description convert long to wide for SQL dataframe from spark database
#' @param data fileterd data from spark database 
#' @param mainCode_ mainCode selection
#' @param subCode_ subCode selection
#' @return a expression matrix with gene row and sample column
#' @export
#' 
convert_long_to_wide <- function(data, mainCode_, subCode_) {
  # 添加樣本名稱列
  data <- data %>%
    mutate(sample_name = paste(!!sym(mainCode_), !!sym(subCode_), sep = "_"))
  
  # 轉換為寬表格，基因作為行，樣本作為列
  rna_expr <- data %>%
    select(Gene, sample_name, value) %>%
    pivot_wider(names_from = sample_name, values_from = value, values_fill = NULL)
  
  return(rna_expr)
}

#' @title sample information table
#' @description generate sample information table for MAE colData 
#' @param data fileterd data from spark database 
#' @param mainCode_ mainCode selection
#' @param subCode_ subCode selection
#' @return a table about sample codes
#' @export
#' 
generate_sample_info <- function(data, mainCode_, subCode_) {
  # 提取 mainCode 和 subCode 並生成樣本代號
  sample_info <- data %>%
    distinct(!!sym(mainCode_), !!sym(subCode_)) %>%
    mutate(sample_id = paste0("S", sprintf("%03d", row_number())))
  
  # 整理樣本資訊表格
  sample_info_table <- sample_info %>%
    select(sample_id, !!sym(mainCode_), !!sym(subCode_))
  
  return(sample_info_table)
}




#' Trigger a Spark Query by Table Name Pattern
#'
#' This function connects to a Spark database and retrieves data from a table that matches a specified
#' regular expression pattern. Instead of relying on a fixed table index, it filters table names using
#' the provided pattern, and if multiple tables match, selects the first one. It then executes a query
#' to retrieve all data from the matched table.
#'
#' @param master Character. The Spark master connection string (e.g., "sc://localhost:7077").
#' @param method Character. The method used to connect to Spark (e.g., "spark_connect").
#' @param version Character. The Spark version (e.g., "2.4.3").
#' @param pattern Character. A regular expression used to filter table names, for example "^normcounts",
#'   "^exacttest", or "^coldata".
#' @param output_label Character or NULL. An optional label used to name the returned result. If provided,
#'   the function will return a named list. Default is \code{NULL}.
#' @param verbose Logical. If \code{TRUE}, prints detailed debugging messages. Default is \code{TRUE}.
#'
#' @return A promise that resolves to the query result. If \code{output_label} is provided, the result is
#'   returned as a named list; otherwise, it is returned directly.
#'
#' @importFrom sparklyr spark_connect spark_disconnect
#' @importFrom DBI dbGetQuery dbExecute
#' @importFrom stringr str_equal
#' @importFrom promises future_promise
#'
#' @examples
#' \dontrun{
#'   # Retrieve data from a table whose name starts with "normcounts"
#'   normcount_future <- trigger_cluster_query_by_pattern(
#'     master = "sc://localhost:7077",
#'     method = "spark_connect",
#'     version = "2.4.3",
#'     pattern = "^normcounts",
#'     output_label = "init_tbl_normcount"
#'   )
#'
#'   normcount_future %...>% (function(result) {
#'     print(result)
#'   })
#' }
#'
#' @export
trigger_cluster_query_by_pattern <- function(master, method, version, pattern, output_label = NULL, verbose = TRUE) {
  future_promise({
    # Establish Spark connection
    sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
    on.exit(sparklyr::spark_disconnect(sc_conn))
    
    # Obtain SPARK_USER for connection filtering
    org <- tolower(Sys.getenv("SPARK_USER"))
    if (verbose) print(paste("SPARK_USER:", org))
    
    cond <- ifelse(stringr::str_equal(org, ""), "", sprintf("LIKE '*_%s'", org))
    if (verbose) print(paste("Condition:", cond))
    
    # Retrieve database list and select the first database
    db_list <- DBI::dbGetQuery(sc_conn, sprintf("SHOW DATABASES %s", cond))
    db_name <- db_list[["namespace"]][1]
    if (verbose) print(paste("Using database:", db_name))
    
    # Retrieve all tables in the selected database
    show_init_tbls <- DBI::dbGetQuery(sc_conn, paste0("SHOW TABLES IN ", db_name))
    all_tables <- show_init_tbls$tableName
    
    # Filter tables using the provided pattern
    matched_tables <- all_tables[grepl(pattern, all_tables)]
    if (length(matched_tables) == 0) {
      stop("No table matches the given pattern: ", pattern)
    }
    
    # Select the first matched table
    chosen_table <- matched_tables[1]
    if (verbose) print(paste("Chosen table matching pattern", pattern, ":", chosen_table))
    
    # Switch to the selected database
    DBI::dbExecute(sc_conn, paste0("USE ", db_name))
    
    # Build and execute SQL query
    query <- paste0("SELECT * FROM ", chosen_table)
    result <- DBI::dbGetQuery(sc_conn, query)
    if (verbose) print(paste("Query completed using table", chosen_table))
    
    # Return result with output label if provided
    if (!is.null(output_label)) {
      setNames(list(result), output_label)
    } else {
      result
    }
    
  }, globals = list(master = master, method = method, version = version, verbose = verbose, pattern = pattern, output_label = output_label), seed = TRUE)
}
