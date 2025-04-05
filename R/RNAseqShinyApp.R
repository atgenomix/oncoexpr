#' @import edgeR
#' @import limma
#' @import MultiAssayExperiment
#' @import SummarizedExperiment
#' @import pcaPP
#' @import reshape2
#' @import stringr
#' @import readxl
#' @import ggplot2
#' @import tidyr
#' @import tibble
#' @import viridis
#' @import RColorBrewer
#' @import pheatmap
#' @import ggrepel
#' @import readr
#' @import dplyr
#' @import sparklyr
#' @import org.Hs.eg.db
#' @import enrichplot
#' @import clusterProfiler
#' @import DBI
#' @import shinybusy
#' @import shiny
#' @import bslib
#' @import DT
#' @import InteractiveComplexHeatmap
#' @import ggiraph
#' @import shinydashboard
#' @import promises
#' @import future
#' @import purrr
#' @import shinycssloaders
#' @importFrom ggpubr color_palette
#' @importFrom enrichplot color_palette
#' @importFrom DT dataTableOutput renderDataTable

#'
NULL
#' @keywords internal
#' RNAseqShinyAppSpark
#'
#' This function creates a Shiny application for RNA sequencing data analysis using Spark.
#' It provides an interactive interface for users to explore and analyze RNAseq data.
#'

#' @title shiny app for RNAseq for public use
#' @description start the RNAseq shiny app with spark connection for shiny proxy at the seqslab console
#' @return A Shiny application object.
#'
#' @examples
#' \dontrun{
#' library(sparklyr)
#' sc <- spark_connect(master = "local")
#' input_data <- read.csv("path/to/rnaseq_data.csv")
#' RNAseqShinyAppSpark(input_data, sc, "RNAseq Analysis", 8080)
#' }
#'
#' @export


RNAseqShinyAppSpark <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {
  plan(multisession, workers = 4)
  # plan(sequential)
  print(future::plan())
  ui <- fluidPage(
    navbarPage(
      title = "RNAseq App",
      header = tagList(
        tags$style(".shinybusy-overlay {opacity: 0.7; background-color: #7c7c7c;}"),
        add_busy_spinner(
          spin = "fading-circle",
          position = "full-page",
          timeout = 1000
        )
      ),
      tabPanel(
        title = "Gene Expression Profile",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            h4("Analysis Runs"),
            dbBrowserUI("dbBrowser1")
          ),
          mainPanel(
            fluidRow(
              column(
                width = 12,
                tabsetPanel(
                  tabPanel(
                    "normCount Table",
                    withSpinner(DT::dataTableOutput("wide_table_dt", width = "100%"))
                  ),
                  tabPanel(
                    "DEG Table",
                    downloadButton("download_DEG", "Download DEG CSV"),
                    withSpinner((DT::dataTableOutput("DEG_table", width = "100%")))
                  )
                )
              )
            )
          )
        )
      ),
      tabPanel(
        title = "PCA",
        mainPanel(
          pcaModuleUI("pca1")
        )
      ),
      tabPanel(
        title = "Differential Expression Analysis",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            sliderInput("lfc_cut", "Fold Change Threshold (log2):",
              min = 0, max = 10, value = 1, step = 0.1
            ),
            sliderInput("pval_cut", "p-value Threshold:",
              min = 0.001, max = 1, value = 0.05, step = 0.001
            ),
            sliderInput("pointSize", "Point Size:",
              min = 1, max = 5, value = 2, step = 0.5
            ),
            sliderInput("ptAlpha", "Transparent:",
              min = 0.1, max = 1, value = 0.6, step = 0.1
            ),
            sliderInput("labelSize", "Gene Label Size:",
              min = 1, max = 6, value = 3, step = 0.5
            ),
            numericInput("topN", "Label N number of Genes (0 is no labeling):",
              value = 100, min = 0, max = 1000
            ),
            actionButton("run_DEG", "Update DEG")
          ),
          mainPanel(
            width = 12,
            tabsetPanel(
              tabPanel(
                "Volcano Plot",
                fluidRow(
                  column(
                    width = 12,
                    interactivePlotsUI("volcano_plots")
                  ),
                  column(
                    width = 12,
                    mod_geneSelector_ui("gene_selector")
                  )
                )
              ),
              tabPanel(
                "Heatmap",
                fluidRow(
                  column(
                    width = 6,
                    box(
                      title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                      originalHeatmapOutput("ht", height = 1000, containment = TRUE)
                    )
                  ),
                  column(
                    width = 6,
                    box(
                      title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                      subHeatmapOutput("ht", title = NULL, containment = TRUE)
                    )
                  )
                )
              ),
              tabPanel(
                "Gene Set Enrichment",
                fluidRow(
                  column(
                    12,
                    h4("UPregulated DEGs"),
                    tabsetPanel(
                      tabPanel("KEGG", withSpinner(plotOutput("G1_KEGG"))),
                      tabPanel("MF", withSpinner(plotOutput("G1_MF"))),
                      tabPanel("BP", withSpinner(plotOutput("G1_BP"))),
                      tabPanel("CC", withSpinner(plotOutput("G1_CC"))),
                      tabPanel("GSEA(KEGG)", withSpinner(gseaFCModuleUI("gsea_up")))
                    )
                  )
                ),
                fluidRow(
                  column(
                    12,
                    h4("DOWNregulated DEGs"),
                    tabsetPanel(
                      tabPanel("KEGG", withSpinner(plotOutput("G2_KEGG"))),
                      tabPanel("MF", withSpinner(plotOutput("G2_MF"))),
                      tabPanel("BP", withSpinner(plotOutput("G2_BP"))),
                      tabPanel("CC", withSpinner(plotOutput("G2_CC"))),
                      tabPanel("GSEA(KEGG)", withSpinner(gseaFCModuleUI("gsea_down")))
                    )
                  )
                )
              )
            )
          )
        )
      )
      
    )
  )

  server <- function(input, output, session) {

    results <- reactiveValues(
      db_info = NULL,
      table_list = NULL,
      grouplist = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      coldata = NULL
    )
    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    DEG_table <- reactiveVal(NULL)
    DEG_summary <- reactiveVal(NULL)
    wide_data <- reactiveVal(NULL)
    maeColData <- reactiveVal(NULL)

    sc <- sparklyr::spark_connect(master = master, method = method, version = version)

    session$onSessionEnded(function() {
      if (!is.null(sc)) {
        sparklyr::spark_disconnect(sc)
        message("Spark connection disconnected.")
      }
    })

    results$db_info <- reactive({
      req(sc)
      start_db <- Sys.time()
      message(sprintf("[db_info] Start: %s", start_db))
      info <- dbBrowserServer("dbBrowser1", sc)
      end_db <- Sys.time()
      message(sprintf("[db_info] End: %s (Duration: %.2f seconds)", end_db, as.numeric(difftime(end_db, start_db, units = "secs"))))
      info
    })

    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      selected_db_name <- results$db_info$selected_db()
      message(sprintf("[DB Selected] %s at %s", selected_db_name, Sys.time()))
      
      start_tbl <- Sys.time()
      DBI::dbExecute(sc, paste0("USE ", selected_db_name))
      tbl_list_query <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", selected_db_name))
      tbls <- tbl_list_query$tableName
      message(sprintf("[Tables] Retrieved at %s: %s", Sys.time(), paste(tbls, collapse = ", ")))
      
      prefix <- c("^normcounts|^exacttest|^coldata")
      tbl_list_query_prefix <- tbl_list_query[grepl(prefix, tbls), ]
      message(sprintf("[Filter] Tables after prefix filter at %s: %s", Sys.time(), paste(tbl_list_query_prefix$tableName, collapse = ", ")))
      
      tbls_with_prefix <- tbl_list_query_prefix$tableName
      tbls_with_time_filter <- get_latest_file_group_df(tbls_with_prefix)
      message("[Time Filter] Result:")
      print(tbls_with_time_filter)
      
      if (sum(tbls_with_time_filter$is_latest) == 0) {
        message(sprintf("[Latest] No latest table found at %s", Sys.time()))
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest == FALSE, ]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest == FALSE, ]
      } else {
        message(sprintf("[Latest] Latest table found at %s", Sys.time()))
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest == TRUE, ]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest == TRUE, ]
      }
      
      tbls_with_prefix_time <- summary_table[["file"]]
      message("[Summary Table] ")
      print(summary_table)
      results$table_list <- tbl_list_query_prefix_time
      
      normcount_tbls <- tbl_list_query_prefix_time[grepl("^normcounts", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]
      exacttest_tbls <- tbl_list_query_prefix_time[grepl("^exacttest", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]
      coldata_tbls <- tbl_list_query_prefix_time[grepl("^coldata", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]
      
      message(sprintf("[Tables] normcount: %s; exacttest: %s; coldata: %s", 
                      paste(normcount_tbls, collapse = ", "),
                      paste(exacttest_tbls, collapse = ", "),
                      paste(coldata_tbls, collapse = ", ")))
      
      req(normcount_tbls, exacttest_tbls, coldata_tbls)
      
      normcount_promise <- future_promise(
        {
          t0 <- Sys.time()
          message(sprintf("[normcount] Start at %s", t0))
          
          sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          normcount_tbl <- tbl(sc_conn, normcount_tbls[1]) %>%
            dplyr::rename(GeneSymbol = genes) %>%
            dplyr::select(-`_c0`) %>%
            dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
          
          normcount <- collect(normcount_tbl)
          t1 <- Sys.time()
          message(sprintf("[normcount] End at %s (Duration: %.2f seconds)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
          normcount
        },
        globals = list(master = master, method = method, version = version,
                      normcount_tbls = normcount_tbls, selected_db_name = selected_db_name),
        seed = TRUE
      )
      
      exacttest_promise <- future_promise(
        {
          t0 <- Sys.time()
          message(sprintf("[exacttest] Start at %s", t0))
          
          sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          exacttest_tbl <- tbl(sc_conn, exacttest_tbls[1]) %>%
            dplyr::rename(GeneSymbol = genes) %>%
            dplyr::select(-`_c0`) %>%
            dplyr::mutate(
              logFC = round(logFC, 4),
              logCPM = round(logCPM, 4)
            ) %>%
            sparklyr::spark_apply(function(df) {
              if ("pvalue" %in% colnames(df)) {
                df$pvalue <- sapply(df$pvalue, function(x) {
                  if (is.na(x)) return(NA_character_)
                  s <- sprintf("%.4e", x)
                  parts <- strsplit(s, "e", fixed = TRUE)[[1]]
                  mantissa <- as.numeric(parts[1])
                  exponent <- as.numeric(parts[2])
                  new_mantissa <- mantissa * 0.1
                  new_exponent <- exponent + 1
                  paste0(sprintf("%.4f", new_mantissa), "e", new_exponent)
                })
              }
              df
            })
          
          exacttest <- collect(exacttest_tbl)
          t1 <- Sys.time()
          message(sprintf("[exacttest] End at %s (Duration: %.2f seconds)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
          exacttest
        },
        globals = list(master = master, method = method, version = version,
                      exacttest_tbls = exacttest_tbls, selected_db_name = selected_db_name),
        seed = TRUE
      )
      
      coldata_promise <- if (length(coldata_tbls) > 0) {
        future_promise(
          {
            t0 <- Sys.time()
            message(sprintf("[coldata] Start at %s", t0))
            
            sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
            on.exit(sparklyr::spark_disconnect(sc_conn))
            
            DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
            query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
            coldata <- DBI::dbGetQuery(sc_conn, query_coldata)
            
            t1 <- Sys.time()
            message(sprintf("[coldata] End at %s (Duration: %.2f seconds)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
            coldata
          },
          globals = list(master = master, method = method, version = version,
                        coldata_tbls = coldata_tbls, selected_db_name = selected_db_name),
          seed = TRUE
        )
      } else {
        normcount_promise %...>% (function(normcount) {
          future_promise(
            {
              message(sprintf("[coldata] Generating random coldata at %s", Sys.time()))
              generate_colData_random(normcount, genecol = "GeneSymbol")
            },
            seed = TRUE
          )
        })
      }
      
      withProgress(message = "Processing data...", value = 0, {
        overall_start <- Sys.time()
        incProgress(0.1, detail = sprintf("Starting normcount query at %s", overall_start))
        normcount_data <- value(normcount_promise)
        t_norm <- Sys.time()
        incProgress(0.3, detail = sprintf("Normcount done (Duration: %.2f sec)", as.numeric(difftime(t_norm, overall_start, units = "secs"))))
        
        incProgress(0.1, detail = sprintf("Starting exacttest query at %s", Sys.time()))
        exacttest_data <- value(exacttest_promise)
        t_exact <- Sys.time()
        incProgress(0.3, detail = sprintf("Exacttest done (Duration: %.2f sec)", as.numeric(difftime(t_exact, t_norm, units = "secs"))))
        
        incProgress(0.1, detail = sprintf("Starting coldata query at %s", Sys.time()))
        coldata <- value(coldata_promise)
        t_col <- Sys.time()
        incProgress(0.1, detail = sprintf("Coldata done (Duration: %.2f sec)", as.numeric(difftime(t_col, t_exact, units = "secs"))))
        
        overall_end <- Sys.time()
        message(sprintf("[Overall] All queries completed (Total Duration: %.2f seconds)",
                        as.numeric(difftime(overall_end, overall_start, units = "secs"))))
        
        results$normcount_data <- normcount_data
        results$exacttest_data <- exacttest_data
        results$coldata <- coldata
        
        message("=== normcount_data ===")
        print(head(results$normcount_data))
        message("=== exacttest_data ===")
        print(head(results$exacttest_data))
        message("=== coldata ===")
        print(head(results$coldata))
      })
    })
      
    observe({
      withProgress(message = "Processing Experiment Data...", value = 0, {
        t0 <- Sys.time()
        incProgress(0.1, detail = "Loading reactive values")
        
        req(DEG_table(), wide_data(), maeColData())
        wide_df <- wide_data()
        deg_df  <- DEG_table()
        sample_info <- maeColData()
        
        incProgress(0.2, detail = "Processing assay_data")
        t1 <- Sys.time()
        message(sprintf("[Assay Data] wide_data loaded at %s (Duration: %.2f sec)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
        
        # 將 wide_data 除去 "GeneSymbol" 欄位後轉成矩陣，並以 GeneSymbol 設定 rownames
        assay_data <- as.matrix(wide_df[, setdiff(colnames(wide_df), "GeneSymbol"), drop = FALSE])
        if ("GeneSymbol" %in% colnames(wide_df)) {
          rownames(assay_data) <- wide_df$GeneSymbol
        }
        
        incProgress(0.3, detail = "Processing DEG data")
        if ("GeneSymbol" %in% colnames(deg_df)) {
          rownames(deg_df) <- deg_df$GeneSymbol
        }
        common_genes <- intersect(rownames(assay_data), rownames(deg_df))
        assay_data <- assay_data[common_genes, , drop = FALSE]
        deg_data_sub <- deg_df[common_genes, , drop = FALSE]
        
        incProgress(0.5, detail = "Processing sample info")
        rownames(sample_info) <- colnames(assay_data)
        
        # 建立 SummarizedExperiment 與 MultiAssayExperiment 物件
        se_expression_matrix <- SummarizedExperiment(
          assays = list(normCount = assay_data),
          colData = sample_info,
          rowData = S4Vectors::DataFrame(deg_data_sub)
        )
        mae <- MultiAssayExperiment(
          experiments = list(RNAseq = se_expression_matrix),
          colData = sample_info
        )
        settingMAE(mae)
        
        t_end <- Sys.time()
        incProgress(1, detail = "Experiment Data Processed")
        message(sprintf("[Experiment Data] Completed at %s (Total Duration: %.2f sec)", t_end, as.numeric(difftime(t_end, t0, units = "secs"))))
      })
    })

    output$DEG_table <- DT::renderDataTable({
      t0 <- Sys.time()
      req(results$coldata, results$normcount_data, results$exacttest_data)
      message(sprintf("[DEG_table] Start at %s", t0))
      dt <- DT::datatable(
        DEG_table(),
        options = list(pageLength = 20, autoWidth = TRUE)
      )
      t1 <- Sys.time()
      message(sprintf("[DEG_table] End at %s (Duration: %.2f sec)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
      dt
    })

    observe({
      withProgress(message = "Updating Results...", value = 0, {
        t0 <- Sys.time()
        incProgress(0.2, detail = "Assigning reactive values")
        
        req(results$exacttest_data, results$normcount_data, results$coldata)
        DEG_table(results$exacttest_data)
        wide_data(results$normcount_data)
        maeColData(results$coldata)
        message(sprintf("[Assign Reactive] Updated at %s", Sys.time()))
        
        incProgress(0.4, detail = "Processing additional parameters")
        # 使用本地變數以避免重複計算
        normCount <- wide_data()
        volcanoData <- DEG_table()
        colData <- maeColData()
        exprData <- transfExprFormat(normCount, colData)
        params <- list(
          lfc_cut    = input$lfc_cut,
          pval_cut   = input$pval_cut,
          pointSize  = input$pointSize,
          ptAlpha    = input$ptAlpha,
          labelSize  = input$labelSize,
          topN       = input$topN,
          use_adjP   = input$use_adjP
        )
        
        # 過濾 DEG_table 取得上調與下調的基因
        DEG_data <- volcanoData
        topGenes <- DEG_data[DEG_data$PValue < input$pval_cut & DEG_data$logFC > input$lfc_cut, "GeneSymbol"]
        downGenes <- DEG_data[DEG_data$PValue < input$pval_cut & DEG_data$logFC < -input$lfc_cut, "GeneSymbol"]
        geneListVec <- c(topGenes, downGenes)
        
        incProgress(0.7, detail = "Updating interactive plots")
        selected_gene <- if (!is.null(geneListReactive)) {
          mod_geneSelector_server("gene_selector", volcanoData, geneListVec)
        } else {
          NULL
        }
        interactivePlotsServer("volcano_plots", volcanoData = volcanoData, exprData = exprData, params = params, selectedGene = selected_gene)
        
        t_end <- Sys.time()
        incProgress(1, detail = "Results updated")
        message(sprintf("[Update Reactive] Completed at %s (Duration: %.2f sec)", t_end, as.numeric(difftime(t_end, t0, units = "secs"))))
      })
    })

    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL)

    geneListReactive <- eventReactive(input$run_DEG, {
      t0 <- Sys.time()
      req(DEG_table(), maeColData(), wide_data())
      DEG_data <- DEG_table()
      topGenes <- DEG_data[DEG_data$PValue < input$pval_cut & DEG_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_data[DEG_data$PValue < input$pval_cut & DEG_data$logFC < -input$lfc_cut, "GeneSymbol"]
      
      topGeneList(topGenes)
      downGeneList(downGenes)
      message(sprintf("[geneListReactive] Top gene count: %s; Down gene count: %s at %s", 
                      length(topGenes), length(downGenes), Sys.time()))
      gene_list <- paste(c(topGenes, downGenes), collapse = ",")
      t1 <- Sys.time()
      message(sprintf("[geneListReactive] Completed at %s (Duration: %.2f sec)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
      gene_list
    })

    observeEvent(geneListReactive(), {
      withProgress(message = "Generating Heatmap...", value = 0, {
        t0 <- Sys.time()
        req(geneListReactive(), settingMAE())
        mae_obj <- settingMAE()
        
        geneListVec <- trimws(unlist(strsplit(geneListReactive(), ",")))
        incProgress(0.5, detail = "Creating heatmap")
        ht <- make_heatmap_mae(mae_obj, geneListVec)
        
        if (!is.null(ht)) {
          makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
        } else {
          output$ht_heatmap <- renderPlot({
            grid::grid.newpage()
            grid::grid.text("No data available.")
          })
        }
        
        t1 <- Sys.time()
        incProgress(1, detail = "Heatmap generated")
        message(sprintf("[Heatmap] Completed at %s (Duration: %.2f sec)", t1, as.numeric(difftime(t1, t0, units = "secs"))))
      })
    })

    # observeEvent(geneListReactive(), {
    #   req(geneListReactive(), settingMAE())
    #   mae <- settingMAE()
      
    #   geneListVec <- unlist(strsplit(geneListReactive(), ","))
    #   geneListVec <- trimws(geneListVec)
      
    #   future_promise({
    #     
    #     make_heatmap_mae(mae, geneListVec)
    #   }) %...>% (function(ht) {
    #     
    #     if (!is.null(ht)) {
    #       makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
    #     } else {
    #       output$ht_heatmap <- renderPlot({
    #         grid::grid.newpage()
    #         grid::grid.text("No data available.")
    #       })
    #     }
    #   }) %...!% (function(e) {
    #   
    #     output$ht_heatmap <- renderPlot({
    #       grid::grid.newpage()
    #       grid::grid.text(paste("An error occurred:", e$message))
    #     })
    #   })
    # })

    result_G1_CC <- reactiveVal(NULL)
    result_G1_BP <- reactiveVal(NULL)
    result_G1_MF <- reactiveVal(NULL)
    result_G1_KEGG <- reactiveVal(NULL)
    result_G2_CC <- reactiveVal(NULL)
    result_G2_BP <- reactiveVal(NULL)
    result_G2_MF <- reactiveVal(NULL)
    result_G2_KEGG <- reactiveVal(NULL)

    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())

      result_G1_CC(NULL)
      result_G1_BP(NULL)
      result_G1_MF(NULL)
      result_G2_CC(NULL)
      result_G2_BP(NULL)
      result_G2_MF(NULL)

      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()

      for (n in seq_along(groups_list)) {
        col <- groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        for (mode in c("CC", "BP", "MF")) {
          future_promise({
            start_time <- Sys.time()
            result <- go_enrich_dotplot(
              gene_list_ = unique(gene_list),
              save_path_ = NULL,
              save_filename_ = NULL,
              mode_ = mode,
              showCategory_ = 10
            )
            end_time <- Sys.time()
            list(
              r = result, c = col, m = mode,
              start_time = start_time,
              end_time = end_time,
              elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
            )
          }) %...>% {
            if (!is.null(.)) {
              var <- paste0(.$c, "_", .$m)
              cat(var, " started at:", as.character(.$start_time), "\n")
              cat(var, " ended at:", as.character(.$end_time), "\n")
              cat(var, " elapsed:", as.character(.$elapsed), " seconds\n")


              if (var == "G1_CC") result_G1_CC(.$r)
              if (var == "G1_BP") result_G1_BP(.$r)
              if (var == "G1_MF") result_G1_MF(.$r)
              if (var == "G2_CC") result_G2_CC(.$r)
              if (var == "G2_BP") result_G2_BP(.$r)
              if (var == "G2_MF") result_G2_MF(.$r)
            }
          }
        }
      }
    })


    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())

      result_G1_KEGG(NULL)
      result_G2_KEGG(NULL)

      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()

      for (n in seq_along(groups_list)) {
        col <- groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        future_promise({
          start_time <- Sys.time()
          result <- kegg_enrich_dotplot(
            gene_list_ = unique(gene_list),
            save_path_ = NULL,
            save_filename_ = NULL,
            showCategory_ = 10
          )
          end_time <- Sys.time()
          list(
            r = result, c = col,
            start_time = start_time,
            end_time = end_time,
            elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
          )
        }) %...>% {
          if (!is.null(.)) {
            var <- paste0(.$c, "_", "KEGG") # 例如 "G1_KEGG"
            cat(var, " started at:", as.character(.$start_time), "\n")
            cat(var, " ended at:", as.character(.$end_time), "\n")
            cat(var, " elapsed:", as.character(.$elapsed), " seconds\n")

            if (var == "G1_KEGG") result_G1_KEGG(.$r)
            if (var == "G2_KEGG") result_G2_KEGG(.$r)
          }
        }
      }
    })


    output$G1_CC <- renderPlot({
      validate(need(result_G1_CC(), "Please update the DEG list, then wait while it loads..."))
      result_G1_CC()
    })
    output$G1_BP <- renderPlot({
      validate(need(result_G1_BP(), "Please update the DEG list, then wait while it loads..."))
      result_G1_BP()
    })
    output$G1_MF <- renderPlot({
      validate(need(result_G1_MF(), "Please update the DEG list, then wait while it loads..."))
      result_G1_MF()
    })
    output$G2_CC <- renderPlot({
      validate(need(result_G2_CC(), "Please update the DEG list, then wait while it loads..."))
      result_G2_CC()
    })
    output$G2_BP <- renderPlot({
      validate(need(result_G2_BP(), "Please update the DEG list, then wait while it loads..."))
      result_G2_BP()
    })
    output$G2_MF <- renderPlot({
      validate(need(result_G2_MF(), "Please update the DEG list, then wait while it loads..."))
      result_G2_MF()
    })
    output$G1_KEGG <- renderPlot({
      validate(need(result_G1_KEGG(), "Please update the DEG list, then wait while it loads..."))
      result_G1_KEGG()
    })
    output$G2_KEGG <- renderPlot({
      validate(need(result_G2_KEGG(), "Please update the DEG list, then wait while it loads..."))
      result_G2_KEGG()
    })



    observeEvent(geneListReactive(), {
      req(DEG_table())
      gseaFCModuleServer("gsea_up", DEG_table = DEG_table, direction = "up")
      gseaFCModuleServer("gsea_down", DEG_table = DEG_table, direction = "down")
    })
    observe({
      req(DEG_table(), wide_data(), maeColData())
      print(head(maeColData()))
      print(head(wide_data()))
      pcaModuleServer("pca1", wide_data(), maeColData())
    })
  }

  for_run <- shinyApp(ui = ui, server = server)
  runapp <- runApp(for_run)
  
  return(runapp)
}
