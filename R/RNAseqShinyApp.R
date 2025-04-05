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
      message(sprintf("[%s] Entering dbBrowser", Sys.time()))
      dbBrowserServer("dbBrowser1", sc)
    })

    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      selected_db_name <- results$db_info$selected_db()
      message(sprintf("[%s] Selected DB: %s", Sys.time(), selected_db_name))
      DBI::dbExecute(sc, paste0("USE ", selected_db_name))
      tbl_list_query <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", selected_db_name))
      tbls <- tbl_list_query$tableName
      message(sprintf("[%s] Retrieved tables: %s", Sys.time(), paste(tbls, collapse = ", ")))

      prefix <- c("^normcounts|^exacttest|^coldata")

      tbl_list_query_prefix <- tbl_list_query[grepl(prefix, tbls), ]
      message(sprintf("[%s] Filtered tables: %s", Sys.time(), paste(tbl_list_query_prefix$tableName, collapse = ", ")))

      tbls_with_prefix <- tbl_list_query_prefix$tableName
      tbls_with_time_filter <- get_latest_file_group_df(tbls_with_prefix)
      message(sprintf("[%s] Time-filtered tables:", Sys.time()))
      print(tbls_with_time_filter)

      if (sum(tbls_with_time_filter$is_latest) == 0) {
        message(sprintf("[%s] No latest table found.", Sys.time()))
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest == FALSE, ]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest == FALSE, ]
      } else {
        message(sprintf("[%s] Latest table found.", Sys.time()))
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest == TRUE, ]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest == TRUE, ]
      }

      tbls_with_prefix_time <- summary_table$"file"
      print("====summary_table====")
      print(summary_table)
      results$table_list <- tbl_list_query_prefix_time

      normcount_tbls <- tbl_list_query_prefix_time[grepl("^normcounts", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]
      exacttest_tbls <- tbl_list_query_prefix_time[grepl("^exacttest", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]
      coldata_tbls <- tbl_list_query_prefix_time[grepl("^coldata", tbls_with_prefix_time, ignore.case = TRUE), "tableName"]

      print("====normcount_tbls====")
      print(normcount_tbls)
      print("====exacttest_tbls====")
      print(exacttest_tbls)
      print("====coldata_tbls====")
      print(coldata_tbls)
      req(normcount_tbls, exacttest_tbls, coldata_tbls)
      normcount_promise <- future_promise(
        {
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying normcounts table", start_time))
          sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
          on.exit(sparklyr::spark_disconnect(sc_conn))
          DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
          normcount_tbl <- tbl(sc_conn, normcount_tbls[1])
          normcount_tbl <- normcount_tbl %>%
                            dplyr::rename(GeneSymbol = genes) %>%
                            dplyr::select(-`_c0`) %>%
                            dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
          normcount <- collect(normcount_tbl)
          end_time <- Sys.time()

          message(sprintf("[%s] Completed normcounts query (Duration: %.2f seconds)", end_time, as.numeric(difftime(end_time, start_time, units = "secs"))))

          normcount
        },
        globals = list(master = master, method = method, version = version, normcount_tbls = normcount_tbls, selected_db_name = selected_db_name),
        seed = TRUE
      )

      exacttest_promise <- future_promise(
        {
          start_time <- Sys.time()
          message(sprintf("[%s] Start querying exacttest table", start_time))
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
                sprintf("%.4f", new_mantissa) %>% paste0("e", new_exponent)
              })
            }
            df
          })
          exacttest <- collect(exacttest_tbl)
          end_time <- Sys.time()
          message(sprintf("[%s] Completed exacttest query (Duration: %.2f seconds)", end_time, as.numeric(difftime(end_time, start_time, units = "secs"))))

          exacttest
        },
        globals = list(master = master, method = method, version = version, exacttest_tbls = exacttest_tbls, selected_db_name = selected_db_name),
        seed = TRUE
      )

      coldata_promise <-
        if (length(coldata_tbls) > 0) {
          future_promise(
            {
              start_time <- Sys.time()
              message(sprintf("[%s] Start querying coldata table", start_time))
              sc_conn <- sparklyr::spark_connect(master = master, method = method, version = version)
              on.exit(sparklyr::spark_disconnect(sc_conn))
              DBI::dbExecute(sc_conn, paste0("USE ", selected_db_name))
              query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
              coldata <- DBI::dbGetQuery(sc_conn, query_coldata)

              end_time <- Sys.time()
              message(sprintf("[%s] Completed coldata query (Duration: %.2f seconds)", end_time, as.numeric(difftime(end_time, start_time, units = "secs"))))

              coldata
            },
            globals = list(master = master, method = method, version = version, coldata_tbls = coldata_tbls, selected_db_name = selected_db_name),
            seed = TRUE
          )
        } else {
          normcount_promise %...>% (function(normcount) {
            future_promise(
              {
                generate_colData_random(normcount, genecol = "GeneSymbol")
              },
              seed = TRUE
            )
          })
        }

       withProgress(message = "Processing data...", value = 0, {
        overall_start <- Sys.time()
        incProgress(0.1, detail = "Starting normcount query")
        normcount_data <- value(normcount_promise)
        incProgress(0.3, detail = sprintf("Completed normcount query at %s", Sys.time()))
        
        incProgress(0.1, detail = "Starting exacttest query")
        exacttest_data <- value(exacttest_promise)
        incProgress(0.3, detail = sprintf("Completed exacttest query at %s", Sys.time()))
        
        incProgress(0.1, detail = "Starting coldata query")
        coldata <- value(coldata_promise)
        incProgress(0.1, detail = sprintf("Completed coldata query at %s", Sys.time()))
        
        overall_end <- Sys.time()
        message(sprintf("[Overall] All queries completed (Total Duration: %.2f seconds)",
                        as.numeric(difftime(overall_end, overall_start, units = "secs"))))
        
        # 更新 reactive values 給 UI 使用
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
      
    output$wide_table_dt <- DT::renderDataTable({
      req(wide_data())
      print("send wide data to UI")
      DT::datatable(
        wide_data(),
        options = list(pageLength = 20, autoWidth = TRUE)
      )
    })

    observe({
      req(DEG_table(), wide_data(), maeColData())

      assay_data <- as.matrix(wide_data()[, -which(colnames(wide_data()) == "GeneSymbol")])
      if ("GeneSymbol" %in% colnames(wide_data())) {
        rownames(assay_data) <- wide_data()[, "GeneSymbol"]
      }

      deg_data <- DEG_table()
      if ("GeneSymbol" %in% colnames(deg_data)) {
        rownames(deg_data) <- deg_data$GeneSymbol
      }

      output$download_DEG <- downloadHandler(
        filename = function() {
          paste("DEG_table_", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(DEG_table(), file, row.names = FALSE)
        }
      )
      common_genes <- intersect(rownames(assay_data), rownames(deg_data))
      assay_data <- assay_data[common_genes, , drop = FALSE]
      deg_data_sub <- deg_data[common_genes, , drop = FALSE]
      sample_info_table <- maeColData()
      rownames(sample_info_table) <- colnames(assay_data)

      se_expression_matrix <- SummarizedExperiment(
        assays = list(normCount = assay_data),
        colData = sample_info_table,
        rowData = S4Vectors::DataFrame(deg_data_sub)
      )

      mae <- MultiAssayExperiment(
        experiments = list(
          RNAseq = se_expression_matrix
        ),
        colData = sample_info_table
      )
      settingMAE(mae)
    })

    output$DEG_table <- renderDT(
      req(results$coldata, results$normcount_data, results$exacttest_data)
        {
          datatable(
            DEG_table(),
            options = list(pageLength = 20, autoWidth = TRUE)
          )
        }
      )

    observe({
      req(results$exacttest_data, results$normcount_data, results$coldata)
      DEG_table(results$exacttest_data)
      wide_data(results$normcount_data)
      maeColData(results$coldata)
      message("assign reactiveVal: DEG_table, wide_data, maeColData")
      params <- reactive({
        list(
          lfc_cut    = input$lfc_cut,
          pval_cut   = input$pval_cut,
          pointSize  = input$pointSize,
          ptAlpha    = input$ptAlpha,
          labelSize  = input$labelSize,
          topN       = input$topN,
          use_adjP   = input$use_adjP
        )
      })


      normCount <- wide_data()
      volcanoData <- DEG_table()
      colData <- maeColData()

      exprData <- transfExprFormat(normCount, colData)


      DEG_table_data <- DEG_table()
      topGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC < -input$lfc_cut, "GeneSymbol"]

      geneListVec <- c(topGenes, downGenes)
      if (!is.null(geneListReactive)) {
        selected_gene <- mod_geneSelector_server("gene_selector", volcanoData, geneListVec)
      } else {
        (selected_gene <- NULL)
      }
      interactivePlotsServer("volcano_plots", volcanoData = volcanoData, exprData = exprData, params = params, selectedGene = selected_gene)
    })

    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL)


    geneListReactive <- eventReactive(input$run_DEG, {
      req(DEG_table(), maeColData(), wide_data())

      DEG_table_data <- DEG_table()
      topGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC < -input$lfc_cut, "GeneSymbol"]

      topGeneList(topGenes)
      downGeneList(downGenes)
      print("=====Case DEG List=====")
      print(length(topGeneList()))
      print("=====Control DEG List=====")
      print(length(downGeneList()))

      gene_list <- paste(c(topGenes, downGenes), collapse = ",")
      gene_list
    })

    observeEvent(geneListReactive(), {
      req(geneListReactive(), settingMAE())
      mae <- settingMAE()

      geneListVec <- unlist(strsplit(geneListReactive(), ","))
      geneListVec <- trimws(geneListVec)
      ht <- make_heatmap_mae(mae, geneListVec)

      if (!is.null(ht)) {
        makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
      } else {
        output$ht_heatmap <- renderPlot({
          grid::grid.newpage()
          grid::grid.text("No data available.")
        })
      }
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
