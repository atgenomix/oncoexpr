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
#'   library(sparklyr)
#'   sc <- spark_connect(master = "local")
#'   input_data <- read.csv("path/to/rnaseq_data.csv")
#'   RNAseqShinyAppSpark(input_data, sc, "RNAseq Analysis", 8080)
#' }
#'
#' @export


RNAseqShinyAppSpark <- function(master = "sc://172.18.0.1:15002", method = "spark_connect", version = "3.5") {

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
                    tabPanel("normCount Table",
                            DT::dataTableOutput("wide_table_dt", width = "100%")
                    ),
                    tabPanel("DEG Table",
                             DT::dataTableOutput('DEG_table', width = "100%")
                    )

                )

              )
            )
          )
        )
      ),
      tabPanel(
        title = "Differential Expression Analysis",
        layout_sidebar(
          full_screen = TRUE,
          sidebar = sidebar(
            style = "min-height: 600px; overflow-y: auto;",
            sliderInput("lfc_cut", "Fold Change Threshold (log2):",
                        min = 0, max = 10, value = 1, step = 0.1),
            sliderInput("pval_cut", "p-value Threshold:",
                        min = 0.001, max = 1, value = 0.05, step = 0.001),
            sliderInput("pointSize", "Point Size:",
                        min = 1, max = 5, value = 2, step = 0.5),
            sliderInput("ptAlpha", "Transparent:",
                        min = 0.1, max = 1, value = 0.6, step = 0.1),
            sliderInput("labelSize", "Gene Label Size:",
                        min = 1, max = 6, value = 3, step = 0.5),
            numericInput("topN", "Label N number of Genes (0 is no labeling):",
                         value = 100, min = 0, max = 1000),

            actionButton("run_DEG", "Update DEG")
          ),
          mainPanel(
            width = 12,
            tabsetPanel(
              tabPanel("Volcano Plot", 
                       interactivePlotsUI("plotVolcano")
              ),

              tabPanel("Gene Set Enrichment",
                        fluidRow(
                          column(2, 
                                textInput("geneLisEnrichment", "DEG Gene List",  value = "EGFR,ESR1,KRAS,ERBB2,AKT1")
                       
                          )
                        ),
                        fluidRow(
                          column(12,
                                h4("UPregulated DEGs"),
                                tabsetPanel(
                                  tabPanel("KEGG", plotOutput("G1_KEGG")),
                                  tabPanel("MF", plotOutput("G1_MF")),
                                  tabPanel("BP", plotOutput("G1_BP")),
                                  tabPanel("CC", plotOutput("G1_CC"))
                                )
                          )
                        ),
                        fluidRow(
                          column(12,
                                h4("DOWNregulated DEGs"),
                                tabsetPanel(
                                  tabPanel("KEGG", plotOutput("G2_KEGG")),
                                  tabPanel("MF", plotOutput("G2_MF")),
                                  tabPanel("BP", plotOutput("G2_BP")),
                                  tabPanel("CC", plotOutput("G2_CC"))
                                )
                          )
                        )
              ),
              tabPanel("Heatmap",
                        fluidRow(
                          column(2,
                                textInput("geneListheatmap", "DEG Gene List", value = "EGFR,ESR1,KRAS,ERBB2,AKT1")
           
                          )
                        ),
                        fluidRow(
                          column(width = 6,
                                box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                                    originalHeatmapOutput("ht", height = 1000, containment = TRUE)
                                )
                          ),
                          column(width = 6,
                                box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                                    subHeatmapOutput("ht", title = NULL, containment = TRUE)
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

    sc <- reactiveVal(NULL)
    results <- reactiveValues(
      db_info = NULL,
      table_list = NULL,
      grouplist = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      coldata = NULL
    )

    observe({
      sc(sparklyr::spark_connect(master = master, method = method, version = version))
    })

    observe({
      req(sc())
      print("dbbrowser")
      results$db_info <- dbBrowserServer("dbBrowser1", sc())
    })

    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      selected_db_name <- results$db_info$selected_db()

      DBI::dbExecute(sc(), paste0("USE ", selected_db_name))
      tbl_list_query <- DBI::dbGetQuery(sc(), paste0("SHOW TABLES IN ", selected_db_name))
      tbls <- tbl_list_query$tableName
      print("====tbls====")
      print(tbls)

      prefix <- c("^normcounts|^exacttest|^coldata")

      tbl_list_query_prefix <- tbl_list_query[grepl(prefix, tbls),]
      print("====tbl_list_query_prefix====")
      print(tbl_list_query_prefix)
      tbls_with_prefix <- tbl_list_query_prefix$tableName
      print("====tbls_with_prefix====")
      print(tbls_with_prefix)

      tbls_with_time_filter <- get_latest_file_group_df(tbls_with_prefix)
      print("====tbls_with_time_filter====")
      print(tbls_with_time_filter)

      if(sum(tbls_with_time_filter$is_latest)==0){
        print("no latest table")
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest==FALSE,]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest==FALSE, ]
      } else {
        print("latest table")
        tbl_list_query_prefix_time <- tbl_list_query_prefix[tbls_with_time_filter$is_latest==TRUE,]
        summary_table <- tbls_with_time_filter[tbls_with_time_filter$is_latest==TRUE,]
      }
      
      tbls_with_prefix_time <- summary_table$"file"
      print("====tbls_with_prefix_time====")
      print(tbls_with_prefix_time)
      print("====summary_table====")
      print(summary_table)
      print("====tbl_list_query_prefix_time====")  
      print(tbl_list_query_prefix_time)
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

      if (length(normcount_tbls) > 0) {
        query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
        results$normcount_data <- DBI::dbGetQuery(sc(), query_normcount)
      }

      if (length(exacttest_tbls) > 0) {
        query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
        results$exacttest_data <- DBI::dbGetQuery(sc(), query_exacttest)
      }

      if (length(coldata_tbls) > 0) {
        query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
        results$coldata <- DBI::dbGetQuery(sc(), query_coldata)
      } else {
        colData <- generate_colData_random(results$normcount_data, genecol = "GeneSymbol")
        results$coldata <- colData
      }
      colnames(results$exacttest_data)[colnames(results$exacttest_data) == "genes"] <- "GeneSymbol"
      colnames(results$normcount_data)[colnames(results$normcount_data) == "genes"] <- "GeneSymbol"
      results$normcount_data <- results$normcount_data[,colnames(results$normcount_data)!="_c0"]
      results$exacttest_data <- results$exacttest_data[,colnames(results$exacttest_data)!="_c0"]
    })

    output$normcount_table <- DT::renderDataTable({
      req(results$normcount_data)
      DT::datatable(results$normcount_data)
    })

    output$exacttest_table <- DT::renderDataTable({
      req(results$exacttest_data)
      DT::datatable(results$exacttest_data)
    })

    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    DEG_table <- reactiveVal(NULL)
    DEG_summary <- reactiveVal(NULL)
    wide_data <- reactiveVal(NULL)
    maeColData <- reactiveVal(NULL)


    output$wide_table_dt <- DT::renderDataTable({
      req(wide_data())
      print("send wide data to UI")
      DT::datatable(
        wide_data(),
        options = list(pageLength = 20, autoWidth = TRUE)
      )
    })

    observeEvent(results$db_info$selected_db(), {
      req(results$coldata, results$normcount_data, results$exacttest_data)
      DEG_table(results$exacttest_data)
      wide_data(results$normcount_data)
      maeColData(results$coldata)
      assay_data <- as.matrix(wide_data()[, -which(colnames(wide_data()) == "GeneSymbol")])
      if ("GeneSymbol" %in% colnames(wide_data())) {
        rownames(assay_data) <- wide_data()[, "GeneSymbol"]
      }

      deg_data <- results$exacttest_data
      if ("GeneSymbol" %in% colnames(deg_data)) {
        rownames(deg_data) <- deg_data$GeneSymbol
      }

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

    observe({
      req(results$coldata, results$normcount_data, results$exacttest_data)
      mae <- settingMAE()
      output$DEG_table <- renderDT({
        datatable(
          DEG_table(),
          options = list(pageLength = 20, autoWidth = TRUE)
        )
      }, server = FALSE)
    })

    observe({
          req(DEG_table(), maeColData(), wide_data())
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
          print("params")
          print(params)

          normCount <- wide_data()
          volcanoData <- DEG_table()
          colData <- maeColData()
          print(str(normCount))
          print(str(volcanoData))
          print(str(colData))
          exprData <- transfExprFormat(normCount, colData)
          interactivePlotsServer("plotVolcano",
                                volcanoData = volcanoData,
                                exprData = exprData, params)
          print(str(exprData))
    })


    
    observe({
        req(results$exacttest_data, results$normcount_data, results$coldata)
        DEG_table(results$exacttest_data)
        wide_data(results$normcount_data)
        maeColData(results$coldata)

        message("run_DEG pressed: reactive values updated.")
    })

    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL)


    geneListReactive <- eventReactive(input$run_DEG, {
      print("geneListReactive section before req")
      req(DEG_table(), maeColData(), wide_data())
      DEG_table_data <- DEG_table()
      print("geneListReactive section 1")
      print(topGeneList())
      print(downGeneList())

      topGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC > input$lfc_cut, "GeneSymbol"]
      downGenes <- DEG_table_data[DEG_table_data$PValue < input$pval_cut & DEG_table_data$logFC < -input$lfc_cut, "GeneSymbol"]
      topGeneList(topGenes)
      downGeneList(downGenes)
      print("geneListReactive section 2")
      print(topGeneList())
      print(downGeneList())

      gene_list <- paste(c(topGenes, downGenes), collapse = ",")
      gene_list
      print("geneListReactive section 3")
      print(gene_list)
    })


    observe({
      new_gene_list <- geneListReactive()
      updateTextInput(session, "geneListheatmap", value = new_gene_list)
      updateTextInput(session, "geneLisEnrichment", value = new_gene_list)
      print("update geneListheatmap and geneLisEnrichment")
      print(new_gene_list)
    })



    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())
      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      

      for (n in seq_len(length(groups_list))) {
        col <- groups_list[n]
        gene_list <- list(topGeneList(), downGeneList())[[n]]
        for (mode in c("CC", "BP", "MF")) {
          VAR <- paste0(col, "_", mode, "GO")
          result <- tryCatch({
            go_enrich_dotplot(
              gene_list_ = unique(gene_list),
              save_path_ = NULL,
              save_filename_ = NULL,
              mode_ = mode,
              showCategory_ = 10
            )
          }, error = function(e) {

            shiny::showNotification(
              paste("GO enrichment error in", col, "mode", mode, ":", e$message),
              type = "error",
              duration = 5
            )
            return(NULL)
          })
          assign(VAR, result, envir = .GlobalEnv)
        }
      }
      output$G1_MF <- renderPlot({ G1_MFGO })
      output$G1_BP <- renderPlot({ G1_BPGO })
      output$G1_CC <- renderPlot({ G1_CCGO })
      output$G2_MF <- renderPlot({ G2_MFGO })
      output$G2_BP <- renderPlot({ G2_BPGO })
      output$G2_CC <- renderPlot({ G2_CCGO })
    })


    observeEvent(geneListReactive(), {
      req(topGeneList(), downGeneList(), settingMAE())
      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- c("G1", "G2")
      
      for (n in seq_len(length(groups_list))) {
        col <- groups_list[n]
        gene_list <- list(topGeneList(), downGeneList())[[n]]
        VAR <- paste0(col, "_", "KEGG")
        result <- tryCatch({
          kegg_enrich_dotplot(
            gene_list_ = unique(gene_list),
            save_path_ = NULL,
            save_filename_ = NULL,
            showCategory_ = 10
          )
        }, error = function(e) {
          shiny::showNotification(
            paste("KEGG enrichment error in", col, ":", e$message),
            type = "error",
            duration = 5
          )
          return(NULL)
        })
        assign(VAR, result, envir = .GlobalEnv)
      }
      output$G1_KEGG <- renderPlot({ G1_KEGG })
      output$G2_KEGG <- renderPlot({ G2_KEGG })
    })


    observeEvent(geneListReactive(), {
      req(geneListReactive(), settingMAE())
      mae <- settingMAE()
      print("Heatmap section 4")
      print(geneListReactive())
      
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
  }
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run)

}



