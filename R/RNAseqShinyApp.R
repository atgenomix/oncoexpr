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
                DT::dataTableOutput("wide_table_dt", width = "100%")
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

            checkboxInput("use_adjP", "Use Adjusted P-value?", value = FALSE),

            actionButton("run_DEG", "Run DEG analysis")
          ),
          mainPanel(
            width = 12,
            tabsetPanel(

              tabPanel("Volcano Plot interaction",


                       interactivePlotsUI("plotVolcano")


              ),

              tabPanel("DEG Table",

                       DT::dataTableOutput('DEG_table', width = "100%")

              )
            )


          )
        )
      ),
      tabPanel(
        title = "Gene-enriched Analysis",
        sidebarPanel(
          actionButton(inputId = "generate_go", label = "Analysis"),
          width = 2
        ),
        mainPanel(
          fluidRow(
            column(
              12,
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
            column(
              12,
              h4("DOWNregulated DEGs"),
              tabsetPanel(
                tabPanel("KEGG", plotOutput("G2_KEGG")),
                tabPanel("MF", plotOutput("G2_MF")),
                tabPanel("BP", plotOutput("G2_BP")),
                tabPanel("CC", plotOutput("G2_CC"))
              )
            )
          )
        )
      ),
      tabPanel(
        title = "Heatmap",
        sidebarLayout(
          sidebarPanel(
            textInput("geneListheatmap", "Heatmap Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2,AKT1"),
            actionButton(inputId = "targetGeneID", label = "Confirm"),
            width = 2
          ),
          mainPanel(
            fluidRow(
              column(width = 6,
                     box(title = "Differential heatmap", width = NULL, solidHeader = TRUE, status = "primary",
                         originalHeatmapOutput("ht", height = 1000, containment = TRUE)
                     )
              ),
              column(width = 6,
                     id = "column2",
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
      #master <- "sc://172.18.0.1:15002"
      #method <- "spark_connect"
      #version <- "3.5"
      sc(sparklyr::spark_connect(master = master, method = method, version = version))
    })

    observe({
      req(sc())
      print("dbbrowser")
      results$db_info <- dbBrowserServer("dbBrowser1", sc())
    })

    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    DEG_table <- reactiveVal(NULL)
    DEG_summary <- reactiveVal(NULL)
    wide_data <- reactiveVal(NULL)
    maeColData <- reactiveVal(NULL)

    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      selected_db_name <- results$db_info$selected_db()

      tbl_list_promise <- future_promise({
        DBI::dbExecute(sc(), paste0("USE ", selected_db_name))
        DBI::dbGetQuery(sc(), paste0("SHOW TABLES IN ", selected_db_name))
      })

      tbl_list_promise %...>% (function(tbl_list) {
        tbls <- tbl_list$tableName

        prefix <- c("^normcounts|^exacttest|^coldata")
        tbls_with_prefix <- tbl_list[grepl(prefix , tbls),]
        results$table_list <- tbls_with_prefix

        normcount_tbls <- tbls_with_prefix[grepl("^normcounts", tbls, ignore.case = TRUE), "tableName"]
        exacttest_tbls <- tbls_with_prefix[grepl("^exacttest", tbls, ignore.case = TRUE), "tableName"]
        coldata_tbls <- tbls_with_prefix[grepl("^coldata", tbls, ignore.case = TRUE), "tableName"]

        normcount_promise <- future_promise({
          query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
          DBI::dbGetQuery(sc(), query_normcount)
        })

        exacttest_promise <- future_promise({
          query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
          DBI::dbGetQuery(sc(), query_exacttest)
        })

        coldata_promise <- future_promise({
          if (length(coldata_tbls) > 0) {
            query_coldata <- paste0("SELECT * FROM ", coldata_tbls[1])
            DBI::dbGetQuery(sc(), query_coldata)
          } else {
            normcount_promise %...>% (function(normcount) {
              generate_colData_random(normcount_data, genecol = "GeneSymbol")
            })
          }
        })

        output$normcount_table <- DT::renderDataTable({
          normcount_promise %...>% {
            colnames(.)[colnames(.) == "genes"] <- "GeneSymbol"
            normcount_data <- .[,colnames(.)!="_c0"]
            DT::datatable(normcount_data)
          }
        })

        output$exacttest_table <- DT::renderDataTable({
          exacttest_promise %...>% {
            colnames(.)[colnames(.) == "genes"] <- "GeneSymbol"
            exacttest_data <- .[,colnames(.)!="_c0"]
            DT::datatable(exacttest_data)
          }
        })

        observeEvent(results$db_info$selected_db(),
                     promise_all(normcount_data = normcount_promise, exacttest_data = exacttest_promise, coldata = coldata_promise) %...>% {
                       req(coldata, normcount_data, exacttest_data)
                       DEG_table(exacttest_data)
                       wide_data(normcount_data)
                       maeColData(coldata)
                       assay_data <- as.matrix(wide_data()[, -which(colnames(wide_data()) == "GeneSymbol")])
                       if ("GeneSymbol" %in% colnames(wide_data())) {
                         rownames(assay_data) <- wide_data()[, "GeneSymbol"]
                       }

                       deg_data <- exacttest_data
                       if ("GeneSymbol" %in% colnames(deg_data)) {
                         rownames(deg_data) <- deg_data$GeneSymbol
                       }

                       common_genes <- intersect(rownames(assay_data), rownames(deg_data))
                       assay_data <- assay_data[common_genes, , drop = FALSE]
                       deg_data_sub <- deg_data[common_genes, , drop = FALSE]
                       sample_info_table <- maeColData()
                       rownames(sample_info_table) <- colnames(assay_data) # The rownames of colData must match the colnames of assay_data

                       se_expression_matrix <- SummarizedExperiment(
                         assays = list(normCount = assay_data), # read count, TPM, COV, FPKM
                         colData = sample_info_table,
                         rowData = S4Vectors::DataFrame(deg_data_sub)
                       )

                       mae <- MultiAssayExperiment(
                         experiments = list(
                           RNAseq = se_expression_matrix   # normCoun
                         ),
                         colData = sample_info_table
                       )
                       settingMAE(mae)
                     })

        observeEvent(input$run_DEG,
                     promise_all(normcount_data = normcount_promise, exacttest_data = exacttest_promise, coldata = coldata_promise) %...>% {
                       req(exacttest_data, normcount_data, coldata)
                       DEG_table(exacttest_data)
                       wide_data(normcount_data)
                       maeColData(coldata)

                       message("run_DEG pressed: reactive values updated.")
                     })

        output$wide_table_dt <- DT::renderDataTable({
          req(wide_data())
          print("send wide data to UI")
          DT::datatable(
            wide_data(),
            options = list(pageLength = 20, autoWidth = TRUE)
          )
        })

        observeEvent(input$run_DEG, {
          mae <- settingMAE()

          output$DEG_table <- renderDT({
            datatable(
              DEG_table(),
              options = list(pageLength = 10, autoWidth = TRUE)
            )
          }, server = FALSE)
        })

        observeEvent(input$run_DEG, {
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
          normCount <- wide_data()
          volcanoData <- DEG_table()
          colData <- maeColData()
          print(str(normCount))
          print(str(volcanoData))
          print(str(colData))
          exprData <- transfExprFormat( normCount, colData)
          interactivePlotsServer("plotVolcano",
                                 volcanoData = volcanoData,
                                 exprData = exprData, params)
          print(str(exprData))
        })

        topGeneList <- reactiveVal(NULL)
        downGeneList <- reactiveVal(NULL)

        observeEvent(input$run_DEG, {
          req(DEG_table())
          DEG_table <- DEG_table()

          DEG_table_filtered <- DEG_table[DEG_table$PValue < input$pval_cut & abs(DEG_table$logFC) > input$lfc_cut, ]
          sorted_DEG <- DEG_table_filtered[order(DEG_table_filtered$logFC, decreasing = TRUE),]
          gene_list_symbol <- sorted_DEG$GeneSymbol

          topGeneList(DEG_table[DEG_table$PValue < input$pval_cut & sign(DEG_table$logFC) == 1, "GeneSymbol"])
          downGeneList(DEG_table[DEG_table$PValue < input$pval_cut & sign(DEG_table$logFC) == -1, "GeneSymbol"])

          gene_list_string <- paste(c(topGeneList(), downGeneList()), collapse = ",")
          updateTextInput(session, "geneListheatmap", value = gene_list_string)

        })

        observeEvent(input$generate_go, {
          req(topGeneList(), downGeneList(), settingMAE())
          mae <- settingMAE()
          sample_info <- colData(mae[["RNAseq"]])
          groups_list <- c("G1", "G2")
          group1_fc_gene_profile <- topGeneList()
          group2_fc_gene_profile <- downGeneList()
          for (n in seq_len(length(groups_list))) {
            col <- groups_list[n]
            gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
            for (mode in c("CC", "BP", "MF")) {
              VAR <- paste0(col, "_", mode, "GO")
              result <- go_enrich_dotplot(
                gene_list_ = unique(gene_list),
                save_path_ = NULL,
                save_filename_ = NULL,
                mode_ = mode,
                showCategory_ = 10
              )
              assign(VAR, result, envir = .GlobalEnv)
            }
          }
          output$G1_MF <- renderPlot({G1_MFGO})
          output$G1_BP <- renderPlot({G1_BPGO})
          output$G1_CC <- renderPlot({G1_CCGO})
          output$G2_MF <- renderPlot({G2_MFGO})
          output$G2_BP <- renderPlot({G2_BPGO})
          output$G2_CC <- renderPlot({G2_CCGO})
        })

        observeEvent(input$generate_go, {
          req(topGeneList(), downGeneList(), settingMAE())
          mae <- settingMAE()
          sample_info <- colData(mae[["RNAseq"]])
          groups_list <- c("G1", "G2")
          group1_fc_gene_profile <- topGeneList()
          group2_fc_gene_profile <- downGeneList()
          for (n in seq_len(length(groups_list))) {
            col <- groups_list[n]
            gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
            VAR <- paste0(col, "_", "KEGG")
            result <- kegg_enrich_dotplot(
              gene_list_ = unique(gene_list),
              save_path_ = NULL,
              save_filename_ = NULL,
              showCategory_ = 10
            )
            assign(VAR, result, envir = .GlobalEnv)
          }
          output$G1_KEGG <- renderPlot({G1_KEGG})
          output$G2_KEGG <- renderPlot({G2_KEGG})
        })

        observeEvent(input$targetGeneID, {
          req(settingMAE())
          mae <- settingMAE()
          geneList <- unlist(strsplit(input$geneListheatmap, ","))
          geneList <- trimws(geneList)
          ht <- make_heatmap_mae(mae, geneList)

          if (!is.null(ht)) {
            makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
          } else {
            output$ht_heatmap <- renderPlot({
              grid::grid.newpage()
              grid::grid.text("No data available.")
            })
          }
        })
      })
    })
  }

  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run)
}

# library(shiny)
# library(DBI)
# library(shinybusy)
# library(bslib)
# library(pheatmap)
# library(oncoexpr)
# library(sparklyr)
# library(InteractiveComplexHeatmap)
# library(promises)
# library(future)
# plan(multisession)
# RNAseqShinyAppSpark()
