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
#' @importFrom ggpubr color_palette
#' @importFrom enrichplot color_palette
#' @importFrom DT dataTableOutput renderDataTable


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

RNAseqShinyAppSpark <- function() {

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
        sidebarLayout(
          sidebarPanel(
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
            tabsetPanel(
              tabPanel("Volcano Plot interaction",
                       #plotOutput("volcano_plot", height = "600px")
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
                tabPanel("MF", plotOutput("G1_MF")),
                tabPanel("BP", plotOutput("G1_BP")),
                tabPanel("CC", plotOutput("G1_CC")),
                tabPanel("KEGG", plotOutput("G1_KEGG"))
              )
            )
          ),
          fluidRow(
            column(
              12,
              h4("DOWNregulated DEGs"),
              tabsetPanel(
                tabPanel("MF", plotOutput("G2_MF")),
                tabPanel("BP", plotOutput("G2_BP")),
                tabPanel("CC", plotOutput("G2_CC")),
                tabPanel("KEGG", plotOutput("G2_KEGG"))
              )
            )
          )
        )
      ),
      tabPanel(
        title = "Target Gene Expression",
        sidebarPanel(
          textInput("geneList", "Target Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2"),
          actionButton(inputId = "targetGeneID", label = "Confirm"),
          width = 2
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Target Gene Expr. Table", 
                     DT::dataTableOutput('target_gene_table', width = "100%", height = "600px")
            ),
            tabPanel(
              title = "Heatmap",
              sidebarPanel(
                textInput("geneListheatmap", "Heatmap Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2,AKT1"),
                width = 2
              ),
              mainPanel(
                originalHeatmapOutput("ht", height = 800, containment = TRUE)
              )
            )
          )
        )
       )

    )
  )

  server <- function(input, output, session) {
    
    sc <- reactiveVal(NULL)
    
    observe({
      #master <- "sc://172.18.0.1:15002"
      #method <- "spark_connect"
      master <- "local"
      method <- "shell"
      version <- "3.5"
      connection <- sparklyr::spark_connect(master = master, method = method, version = version)
      sc(connection)
    })
    
    results <- reactiveValues(
      db_info = NULL,
      table_list = NULL,
      grouplist = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      coldata = NULL
    )
    
    observe({
      req(sc())
      print("test null sc")
      results$db_info <- dbBrowserServer("dbBrowser1", sc())
    })
    
    observeEvent(results$db_info$selected_db(), {
      req(results$db_info$selected_db())
      selected_db_name <- results$db_info$selected_db()
      
      DBI::dbExecute(sc(), paste0("USE ", selected_db_name))
      tbl_list_query <- DBI::dbGetQuery(sc(), paste0("SHOW TABLES IN ", selected_db_name))
      tbls <- tbl_list_query$tableName      

      prefix <- c("^normcounts|^exacttest|^coldata")
      
      tbls_with_prefix <- tbl_list_query[grepl(prefix , tbls),]
      results$table_list <- tbls_with_prefix
      
      normcount_tbls <- tbls_with_prefix[grepl("^normcounts", tbls, ignore.case = TRUE), "tableName"]
      exacttest_tbls <- tbls_with_prefix[grepl("^exacttest", tbls, ignore.case = TRUE), "tableName"]
      coldata_tbls <- tbls_with_prefix[grepl("^coldata", tbls, ignore.case = TRUE), "tableName"]

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

      }else{
        colData <- generate_colData_random(results$normcount_data, genecol = "GeneSymbol") #pseudo coldata
        results$coldata <- colData
      }
      print(str(results$normcount_data))
      print(str(results$exacttest_data))
      colnames(results$exacttest_data)[colnames(results$exacttest_data) == "genes"] <- "GeneSymbol"
      colnames(results$normcount_data)[colnames(results$normcount_data) == "genes"] <- "GeneSymbol"
      results$normcount_data <- results$normcount_data[,colnames(results$normcount_data)!="_c0"]
      print(str(results$exacttest_data))

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
      rownames(sample_info_table) <- colnames(assay_data) # The rownames of colData must match the colnames of assay_data

      se_expression_matrix <- SummarizedExperiment(
        assays = list(normCount = assay_data), #read count, TPM, COV, FPKM
        colData = sample_info_table,
        rowData = S4Vectors::DataFrame(deg_data_sub)
      )

      mae <- MultiAssayExperiment(
        experiments = list(
          RNAseq = se_expression_matrix   #  normCoun
        ),
        colData = sample_info_table
      )
      settingMAE(mae)
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
          exprData <- computeExprData( normCount, colData)
          interactivePlotsServer("plotVolcano", 
                                volcanoData = volcanoData, 
                                exprData = exprData, params)
          print(str(exprData))
    })
    observeEvent(input$run_DEG, {
        req(results$exacttest_data, results$normcount_data, results$coldata)
        DEG_table(results$exacttest_data)
        wide_data(results$normcount_data)
        maeColData(results$coldata)

        message("run_DEG pressed: reactive values updated.")

    })
    # reactive_volcano_plot <- eventReactive(input$run_DEG, {
    #   req(DEG_table(), maeColData(), wide_data())

    #   ggvolcano_custom(
    #     df = DEG_table(),
    #     geneName = DEG_table()$GeneSymbol,
    #     pValCol = "PValue",
    #     logFCCol = "logFC",
    #     coef = 2,
    #     lfc_cut = input$lfc_cut,
    #     pval_cut = input$pval_cut,
    #     useAdjP = FALSE,
    #     title = "Volcano Plot",
    #     topN = input$topN,
    #     geneCol = NULL,
    #     pointSize = input$pointSize, 
    #     ptAlpha = input$ptAlpha,
    #     labelSize = input$labelSize 
    #   )
    # })
    
    # output$volcano_plot <- renderPlot({
    #   reactive_volcano_plot()
    # })
    
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
      updateTextInput(session, "geneList", value = gene_list_string)
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
      req(settingMAE(), wide_data())
      mae <- settingMAE()
      geneList <- unlist(strsplit(input$geneList, ","))
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- rownames(sample_info) 
      expr_profile <- as.data.frame(wide_data())
      if ("GeneSymbol" %in% colnames(expr_profile)) {
        rownames(expr_profile) <- expr_profile[,"GeneSymbol"]
        expr_profile <- expr_profile[,-1]
        print("wide table has GeneSymbol column as rownames")
      }
      targetGeneExpr <- target_exprofile( 
        geneList_ = geneList, 
        groups_list_ = groups_list,
        expr_profile_ = expr_profile
      )
      output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
    })
    
    observeEvent(input$targetGeneID, {
      req(settingMAE())
      mae <- settingMAE()
      geneList <- unlist(strsplit(input$geneList, ","))
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

  }
  
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run)
}
