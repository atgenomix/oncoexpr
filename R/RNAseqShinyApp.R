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
#' @importFrom ggpubr color_palette
#' @importFrom enrichplot color_palette
#' @importFrom DT dataTableOutput renderDataTable



#' @title shiny app for this analysis 
#' @description start the shiny app
#' @param host_ host name
#' @param port_ port number
#' @return start the UI and Server for analysis
#' @export
oncoExprApp <- function(host_ = NULL, port_ = NULL){
  
  options(shiny.maxRequestSize = 500*1024^2)  # 增加上傳限制為 100 MB
  ui <- fluidPage(
    navbarPage("Acrocyte RNAseq App (Beta)",
               
               tabPanel("Gene Expression (RNA-seq)",
                        sidebarPanel(
                          textInput("prefixSample", "Sample Code", value = "BC021"),
                          hr(),
                          textInput(inputId ="groupNameS1", "Sample 1 ", value = "tissue"),
                          textInput(inputId ="groupNameS2", "Sample 2", value = "ctc"),
                          fileInput(inputId = "gtfFileS1", "Upload Sample 1 file", buttonLabel = "Upload..."),
                          fileInput(inputId = "gtfFileS2", "Upload Sample 2 file", buttonLabel = "Upload..."),
                          actionButton("prepare_data", "準備資料"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Sample 1 Expr", DT::dataTableOutput('ExprTableS1',width="100%")),
                            tabPanel("Sample 2 Expr", DT::dataTableOutput('ExprTableS2',width="100%"))
                          )
                        )
                        
               ),
               tabPanel("Target Gene Expression",
                        sidebarPanel(
                          textInput("geneList", "Target Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2,AKT1,PIK3CA,ERBB3,CCND1,SF3B1,FGFR1,FBXW7,TP53,BRCA1,BRCA2"),
                          actionButton(inputId = "targetGeneID", label = "Confirm"),width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Target Gene Expr. Table", DT::dataTableOutput('target_gene_table', width="100%", height = "600px")),
                          )
                        )
               ),
               tabPanel("Expression Pattern", 
                        sidebarPanel(
                          actionButton(inputId = "generate_scatter", label = "Expression Analysis"),
                          width=2
                        ),
                        mainPanel(
                          width = 10, 
                          tabsetPanel(
                            tabPanel("Comparasion of Expr. Pattern", plotOutput("scatter_plot", width = "1000px", height = "800px"))
                          )
                        )
               ),
               tabPanel("Gene-enriched Analysis", 
                        sidebarPanel(
                          actionButton(inputId = "generate_go", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            id = "mainTabs",
                            
                            # Sub1 tab
                            tabPanel("Enrichment",
                                     fluidRow(
                                       column(12,
                                              # Upper section tabs for Sub1
                                              h4("Tissue"),
                                              tabsetPanel(
                                                id = "sub1Upper",
                                                tabPanel("MF", plotOutput("tissueMF")),
                                                tabPanel("BP", plotOutput("tissueBP")),
                                                tabPanel("CC", plotOutput("tissueCC")),
                                                tabPanel("KEGG", plotOutput("tissueKEGG"))
                                              )
                                       )
                                     ),
                                     fluidRow(
                                       column(12,
                                              # Lower section tabs for Sub1
                                              h4("CTC"),
                                              tabsetPanel(
                                                id = "sub1Lower",
                                                tabPanel("MF", plotOutput("ctcMF")),
                                                tabPanel("BP", plotOutput("ctcBP")),
                                                tabPanel("CC", plotOutput("ctcCC")),
                                                tabPanel("KEGG", plotOutput("ctcKEGG"))
                                                
                                              )
                                       )
                                     )
                            )
                            
                            
                          )
                        )
               ),
               tabPanel("TissueTyper (GTEx)", 
                        sidebarPanel(
                          actionButton(inputId = "generate_tissuetyper", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Normal Tissue Analysis", plotOutput("gtexCorPlot", width = "800px", height = "800px")),
                            tabPanel("Bulk sample cor. rank", DT::dataTableOutput('rankTable1',width="100%")),
                            tabPanel("CTC sample cor. rank", DT::dataTableOutput('rankTable2',width="100%"))
                          )
                        )
               ),
               tabPanel("Immune Fraction", 
                        sidebarPanel(
                          actionButton(inputId = "generate_immune", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Immune Cell Type Prediction", plotOutput("immune_plot", width = "1000px", height = "800px"))
                          )
                        )
               )
               
    )
  )
  
  server <- function(input, output, session) {
    #設定所需檔案之路徑
    #disease_list_path <- system.file("extdata", "0_疾病清單.xlsx", package = "alprogar")
    #rsID_and_Gene_name_path <- file.path(.libPaths()[1],'alprogar', 'extdata',"/")                
    #download_path <- file.path(path.expand('~') , "alprogar")
    required_pkgs <- c("pcaPP", "reshape2", "stringr", "readxl","ggplot2","tidyr","tibble", "viridis","RColorBrewer","pheatmap","ggpubr","ggrepel","readr","dplyr","sparklyr","org.Hs.eg.db","enrichplot","clusterProfiler")
    invisible(lapply(required_pkgs, library, character.only = TRUE))
    gtex_expr_path <- system.file("extdata", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", package = "oncoexpr")
    
    results <- reactiveValues(
      geneList = NULL,
      date = NULL,
      groups_list = NULL,
      group1_ = NULL,
      group2_ = NULL,
      expr_data = NULL,
      scatter_plot = NULL,
      immune_plot = NULL,
      go_dotplot = list(),
      kegg_dotplot = NULL,
      tissue_heatmap = NULL,
      special_gene_list = NULL,
      group1_fc_gene_profile = NULL,
      group2_fc_gene_profile = NULL
    )
    #showtext_auto() 
    
    observeEvent(input$prepare_data, {
      req(input$gtfFileS1, input$gtfFileS2)
      date <- Sys.Date() 
      prefixSample <- input$prefixSample
      save_path <- input$save_path
      gtexFilePath <- input$gtexFilePath
      
      geneList <- unlist(strsplit(input$geneList, ","))
      print(input$gtfFileS1$datapath)
      print(input$gtfFileS2$datapath)
      
      basenameS1 <- basename(input$gtfFileS1$datapath)
      basenameS2 <- basename(input$gtfFileS2$datapath)
      dirnameS1 <- dirname(input$gtfFileS1$datapath)
      dirnameS2 <- dirname(input$gtfFileS2$datapath)
      
      groups_list <- c(input$groupNameS1, input$groupNameS2)
      files_list <- c(basenameS1, basenameS2)
      inputpath <- c(dirnameS1, dirnameS2)
      input_and_prep_file_directly(c(input$gtfFileS1$datapath , input$gtfFileS2$datapath), inputpath, files_list, groups_list)
      
      group1_ <- get(paste0(groups_list[1], "_rna_expr"))
      group2_ <- get(paste0(groups_list[2], "_rna_expr"))
      output$ExprTableS1 <- DT::renderDataTable({group1_ })
      output$ExprTableS2 <- DT::renderDataTable({group2_ })
      # 分析數據
      Prep4ExprComp_result <- Prep4ExprComp(group1 = group1_, group2 = group2_, groups_name = groups_list)
      data_v2 <- Prep4ExprComp_result$result
      results$special_gene_list <- Prep4ExprComp_result$special_gene_list
      
      # 篩選顯著基因表現譜
      results$group1_fc_gene_profile <- data_v2[((data_v2[[groups_list[1]]] + 1) / (data_v2[[groups_list[2]]] + 1)) > 1,]
      results$group2_fc_gene_profile <- data_v2[((data_v2[[groups_list[2]]] + 1) / (data_v2[[groups_list[1]]] + 1)) > 1,]
      results$group1_ <- group1_
      results$group2_ <- group2_
      # 保存資料
      results$expr_data <- data_v2
      results$geneList <- geneList
      results$date_ <- date 
      results$groups_list <- groups_list
      
    })
    
    observeEvent(input$targetGeneID, {
      req(results$expr_data, results$date_, results$groups_list)
      saveFileName <- paste0(input$prefixSample, "_target_gene_expr_", results$date_ ,".csv")
      geneList <- unlist(strsplit(input$geneList, ","))
      targetGeneExpr <- target_exprofile( 
        geneList_ = geneList, 
        groups_list_ = results$groups_list, 
        save_path_ = input$save_path, 
        save_filename_ = saveFileName, 
        expr_profile_ = results$expr_data
      )
      output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
    })
    observeEvent(input$generate_scatter, {
      req(results$expr_data, results$special_gene_list)
      
      # 調用 expr_pattern 函數，直接返回 ggplot 對象
      output$scatter_plot <- renderPlot({
        expr_pattern(
          data_v2 = results$expr_data, 
          special_gene_list = results$special_gene_list, 
          save_path_ = input$save_path, 
          save_filename_ = paste0(input$prefixSample, "_scatter_plot.png")
        )
      })
    })
    
    
    observeEvent(input$generate_tissuetyper, {
      req(results$group1_, results$group2_, results$groups_list) 
      output$gtexCorPlot <-  renderPlot({
        tissueTyperResult <- TissueTyper( 
          group1_ = results$group1_,
          group2_ = results$group2_,
          save_path = input$save_path, 
          gtexFilePath = gtex_expr_path, 
          groups_list = results$groups_list, 
          sample_id = input$prefixSample 
        )
        tissueTyperResult$plot
        output$rankTable1 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank1)})
        output$rankTable2 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank2)})
        
      })
    })
    
    # 生成免疫分量分析圖
    observeEvent(input$generate_immune, {
      req(results$expr_data, results$group1_, results$group2_, results$groups_list) 
      output$immune_plot <-  renderPlot({
        ImmuneFractionPlot(
          
          group1_ = results$group1_,
          group2_ = results$group2_,
          groups_list = results$groups_list,
          save_path = input$save_path
          
        )
        
      })
    })
    
    observeEvent(input$generate_go, {
      req(results$groups_list, results$group1_fc_gene_profile, results$group2_fc_gene_profile, results$date_)
      group1_fc_gene_profile <- results$group1_fc_gene_profile
      group2_fc_gene_profile <- results$group2_fc_gene_profile
      prefixSample <- input$prefixSample
      for(n in seq_len(length(results$groups_list)) ){
        print(col)
        col <- results$groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        save_filename <- paste0(paste(prefixSample, col, sep="_"), "_go_enrich_plot_", results$date_ ,".png")
        print("setting ok")
        print(save_filename)
        for( mode in c("CC", "BP", "MF")){
          print(mode)
          VAR <- paste0(col, mode, "GO")
          assign(VAR, go_enrich_dotplot( 
            
            gene_list_= gene_list, 
            save_path_ = input$save_path,
            save_filename_ = save_filename, 
            mode_ = mode, 
            showCategory_ = 10), envir=.GlobalEnv
            
          )
          
          
        }
        
      }
      output$tissueMF <- renderPlot({tissueMFGO})
      output$tissueBP <- renderPlot({tissueBPGO})
      output$tissueCC <- renderPlot({tissueCCGO})
      output$ctcMF <- renderPlot({ctcMFGO})
      output$ctcBP <- renderPlot({ctcBPGO})
      output$ctcCC <- renderPlot({ctcCCGO})
      
    })
    
    observeEvent(input$generate_go, {
      req(results$groups_list, results$group1_fc_gene_profile, results$group2_fc_gene_profile, results$date_)
      group1_fc_gene_profile <- results$group1_fc_gene_profile
      group2_fc_gene_profile <- results$group2_fc_gene_profile
      prefixSample <- input$prefixSample
      for(n in seq_len(length(results$groups_list)) ){
        print(col)
        col <- results$groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        save_filename <- paste0(paste(prefixSample, col, sep="_"), "_kegg_enrich_plot_", results$date_ ,".png")
        print(col)
        VAR <- paste0(col,"KEGG")
        assign(VAR, kegg_enrich_dotplot ( 
          gene_list_ = gene_list, 
          save_path_ = input$save_path,
          save_filename_ = save_filename, 
          showCategory_ = 10
        )
        )
      }
      
      output$tissueKEGG <- renderPlot({tissueKEGG})
      output$ctcKEGG <- renderPlot({ctcKEGG})
      
    })
    
    
  }
  
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run, host = host_, port = port_)
  
}



#' @title shiny app for RNAseq  
#' @description start the shiny app with spark connection
#' @param host_ host name
#' @param port_ port number
#' @return start the UI and Server for analysis
#' @export
oncoExprAppSpark <- function(host_ = NULL, port_ = NULL){
  options(shiny.maxRequestSize = 10 * 1024^5) # 10 MB
  ui <- fluidPage(
    navbarPage("Acrocyte RNAseq App (Beta)",
               tabPanel("Spark DB browser",
                        tags$style(".shinybusy-overlay {opacity: 0.7; background-color: #7c7c7c;}"),
                        add_busy_spinner(
                          spin = "fading-circle",
                          position = "full-page",
                          timeout = 1000
                        ),
                        sidebarLayout(
                          sidebarPanel(
                            # 第一塊：Spark DB 資料庫操作
                            wellPanel(
                              h4("Spark DB Data"),
                              #selectInput("selected_db", "Choose Database:",
                              #            choices = character(0),    # 先給空
                              #            selected = NULL),
                              #uiOutput("table_selector"),
                              #actionButton("preview_btn", "Show Table"),
                              #hr(),
                              #以shinymodule取代
                              #textInput("mainCode_input", "Main Code:", value = "BC026"),
                              #textInput("subCode_input", "Sub Code:", value = "PBMC"),
                              #selectInput("normMethod_input", "Normalization Method:",
                              #            choices = c("max_Cov", "max_FPKM", "max_TPM"),
                              #            selected = "max_TPM"),
                              #actionButton("preview_2_btn", "Query & Preview"),
                              dbBrowserUI("dbBrowser1"),
                              sampleSelectionUI("sampleSelector"),
                              actionButton("mae_start", "build MAE")
                            ),
                            
                            # 第二塊：原本 "Prepare Samples" 的功能，放在同一 sidebarPanel 裡
                            wellPanel(
                              h4("Local Upload / Prepare Samples"),
                              fileInput(
                                "files",
                                "選擇多個檔案上傳",
                                multiple = TRUE,
                                accept = c(".csv", ".txt", ".gtf")
                              ),
                              actionButton("enable_edit", "啟用編輯"),
                              actionButton("add_row", "新增空白列"),
                              actionButton("add_column", "新增空白欄位"),
                              actionButton("rename_column", "修改欄位名稱"),
                              actionButton("prepare_data", "準備資料")
                            )
                          ),
                          
                          
                          mainPanel(
                            tabsetPanel(
                              tabPanel("Sample Information",
                                       # 其他想顯示的元件放這裡...
                                       h3("樣本資訊"),
                                       DTOutput("file_table"),
                                       downloadButton("download_table", "下載註記資料"),
                                       verbatimTextOutput("error_message") # 顯示錯誤訊息
                              ),
                              tabPanel("Entire Table",
                                       DT::dataTableOutput("table_preview_entire", width = "100%")
                              ),
                              tabPanel("Filtered Table",
                                       #DTOutput
                                       DT::dataTableOutput("table_preview_filtered", width = "100%")
                              )
                              
                            )
                          )
                          
                        )
               ),
               tabPanel("Gene Expression (RNA-seq)",
                        tabsetPanel(
                          tabPanel("Selecting uploaded samples",
                                   sidebarPanel(
                                     textInput("prefixSample", "Sample Code", value = "BC021"),
                                     hr(),
                                     textInput(inputId ="groupNameS1", "Sample 1 ", value = "tissue"),
                                     textInput(inputId ="groupNameS2", "Sample 2", value = "ctc"),
                                     #fileInput(inputId = "gtfFileS1", "Upload Sample 1 file", buttonLabel = "Upload..."),
                                     #fileInput(inputId = "gtfFileS2", "Upload Sample 2 file", buttonLabel = "Upload..."),
                                     checkboxGroupInput(
                                       "file_selection",
                                       "勾選已上傳檔案進行操作：",
                                       choices = NULL # 由 server 動態更新
                                     ),
                                     actionButton("prepare_data", "準備資料"),
                                     width=2
                                   ),
                                   mainPanel(
                                     verbatimTextOutput("selected_files"), # 顯示選取的檔案名稱
                                     tabsetPanel(
                                       tabPanel("Sample 1 Expr", DT::dataTableOutput('ExprTableS1',width="100%")),
                                       tabPanel("Sample 2 Expr", DT::dataTableOutput('ExprTableS2',width="100%"))
                                     )
                                   )
                          ),
                          tabPanel("Sample Information",
                                   
                                   sidebarPanel(
                                     textInput("mainCode_input", "Main Code:", value = "BC026"),
                                     textInput("subCode_input", "Sub Code:", value = "PBMC"),
                                     selectInput("normMethod_input", "Normalization Method:",
                                                 choices = c("max_Cov", "max_FPKM", "max_TPM"),
                                                 selected = "max_TPM"),
                                     actionButton("preview_3_btn", "Query & Preview")
                                   ),
                                   mainPanel(
                                     h3("樣本資訊"),
                                     DTOutput("file_table_2"),
                                     downloadButton("download_table", "下載註記資料"),
                                     verbatimTextOutput("error_message") # 顯示錯誤訊息
                                   ),
                                   
                          ),
                          # 子頁 1: Wide Table --------------------------------
                          tabPanel("Gene Expression Profile",
                                   fluidRow(
                                     column(
                                       width = 12,
                                       br(),
                                       actionButton("run_wide_conversion", "Convert to Wide Table"),
                                       br(), br(),
                                       DT::dataTableOutput("wide_table_dt", width = "100%")
                                     )
                                   )
                          ),
                          
                          tabPanel("Differential Expression Analysis",
                                   sidebarLayout(
                                     sidebarPanel(
                                       sliderInput("lfc_cut", "Fold Change 閾值 (log2):", 
                                                   min = 0, max = 3, value = 1, step = 0.1),
                                       sliderInput("pval_cut", "p-value 閾值:", 
                                                   min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                                       sliderInput("pointSize", "點大小:", 
                                                   min = 1, max = 5, value = 2, step = 0.5),
                                       sliderInput("ptAlpha", "點透明度:", 
                                                   min = 0.1, max = 1, value = 0.6, step = 0.1),
                                       sliderInput("labelSize", "基因標籤字型大小:", 
                                                   min = 1, max = 6, value = 3, step = 0.5),
                                       numericInput("topN", "標記前 N 個顯著基因 (0 表示不標記):", 
                                                    value = 0, min = 0, max = 100),
                                       checkboxInput("use_adjP", "Use Adjusted P-value?", value = FALSE),
                                       actionButton("run_limma", "Run limma analysis")
                                     ),
                                     
                                     mainPanel(
                                       tabsetPanel(
                                         tabPanel("Volcano Plot interaction",
                                                  plotOutput("volcano_plot", height = "600px")
                                         ),
                                         tabPanel("DEG Table", 
                                                  DT::dataTableOutput('limma_table',width="100%")
                                         ),
                                         tabPanel("Up & Dwown Summary", 
                                                  DT::dataTableOutput('limma_summary',width="100%")
                                         )
                                       )
                                     )
                                   )
                                   
                                   
                                   
                                   
                          ),
                          tabPanel("Target Gene Expression",
                                   sidebarPanel(
                                     textInput("geneList", "Target Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2,AKT1,PIK3CA,ERBB3,CCND1,SF3B1,FGFR1,FBXW7,TP53,BRCA1,BRCA2"),
                                     actionButton(inputId = "targetGeneID", label = "Confirm"),width=2
                                   ),
                                   mainPanel(
                                     tabsetPanel(
                                       tabPanel("Target Gene Expr. Table", DT::dataTableOutput('target_gene_table', width="100%", height = "600px")),
                                     )
                                   )
                          ),
                          tabPanel("Expression Pattern", 
                                   sidebarPanel(
                                     actionButton(inputId = "generate_scatter", label = "Expression Analysis"),
                                     width=2
                                   ),
                                   mainPanel(
                                     width = 10, 
                                     tabsetPanel(
                                       tabPanel("Comparasion of Expr. Pattern", plotOutput("scatter_plot", width = "1000px", height = "800px"))
                                     )
                                   )
                          )
                        )
               ),
               
               tabPanel("Gene-enriched Analysis", 
                        sidebarPanel(
                          sliderInput("top_gene_number", "Gene number", 
                                      min = 3, max = 20000, value = 500, step = 1),
                          actionButton(inputId = "generate_genelist", label = "Generate Gene List"),
                          actionButton(inputId = "generate_go", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            id = "mainTabs",
                            
                            # Sub1 tab
                            tabPanel("Enrichment",
                                     fluidRow(
                                       column(12,
                                              # Upper section tabs for Sub1
                                              h4("Tissue"),
                                              tabsetPanel(
                                                id = "sub1Upper",
                                                tabPanel("MF", plotOutput("tissueMF")),
                                                tabPanel("BP", plotOutput("tissueBP")),
                                                tabPanel("CC", plotOutput("tissueCC")),
                                                tabPanel("KEGG", plotOutput("tissueKEGG"))
                                              )
                                       )
                                     ),
                                     fluidRow(
                                       column(12,
                                              # Lower section tabs for Sub1
                                              h4("CTC"),
                                              tabsetPanel(
                                                id = "sub1Lower",
                                                tabPanel("MF", plotOutput("ctcMF")),
                                                tabPanel("BP", plotOutput("ctcBP")),
                                                tabPanel("CC", plotOutput("ctcCC")),
                                                tabPanel("KEGG", plotOutput("ctcKEGG"))
                                                
                                              )
                                       )
                                     ),
                                     fluidRow(
                                       column(12,
                                              # Lower section tabs for Sub1
                                              h4("PBMC"),
                                              tabsetPanel(
                                                id = "sub2Lower",
                                                tabPanel("MF", plotOutput("ctcMF")),
                                                tabPanel("BP", plotOutput("ctcBP")),
                                                tabPanel("CC", plotOutput("ctcCC")),
                                                tabPanel("KEGG", plotOutput("ctcKEGG"))
                                                
                                              )
                                       )
                                     )
                            )
                            
                            
                          )
                        )
               ),
               tabPanel("TissueTyper (GTEx)", 
                        sidebarPanel(
                          actionButton(inputId = "generate_tissuetyper", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Normal Tissue Analysis", plotOutput("gtexCorPlot", width = "800px", height = "800px")),
                            tabPanel("Bulk sample cor. rank", DT::dataTableOutput('rankTable1',width="100%")),
                            tabPanel("CTC sample cor. rank", DT::dataTableOutput('rankTable2',width="100%"))
                          )
                        )
               ),
               tabPanel("Immune Fraction", 
                        sidebarPanel(
                          actionButton(inputId = "generate_immune", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Immune Cell Type Prediction", plotOutput("immune_plot", width = "1000px", height = "800px"))
                          )
                        )
               )
               
    )
  )
  
  server <- function(input, output, session) {
    #設定所需檔案之路徑
    #disease_list_path <- system.file("extdata", "0_疾病清單.xlsx", package = "alprogar")
    #rsID_and_Gene_name_path <- file.path(.libPaths()[1],'alprogar', 'extdata',"/")                
    #download_path <- file.path(path.expand('~') , "alprogar")
    #--------------------------------
    #以下用dbBrowserServer取代
    selected_db <- ""
    sc <- spark_connect(
      master = "sc://172.18.0.1:15002", #localhost 改成172.18.0.1
      method = "spark_connect",
      version = "3.5"
    )
    #--------------------------------
    
    db_info <- dbBrowserServer("dbBrowser1", sc)
    required_pkgs <- c("edgeR", "limma", "MultiAssayExperiment", "SummarizedExperiment","pcaPP", "reshape2", "stringr", "readxl","ggplot2","tidyr","tibble", "viridis","RColorBrewer","pheatmap","ggpubr","ggrepel","readr","dplyr","sparklyr","org.Hs.eg.db","enrichplot","clusterProfiler")
    invisible(lapply(required_pkgs, library, character.only = TRUE))
    gtex_expr_path <- system.file("extdata", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", package = "oncoexpr")
    
    results <- reactiveValues(
      geneList = NULL,
      date = NULL,
      groups_list = NULL,
      group1_ = NULL,
      group2_ = NULL,
      expr_data = NULL,
      scatter_plot = NULL,
      immune_plot = NULL,
      go_dotplot = list(),
      kegg_dotplot = NULL,
      tissue_heatmap = NULL,
      special_gene_list = NULL,
      tissue_fc_gene_profile = NULL,
      ctc_fc_gene_profile = NULL
    )
    
    initial_colData <- data.frame(
      mainCode = character(),
      subCode  = character(),
      #fileName = character(),
      stringsAsFactors = FALSE
    )
    set_colData <- reactiveVal(initial_colData)
    edit_mode <- reactiveVal(FALSE)
    # 編輯模式開關
    observeEvent(input$enable_edit, {
      edit_mode(TRUE)
    })
    
    # 儲存錯誤訊息
    error_message <- reactiveVal(NULL)
    
    # 動態更新檔案選項
    observeEvent(input$files, {
      req(input$files)
      updated_data <- data.frame(
        MainCode = "",
        SubCode = "",
        FileName = input$files$name,
        #FileSize = input$files$size,
        FilePath = input$files$datapath,
        stringsAsFactors = FALSE
      )
      #sample_info_table <- generate_sample_info(dat_long)
      updated_data <- rbind(set_colData(), updated_data)
      set_colData(updated_data)
      
      updateCheckboxGroupInput(
        session,
        "file_selection",
        choices = updated_data$FileName,
        selected = NULL
      )
    })
    # 顯示檔案資訊表格
    output$file_table <- renderDT({
      datatable(
        set_colData(),
        editable = edit_mode(), # 編輯模式根據按鈕控制
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    }, server = FALSE)
    output$file_table_2 <- renderDT({
      datatable(
        set_colData(),
        editable = edit_mode(), # 編輯模式根據按鈕控制
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    }, server = FALSE)
    # 監聽表格編輯事件
    observeEvent(input$file_table_cell_edit, {
      req(edit_mode()) # 確保只有在編輯模式啟用時可編輯
      info <- input$file_table_cell_edit
      df <- set_colData()
      df[info$row, info$col] <- info$value
      set_colData(df)
    })
    
    # 顯示選擇的檔案
    output$selected_files <- renderText({
      if (is.null(input$file_selection)) {
        "尚未選擇檔案"
      } else {
        paste("選擇的檔案：", paste(input$file_selection, collapse = ", "))
      }
    })
    # 使用者選擇的檔案進行處理
    observeEvent(input$prepare_data, {
      req(input$file_selection)
      geneList <- unlist(strsplit(input$geneList, ","))
      # 提取選擇的檔案的完整路徑
      selected_files <- set_colData()[set_colData()$FileName %in% input$file_selection, ]
      selected_paths <- selected_files$FilePath
      print(selected_paths)
      if (length(selected_paths) < 2) {
        showModal(modalDialog(
          title = "錯誤",
          "請至少選擇兩個檔案進行分析。",
          easyClose = TRUE,
          footer = modalButton("關閉")
        ))
        return()
      }
      
      # 動態生成 groups_list 和 entirepath
      groups_list <- c(input$groupNameS1, input$groupNameS2)
      entirepath <- selected_paths  # 完整路徑直接使用
      # 執行數據處理函數
      input_and_prep_file_directly(entirepath = entirepath, groups_list = groups_list)
      
      print(groups_list[1])
      print(groups_list[2])
      print(entirepath)
      # 獲取處理結果
      group1_ <- get(paste0(groups_list[1], "_rna_expr"))
      group2_ <- get(paste0(groups_list[2], "_rna_expr"))
      # 顯示處理結果
      
      output$ExprTableS1 <- DT::renderDataTable({ group1_ })
      output$ExprTableS2 <- DT::renderDataTable({ group2_ })
      
      # 分析數據
      Prep4ExprComp_result <- Prep4ExprComp(group1 = group1_, group2 = group2_, groups_name = groups_list)
      print(3)
      data_v2 <- Prep4ExprComp_result$result
      results$special_gene_list <- Prep4ExprComp_result$special_gene_list
      
      # 篩選顯著基因表現譜
      results$tissue_fc_gene_profile <- data_v2[((data_v2[[groups_list[1]]] + 1) / (data_v2[[groups_list[2]]] + 1)) > 1,]
      results$ctc_fc_gene_profile <- data_v2[((data_v2[[groups_list[2]]] + 1) / (data_v2[[groups_list[1]]] + 1)) > 1,]
      results$group1_ <- group1_
      results$group2_ <- group2_
      # 保存資料
      results$expr_data <- data_v2
      results$geneList <- geneList
      results$date_ <- date 
      results$groups_list <- groups_list
      results$date_ <- Sys.Date()
    })
    
    
    
    
    observeEvent(input$targetGeneID, {
      req(settingMAE(), wide_data())
      #saveFileName <- paste0(input$prefixSample, "_target_gene_expr_", results$date_ ,".csv")
      mae <- settingMAE()
      geneList <- unlist(strsplit(input$geneList, ","))
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- rownames(sample_info) 
      print(groups_list)
      print(head(wide_data()))
      expr_profile <- as.data.frame(wide_data())
      rownames(expr_profile) <- expr_profile[,"Gene"]
      expr_profile <- expr_profile[,-1]
      targetGeneExpr <- target_exprofile( 
        geneList_ = geneList, 
        #groups_list_ = results$groups_list, 
        groups_list_ = groups_list,
        #save_path_ = input$save_path, 
        #save_filename_ = saveFileName, 
        #expr_profile_ = results$expr_data
        expr_profile_ = expr_profile
      )
      print(targetGeneExpr)
      output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
    })
    observeEvent(input$generate_scatter, {
      req(results$expr_data, results$special_gene_list, results$groups_list)
      
      # 調用 expr_pattern 函數，直接返回 ggplot 對象
      output$scatter_plot <- renderPlot({
        expr_pattern(
          data_v2 = results$expr_data, 
          special_gene_list = results$special_gene_list, 
          save_path_ = input$save_path, 
          x_col = results$groups_list[1], y_col = results$groups_list[2],
          save_filename_ = paste0(input$prefixSample, "_scatter_plot.png")
        )
      })
    })
    
    
    observeEvent(input$generate_tissuetyper, {
      req(results$group1_, results$group2_, results$groups_list) 
      output$gtexCorPlot <-  renderPlot({
        tissueTyperResult <- TissueTyper( 
          group1_ = results$group1_,
          group2_ = results$group2_,
          save_path = input$save_path, 
          gtexFilePath = gtex_expr_path, 
          groups_list = results$groups_list, 
          sample_id = input$prefixSample 
        )
        tissueTyperResult$plot
        output$rankTable1 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank1)})
        output$rankTable2 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank2)})
        
      })
    })
    
    # 生成免疫分量分析圖
    observeEvent(input$generate_immune, {
      req(settingMAE()) 
      
      mae <- settingMAE()
      print(mae)
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- unique(sample_info[,"subCode"]) 
      output$immune_plot <-  renderPlot({
        ImmuneFractionPlot_mae( mae = mae,
                                groups_list = groups_list                
        )
        
      })
    })
    
    
    observeEvent(input$generate_go, {
      req(topGeneList(),downGeneList(),settingMAE())                
      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- unique(sample_info[,"subCode"]) 
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()
      #prefixSample <- input$prefixSample
      for(n in seq_len(length(groups_list)) ){
        col <- groups_list[n]
        print(col)
        print(n)
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        #save_filename <- paste0(paste(prefixSample, col, sep="_"), "_go_enrich_plot_", results$date_ ,".png")
        print("setting ok")
        #print(save_filename)
        for( mode in c("CC", "BP", "MF")){
          print(mode)
          VAR <- paste0(col, mode, "GO")
          print(VAR)
          assign(VAR, go_enrich_dotplot( 
            
            gene_list_= gene_list, 
            save_path_ = NULL,
            save_filename_ = NULL, 
            mode_ = mode, 
            showCategory_ = 10), envir=.GlobalEnv
            
          )
          
        }
        
      }
      
      output$tissueMF <- renderPlot({BULKMFGO})
      output$tissueBP <- renderPlot({BULKBPGO})
      output$tissueCC <- renderPlot({BULKCCGO})
      output$ctcMF <- renderPlot({CTCMFGO})
      output$ctcBP <- renderPlot({CTCBPGO})
      output$ctcCC <- renderPlot({CTCCCGO})
      output$pbmcMF <- renderPlot({PBMCMFGO})
      output$pbmcBP <- renderPlot({PBMCBPGO})
      output$pbmcCC <- renderPlot({PBMCCCGO})
      
      
    })
    
    preview_data_entire <- reactiveVal(NULL)
    
    
    maeColData <- reactiveVal(NULL)
    
    db_info <- dbBrowserServer("dbBrowser1", sc)                
    output$table_preview_entire <- DT::renderDataTable({ 
      preview_data_entire() 
    })
    
    mainCodes <- reactive({
      codes <- strsplit(input$mainCode_input, ",")[[1]]
      codes <- trimws(codes)
      codes
    })
    
    subCodes <- reactive({
      codes <- strsplit(input$subCode_input, ",")[[1]]
      codes <- trimws(codes)
      codes
    })
    
    colDataVal <- reactiveVal(data.frame())
    set_colData <- function(newData) {
      colDataVal(newData)
    }
    sample_mod_return <- sampleSelectionServer(
      id            = "sampleSelector",
      sc            = sc,
      selected_db   = db_info$selected_db,
      selected_table= db_info$selected_table,
      set_colData   = set_colData
    )
    
    output$file_table <- renderDT({
      req(colDataVal())
      datatable(colDataVal())
    })
    #sample-selected table
    output$table_preview_filtered  <- DT::renderDT({
      datatable(sample_mod_return$filteredTable())
    })
    
    # 用來存放 wide table 的 reactiveVal
    wide_data <- reactiveVal(NULL)
    
    # 監聽「Convert to Wide Table」按鈕
    observeEvent(input$run_wide_conversion, {
      #req(preview_data_filtered())  
      req(sample_mod_return$filteredTable())
      #dat_long <- preview_data_filtered()
      dat_long <- sample_mod_return$filteredTable()
      
      
      dat_wide <- convert_long_to_wide(dat_long)
      
      wide_data(dat_wide)
    })
    
    output$wide_table_dt <- DT::renderDataTable({
      req(wide_data())
      DT::datatable(
        wide_data(),
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    })
    
    
    
    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    limma_table <- reactiveVal(NULL)
    limma_summary <- reactiveVal(NULL)
    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL) 
    
    observeEvent(input$mae_start,{
      
      #dat_long <- preview_data_filtered()
      dat_long <- sample_mod_return$filteredTable()
      dat_wide <- convert_long_to_wide(dat_long)
      
      sample_info_table <- generate_sample_info(dat_long)
      
      assay_data <- as.matrix(dat_wide[, -1])  
      rownames(assay_data) <- dat_wide$Gene
      colnames(assay_data) <- colnames(dat_wide)[-1]
      rownames(sample_info_table) <- paste0(sample_info_table$mainCode, "_", sample_info_table$subCode)
      
      maeColData(sample_info_table)
      
      
      se_expression_matrix <- SummarizedExperiment(
        assays = list(max_TPM = assay_data), #read count, TPM, COV, FPKM
        colData = sample_info_table
      )
      
      
      mae <- MultiAssayExperiment(
        experiments = list(RNAseq = se_expression_matrix),  #expr_matrix, DEGresult, GSEAresult
        colData = sample_info_table
      )
      settingMAE(mae)
      print(settingMAE())
      print("OK")
      
    })
    
    
    observeEvent(input$run_limma, {
      #req(preview_data_filtered(), settingMAE())
      req(sample_mod_return$filteredTable(), settingMAE())
      mae <- settingMAE()
      
      limma_result <- LimmaMAE(
        mae        = mae,
        assayName  = "RNAseq",  
        subCodeCol = "subCode", 
        coef       = 2,  
        pval_cut   = input$pval_cut,
        lfc_cut    = input$lfc_cut,
        useAdjP    = input$use_adjP
      )
      limma_table(limma_result$topTable)
      limma_summary(limma_result$upDownSummary)
      output$limma_table <- renderDT({
        datatable(
          limma_table(),
          #editable = edit_mode(), # 編輯模式根據按鈕控制
          options = list(pageLength = 10, autoWidth = TRUE)
        )
      }, server = FALSE)
      
      
      output$limma_summary <- renderDT({
        datatable(
          limma_summary(),
          #editable = edit_mode(), # 編輯模式根據按鈕控制
          options = list(pageLength = 10, autoWidth = TRUE)
        )
      }, server = FALSE)
      p <- ggvolcano_limma(
        fit       = limma_result$fit,   # 實際應替換為您的 limma_fit 物件
        coef      = 2,
        lfc_cut   = input$lfc_cut,
        pval_cut  = input$pval_cut,
        useAdjP   = input$use_adjP,
        title     = "Volcano Plot (limma)",
        topN      = input$topN,
        pointSize = input$pointSize,
        ptAlpha   = input$ptAlpha,
        labelSize = input$labelSize
      )
      # 把整個結果放到 volcano_res() 方便後續需要
      volcano_res(list(
        limma_result = limma_result
      ))
    })
    
    reactive_volcano_plot <- reactive({
      req(volcano_res())
      ggvolcano_limma(
        fit       = volcano_res()$limma_result$fit,
        coef      = 2,
        lfc_cut   = input$lfc_cut,
        pval_cut  = input$pval_cut,
        useAdjP   = input$use_adjP,
        title     = "Volcano Plot (limma)",
        topN      = input$topN,
        pointSize = input$pointSize,
        ptAlpha   = input$ptAlpha,
        labelSize = input$labelSize
      )
    })
    
    output$volcano_plot <- renderPlot({
      reactive_volcano_plot()
    })
    
    observeEvent(input$generate_genelist, {
      req(limma_table(),input$top_gene_number)
      DEG_table <- limma_table() 
      sorted_DEG <- DEG_table[order(DEG_table[,"logFC"], decreasing=TRUE),]
      topGeneList(sorted_DEG[1:input$top_gene_number,])
      downGeneList(sorted_DEG[(nrow(sorted_DEG)-input$top_gene_number):nrow(sorted_DEG),])
      print(nrow(topGeneList()))
      print(nrow(downGeneList()))
      print(input$top_gene_number)
      
    })
    
  }
  
  
  
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run, host = host_, port = port_)
  
}






#' @title shiny app for RNAseq  demo version
#' @description start the shiny app with spark connection for shiny proxy at the seqslab console
#' @return start the UI and Server for analysis
#' @export
demo_oncoExprAppSpark <- function(){
  #options(shiny.maxRequestSize = 10 * 1024^5) # 10 MB
  required_pkgs <- c( "edgeR", 
                      "limma", 
                      "MultiAssayExperiment", 
                      "SummarizedExperiment", 
                      "pcaPP", 
                      "reshape2", 
                      "stringr", 
                      "readxl", 
                      "ggplot2", 
                      "tidyr", 
                      "tibble", 
                      "viridis", 
                      "RColorBrewer", 
                      "pheatmap", 
                      "ggpubr", 
                      "ggrepel", 
                      "readr", 
                      "dplyr", 
                      "sparklyr", 
                      "org.Hs.eg.db", 
                      "enrichplot", 
                      "clusterProfiler",
                      "DT",
                      "DBI",
                      "shinybusy",
                      "shiny")
  invisible(lapply(required_pkgs, library, character.only = TRUE))
  ui <- fluidPage(
    navbarPage("RNAseq App (Beta)",
               tabPanel("Spark DB browser",
                        tags$style(".shinybusy-overlay {opacity: 0.7; background-color: #7c7c7c;}"),
                        add_busy_spinner(
                          spin = "fading-circle",
                          position = "full-page",
                          timeout = 1000
                        ),
                        sidebarLayout(
                          sidebarPanel(
                            # 第一塊：Spark DB 資料庫操作
                            wellPanel(
                              h4("Spark DB Data"),
                              sparkConnectionUI("spark_mod"),
                              dbBrowserUI("dbBrowser1"),
                              #sampleSelectionUI("sampleSelector"),
                              mod_filter_sidebar_ui("filter1"),
                              filterModuleUI("filter"),
                              actionButton("test_btn", "Get"),
                              actionButton("mae_start", "Send")
                            )
                            
                            # 第二塊：原本 "Prepare Samples" 的功能，放在同一 sidebarPanel 裡
                            # wellPanel(
                            #     h4("Local Upload / Prepare Samples"),
                            #     fileInput(
                            #         "files",
                            #         "選擇多個檔案上傳",
                            #         multiple = TRUE,
                            #         accept = c(".csv", ".txt", ".gtf")
                            #     ),
                            #     actionButton("enable_edit", "啟用編輯"),
                            #     actionButton("add_row", "新增空白列"),
                            #     actionButton("add_column", "新增空白欄位"),
                            #     actionButton("rename_column", "修改欄位名稱"),
                            #     actionButton("prepare_data", "準備資料")
                            # )
                          ),
                          
                          
                          mainPanel(
                            tabsetPanel(
                              tabPanel("SparkConnection",
                                       verbatimTextOutput("conn_status")
                              ),
                              # tabPanel("Sample Information (colData)",
                              #     # 其他想顯示的元件放這裡...
                              #     h3("樣本資訊"),
                              #     DTOutput("file_table"),
                              #     downloadButton("download_table", "下載註記資料"),
                              #     verbatimTextOutput("error_message") # 顯示錯誤訊息
                              # ),
                              # tabPanel("Entire Table",
                              #     DT::dataTableOutput("table_preview_entire", width = "100%")
                              # ),
                              # tabPanel("Filtered Table (longData)",
                              #     #DTOutput
                              #     DT::dataTableOutput("table_preview_filtered", width = "100%")
                              # ),
                              tabPanel("Sample Information (colData)", DTOutput(NS("filter1", "colData_table"))),
                              tabPanel("Filtered Data (pivot long)", DTOutput(NS("filter1", "filtered_data_table"))),
                              # tabPanel("A",
                              
                              #     DT::dataTableOutput("filtered_table")
                              # ),
                              tabPanel("Group Check (subCode)",DT::dataTableOutput("filtered_colData_table"))
                              
                              
                            ),
                            
                          )
                          
                        )
               ),
               tabPanel("Gene Expression (RNA-seq)",
                        tabsetPanel(
                          # tabPanel("Selecting uploaded samples",
                          #     sidebarPanel(
                          #         textInput("prefixSample", "Sample Code", value = "BC021"),
                          #         hr(),
                          #         textInput(inputId ="groupNameS1", "Sample 1 ", value = "tissue"),
                          #         textInput(inputId ="groupNameS2", "Sample 2", value = "ctc"),
                          #         #fileInput(inputId = "gtfFileS1", "Upload Sample 1 file", buttonLabel = "Upload..."),
                          #         #fileInput(inputId = "gtfFileS2", "Upload Sample 2 file", buttonLabel = "Upload..."),
                          #         checkboxGroupInput(
                          #             "file_selection",
                          #             "勾選已上傳檔案進行操作：",
                          #             choices = NULL # 由 server 動態更新
                          #         ),
                          #         actionButton("prepare_data", "準備資料"),
                          #         width=2
                          #         ),
                          #     mainPanel(
                          #         verbatimTextOutput("selected_files"), # 顯示選取的檔案名稱
                          #         tabsetPanel(
                          #                 tabPanel("Sample 1 Expr", DT::dataTableOutput('ExprTableS1',width="100%")),
                          #                 tabPanel("Sample 2 Expr", DT::dataTableOutput('ExprTableS2',width="100%"))
                          #             )
                          #         )
                          # ),
                          # tabPanel("Sample Information",
                          #     tabPanel("Selecting uploaded samples",
                          #         sidebarPanel(
                          #             sampleSelectionUI("sampleSelector_2")
                          #         )
                          #     )
                          
                          # ),
                          # 子頁 1: Wide Table --------------------------------
                          tabPanel("Gene Expression Profile",
                                   fluidRow(
                                     column(
                                       width = 12,
                                       br(),
                                       actionButton("run_wide_conversion", "Convert to Wide Table"),
                                       br(), br(),
                                       DT::dataTableOutput("wide_table_dt", width = "100%")
                                     )
                                   )
                          ),
                          
                          tabPanel("Differential Expression Analysis",
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
                                                    value = 10, min = 0, max = 1000),
                                       checkboxInput("use_adjP", "Use Adjusted P-value?", value = FALSE),
                                       actionButton("run_limma", "Run limma analysis")
                                     ),
                                     
                                     mainPanel(
                                       tabsetPanel(
                                         tabPanel("Volcano Plot interaction",
                                                  plotOutput("volcano_plot", height = "600px")
                                         ),
                                         tabPanel("DEG Table", 
                                                  DT::dataTableOutput('limma_table',width="100%")
                                         ),
                                         tabPanel("Up & Down Summary", 
                                                  DT::dataTableOutput('limma_summary',width="100%")
                                         )
                                       )
                                     )
                                   )
                                   
                                   
                                   
                                   
                          ),
                          tabPanel("Target Gene Expression",
                                   sidebarPanel(
                                     textInput("geneList", "Target Gene List (sep by comma without space)", value = "EGFR,ESR1,KRAS,ERBB2,AKT1,PIK3CA,ERBB3,CCND1,SF3B1,FGFR1,FBXW7,TP53,BRCA1,BRCA2"),
                                     actionButton(inputId = "targetGeneID", label = "Confirm"),width=2
                                   ),
                                   mainPanel(
                                     tabsetPanel(
                                       tabPanel("Target Gene Expr. Table", DT::dataTableOutput('target_gene_table', width="100%", height = "600px")),
                                     )
                                   )
                          )#,
                          # tabPanel("Expression Pattern", 
                          #     sidebarPanel(
                          #         actionButton(inputId = "generate_scatter", label = "Expression Analysis"),
                          #         width=2
                          #         ),
                          #     mainPanel(
                          #         width = 10, 
                          #         tabsetPanel(
                          #             tabPanel("Comparasion of Expr. Pattern", plotOutput("scatter_plot", width = "1000px", height = "800px"))
                          #             )
                          #         )
                          # )
                        )
               ),
               
               tabPanel("Gene-enriched Analysis", 
                        sidebarPanel(
                          sliderInput("top_gene_number", "Gene number", 
                                      min = 3, max = 20000, value = 500, step = 1),
                          actionButton(inputId = "generate_genelist", label = "Generate Gene List"),
                          actionButton(inputId = "generate_go", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            id = "mainTabs",
                            
                            # Sub1 tab
                            tabPanel("Enrichment",
                                     fluidRow(
                                       column(12,
                                              # Upper section tabs for Sub1
                                              h4("Tumor"),
                                              tabsetPanel(
                                                id = "sub1Upper",
                                                tabPanel("MF", plotOutput("tissueMF")),
                                                tabPanel("BP", plotOutput("tissueBP")),
                                                tabPanel("CC", plotOutput("tissueCC")),
                                                tabPanel("KEGG", plotOutput("tissueKEGG"))
                                              )
                                       )
                                     ),
                                     fluidRow(
                                       column(12,
                                              # Lower section tabs for Sub1
                                              h4("Normal"),
                                              tabsetPanel(
                                                id = "sub1Lower",
                                                tabPanel("MF", plotOutput("ctcMF")),
                                                tabPanel("BP", plotOutput("ctcBP")),
                                                tabPanel("CC", plotOutput("ctcCC")),
                                                tabPanel("KEGG", plotOutput("ctcKEGG"))
                                                
                                              )
                                       )
                                     )#,
                                     # fluidRow(
                                     #     column(12,
                                     #         # Lower section tabs for Sub1
                                     #         h4("PBMC"),
                                     #         tabsetPanel(
                                     #             id = "sub2Lower",
                                     #             tabPanel("MF", plotOutput("pbmcMF")),
                                     #             tabPanel("BP", plotOutput("pbmcBP")),
                                     #             tabPanel("CC", plotOutput("pbmcCC")),
                                     #             tabPanel("KEGG", plotOutput("pbmcKEGG"))
                                     
                                     #         )
                                     #     )
                                     # )
                            )
                          )
                        )
               ),
               # tabPanel("TissueTyper (GTEx)", 
               #     sidebarPanel(
               #         actionButton(inputId = "generate_tissuetyper", label = "Analysis"),
               #         width=2
               #         ),
               #     mainPanel(
               #         tabsetPanel(
               #             tabPanel("Normal Tissue Analysis", plotOutput("gtexCorPlot", width = "800px", height = "800px")),
               #             tabPanel("Bulk sample cor. rank", DT::dataTableOutput('rankTable1',width="100%")),
               #             tabPanel("CTC sample cor. rank", DT::dataTableOutput('rankTable2',width="100%"))
               #             )
               #         )
               # ),
               tabPanel("Immune Fraction", 
                        sidebarPanel(
                          actionButton(inputId = "generate_immune", label = "Analysis"),
                          width=2
                        ),
                        mainPanel(
                          tabsetPanel(
                            tabPanel("Immune Cell Type Prediction", plotOutput("immune_plot", width = "1000px", height = "800px"))
                          )
                        )
               )
               
    )
  )
  
  server <- function(input, output, session) {
    
    selected_db <- ""
    demo_tcga_brca_path <- system.file("extdata", "demo_tcga_brca_example.csv", package = "oncoexpr")
    test_se2lf <- read.csv(demo_tcga_brca_path, sep = ",", header=TRUE) 
    sc <- sparkConnectionServer("spark_mod")
    # 示範使用 sc 連線物件，例如顯示連線狀態
    output$conn_status <- renderPrint({
      if (is.null(sc())) {
        "尚未建立 Spark 連線"
      } else {
        # 可以顯示連線資訊，例如：
        cat("Spark 連線資訊：\n")
        print(sc())
      }
    })
    
    #sc <- spark_connect(master = "local")
    # sc <- spark_connect(
    #     master = "sc://localhost:15002",
    #     method = "spark_connect",
    #     version = "3.5"
    # )
    print("TCGA")
    
    observe({   print(sc()$"master")
      req(sc()$"master"=="local[8]", test_se2lf)  # 確保 Spark 連線存在
      spark_df <<- copy_to(sc(), as.data.frame(test_se2lf), "tcga_brca_long", overwrite = TRUE)
      #spark_colData <<- copy_to(sc(), as.data.frame(brca_colData), "tcga_brca_colData", overwrite = TRUE)
      
    })
    
    #db_info <- reactiveValues(selected_table = NULL)
    results <- reactiveValues(
      geneList = NULL,
      date = NULL,
      groups_list = NULL,
      group1_ = NULL,
      group2_ = NULL,
      expr_data = NULL,
      scatter_plot = NULL,
      immune_plot = NULL,
      go_dotplot = list(),
      kegg_dotplot = NULL,
      tissue_heatmap = NULL,
      special_gene_list = NULL,
      tissue_fc_gene_profile = NULL,
      ctc_fc_gene_profile = NULL,
      db_info =NULL
    )
    observe({
      req(sc()) 
      results$db_info <- dbBrowserServer("dbBrowser1", sc())
      
    })
    
    gtex_expr_path <- system.file("extdata", "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", package = "oncoexpr")
    spark_connection <- reactiveValues(sc =NULL)
    
    
    initial_colData <- data.frame(
      mainCode = character(),
      subCode  = character(),
      #fileName = character(),
      stringsAsFactors = FALSE
    )
    set_colData <- reactiveVal(initial_colData)
    edit_mode <- reactiveVal(FALSE)
    # 編輯模式開關
    observeEvent(input$enable_edit, {
      edit_mode(TRUE)
    })
    
    # 儲存錯誤訊息
    error_message <- reactiveVal(NULL)
    
    # 動態更新檔案選項
    observeEvent(input$files, {
      req(input$files)
      updated_data <- data.frame(
        MainCode = "",
        SubCode = "",
        FileName = input$files$name,
        #FileSize = input$files$size,
        FilePath = input$files$datapath,
        stringsAsFactors = FALSE
      )
      #sample_info_table <- generate_sample_info(dat_long)
      updated_data <- rbind(set_colData(), updated_data)
      set_colData(updated_data)
      
      updateCheckboxGroupInput(
        session,
        "file_selection",
        choices = updated_data$FileName,
        selected = NULL
      )
    })
    # 顯示檔案資訊表格
    output$file_table <- renderDT({
      datatable(
        set_colData(),
        editable = edit_mode(), # 編輯模式根據按鈕控制
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    }, server = FALSE)
    output$file_table_2 <- renderDT({
      datatable(
        set_colData(),
        editable = edit_mode(), # 編輯模式根據按鈕控制
        options = list(pageLength = 10, autoWidth = TRUE)
      )
    }, server = FALSE)
    # 監聽表格編輯事件
    observeEvent(input$file_table_cell_edit, {
      req(edit_mode()) # 確保只有在編輯模式啟用時可編輯
      info <- input$file_table_cell_edit
      df <- set_colData()
      df[info$row, info$col] <- info$value
      set_colData(df)
    })
    
    # 顯示選擇的檔案
    output$selected_files <- renderText({
      if (is.null(input$file_selection)) {
        "尚未選擇檔案"
      } else {
        paste("選擇的檔案：", paste(input$file_selection, collapse = ", "))
      }
    })
    # 使用者選擇的檔案進行處理
    observeEvent(input$prepare_data, {
      req(input$file_selection)
      geneList <- unlist(strsplit(input$geneList, ","))
      # 提取選擇的檔案的完整路徑
      selected_files <- set_colData()[set_colData()$FileName %in% input$file_selection, ]
      selected_paths <- selected_files$FilePath
      print(selected_paths)
      if (length(selected_paths) < 2) {
        showModal(modalDialog(
          title = "錯誤",
          "請至少選擇兩個檔案進行分析。",
          easyClose = TRUE,
          footer = modalButton("關閉")
        ))
        return()
      }
      
      # 動態生成 groups_list 和 entirepath
      groups_list <- c(input$groupNameS1, input$groupNameS2)
      entirepath <- selected_paths  # 完整路徑直接使用
      # 執行數據處理函數
      input_and_prep_file_directly(entirepath = entirepath, groups_list = groups_list)
      
      
      print(groups_list[1])
      print(groups_list[2])
      print(entirepath)
      # 獲取處理結果
      group1_ <- get(paste0(groups_list[1], "_rna_expr"))
      group2_ <- get(paste0(groups_list[2], "_rna_expr"))
      print(head(group2_ ))
      # 顯示處理結果
      
      output$ExprTableS1 <- DT::renderDataTable({ group1_ })
      output$ExprTableS2 <- DT::renderDataTable({ group2_ })
      
      # 分析數據
      Prep4ExprComp_result <- Prep4ExprComp(group1 = group1_, group2 = group2_, groups_name = groups_list)
      data_v2 <- Prep4ExprComp_result$result
      results$special_gene_list <- Prep4ExprComp_result$special_gene_list
      
      # 篩選顯著基因表現譜
      results$tissue_fc_gene_profile <- data_v2[((data_v2[[groups_list[1]]] + 1) / (data_v2[[groups_list[2]]] + 1)) > 1,]
      results$ctc_fc_gene_profile <- data_v2[((data_v2[[groups_list[2]]] + 1) / (data_v2[[groups_list[1]]] + 1)) > 1,]
      results$group1_ <- group1_
      results$group2_ <- group2_
      # 保存資料
      results$expr_data <- data_v2
      results$geneList <- geneList
      results$date_ <- date 
      results$groups_list <- groups_list
      results$date_ <- Sys.Date()
    })
    
    
    
    
    observeEvent(input$targetGeneID, {
      req(settingMAE(), wide_data())
      #saveFileName <- paste0(input$prefixSample, "_target_gene_expr_", results$date_ ,".csv")
      mae <- settingMAE()
      geneList <- unlist(strsplit(input$geneList, ","))
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- rownames(sample_info) 
      print(groups_list)
      print(head(wide_data()))
      expr_profile <- as.data.frame(wide_data())
      rownames(expr_profile) <- expr_profile[,"Gene"]
      expr_profile <- expr_profile[,-1]
      targetGeneExpr <- target_exprofile( 
        geneList_ = geneList, 
        #groups_list_ = results$groups_list, 
        groups_list_ = groups_list,
        #save_path_ = input$save_path, 
        #save_filename_ = saveFileName, 
        #expr_profile_ = results$expr_data
        expr_profile_ = expr_profile
      )
      print(targetGeneExpr)
      output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
    })
    observeEvent(input$generate_scatter, {
      req(results$expr_data, results$special_gene_list, results$groups_list)
      
      # 調用 expr_pattern 函數，直接返回 ggplot 對象
      output$scatter_plot <- renderPlot({
        expr_pattern(
          data_v2 = results$expr_data, 
          special_gene_list = results$special_gene_list, 
          save_path_ = input$save_path, 
          x_col = results$groups_list[1], y_col = results$groups_list[2],
          save_filename_ = paste0(input$prefixSample, "_scatter_plot.png")
        )
      })
    })
    
    
    observeEvent(input$generate_tissuetyper, {
      req(results$group1_, results$group2_, results$groups_list) 
      output$gtexCorPlot <-  renderPlot({
        tissueTyperResult <- TissueTyper( 
          group1_ = results$group1_,
          group2_ = results$group2_,
          save_path = input$save_path, 
          gtexFilePath = gtex_expr_path, 
          groups_list = results$groups_list, 
          sample_id = input$prefixSample 
        )
        tissueTyperResult$plot
        output$rankTable1 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank1)})
        output$rankTable2 <- DT::renderDataTable({DT::datatable(tissueTyperResult$rank2)})
        
      })
    })
    
    # 生成免疫分量分析圖
    observeEvent(input$generate_immune, {
      req(settingMAE()) 
      
      mae <- settingMAE()
      print(mae)
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- unique(sample_info[,input[["filter-subCode_col"]]]) 
      output$immune_plot <-  renderPlot({
        ImmuneFractionPlot_mae( mae = mae,
                                groups_list = groups_list                
        )
        
      })
    })
    
    
    observeEvent(input$generate_go, {
      req(topGeneList(),downGeneList(),settingMAE())                
      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- unique(sample_info[,input[["filter-subCode_col"]]]) 
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()
      #prefixSample <- input$prefixSample
      for(n in seq_len(length(groups_list)) ){
        col <- groups_list[n]
        print(col)
        print(n)
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        #save_filename <- paste0(paste(prefixSample, col, sep="_"), "_go_enrich_plot_", results$date_ ,".png")
        print("setting ok")
        #print(save_filename)
        for( mode in c("CC", "BP", "MF")){
          print(mode)
          VAR <- paste0(col, mode, "GO")
          print(VAR)
          assign(VAR, go_enrich_dotplot( 
            
            gene_list_= rownames(gene_list), 
            save_path_ = NULL,
            save_filename_ = NULL, 
            mode_ = mode, 
            showCategory_ = 10), envir=.GlobalEnv
            
          )
          
        }
        
      }
      
      # output$tissueMF <- renderPlot({BULKMFGO})
      # output$tissueBP <- renderPlot({BULKBPGO})
      # output$tissueCC <- renderPlot({BULKCCGO})
      # output$ctcMF <- renderPlot({CTCMFGO})
      # output$ctcBP <- renderPlot({CTCBPGO})
      # output$ctcCC <- renderPlot({CTCCCGO})
      # output$pbmcMF <- renderPlot({PBMCMFGO})
      # output$pbmcBP <- renderPlot({PBMCBPGO})
      # output$pbmcCC <- renderPlot({PBMCCCGO})
      
      output$tissueMF <- renderPlot({TumorMFGO})
      output$tissueBP <- renderPlot({TumorBPGO})
      output$tissueCC <- renderPlot({TumorCCGO})
      output$ctcMF <- renderPlot({NormalMFGO})
      output$ctcBP <- renderPlot({NormalBPGO})
      output$ctcCC <- renderPlot({NormalCCGO})
      
      
    })
    
    observeEvent(input$generate_go, {
      req(topGeneList(),downGeneList(),settingMAE())
      mae <- settingMAE()
      sample_info <- colData(mae[["RNAseq"]])
      groups_list <- unique(sample_info[,input[["filter-subCode_col"]]]) 
      group1_fc_gene_profile <- topGeneList()
      group2_fc_gene_profile <- downGeneList()
      #prefixSample <- input$prefixSample
      for(n in seq_len(length(groups_list)) ){
        print(col)
        col <- groups_list[n]
        gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
        #save_filename <- paste0(paste(prefixSample, col, sep="_"), "_kegg_enrich_plot_", results$date_ ,".png")
        print(col)
        VAR <- paste0(col,"KEGG")
        assign(VAR, kegg_enrich_dotplot ( 
          gene_list_ = rownames(gene_list), 
          save_path_ = NULL,
          save_filename_ = NULL, 
          showCategory_ = 10
        )
        )
      }
      
      output$tissueKEGG <- renderPlot({TumorKEGG})
      output$ctcKEGG <- renderPlot({NormalKEGG})
      
    })
    print("建立變數空值")
    preview_data_entire <- reactiveVal(NULL)
    maeColData <- reactiveVal(NULL)
    print("取得spark deltaTable")
    #db_info <- dbBrowserServer("dbBrowser1", sc)                
    # output$table_preview_entire <- DT::renderDataTable({ 
    #     preview_data_entire() 
    # })
    print("輸入mainCode欄位的特定內容，篩選資料用")
    mainCodes <- reactive({
      codes <- strsplit(input$mainCode_input, ",")[[1]]
      codes <- trimws(codes)
      codes
    })
    print("輸入subCode欄位的特定內容，篩選資料用")
    subCodes <- reactive({
      codes <- strsplit(input$subCode_input, ",")[[1]]
      codes <- trimws(codes)
      codes
    })
    print("建立colDatad空值，也就是樣本資訊")
    colDataVal <- reactiveVal(data.frame())
    set_colData <- function(newData) {
      colDataVal(newData)
    }
    
    print("針對spark deltaTable篩選樣本")
    reactive({
      req(sc())
      sample_mod_return <- sampleSelectionServer(
        id            = "sampleSelector",
        sc            = sc(),
        selected_db   = db_info$selected_db,
        selected_table= db_info$selected_table,
        set_colData   = set_colData
      )
    })
    #sample_mod_return
    #subColData = set_colData,
    #filteredTable  = preview_data_filtered
    filteredTable <- reactiveVal(NULL)
    
    #"--------------"
    observeEvent(input$test_btn,{
      req(results$db_info$selected_table)
      print("測試，新增自動選擇mainCode subCode功能")
      use_table <- results$db_info$selected_table #基於選擇的表單來取得資料
      filtered_table_test <- tbl(sc(), use_table()) #test_se2lf
      long_df_reactive <- reactive({ req(filtered_table_test); as.data.frame(filtered_table_test) })
      
      
      filtered_data_longformat <- mod_filter_server("filter1", long_df = long_df_reactive)
      
      print(head(filtered_data_longformat))
      print("準備進行main cod sub code選擇")
      filtered_data_main_and_sub <- filterModuleServer("filter", filtered_data_longformat$filtered_data)
      print(head(filtered_data_main_and_sub))
      
      print("傳到外部")
      filteredTable(as.data.frame(filtered_data_main_and_sub()))
      print(head(filteredTable()))
      
      
      output$filtered_table <- DT::renderDataTable({
        req(filtered_data_main_and_sub())
        DT::datatable(filtered_data_main_and_sub(), options = list(pageLength = 10, autoWidth = TRUE))
      })
      
      colData_reactive <- reactive({
        req(filtered_data_main_and_sub(), input[["filter-mainCode_col"]], input[["filter-subCode_col"]])
        
        df <- generate_coldata_for_local_tcga_data(
          filtered_data_main_and_sub(), 
          input[["filter-mainCode_col"]], 
          input[["filter-subCode_col"]]
        )
        
        validate(need(is.data.frame(df), "No valid data available"))
        df
      })
      
      
      output$filtered_colData_table <- DT::renderDataTable({
        req(colData_reactive())  
        df <- colData_reactive()  # 避免重複計算
        DT::datatable(df, options = list(pageLength = 10, autoWidth = TRUE))
      })
    })
    
    #"--------------"
    
    print("回傳至UI sample information")
    output$file_table <- renderDT({
      req(colDataVal())
      datatable(colDataVal())
    })
    
    #return the tabel to UI 
    # output$table_preview_filtered  <- DT::renderDT({
    #         datatable(sample_mod_return$filteredTable())
    # })
    # 用來存放 wide table 的 reactiveVal
    wide_data <- reactiveVal(NULL)
    
    # 監聽「Convert to Wide Table」按鈕
    observeEvent(input$run_wide_conversion, {
      #req(preview_data_filtered())  
      #req(sample_mod_return$filteredTable())
      req(filteredTable(), input[["filter-mainCode_col"]], input[["filter-subCode_col"]])
      #dat_long <- preview_data_filtered()
      #dat_long <- sample_mod_return$filteredTable()
      dat_long <- filteredTable()
      dat_wide <- convert_long_to_wide(dat_long, input[["filter-mainCode_col"]] ,input[["filter-subCode_col"]] )
      wide_data(dat_wide)
      print("生成wide data")
    })
    
    output$wide_table_dt <- DT::renderDataTable({
      req(wide_data())
      print("傳送wide data至UI")
      DT::datatable(
        wide_data(),
        options = list(pageLength = 20, autoWidth = TRUE)
      )
    })
    
    
    
    volcano_res <- reactiveVal(NULL)
    settingMAE <- reactiveVal(NULL)
    limma_table <- reactiveVal(NULL)
    limma_summary <- reactiveVal(NULL)
    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL) 
    # observeEvent(input$mae_start, {
    #     print("開始 MAE 建立")
    #     req(filteredTable())
    #     print("filteredTable() 存在")
    #     dat_long <- filteredTable()
    #     print(class(dat_long))  # 應該是 "data.frame" 或 "tibble"
    #     print(dim(dat_long))    # 應該顯示行數與列數
    # })
    
    observeEvent(input$mae_start,{
      req(filteredTable(),input[["filter-mainCode_col"]] ,input[["filter-subCode_col"]])
      #dat_long <- preview_data_filtered()
      dat_long <- filteredTable()
      print(head(dat_long))
      #dat_long <- sample_mod_return$filteredTable()
      dat_wide <- convert_long_to_wide(dat_long, input[["filter-mainCode_col"]] ,input[["filter-subCode_col"]] )
      #dat_wide <- convert_long_to_wide(dat_long)
      print(dat_wide)
      print("生成sample info")
      sample_info_table <- generate_sample_info(dat_long, input[["filter-mainCode_col"]] ,input[["filter-subCode_col"]])
      print(head(sample_info_table))
      assay_data <- as.matrix(dat_wide[, -1])  
      rownames(assay_data) <- dat_wide$Gene
      colnames(assay_data) <- colnames(dat_wide)[-1]
      print(rownames(sample_info_table))
      print(paste0(sample_info_table[,input[["filter-mainCode_col"]]], "_", sample_info_table[,input[["filter-subCode_col"]]]))
      print("testttt")
      rownames(sample_info_table) <- paste0(sample_info_table[,input[["filter-mainCode_col"]]], "_", sample_info_table[,input[["filter-subCode_col"]]])
      
      print("ddddd")
      maeColData(sample_info_table)
      
      
      se_expression_matrix <- SummarizedExperiment(
        assays = list(max_TPM = assay_data), #read count, TPM, COV, FPKM
        colData = sample_info_table
      )
      
      mae <- MultiAssayExperiment(
        experiments = list(RNAseq = se_expression_matrix),  #expr_matrix, DEGresult, GSEAresult
        colData = sample_info_table
      )
      settingMAE(mae)
      print(settingMAE())
      print("OK")
      
    })
    
    
    observeEvent(input$run_limma, {
      #req(preview_data_filtered(), settingMAE())
      #req(sample_mod_return$filteredTable(), settingMAE())
      req(settingMAE())
      mae <- settingMAE()
      
      limma_result <- LimmaMAE(
        mae        = mae,
        assayName  = "RNAseq",  
        #subCodeCol = "subCode",
        subCodeCol = input[["filter-subCode_col"]], 
        coef       = 2,  
        pval_cut   = input$pval_cut,
        lfc_cut    = input$lfc_cut,
        useAdjP    = input$use_adjP
      )
      limma_table(limma_result$topTable)
      limma_summary(limma_result$upDownSummary)
      output$limma_table <- renderDT({
        datatable(
          limma_table(),
          #editable = edit_mode(), # 編輯模式根據按鈕控制
          options = list(pageLength = 10, autoWidth = TRUE)
        )
      }, server = FALSE)
      
      
      output$limma_summary <- renderDT({
        datatable(
          limma_summary(),
          #editable = edit_mode(), # 編輯模式根據按鈕控制
          options = list(pageLength = 10, autoWidth = TRUE)
        )
      }, server = FALSE)
      p <- ggvolcano_limma(
        fit       = limma_result$fit,   # 實際應替換為您的 limma_fit 物件
        coef      = 2,
        lfc_cut   = input$lfc_cut,
        pval_cut  = input$pval_cut,
        useAdjP   = input$use_adjP,
        title     = "Volcano Plot (limma)",
        topN      = input$topN,
        pointSize = input$pointSize,
        ptAlpha   = input$ptAlpha,
        labelSize = input$labelSize
      )
      # 把整個結果放到 volcano_res() 方便後續需要
      volcano_res(list(
        limma_result = limma_result
      ))
    })
    
    reactive_volcano_plot <- reactive({
      req(volcano_res())
      ggvolcano_limma(
        fit       = volcano_res()$limma_result$fit,
        coef      = 2,
        lfc_cut   = input$lfc_cut,
        pval_cut  = input$pval_cut,
        useAdjP   = input$use_adjP,
        title     = "Volcano Plot (limma)",
        topN      = input$topN,
        pointSize = input$pointSize,
        ptAlpha   = input$ptAlpha,
        labelSize = input$labelSize
      )
    })
    
    output$volcano_plot <- renderPlot({
      reactive_volcano_plot()
    })
    
    observeEvent(input$generate_genelist, {
      req(limma_table(),input$top_gene_number)
      DEG_table <- limma_table() 
      sorted_DEG <- DEG_table[order(DEG_table[,"logFC"], decreasing=TRUE),]
      topGeneList(sorted_DEG[1:input$top_gene_number,])
      downGeneList(sorted_DEG[(nrow(sorted_DEG)-input$top_gene_number):nrow(sorted_DEG),])
      print(nrow(topGeneList()))
      print(nrow(downGeneList()))
      print(input$top_gene_number)
      
    })
    
  }
  
  
  
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run)
  
}



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
      title = "RNAseq App (Beta)",
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
            #actionButton("get_tbl", "Send")
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
                       plotOutput("volcano_plot", height = "600px")
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
          textInput("geneList", "Target Gene List (sep by comma without space)",
                    value = "EGFR,ESR1,KRAS,ERBB2,AKT1,PIK3CA,ERBB3,CCND1,SF3B1,FGFR1,FBXW7,TP53,BRCA1,BRCA2"),
          actionButton(inputId = "targetGeneID", label = "Confirm"),
          width = 2
        ),
        mainPanel(
          tabsetPanel(
            tabPanel("Target Gene Expr. Table", 
                     DT::dataTableOutput('target_gene_table', width = "100%", height = "600px")
            )
          )
        )
      )
    )
  )

  server <- function(input, output, session) {
    sc <- reactiveVal(NULL)
    observe({
      master <- "sc://172.18.0.1:15002"
      #master <- "sc://localhost:15002"
      method <- "spark_connect"
      version <- "3.5"
      connection <- sparklyr::spark_connect(master = master, method = method, version = version)
      sc(connection)
      print("connected")
    })
    
    results <- reactiveValues(
      db_info = NULL,
      normcount_data = NULL,
      exacttest_data = NULL,
      table_list = NULL,
      coldata = NULL
    )
    
    observe({
      req(sc())
      results$db_info <- dbBrowserServer("dbBrowser1", sc())
      print("db info 更新")
    })
    
    observeEvent(input$get_tbl, {
        req(results$db_info$selected_db())
        print("成功取得資料庫名稱")
        selected_db_name <- results$db_info$selected_db()
        print(paste0("使用的資料庫：", selected_db_name))
        
        DBI::dbExecute(sc(), paste0("USE ", selected_db_name))
        tbl_list_query <- DBI::dbGetQuery(sc(), paste0("SHOW TABLES IN ", selected_db_name))
        tbls <- tbl_list_query$tableName      
        print(tbls)
        
        # 篩選出符合條件的 table
        #prefix <- c("normcount", "exacttest", "coldata")
        tbls_with_prefix <- tbl_list_query[grepl("^normcount|^exacttest|^coldata" , tbl_list_query$"tableName"),]



        print(tbls_with_prefix)
        results$table_list <- tbls_with_prefix
        
        # 取得 normcount 資料表 (假設取第一個符合條件的)
        normcount_tbls <- tbls_with_prefix[grepl("^normcount", tbls_with_prefix$"tableName", ignore.case = TRUE),"tableName"]
        exacttest_tbls <- tbls_with_prefix[grepl("^exacttest", tbls_with_prefix$"tableName", ignore.case = TRUE),"tableName"]


        if(length(normcount_tbls) > 0){
            query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
            results$normcount_data <- DBI::dbGetQuery(sc(), query_normcount)
        }
        colnames(results$normcount_data)[colnames(results$normcount_data)=="genes"] <- "GeneSymbol"

        if(length(exacttest_tbls) > 0){
            query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
            results$exacttest_data <- DBI::dbGetQuery(sc(), query_exacttest)
        }
        colnames(results$exacttest_data)[colnames(results$exacttest_data)=="genes"] <- "GeneSymbol"
        
        colData <- generate_colData_random(results$normcount_data, genecol = "GeneSymbol")
        results$coldata <- colData
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
      sample_info_table <- maeColData()
      rownames(sample_info_table) <- colnames(assay_data) # The rownames of colData must match the colnames of assay_data
      se_expression_matrix <- SummarizedExperiment(
        assays = list(normCount = assay_data), #read count, TPM, COV, FPKM
        colData = sample_info_table
      )
      mae <- MultiAssayExperiment(
        experiments = list(RNAseq = se_expression_matrix),  #expr_matrix, DEGresult, GSEAresult
        colData = sample_info_table
      )
      settingMAE(mae)
      print(settingMAE())
      print("OK")
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
    
    reactive_volcano_plot <- eventReactive(input$run_DEG, {
      req(DEG_table())
      ggvolcano_custom(
        df = DEG_table(),
        geneName = DEG_table()$GeneSymbol,
        pValCol = "PValue",
        logFCCol = "logFC",
        coef = 2,
        lfc_cut = input$lfc_cut,
        pval_cut = input$pval_cut,
        useAdjP = FALSE,
        title = "Volcano Plot",
        topN = input$topN,
        geneCol = NULL,
        pointSize = input$pointSize, 
        ptAlpha = input$ptAlpha,
        labelSize = input$labelSize 
      )
    })
    
    output$volcano_plot <- renderPlot({
      reactive_volcano_plot()
    })
    
    topGeneList <- reactiveVal(NULL)
    downGeneList <- reactiveVal(NULL)
    
    observeEvent(input$generate_go, {
      req(DEG_table())
      DEG_table <- DEG_table()
      
      DEG_table_filtered <- DEG_table[DEG_table$PValue < input$pval_cut & abs(DEG_table$logFC) > input$lfc_cut, ]
      sorted_DEG <- DEG_table_filtered[order(DEG_table_filtered$logFC, decreasing = TRUE),]
      gene_list_symbol <- sorted_DEG$GeneSymbol
      print(gene_list_symbol)
      topGeneList(DEG_table[DEG_table$PValue < input$pval_cut & sign(DEG_table_filtered$logFC) == 1, "GeneSymbol"])
      downGeneList(DEG_table[DEG_table$PValue < input$pval_cut & sign(DEG_table_filtered$logFC) == -1, "GeneSymbol"]) 
      print(head(topGeneList()))
      print(head(downGeneList()))
      gene_list_string <- paste(c(topGeneList(), downGeneList()), collapse = ",")
      updateTextInput(session, "geneList", value = gene_list_string)
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
      print(targetGeneExpr)
      output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
    })
  }
  
  for_run <- shinyApp(ui = ui, server = server)
  runApp(for_run)
}
