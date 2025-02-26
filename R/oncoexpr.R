
#' @title filtering table for analysis (UI)
#' @description filtering table for analysis
#' @param id id
#' @return filtered table
#' @export
#' 
filterModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Sample Grouping"),
    selectInput(ns("mainCode_col"), "mainCode (Patient ID, Sample ID...)：", choices = NULL),
    uiOutput(ns("mainCode_values")),
    selectInput(ns("subCode_col"), "subCode (tissue type, drug, dose ...)：", choices = NULL),
    uiOutput(ns("subCode_values"))
  )
}
#' @title filtering table for analysis (server)
#' @description filtering table for analysis
#' @param id id
#' @param data data
#' @return filtered table
#' @export
filterModuleServer <- function(id, data) {
  moduleServer(id, function(input, output, session) {
    # 更新欄位選單：讓使用者從 data 中選擇篩選的欄位
    observe({
      req(data()) 
      print("獲取data()")
      updateSelectInput(session, "mainCode_col", choices = names(data()))
      updateSelectInput(session, "subCode_col", choices = names(data()))
    })
    
    # 根據選定欄位取得獨特值，並加入 "all" 選項（"all" 表示不進行進一步篩選）
    output$mainCode_values <- renderUI({
      req(input$mainCode_col)
      vals <- data() %>%
        select(!!sym(input$mainCode_col)) %>%
        distinct() %>%
        .[[1]] %>% as.character() %>% trimws()
      selectInput(session$ns("mainCode_select"), "選擇 mainCode 值：", 
                  choices = c("all", vals), multiple = TRUE, selected = "all")
    })
    
    output$subCode_values <- renderUI({
      req(input$subCode_col)
      vals <- data() %>%
        select(!!sym(input$subCode_col)) %>%
        distinct() %>%
        .[[1]] %>% as.character() %>% trimws()
      selectInput(session$ns("subCode_select"), "選擇 subCode 值：", 
                  choices = c("all", vals), multiple = TRUE, selected = "all")
    })
    
    # 篩選資料：當選擇 "all" 或沒有進行細部選擇時，就不進行該欄位的過濾
    filtered_data <- reactive({
      req(data())
      df <- data()
      if (!is.null(input$mainCode_select) && length(input$mainCode_select) > 0 &&
          !("all" %in% input$mainCode_select)) {
        df <- df %>% filter(as.character(.data[[input$mainCode_col]]) %in% input$mainCode_select)
      }
      if (!is.null(input$subCode_select) && length(input$subCode_select) > 0 &&
          !("all" %in% input$subCode_select)) {
        df <- df %>% filter(as.character(.data[[input$subCode_col]]) %in% input$subCode_select)
      }
      df
    })
    
    # 回傳篩選後的資料
    return(filtered_data)
  })
}



#' @title filtering table for viewing data (UI)
#' @description filtering table for viewing data
#' @param id id
#' @return filtered table
#' @export

mod_filter_sidebar_ui <- function(id) {
  ns <- NS(id)
  tagList(
    # 讓使用者選擇要篩選的欄位（例如：Sample, Tissue, Gene）
    selectizeInput(ns("filter_columns"),
                   "column to filter",
                   choices = NULL,
                   multiple = TRUE),
    # 根據選取的欄位，動態產生相應的多選下拉選單
    uiOutput(ns("value_select_ui"))
  )
}

#' @title  filtering table for viewing data (server)
#' @description filtering data with long format
#' @param id id
#' @param long_df data with long format
#' @return filtered table
#' @export


mod_filter_server <- function(id, long_df) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 使用 reactiveValues 儲存每個欄位目前的選擇值
    rv_filter <- reactiveValues()
    
    # 取得 colData 資料：排除掉 "Value" 欄位 (以 dplyr::select 實作)
    colData_df <- reactive({
      req(long_df())
      long_df() %>% select(-value)
    })
    
    # 更新可供篩選的欄位選項（依據 colData 的欄位名稱）
    observe({
      req(colData_df())
      updateSelectizeInput(session, "filter_columns",
                           choices = names(colData_df()),
                           server = TRUE)
    })
    
    # 當使用者變更選擇的欄位時，移除那些已儲存但不再選擇的欄位
    observe({
      req(input$filter_columns)
      stored <- names(rv_filter)
      for (col in stored) {
        if (!(col %in% input$filter_columns)) {
          rv_filter[[col]] <- NULL
        }
      }
    })
    
    # 根據使用者選擇的欄位動態產生多選下拉選單
    output$value_select_ui <- renderUI({
      req(colData_df())
      selected_cols <- input$filter_columns
      if (is.null(selected_cols) || length(selected_cols) == 0) {
        return(NULL)
      }
      ui_list <- lapply(selected_cols, function(col) {
        # 以 dplyr::distinct 取得該欄位的不重複值，並排序
        choices <- colData_df() %>% 
          distinct(.data[[col]]) %>% 
          pull() %>% 
          sort()
        choices <- c("All", choices)
        selected_val <- rv_filter[[col]]
        if (is.null(selected_val)) {
          selected_val <- "All"
        }
        selectizeInput(ns(paste0("filter_", col)),
                       label = paste("Select", col, "items"),
                       choices = choices,
                       selected = selected_val,
                       multiple = TRUE)
      })
      do.call(tagList, ui_list)
    })
    
    # 監控各個動態下拉選單，更新 rv_filter，同時若同時選到 "All" 與其他選項時，移除 "All"
    observe({
      req(input$filter_columns)
      for (col in input$filter_columns) {
        input_id <- paste0("filter_", col)
        cur_val <- tryCatch(input[[input_id]], error = function(e) NULL)
        if (!is.null(cur_val)) {
          if ("All" %in% cur_val && length(cur_val) > 1) {
            cur_val <- setdiff(cur_val, "All")
            updateSelectizeInput(session, input_id, selected = cur_val)
          }
          rv_filter[[col]] <- cur_val
        }
      }
    })
    
    # 根據儲存的篩選值，使用 dplyr::filter 動態過濾原始 long format 資料
    filtered_df <- reactive({
      req(long_df())
      df <- long_df()
      selected_cols <- input$filter_columns
      if (!is.null(selected_cols) && length(selected_cols) > 0) {
        for (col in selected_cols) {
          filter_val <- rv_filter[[col]]
          # 若該欄位的選項為 NULL 或包含 "All"，則不做篩選
          if (!is.null(filter_val) && !("All" %in% filter_val)) {
            df <- df %>% filter(.data[[col]] %in% filter_val)
          }
        }
      }
      df
    })
    
    # 產生 colData 表格輸出 (使用 dplyr::distinct 過濾重複)
    output$colData_table <- renderDT({
      if (!is.null(input$filter_columns) && length(input$filter_columns) > 0) {
        df <- colData_df() %>% 
          select(any_of(input$filter_columns)) %>% 
          distinct()
      } else {
        df <- colData_df() %>% distinct()
      }
      df
    }, options = list(pageLength = 5))
    
    # 產生依篩選條件過濾後的 long format 資料輸出
    output$filtered_data_table <- renderDT({
      filtered_df()
    }, options = list(pageLength = 5))
    
    return(list(filtered_data = filtered_df))
  })
}


#' @title colData for grouping to analysis
#' @description colData for grouping to analysis
#' @param filtered_table table with selecting column
#' @param mainCode_ main column for patient ID
#' @param subCode_ sub column for grouping or classification
#' @return filtered table
#' @export

generate_coldata_for_local_tcga_data <- function(filtered_table, mainCode_, subCode_){
  selected_df <- filtered_table %>% select(mainCode_, subCode_) %>% distinct() 
  print(dim(selected_df))
  print(mainCode_)
  print(subCode_)
  
  # 直接修改 data frame 的欄位名稱
  names(selected_df)[names(selected_df) == mainCode_] <- "mainCode"
  names(selected_df)[names(selected_df) == subCode_] <- "subCode"
  
  return(selected_df)
}




#' @title UI for connecting to spark 
#' @description spark connection shiny UI
#' @param id id
#' @return pivot long format from se object
#' @export


sparkConnectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    radioButtons(ns("conn_option"), "Connection Setting",
                 choices = c("Default" = "default", "Custom" = "custom"),
                 selected = "default", inline = TRUE),
    # 當選擇自定義連線時，顯示輸入參數的區塊
    conditionalPanel(
      condition = sprintf("input['%s'] == 'custom'", ns("conn_option")),
      textInput(ns("spark_master"), "Spark Master", value = "sc://localhost:15002"),
      textInput(ns("spark_method"), "Spark Method", value = "spark_connect"),
      textInput(ns("spark_version"), "Spark Version", value = "3.5")
      # 若有其他參數（例如 spark_home、其他 config 設定）也可在此新增
    ),
    actionButton(ns("connect"), "Connect")
  )
}


#' @title Server for connecting to spark 
#' @description spark connection shiny server
#' @param id id  
#' @return pivot long format from se object
#' @export

# Server module: 根據使用者輸入參數建立 spark 連線，並以 reactive 形式回傳 sc
sparkConnectionServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # 以 reactiveVal 儲存連線物件
    sc <- reactiveVal(NULL)
    
    observeEvent(input$connect, {
      # 根據使用者選擇決定參數設定
      if (input$conn_option == "default") {
        master <- "local"
        print(input$conn_option)
        #config <- sparklyr::spark_config()  # 使用預設 config
        connection <<- sparklyr::spark_connect(master = master)
        sc(connection)
      } else {
        master <- input$spark_master
        method <- input$spark_method
        version <- input$spark_version
        #config <- sparklyr::spark_config()
        # 若需要根據其他使用者輸入設定 config，可在此處調整
        connection <<- sparklyr::spark_connect(master = master, method = method, version = version )
        sc(connection)
      }
      
      # 嘗試建立 spark 連線
      
      
      
    })
    
    # 回傳 sc 連線物件（為 reactive）
    return(sc)
  })
}



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

#' @title limma analysis with MAE input
#' @description differential expression by limma 
#' @param mae multiple assays experiment
#' @param assayName assay 
#' @param subCodeCol which col for comparasion
#' @param coef coefficient for limma 
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @return limma analysis table 
#' @export

LimmaMAE <- function(mae,
                     assayName = "RNAseq",
                     # 假設分組資訊在 colData 中的 subCode 欄位
                     subCodeCol = "subCode",
                     coef = 2,            # 要比較的係數 (group 的哪一層)
                     pval_cut = 0.05,
                     lfc_cut  = 1,
                     useAdjP  = FALSE) {
  
  # 從 mae 中擷取資料：這裡示範從名為 "RNAseq" 的 experiment 拿資料
  # 且假設其中的 assays 包含一個名為 "max_TPM" 的矩陣
  rna_counts  <- assays(mae[[assayName]])$max_TPM
  sample_info <- colData(mae[[assayName]])
  print(sample_info)
  # 構建分組因子 (這裡以 subCode 作為分組資訊)
  group <- factor(sample_info[[subCodeCol]])
  print(group)
  # 建立 DGEList，並進行標準化
  dge <- DGEList(counts = rna_counts, group = group)
  dge <- calcNormFactors(dge)  # TMM 等方法
  
  # 建置設計矩陣 (含截距項)
  design <- model.matrix(~ group)
  
  # voom 轉換
  v <- voom(dge, design)
  
  # 線性模型擬合 + Bayes校正
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  
  # 使用 decideTests() 判定哪些基因顯著差異
  # 在此設定 p.value 與 lfc 閾值，作為上/下調判定標準
  dt <- decideTests(fit, p.value = pval_cut, lfc = lfc_cut)
  sumDT <- summary(dt)  # 統計上調/下調/不顯著的數目
  
  # 用 topTable() 匯出完整結果 (logFC, P.Value, adj.P.Val 等)
  # coef = 2 代表第二個參數 (group 的對比)，請根據實際情況調整
  topTab <- topTable(fit, coef = coef, number = Inf)
  
  # 最後以 list 形式回傳，以便後續使用
  return(list(
    fit            = fit,      # 可給 ggvolcano_limma() 繪圖
    topTable       = topTab,   # 所有基因的差異分析統計
    decideTests    = dt,       # limma 對基因上調/下調/不顯著的判定
    upDownSummary  = sumDT     # 上下調總結
  ))
}

#' @title limma analysis with MAE input
#' @description differential expression by limma 
#' @param fit limma fit result
#' @param assayName assay 
#' @param subCodeCol which col for comparasion
#' @param coef coefficient for limma 
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @param ptAlpha transparent
#' @param labelSize gene symbol text size
#' @param pointSize point plot size
#' @param topN show a number of genes symbol with high log fold-change
#' @return limma analysis table 
#' @export


ggvolcano_limma <- function(fit,
                            coef      = 2,
                            lfc_cut   = 1,
                            pval_cut  = 0.05,
                            useAdjP   = FALSE,
                            title     = "Volcano Plot",
                            topN      = 0,
                            geneCol   = NULL,
                            pointSize = 2,       # 加入點大小參數
                            ptAlpha   = 0.6,     # 加入點透明度參數
                            labelSize = 3,       # 加入基因標籤字型大小參數
                            ...) {
  
  # 1) 從 limma 的 fit 物件取出指定係數的結果
  tt <- topTable(
    fit      = fit,
    coef     = coef,
    number   = Inf,
    ...
  )
  
  # 如果沒有 geneCol，則使用 rownames 作為 gene 名稱
  if (is.null(geneCol)) {
    tt$gene <- rownames(tt)
  } else {
    if (!geneCol %in% colnames(tt)) {
      stop("您指定的 geneCol 不存在於 topTable() 輸出的欄位中，請確認。")
    }
    tt$gene <- tt[[geneCol]]
  }
  
  # 2) 選擇 p-value 或 adj.P.Val
  pvalCol <- if (useAdjP) "adj.P.Val" else "P.Value"
  if (!pvalCol %in% colnames(tt)) {
    stop(paste0("無法在 topTable 結果中找到欄位: ", pvalCol))
  }
  
  # 3) 建立繪圖用的資料框
  plotData <- data.frame(
    gene  = tt$gene,
    logFC = tt$logFC,
    pval  = tt[[pvalCol]]
  )
  
  # 4) 根據閾值新增顏色欄位（上調：red，下調：blue，其餘 grey）
  plotData$color <- "grey"
  plotData$color[plotData$logFC >=  lfc_cut & plotData$pval <= pval_cut] <- "red"
  plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"
  
  # 5) 計算 -log10(p-value)
  plotData$negLogP <- -log10(plotData$pval)
  
  # 6) 繪製基本圖形，將點大小、透明度參數納入
  p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) +
    geom_point(size = pointSize, alpha = ptAlpha) +
    scale_color_identity() +
    geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
    theme_classic() +
    labs(
      title = title,
      x     = expression(log[2]~"Fold Change"),
      y     = bquote(-log[10]~.(pvalCol))
    )
  
  # 7) 標註前 topN 個最顯著基因（若有指定 topN > 0）
  if (topN > 0) {
    topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
    p <- p + 
      ggrepel::geom_text_repel(
        data = topGenesDF,
        aes(label = gene),
        size = labelSize,
        max.overlaps = Inf
      )
  }
  
  return(p)
}


#' @title gtf to expr profile
#' @description 輸入RNA-seq pipline生成的 .gtf file
#' @param gtf_data gtf_data的絕對路徑
#' @return 經過處理後的表現量表格，包括Cov, FPKM, TPM三種量化指標
#' @export

gtf2expr <- function(gtf_data) {
  colnames(gtf_data) <- c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")
  
  gtf_data <- gtf_data %>%
    filter(feature == "transcript") %>%
    mutate(
      gene_id = ifelse(grepl('gene_id', attribute), sub('.*gene_id "([^"]+)".*', '\\1', attribute), NA),
      Gene = ifelse(grepl('gene_name', attribute), sub('.*gene_name "([^"]+)".*', '\\1', attribute), NA),
      cov = ifelse(grepl('cov', attribute), as.numeric(sub('.*cov "([^"]+)".*', '\\1', attribute)), NA),
      FPKM = ifelse(grepl('FPKM', attribute), as.numeric(sub('.*FPKM "([^"]+)".*', '\\1', attribute)), NA),
      TPM = ifelse(grepl('TPM', attribute), as.numeric(sub('.*TPM "([^"]+)".*', '\\1', attribute)), NA)
    )
  
  gene_expression <- gtf_data %>%
    group_by(Gene) %>%
    summarise(
      max_cov = ifelse(all(is.na(cov)), NA, max(cov, na.rm = TRUE)),
      max_FPKM = ifelse(all(is.na(FPKM)), NA, max(FPKM, na.rm = TRUE)),
      max_TPM = ifelse(all(is.na(TPM)), NA, max(TPM, na.rm = TRUE))
    )
  
  return(gene_expression)
}

#' @title prepare for expression pattern
#' @description compare two samples' expression pattern
#' @param group1 first sample 
#' @param group2 sencond sample
#' @param groups_name sample type name e.g. tissue & ctc
#' @return a list with dataframe and specific gene list 
#' @export

Prep4ExprComp <- function(group1, group2, groups_name = c("tissue", "ctc")){
  group2_matched <- group2[match(group1$"Gene", group2$"Gene"),]
  data_ <- data.frame(group1$max_TPM, group2_matched$max_TPM)
  rownames(data_) <- group1$Gene
  print(colnames(data_))
  rm_noise_ <- data_[,1] > 1 | data_[,2] > 1
  
  only_group1_high_gene_profile <- data_[data_[,"group1.max_TPM"] < 100 & data_[,"group2_matched.max_TPM"] < 0.001,]
  only_group2_high_gene_profile <- data_[data_[,"group1.max_TPM"] < 0.001 & data_[,"group2_matched.max_TPM"] < 100,]
  
  special_gene_list <- c(rownames(only_group1_high_gene_profile), rownames(only_group2_high_gene_profile))
  
  data_v2 <- as.data.frame(data_[rm_noise_,])
  data_v2 <- data_v2 + 0.0001
  colnames(data_v2) <- groups_name
  list_ <- list("special_gene_list" = special_gene_list, "result" = data_v2)
  return(list_)
  
}

#' @title target gene expression
#' @description 針對目標基因回傳兩個樣本的基因表現量
#' @param selectMode select or exclude to gene list
#' @param geneList_ interesting genes
#' @param groups_list_ term of group
#' @param save_path_ save path
#' @param save_filename_  save name
#' @param expr_profile_ target gene expression profile
#' @return the expression of genes which you are interesting
#' @export

target_exprofile <- function(selectMode = TRUE, geneList_, groups_list_ , save_path_ = NULL, save_filename_ = "target_gene_expr.csv", expr_profile_){
  
  target_exprprfile_ <-  expr_profile_ %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  
  colnames(target_exprprfile_) <- groups_list_
  
  if(!is.null(save_path_)){
    #collect from instance to local 
    write.csv(target_exprprfile_, paste0(save_path_, save_filename_))
    print(paste("save at the path:", paste0(save_path_, save_filename_)))
    
  }else(
    
    print("do not save file")
    
  )
  
  return(target_exprprfile_ )
  
}

#' @title expression pattern analysis
#' @description expression pattern analysis for two samples by scatter plot 
#' @param data_v2 dataframe derived from Prep4ExprComp function
#' @param special_gene_list gene list derived from Prep4ExprComp function
#' @param save_path_ where you want to save file
#' @param save_filename_ plot filename
#' @param cutoff_tpm what target genes would be shown in red color
#' @param x_col x axis tile of scatter plot
#' @param y_col y axis tile of scatter plot
#' @return scatter plot for expression pattern analysis
#' @export

expr_pattern <- function(data_v2, special_gene_list, save_path_ = NULL, save_filename_ = "acrocyte_scatter_plot.png", cutoff_tpm = 7, x_col = "tissue", y_col = "ctc"){
  data_v2$distance_to_diagonal <- abs(log2(data_v2[,x_col]) - log2(data_v2[,y_col]))
  color_ <- ifelse(data_v2$distance_to_diagonal>cutoff_tpm, "outlier", "consistent")
  data_v2$color_ <- color_
  genes_beyond <- data_v2[data_v2$distance_to_diagonal>cutoff_tpm, ]
  genes_beyond_v2 <- genes_beyond[!is.element(rownames(genes_beyond), special_gene_list), ]
  genes_beyond_v2$label_ <- rownames(genes_beyond_v2)
  scatter_plot <- ggplot(data_v2, aes(x = .data[[x_col]], y = .data[[y_col]]))+
    geom_point(aes(color = color_), size=1.5, alpha=0.4) +  
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") + 
    geom_text_repel(data = genes_beyond_v2, aes(label = .data[["label_"]]), size = 1.5, force = 1, max.overlaps = Inf, segment.size=0.1) +
    labs(x = x_col, y =  y_col ) +
    scale_x_continuous(limits = c(NA, 200000),expand = c(0, 0), trans = "log10", 
                       labels = scales::number_format(accuracy = 0.0001)) +
    scale_y_continuous(limits = c(NA, 200000),expand = c(0, 0), trans = "log10", 
                       labels = scales::number_format(accuracy = 0.0001)) +
    coord_fixed(ratio = 1) +
    scale_color_manual(values = c("blue", "red")) 
  save_parameter <- paste0(save_path_, save_filename_)
  ggsave(save_parameter, plot = scatter_plot , width = 28, height = 28, units="cm",dpi = 300)
  return(scatter_plot)
  
}

#' @title Gene Ontology Enrichment Analysis
#' @description a set of genes for enrichment analysis with Gene Ontology database
#' @param gene_list_ a set of genes with specific expression level
#' @param save_path_ what path you want to save the files
#' @param save_filename_ what name you call the files
#' @param mode_ Gene Ontology mode such as MF, BP, CC
#' @param showCategory_ the number of GO term you want to show in doplot, default: 10
#' @return the dotplot for GO enrichment
#' @export
go_enrich_dotplot <- function(gene_list_, save_path_=NULL, save_filename_=NULL, mode_ = "BP", showCategory_ = 10){
  gene_ids <- bitr(gene_list_, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  go_enrich_result <- enrichGO(gene = gene_ids$ENTREZID,
                               OrgDb = org.Hs.eg.db,
                               keyType = "ENTREZID",
                               ont = mode_,
                               pvalueCutoff = 0.05, 
                               qvalueCutoff = 0.2,
                               readable = TRUE)
  go_dotplot_ <- dotplot(go_enrich_result, showCategory = showCategory_, font.size = 15)
  save_parameter <- paste0(save_path_, paste0(mode_,"_"), save_filename_)
  if(length(save_path_)!=0){
    ggsave(save_parameter, plot(go_dotplot_), width = 32, height = 28,  units="cm", dpi = 300)
  }else{
    print("not saving the figure")
  }
  return(go_dotplot_) 
  
}
#' @title KEGG Enrichment Analysis
#' @description a set of genes for enrichment analysis with KEGG database
#' @param gene_list_ a set of genes with specific expression level
#' @param save_path_ what path you want to save the files
#' @param save_filename_ what name you call the files
#' @param showCategory_ the number of GO term you want to show in doplot, default: 10
#' @return the dotplot for KEGG enrichment
#' @export
#' 
kegg_enrich_dotplot <- function(gene_list_, save_path_=NULL, save_filename_=NULL, showCategory_ = 10){
  save_parameter <- paste0(save_path_, save_filename_)
  gene_ids <- bitr(gene_list_, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  kegg_enrich_result <- enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05)
  kegg_dotplot_ <- dotplot(kegg_enrich_result, showCategory = showCategory_ , font.size = 15)
  if(length(save_path_)!=0){
    ggsave(   save_parameter, 
              plot = kegg_dotplot_, 
              width = 32, 
              height = 28, 
              units="cm", 
              dpi = 300)
  }else{
    print("not saving the figure")
  }
  
  return(kegg_dotplot_)
}

#' @title prepare data
#' @description precessing the gtf file into expression profile
#' @param entirepath entire path for gtf files
#' @param groups_list the group name of samples
#' @return expression profile
#' @export
input_and_prep_file <-  function(entirepath , groups_list ){
  groups_path_variable <- paste0(groups_list, c("_gtf_file"))
  groups_data_variable <- paste0(groups_list, c("_gtf_data"))
  groups_expr_variable <- paste0(groups_list, c("_rna_expr"))
  print(groups_expr_variable)
  
  names(groups_path_variable) <- groups_list
  names(groups_data_variable) <- groups_list
  names(groups_expr_variable) <- groups_list
  
  print(groups_expr_variable)
  lapply(seq_len(length(groups_list)), function(i) assign(groups_data_variable[[groups_list[i]]], read_tsv(get(entirepath[i]), comment = "#", col_names = FALSE), envir=.GlobalEnv) )
  lapply(seq_len(length(groups_list)), function(i) assign(groups_expr_variable[[groups_list[i]]], gtf2expr(get(groups_data_variable[[groups_list[i]]])), envir=.GlobalEnv) )
  
  print(dim(groups_expr_variable[[groups_list[1]]]))
  print(dim(groups_expr_variable[[groups_list[2]]]))
  list_ <- list("path" = groups_path_variable, "gtf" = groups_data_variable, "expr" = groups_expr_variable)
  return(list_)
  
}

#' @title prepare data
#' @description precessing the gtf file into expression profile
#' @param entirepath entire path for gtf files
#' @param groups_list the group name of samples
#' @return expression profile
#' @export

input_and_prep_file_directly <-  function(entirepath, groups_list){
  groups_path_variable <- paste0(groups_list, c("_gtf_file"))
  groups_data_variable <- paste0(groups_list, c("_gtf_data"))
  groups_expr_variable <- paste0(groups_list, c("_rna_expr"))
  print(groups_expr_variable)
  names(groups_path_variable) <- groups_list
  names(groups_data_variable) <- groups_list
  names(groups_expr_variable) <- groups_list
  
  print(groups_expr_variable)
  assign(groups_path_variable[[groups_list[1]]], entirepath[1], envir=.GlobalEnv)
  assign(groups_path_variable[[groups_list[2]]], entirepath[2], envir=.GlobalEnv)
  
  assign(groups_data_variable[[groups_list[1]]], read_tsv(get(groups_path_variable[[groups_list[1]]]), comment = "#", col_names = FALSE), envir=.GlobalEnv)
  assign(groups_data_variable[[groups_list[2]]], read_tsv(get(groups_path_variable[[groups_list[2]]]), comment = "#", col_names = FALSE), envir=.GlobalEnv)
  
  
  assign(groups_expr_variable[[groups_list[1]]], gtf2expr(get(groups_data_variable[[groups_list[1]]])), envir=.GlobalEnv)
  assign(groups_expr_variable[[groups_list[2]]], gtf2expr(get(groups_data_variable[[groups_list[2]]])), envir=.GlobalEnv)
  print(dim(groups_expr_variable[[groups_list[1]]]))
  print(dim(groups_expr_variable[[groups_list[2]]]))
  list_ <- list("path" = groups_path_variable, "gtf" = groups_data_variable, "expr" = groups_expr_variable)
  return(list_)
  
}
#' @title TissueTyper
#' @description predict the original tissue what the sample from 
#' @param group1_ sample 1
#' @param group2_ sample 2
#' @param save_path what the path you want to save the files
#' @param gtexFilePath "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
#' @param groups_list sample types
#' @param sample_id samply id
#' @return heatmap plot for correlation coefficient between multiple tissues of GTEx and your samples
#' @export
TissueTyper <- function(selectMode = FALSE, geneList_ = NULL , group1_ , group2_ , save_path, gtexFilePath, groups_list = c("tissue", "ctc"), sample_id ="BC021-RNA"){
  print("讀取檔案")
  date_ <- Sys.Date() 
  print(date_)
  gtex_tissue_gct_data <- read.delim(gtexFilePath, skip=2)
  head(gtex_tissue_gct_data)
  #gtex_gene_expression <- gtex_tissue_gct_data %>%
  #                          dplyr::group_by(Description) %>%
  #                          dplyr::summarise(across(where(is.numeric), max))%>%
  #                          ungroup()                   
  print("將表現量型別轉成數值並取最大值來合併重複基因")
  
  gtex_gene_expression <- gtex_tissue_gct_data %>%
    group_by(Description) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
    ungroup()
  print("取最大值")                                    
  colnames(gtex_gene_expression)[2] <- "ENSG_id"
  colnames(gtex_gene_expression)[1] <- "Gene"
  print("修改colnames")
  group1_tpm <- group1_[,c("Gene","max_TPM")]
  group2_tpm <- group2_[,c("Gene","max_TPM")]
  print("取得樣本資訊")
  
  colnames(group1_tpm) <- c("Gene", paste0(groups_list[1],"_TPM"))
  colnames(group2_tpm) <- c("Gene", paste0(groups_list[2],"_TPM"))
  print("修改colnames")
  
  merged_df <- gtex_gene_expression %>%
    right_join(group1_tpm, by = "Gene") %>%
    right_join(group2_tpm, by = "Gene")
  print("合併GTEx和兩樣本表現量")
  merged_df <- as.data.frame(na.omit(merged_df))
  print(dim(merged_df))
  #20241114 add new function: select genes
  rownames(merged_df) <- merged_df$"Gene"
  merged_df <-  merged_df %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  print(dim(merged_df))
  cor_matrix<- cor.fk(merged_df[,-c(1:2)]) #"kendall 用cor.fk"
  print("計算kendall's tau")
  cor_data <- melt(cor_matrix)
  colnames(cor_data) <- c("Group1", "Group2", "Correlation")
  
  group1_cor_result <- cor_data[grepl(paste0(groups_list[1],"_TPM"), cor_data$'Group1'),]
  group2_cor_result <- cor_data[grepl(paste0(groups_list[2],"_TPM"), cor_data$'Group1'),]
  
  group1_cor_rank <- group1_cor_result[order(group1_cor_result$'Correlation', decreasing=TRUE),]
  group2_cor_rank <- group2_cor_result[order(group2_cor_result$'Correlation', decreasing=TRUE),]
  print("完成相關係數排序，準備輸出成檔案")
  write.csv(group1_cor_rank, paste0(save_path, paste(c(sample_id, date_ , groups_list[1], "cor_rank.csv"), collapse="_")))
  write.csv(group2_cor_rank, paste0(save_path, paste(c(sample_id, date_ , groups_list[2], "cor_rank.csv"), collapse="_")))
  print("準備繪製成heatmap")
  color_palette <- colorRampPalette( c("white","orange", "red"))(100)
  breaks <- seq( 0, 1, length.out = length(color_palette) + 1)
  acrocyte_heatmap_cor_plot <- pheatmap(
    cor_matrix,
    color = color_palette,
    breaks = breaks,
    na_col = "grey", main = "TissueTyper"
  )
  print("輸出heatmap")
  print(paste0(sample_id, "_", date_, "_acrocyte_GTEx_cor_analysis.png"))
  ggsave(   paste0(sample_id, "_", date_ ,"_acrocyte_GTEx_cor_analysis.png"), 
            plot = acrocyte_heatmap_cor_plot, 
            path = save_path, 
            width = 30, 
            height = 28, 
            units="cm", 
            dpi = 300
  )
  list_ <- list("plot"= acrocyte_heatmap_cor_plot, "rank1"= group1_cor_rank, "rank2"= group2_cor_rank)
  return(list_ )
  
  
  
}



#' @title multiple samples TissueTyper
#' @description predict the original tissue what the sample from 
#' @param selectMode select or excluding
#' @param geneList_ target gene set
#' @param sample_list_ multiple sample id  
#' @param save_path what the path you want to save the files
#' @param gtexFilePath "GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct"
#' @param groups_list sample types
#' @return heatmap plot for correlation coefficient between multiple tissues of GTEx and your samples
#' @export
#' 
#' 

mTissueTyper <- function(selectMode = FALSE, geneList_ = NULL , sample_list_, save_path, gtexFilePath, groups_list = c("tissue", "ctc")){
  print(sample_list_)
  rna_expr <- paste0(unlist(sample_list_), c("_rna_expr"))
  samples_data <- mget(rna_expr, envir = .GlobalEnv)  # 防止變數不存在
  print(samples_data[[1]])
  print("讀取檔案")
  date_ <- Sys.Date() 
  print(date_)
  gtex_tissue_gct_data <- read.delim(gtexFilePath, skip=2)
  head(gtex_tissue_gct_data)
  #gtex_gene_expression <- gtex_tissue_gct_data %>%
  #                          dplyr::group_by(Description) %>%
  #                          dplyr::summarise(across(where(is.numeric), max))%>%
  #                          ungroup()                   
  print("將表現量型別轉成數值並取最大值來合併重複基因")
  
  gtex_gene_expression <- gtex_tissue_gct_data %>%
    group_by(Description) %>%
    summarise(across(where(is.numeric), ~ max(.x, na.rm = TRUE))) %>%
    ungroup()
  print("取最大值")                                    
  colnames(gtex_gene_expression)[2] <- "ENSG_id"
  colnames(gtex_gene_expression)[1] <- "Gene"
  print("修改colnames")
  library(purrr)
  processed_list <- lapply(seq_along(samples_data), function(i) {
    samples_data[[i]] %>%
      select(Gene, max_TPM) %>%
      rename(!!unlist(sample_list)[i] := max_TPM)
  })
  
  
  # 使用 reduce 合併處理後的資料框
  merged_df <- purrr::reduce(processed_list, function(df1, df2) {
    right_join(df1, df2, by = "Gene")
  }, .init = gtex_gene_expression)
  
  write.csv(merged_df, paste0(save_path, paste0(paste(patient_code, collapse="_"), "_", date_, "_acrocyte_GTEx_merged_df.csv")))
  
  print("合併GTEx和兩樣本表現量")
  merged_df <- as.data.frame(na.omit(merged_df))
  print(dim(merged_df))
  #20241114 add new function: select genes
  rownames(merged_df) <- merged_df$"Gene"
  merged_df <-  merged_df %>% filter(  if (selectMode) rownames(.) %in% geneList_ else !rownames(.) %in% geneList_ )
  print(dim(merged_df))
  cor_matrix <- cor.fk(merged_df[,-c(1:2)]) #"kendall 用cor.fk"
  print("計算kendall's tau")
  cor_data <- melt(cor_matrix)
  colnames(cor_data) <- c("Group1", "Group2", "Correlation")
  print("處理數據表格")
  cor_results <- paste0(unlist(sample_list_), "_cor_result")
  print(cor_results)
  cor_results_order <- paste0(unlist(sample_list_), "_cor_order")
  lapply(1:length(unlist(sample_list_)), function(s) assign(cor_results[s], cor_data[grepl(unlist(sample_list_)[s], cor_data$'Group1'),], envir = .GlobalEnv) )
  print(get(cor_results[1]))
  mget_cor_results <- mget(cor_results, envir = .GlobalEnv)
  
  lapply(seq_along(mget_cor_results), function(s) {
    assign(
      cor_results_order[s], 
      mget_cor_results[[s]][order(mget_cor_results[[s]]$Correlation, decreasing = TRUE), ], 
      envir = .GlobalEnv
    )
    write.csv(get(cor_results_order[s]), paste0(save_path, paste(c(unlist(sample_list_)[s], date_ , "cor_rank.csv"), collapse="_")))
  })
  print("完成相關係數排序，準備輸出成檔案")
  
  print("準備繪製成heatmap")
  color_palette <- colorRampPalette( c("white","orange", "red"))(100)
  breaks <- seq( 0, 1, length.out = length(color_palette) + 1)
  acrocyte_heatmap_cor_plot <- pheatmap(
    cor_matrix,
    color = color_palette,
    breaks = breaks,
    na_col = "grey", main = "TissueTyper"
  )
  print("輸出heatmap")
  print(paste0(paste(patient_code, collapse="_"), "_", date_, "_acrocyte_GTEx_cor_analysis.png"))
  ggsave( paste0(paste(patient_code, collapse="_"), "_", date_, "_acrocyte_GTEx_cor_analysis.png"), 
          plot = acrocyte_heatmap_cor_plot, 
          path = save_path, 
          width = 30, 
          height = 28, 
          units="cm", 
          dpi = 300
  )
  list_ <- list("plot"= acrocyte_heatmap_cor_plot, "merged_df" = merged_df, "cor_matrix" = cor_matrix)
  return(list_ )
  
  
  
}


#' @title Immune Cell Fraction Prediction
#' @description predict the composition of immune cell types for your samples
#' @param group1_ sample 1
#' @param group2_ sample 2
#' @param groups_list sample types
#' @param method quantiseq, xCell, CIBERSORT
#' @param tumorStatus if your sample is tumor TRUE or FALSE
#' @param save_path what the path you want to save the files
#' @return plot for immune cell faction of your samples
#' @export
ImmuneFractionPlot <- function(group1_, group2_, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL){
  
  tryCatch({  
    print("ensuring the dataframe be correct")
    if (!all(c("Gene", "max_TPM") %in% colnames(group1_)) ||
        !all(c("Gene", "max_TPM") %in% colnames(group2_))) {
      stop("Input data frames must contain 'Gene' and 'max_TPM' columns")
    }
    # 2. 匹配基因並創建數據框
    print("match genes position")
    group2_matched <- group2_[match(group1_$"Gene", group2_$"Gene"), ]
    
    # 3. 創建並處理合併後的數據框
    print("creat merged dataframe")
    df <- data.frame(
      Gene = group1_$"Gene",
      group1 = group1_$"max_TPM",
      group2 = group2_matched$"max_TPM",
      stringsAsFactors = FALSE
    )
    # 4. 移除 NA 值的行
    print("remove na value")
    df <- df[complete.cases(df), ]
    
    print("set colnames")        
    colnames(df) <- c("Gene", groups_list)
    
    # 6. 準備用於 deconvolution 的矩陣
    print("transfer Gene column to be rownames")
    expression_matrix <- df %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()
    
    print(head(expression_matrix))
    # 7. 運行 deconvolution
    res_quantiseq <- immunedeconv::deconvolute(
      expression_matrix,
      method = method_,
      tumor = tumorStatus 
    )
    
    
  }, error = function(e) {
    message("Error in immuneFractionPlot: ", e$message)
    # 診斷信息
    message("\nDiagnostic information:")
    message("group1 structure:")
    print(str(group1_))
    message("\ngroup2 structure:")
    print(str(group2_))
    message("\ngroups_list:")
    print(groups_list)
    return(NULL)
  })
  my_colors <- c(
    # T細胞
    "T cell CD4+ (non-regulatory)" = "#D62728",     # 深紅色
    "T cell CD8+" = "#FF7F0E",                      # 橙色
    "T cell regulatory (Tregs)" = "#FF9896",        # 浅红色
    
    # Macrophage M1 和 M2
    "Macrophage M1" = "#2CA02C",                    # 綠色
    "Macrophage M2" = "#98DF8A",                    # 淺綠色
    # NK cell 和 B cell
    "NK cell" = "#1F77B4",                          # 藍色
    "B cell" = "#AEC7E8",                           # 淺藍色
    # Neutrophil 和 Monocyte
    "Neutrophil" = "#9467BD",                       # 深紫色
    "Monocyte" = "#8C564B",                         # 棕色
    # Myeloid dendritic cell
    "Myeloid dendritic cell" = "#BC80BD",           # 淺紫色
    # 未分類细胞
    "uncharacterized cell" = "#D0E1F9"              # 淺粉藍色
  )
  
  p <-  res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    #scale_fill_manual(values = my_colors)+
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  return(p)
  ggsave(past0(method_, "_", "immuneFractionPlot.jpg"), path = save_path ,plot = p, device = "png", width = 35, height = 25, units = "cm", dpi = 300)
  
  
}

#' @title Immune Cell Fraction Prediction
#' @description predict the composition of immune cell types for your samples
#' @param mae multiple assays experiments
#' @param groups_list sample types
#' @param method quantiseq, xCell, CIBERSORT
#' @param tumorStatus if your sample is tumor TRUE or FALSE
#' @param save_path what the path you want to save the files
#' @return plot for immune cell faction of your samples
#' @export
ImmuneFractionPlot_mae <- function(mae, groups_list, method_ = "quantiseq", tumorStatus = TRUE, save_path = NULL){
  
  rna_counts  <- assays(mae[["RNAseq"]])$max_TPM
  # 6. 準備用於 deconvolution 的矩陣
  print("transfer Gene column to be rownames")
  
  expression_matrix <- rna_counts %>%
    #tibble::column_to_rownames("Gene") %>% #rownames 要是gene symbol
    as.matrix()
  
  print(head(expression_matrix))
  # 7. 運行 deconvolution
  res_quantiseq <- immunedeconv::deconvolute(
    expression_matrix,
    method = method_,
    tumor = tumorStatus 
  )
  
  
  
  my_colors <- c(
    # T細胞
    "T cell CD4+ (non-regulatory)" = "#D62728",     # 深紅色
    "T cell CD8+" = "#FF7F0E",                      # 橙色
    "T cell regulatory (Tregs)" = "#FF9896",        # 浅红色
    
    # Macrophage M1 和 M2
    "Macrophage M1" = "#2CA02C",                    # 綠色
    "Macrophage M2" = "#98DF8A",                    # 淺綠色
    # NK cell 和 B cell
    "NK cell" = "#1F77B4",                          # 藍色
    "B cell" = "#AEC7E8",                           # 淺藍色
    # Neutrophil 和 Monocyte
    "Neutrophil" = "#9467BD",                       # 深紫色
    "Monocyte" = "#8C564B",                         # 棕色
    # Myeloid dendritic cell
    "Myeloid dendritic cell" = "#BC80BD",           # 淺紫色
    # 未分類细胞
    "uncharacterized cell" = "#D0E1F9"              # 淺粉藍色
  )
  
  p <-  res_quantiseq %>%
    gather(sample, fraction, -cell_type) %>%
    mutate(sample = factor(sample, levels = rev(unique(sample)))) %>%
    ggplot(aes(x = fraction, y = sample, fill = cell_type)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    scale_fill_brewer(palette = "Paired") +
    #scale_fill_manual(values = my_colors)+
    scale_y_discrete(limits = rev(levels(res_quantiseq))) +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14)
    )
  return(p)
  #ggsave(past0(method_, "_", "immuneFractionPlot.jpg"), path = save_path ,plot = p, device = "png", width = 35, height = 25, units = "cm", dpi = 300)
  
  
}




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

#' @title sampleSelect shinmy UI module
#' @description set shiny UI module for sample selection
#' @param id module id for UI and Server
#' @return UI module 
#' @export

sampleSelectionUI <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Sample Selection (Module)"),
    textInput(ns("mainCode_input"), "Main Code:", value = "BC026"),
    textInput(ns("subCode_input"), "Sub Code:", value = "PBMC"),
    selectInput(
      ns("normMethod_input"),
      "Normalization Method:",
      choices = c("max_Cov", "max_FPKM", "max_TPM"),
      selected = "max_TPM"
    ),
    actionButton(ns("preview_btn"), "Query & Preview test")
    # 不要放 output$file_table，或用其他方式在主畫面顯示
  )
}


#' @title sampleSelect shinmy Server module
#' @description set shiny Server module for sample selection
#' @param id module id for UI and Server
#' @param sc spark connection
#' @param selected_db select database from spark 
#' @param selected_table select deltatable from selected database
#' @param set_colData a reactive Value for storage dynamic dataframe
#' @return Sample selection Server module 
#' @export

sampleSelectionServer <- function(
    id,
    sc,
    selected_db,
    selected_table,
    set_colData
) {
  
  
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    #spark_entired_preview <- reactiveVal(NULL) #用來動態儲存spark裡面完整的表格呈現於UI，因為檔案較大，只取前20行，回傳至local
    preview_data_filtered <- reactiveVal(NULL) #用來動態儲存基於篩選的特定mainCode subCode回傳至local，以便後續分析
    
    print(reactive({input$mainCode_input}))
    print(reactive({input$subCode_input}))
    # 將 input$mainCode_input / input$subCode_input -> 向量
    mainCodes <- reactive({
      strsplit(input$mainCode_input, ",")[[1]] |> trimws()
    })
    subCodes <- reactive({
      strsplit(input$subCode_input, ",")[[1]] |> trimws()
    })
    #print(mainCodes())
    #print(subCodes())
    observeEvent(input$preview_btn, {
      # 切換到選定的資料庫
      print("確定按下按鈕了")
      DBI::dbExecute(sc, paste0("USE spark_catalog.", selected_db()))
      print("確定切換資料庫")
      #req(selected_db(), selected_table())  # 確保有 DB 與 Table
      #req(mainCodes(), subCodes(), input$normMethod_input)
      print("確定req的東西都有了")
      
      # 呼叫您自訂的查詢函式 (需事先在全域或其他地方定義好)
      # 注意這裡的 use_table = selected_table() 是取 reactive 的值
      result <- query_spark_rnaseqdata(
        sc         = sc,
        use_table  = selected_table(),
        mainCode   = mainCodes(),
        subCode    = subCodes(),
        normMethod = input$normMethod_input
      )
      print(head(result$data))
      # 檢查結果
      
      df_renamed <- result$data |> dplyr::rename(
        mainCode = main_code,
        subCode  = sub_code
      )
      preview_data_filtered(df_renamed)  # 回傳至reactive儲存
      
      print(dim(preview_data_filtered()))
      print(colnames(df_renamed))
      
      # 如果您想用部分資料展示預覽，可以先取 head()
      #spark_data_preview <- df_renamed %>% head(20) %>% collect()
      #spark_entired_preview(spark_data_preview)
      print(colnames(df_renamed))
      required_cols <- c("mainCode","subCode")
      print(required_cols %in% colnames(df_renamed))
      
      distinct_data <- df_renamed %>%
        select(mainCode, subCode) %>%
        distinct()
      print(dim(distinct_data))
      set_colData(distinct_data) # sample information (MAE colData)
      print(head(set_colData()))
      
      
    })
    
    
    return(list(
      #這裡只能回傳reactiveVal()，不能回傳moduleServer 裡面的 df
      #previewTable = spark_entired_preview,
      subColData = set_colData,
      filteredTable  = preview_data_filtered
      
    ))
  })
}



#' @title dbBrowerUI
#' @description search spark database
#' @param id module id for UI and Server
#' @return database browser UI module 
#' @export

dbBrowserUI <- function(id) {
  ns <- NS(id)
  tagList(
    # 可以是 fluidRow/sidebarPanel/... 依您需要排版
    selectInput(
      inputId = ns("selected_db"),
      label   = "Choose Database:",
      choices = character(0)  # 初始給空
    ),
    # 動態顯示對應資料表選項
    #uiOutput(ns("table_selector"))
  )
}

#' @title dbBrowerSever
#' @description search spark database
#' @param id module id for UI and Server
#' @param sc spark connection
#' @return database browser Server module 
#' @export
dbBrowserServer <- function(id, sc) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # 動態更新資料庫清單
    observe({
      db_list_query <- dbGetQuery(sc, "SHOW DATABASES")
      # 假設第一欄就是資料庫名稱
      db_list <- db_list_query
      
      # 用 updateSelectInput 來更新 UI
      updateSelectInput(
        session,
        "selected_db",
        choices = db_list,
        selected = db_list[1]
      )
    })
    
    # 動態生成 "選擇資料表" 的下拉選單
    # output$table_selector <- renderUI({
    #   req(input$selected_db)
    #   # 切換到使用者選擇的 database
    #   DBI::dbExecute(sc, paste0("USE ", input$selected_db))
      
    #   tbl_list_query <- DBI::dbGetQuery(sc, paste0("SHOW TABLES IN ", input$selected_db))
    #   # 假設欄位名為 tableName
    #   tbls <- tbl_list_query$tableName
    #   if (length(tbls) == 0) {
    #     selectInput(ns("selected_table"), "Choose Table:", choices = NULL)
    #   } else {
    #     selectInput(ns("selected_table"), "Choose Table:", choices = tbls, selected = tbls[1])
    #   }
    # })


    # 為了讓主應用拿到目前選擇的 DB & Table，回傳 reactive
     selected_db    <- reactive({ input$selected_db    })
    # selected_table <- reactive({ input$selected_table })
    
    return(list(
      selected_db    = selected_db
      #selected_table = selected_table
    ))
  })
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
      master = "sc://localhost:15002",
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
                                                   min = 0, max = 3, value = 1, step = 0.1),
                                       sliderInput("pval_cut", "p-value Threshold:", 
                                                   min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                                       sliderInput("pointSize", "Point Size:", 
                                                   min = 1, max = 5, value = 2, step = 0.5),
                                       sliderInput("ptAlpha", "Transparent:", 
                                                   min = 0.1, max = 1, value = 0.6, step = 0.1),
                                       sliderInput("labelSize", "Gene Label Size:", 
                                                   min = 1, max = 6, value = 3, step = 0.5),
                                       numericInput("topN", "Label N number of Genes (0 is no labeling):", 
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

#' @title ggvolcano_custom
#' @description ggvolcano_custom for RNAseq analysis with custom parameters
#' @param df DEG table
#' @param geneName gene symbol
#' @param pValCol which col for p-value
#' @param logFCCol which col for log fold-change
#' @param coef coefficient for limma 
#' @param pval_cut p-value for cut-off label
#' @param lfc_cut log fold-change cut-off label
#' @param useAdjP adjusted p-value
#' @param title plot title
#' @param ptAlpha transparent
#' @param labelSize gene symbol text size
#' @param geneCol gene symbol column
#' @param pointSize point plot size
#' @param topN show a number of genes symbol with high log fold-change
#' @return ggplot object
#' @export
ggvolcano_custom <- function(   df       ,
                                geneName  ,
                                pValCol = "PValue",
                                logFCCol = "logFC",
                                coef      = 2,
                                lfc_cut   = 1,
                                pval_cut  = 0.05,
                                useAdjP   = FALSE,
                                title     = "Volcano Plot",
                                topN      = 20,
                                geneCol   = NULL,
                                pointSize = 2,       # 加入點大小參數
                                ptAlpha   = 0.6,     # 加入點透明度參數
                                labelSize = 3       # 加入基因標籤字型大小參數
                            ){


                                        plotData <- data.frame(
                                            gene  = geneName,
                                            logFC = df[[logFCCol]],
                                            pval  = df[[pValCol]]
                                        )

                                        plotData$color <- "grey"
                                        plotData$color[plotData$logFC >=  lfc_cut & plotData$pval <= pval_cut] <- "red"
                                        plotData$color[plotData$logFC <= -lfc_cut & plotData$pval <= pval_cut] <- "blue"
                                        
                                        plotData$negLogP <- -log10(plotData$pval)
                                        
                                        p <- ggplot(plotData, aes(x = logFC, y = negLogP, color = color)) +
                                            geom_point(size = pointSize, alpha = ptAlpha) +
                                            scale_color_identity() +
                                            geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "black") +
                                            geom_hline(yintercept = -log10(pval_cut), linetype = "dashed", color = "black") +
                                            theme_classic() +
                                            labs(
                                            title = title,
                                            x     = expression(log[2]~"Fold Change"),
                                            y     = bquote(-log[10]~.(pValCol))
                                            )
                                        

                                        if (topN > 0) {
                                            topGenesDF <- head(plotData[order(plotData$pval), ], n = topN)
                                            p <- p + 
                                            ggrepel::geom_text_repel(
                                                data = topGenesDF,
                                                aes(label = gene),
                                                size = labelSize,
                                                max.overlaps = Inf
                                            )
                                        }
                                        
                                        return(p)
                            }


#' @title shiny app for RNAseq for public use
#' @description start the RNAseq shiny app with spark connection for shiny proxy at the seqslab console
#' @return start the UI and Server for analysis
#' @export
RNAseqShinyAppSpark <- function(){
    ui <- fluidPage(
        navbarPage("RNAseq App (Beta)",
            tabPanel("Gene Expression (RNA-seq)",
                tags$style(".shinybusy-overlay {opacity: 0.7; background-color: #7c7c7c;}"),
                     add_busy_spinner(
                         spin = "fading-circle",
                         position = "full-page",
                         timeout = 1000
                 ),

                layout_sidebar(
                    sidebar = sidebar(
                        h4("Spark DB Data"),
                        #actionButton("connect", "Connect"),
                        dbBrowserUI("dbBrowser1"),
                        actionButton("get_tbl", "Get"),
                        actionButton("mae_start", "Send")

                    ),
                     mainPanel(
                            tabsetPanel(
                                tabPanel("Gene Expression Profile",
                                    fluidRow(
                                        column(
                                            width = 12,
                                            #actionButton("run_wide_conversion", "Convert to Wide Table"),
                                            DT::dataTableOutput("wide_table_dt", width = "100%")
                                        )
                                    )
                                ),

                                tabPanel("Differential Expression Analysis",
                                                sidebarLayout(
                                                    sidebarPanel(
                                                        sliderInput("lfc_cut", "Fold Change Threshold (log2):", 
                                                                    min = 0, max = 3, value = 1, step = 0.1),
                                                        sliderInput("pval_cut", "p-value Threshold:", 
                                                                    min = 0.001, max = 0.1, value = 0.05, step = 0.001),
                                                        sliderInput("pointSize", "Point Size:", 
                                                                    min = 1, max = 5, value = 2, step = 0.5),
                                                        sliderInput("ptAlpha", "Transparent:", 
                                                                    min = 0.1, max = 1, value = 0.6, step = 0.1),
                                                        sliderInput("labelSize", "Gene Label Size:", 
                                                                    min = 1, max = 6, value = 3, step = 0.5),
                                                        numericInput("topN", "Label N number of Genes (0 is no labeling):", 
                                                                    value = 0, min = 0, max = 100),
                                                        checkboxInput("use_adjP", "Use Adjusted P-value?", value = FALSE),
                                                        actionButton("run_DEG", "Run DEG analysis")
                                                    ),

                                                    mainPanel(
                                                        tabsetPanel(
                                                            tabPanel("Volcano Plot interaction",
                                                                plotOutput("volcano_plot", height = "600px")
                                                            ),
                                                            tabPanel("DEG Table", 
                                                                DT::dataTableOutput('DEG_table',width="100%")
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
                                tabPanel("Gene-enriched Analysis", 
                                sidebarPanel(
                                    #sliderInput("top_gene_number", "Gene number", min = 3, max = 20000, value = 500, step = 1),
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
                                                    h4("T"),
                                                    tabsetPanel(
                                                        id = "sub1Upper",
                                                        tabPanel("MF", plotOutput("T_MF")),
                                                        tabPanel("BP", plotOutput("T_BP")),
                                                        tabPanel("CC", plotOutput("T_CC")),
                                                        tabPanel("KEGG", plotOutput("T_KEGG"))
                                                    )
                                                )
                                            ),
                                            fluidRow(
                                                column(12,
                                                    # Lower section tabs for Sub1
                                                    h4("NT"),
                                                    tabsetPanel(
                                                        id = "sub1Lower",
                                                        tabPanel("MF", plotOutput("NT_MF")),
                                                        tabPanel("BP", plotOutput("NT_BP")),
                                                        tabPanel("CC", plotOutput("NT_CC")),
                                                        tabPanel("KEGG", plotOutput("NT_KEGG"))
                                                        
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
            )
        )
    )
    
    server <- function(input, output, session) {
        sc <- reactiveVal(NULL)
        
        observe({
            master <- "sc://localhost:15002"
            method <- "spark_connect"
            version <- "3.5"
            connection <<- sparklyr::spark_connect(master = master, method = method, version = version )
            sc(connection)
            print("完成連線")
        })
        
        output$conn_status <- renderPrint({
            if (is.null(sc())) {
                "尚未建立 Spark 連線"
            } else {
                cat("Spark 連線資訊：\n")
                print(sc())
            }
        })
        
        results <- reactiveValues(
            db_info = NULL,
            normcount_data = NULL,
            exacttest_data = NULL,
            table_list = NULL,
            coldata =NULL
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
            prefix <- c("normcount", "exacttest", "coldata")
            tbls_with_prefix <- tbls[sapply(tbls, function(x) {
                any(sapply(prefix, function(p) grepl(paste0("^", p), x, ignore.case = TRUE)))
            })]
            print(tbls_with_prefix)
            results$table_list <- tbls_with_prefix
            
            # 取得 normcount 資料表 (假設取第一個符合條件的)
            normcount_tbls <- tbls_with_prefix[grepl("^normcount", tbls_with_prefix, ignore.case = TRUE)]
            if(length(normcount_tbls) > 0){
                query_normcount <- paste0("SELECT * FROM ", normcount_tbls[1])
                results$normcount_data <- DBI::dbGetQuery(sc(), query_normcount)
            }
            colnames(results$normcount_data)[1] <- "GeneSymbol"

            # 取得 exacttest 資料表 (同理)
            exacttest_tbls <- tbls_with_prefix[grepl("^exacttest", tbls_with_prefix, ignore.case = TRUE)]
            if(length(exacttest_tbls) > 0){
                query_exacttest <- paste0("SELECT * FROM ", exacttest_tbls[1])
                results$exacttest_data <- DBI::dbGetQuery(sc(), query_exacttest)
            }
            colnames(results$exacttest_data)[which(colnames(results$exacttest_data)=="genes")] <- "GeneSymbol"
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
        maeColData<- reactiveVal(NULL)
        
        output$wide_table_dt <- DT::renderDataTable({
            req(wide_data())
            print("傳送wide data至UI")
            DT::datatable(
                wide_data(),
                options = list(pageLength = 20, autoWidth = TRUE)
            )
        })


        observeEvent(input$mae_start,{
            req(results$coldata, results$normcount_data, results$exacttest_data)
            DEG_table(results$exacttest_data)
            wide_data(results$normcount_data)
            maeColData(results$coldata)
            print(head(DEG_table()))
            print(head(wide_data()))
            print(head(results$coldata))

            assay_data <- as.matrix(wide_data()[, -which(colnames(wide_data())=="GeneSymbol")])
            sample_info_table <- maeColData()
            print(dim(assay_data))
            print(dim(sample_info_table))
            rownames(sample_info_table) <- colnames(assay_data) #colData的rownames要和assay_data的conames一致




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


        reactive_volcano_plot <- eventReactive(input$run_DEG , {
            req(DEG_table())
            ggvolcano_custom(
                        df = DEG_table(),
                        geneName = DEG_table()$GeneSymbol,
                        pValCol = "PValue",
                        logFCCol = "logFC",
                        coef      = 2,
                        lfc_cut   = input$lfc_cut,
                        pval_cut  = input$pval_cut,
                        useAdjP   = FALSE,
                        title     = "Volcano Plot",
                        topN      = input$topN,
                        geneCol   = NULL,
                        pointSize = input$pointSize,       # 加入點大小參數
                        ptAlpha   = input$ptAlpha,     # 加入點透明度參數
                        labelSize = input$labelSize       # 加入基因標籤字型大小參數
            )
        })



        output$volcano_plot <- renderPlot({
            reactive_volcano_plot()
            
        })
        observeEvent(input$targetGeneID, {
                    req(settingMAE(), wide_data())
                    mae <- settingMAE()
                    geneList <- unlist(strsplit(input$geneList, ","))
                    sample_info <- colData(mae[["RNAseq"]])
                    groups_list <- rownames(sample_info) 
                    expr_profile <- as.data.frame(wide_data())
                    if("GeneSymbol" %in% colnames(expr_profile)){
                            rownames(expr_profile) <- expr_profile[,"GeneSymbol"]
                            expr_profile <- expr_profile[,-1]
                        print("wide table 的Gene欄位轉成row names")
                    }
                    
                    targetGeneExpr <- target_exprofile( 
                                                        geneList_ = geneList, 
                                                        groups_list_ = groups_list,
                                                        expr_profile_ = expr_profile
                    )
                    print(targetGeneExpr)
                    output$target_gene_table <- DT::renderDataTable({DT::datatable(targetGeneExpr)})
        })


        observeEvent(input$generate_go, {
                req(topGeneList(),downGeneList(),settingMAE())                
                mae <- settingMAE()
                sample_info <- colData(mae[["RNAseq"]])
                groups_list <- unique(sample_info[, "subCode"]) 
                group1_fc_gene_profile <- topGeneList()
                group2_fc_gene_profile <- downGeneList()
                for(n in seq_len(length(groups_list)) ){
                    col <- groups_list[n]
                    print(col)
                    print(n)
                    gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])
                    print(head(gene_list))
                    print("setting ok")
                    for( mode in c("CC", "BP", "MF")){
                        print(mode)
                        VAR <- paste0(col,"_", mode, "GO")
                        print(VAR)
                        assign(VAR, go_enrich_dotplot( 
                        
                                                        gene_list_= unique(gene_list), 
                                                        save_path_ = NULL,
                                                        save_filename_ = NULL, 
                                                        mode_ = mode, 
                                                        showCategory_ = 10), envir=.GlobalEnv
                        
                        )

                    }

                }


                output$T_MF <- renderPlot({T_MFGO})
                output$T_BP <- renderPlot({T_BPGO})
                output$T_CC <- renderPlot({T_CCGO})
                output$NT_MF <- renderPlot({NT_MFGO})
                output$NT_BP <- renderPlot({NT_BPGO})
                output$NT_CC <- renderPlot({NT_CCGO})
                
                
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
                    print(col)
                    col <- groups_list[n]
                    gene_list <- get(c("group1_fc_gene_profile", "group2_fc_gene_profile")[n])

                    print(col)
                    VAR <- paste0(col,"_","KEGG")

                    assign(VAR, kegg_enrich_dotplot ( 
                                                        gene_list_ = unique(gene_list), 
                                                        save_path_ = NULL,
                                                        save_filename_ = NULL, 
                                                        showCategory_ = 10
                                )
                    )
                }

                output$T_KEGG <- renderPlot({T_KEGG})
                output$NT_KEGG <- renderPlot({NT_KEGG})

        })

        topGeneList <- reactiveVal(NULL)
        downGeneList <- reactiveVal(NULL)

        observeEvent(input$generate_genelist, {
                    req(DEG_table())
                    DEG_table <- DEG_table()

                    DEG_table_filtered <- DEG_table[DEG_table$PValue<0.05, ]
                    sorted_DEG <- DEG_table_filtered[order(DEG_table_filtered[,"logFC"], decreasing=TRUE),]
                    gene_list_symbol <- sorted_DEG$GeneSymbol
                    print(dim(DEG_table_filtered))
                    print(gene_list_symbol)

                    #topGeneList(gene_list_symbol[1:input$top_gene_number, "SYMBOL"])
                    #downGeneList(gene_list_symbol[(nrow(gene_list_symbol)-input$top_gene_number):nrow(gene_list_symbol), "SYMBOL"])
                    topGeneList( DEG_table[DEG_table$PValue<0.05 & sign(DEG_table_filtered[,"logFC"]) == 1, "GeneSymbol"]) #選所有顯著差異基因
                    downGeneList( DEG_table[DEG_table$PValue<0.05 & sign(DEG_table_filtered[,"logFC"]) == -1, "GeneSymbol"]) 
                    print(head(topGeneList()))
                    print("topGeneList")
                    print(head(downGeneList()))
                    print("downGeneList")


        })

    }
    
    for_run <- shinyApp(ui = ui, server = server)
    runApp(for_run)
}

#' @title generate pseudo colData with random classification
#' @description generate random classification for colData
#' @param normcount normCount table to get colnames (sample id)
#' @param genecol to remove gene column
#' @return pseudo colData
#' @export

generate_colData_random <- function(normcount = normCount, genecol = "gene_id" ){
    rownames(normcount) <- normcount[, genecol]  # 假設第一列是基因名稱
    normcount <- normcount[, -which(colnames(normcount) == genecol)]  # 去除基因名稱列，確保只有數值矩陣
    classification <- sample(c("T","NT"), ncol(normcount), replace = TRUE)
    colData <- data.frame("mainCode" =colnames(normcount) , "subCode" = classification)
    rownames(colData) <- colnames(normcount)
    return(colData)
}
