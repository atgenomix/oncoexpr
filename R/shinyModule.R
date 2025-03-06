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
      label   = "",
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
     selected_db    <- reactive({ input$selected_db })
    # selected_table <- reactive({ input$selected_table })
    
    return(list(
      selected_db    = selected_db
      #selected_table = selected_table
    ))
  })
}

