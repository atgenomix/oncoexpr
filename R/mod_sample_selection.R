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
    observe({
      req(data()) 
      updateSelectInput(session, "mainCode_col", choices = names(data()))
      updateSelectInput(session, "subCode_col", choices = names(data()))
    })
    
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

    selectizeInput(ns("filter_columns"),
                   "column to filter",
                   choices = NULL,
                   multiple = TRUE),

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
    
    rv_filter <- reactiveValues()
    
    colData_df <- reactive({
      req(long_df())
      long_df() %>% select(-value)
    })

    observe({
      req(colData_df())
      updateSelectizeInput(session, "filter_columns",
                           choices = names(colData_df()),
                           server = TRUE)
    })

    observe({
      req(input$filter_columns)
      stored <- names(rv_filter)
      for (col in stored) {
        if (!(col %in% input$filter_columns)) {
          rv_filter[[col]] <- NULL
        }
      }
    })

    output$value_select_ui <- renderUI({
      req(colData_df())
      selected_cols <- input$filter_columns
      if (is.null(selected_cols) || length(selected_cols) == 0) {
        return(NULL)
      }
      ui_list <- lapply(selected_cols, function(col) {
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

    filtered_df <- reactive({
      req(long_df())
      df <- long_df()
      selected_cols <- input$filter_columns
      if (!is.null(selected_cols) && length(selected_cols) > 0) {
        for (col in selected_cols) {
          filter_val <- rv_filter[[col]]
          if (!is.null(filter_val) && !("All" %in% filter_val)) {
            df <- df %>% filter(.data[[col]] %in% filter_val)
          }
        }
      }
      df
    })
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
  names(selected_df)[names(selected_df) == mainCode_] <- "mainCode"
  names(selected_df)[names(selected_df) == subCode_] <- "subCode"
  return(selected_df)
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
    preview_data_filtered <- reactiveVal(NULL) 
    
    print(reactive({input$mainCode_input}))
    print(reactive({input$subCode_input}))

    mainCodes <- reactive({
      strsplit(input$mainCode_input, ",")[[1]] |> trimws()
    })
    subCodes <- reactive({
      strsplit(input$subCode_input, ",")[[1]] |> trimws()
    })

    observeEvent(input$preview_btn, {


      DBI::dbExecute(sc, paste0("USE spark_catalog.", selected_db()))
      result <- query_spark_rnaseqdata(
        sc         = sc,
        use_table  = selected_table(),
        mainCode   = mainCodes(),
        subCode    = subCodes(),
        normMethod = input$normMethod_input
      )
      print(head(result$data))
      df_renamed <- result$data |> dplyr::rename(
        mainCode = main_code,
        subCode  = sub_code
      )
      preview_data_filtered(df_renamed)
      
      print(dim(preview_data_filtered()))

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
      subColData = set_colData,
      filteredTable  = preview_data_filtered

    ))
  })
}



#' @title DE geneSelector UI
#' @description click gene from DEG table for plot
#' @param id module id for UI and Server
#' @return click gene to render on the volcano plot
#' @export

mod_geneSelector_ui <- function(id) {
  ns <- NS(id)
  tagList(
    h4("Gene Selector"),
    DTOutput(ns("deg_table"))
  )
}

#' @title DE geneSelector server
#' @description click gene from DEG table for plot
#' @param id module id for UI and Server
#' @param deg_table original DEG table
#' @param geneList DEG list
#' @return click gene to render on the volcano plot
#' @export
#' 
mod_geneSelector_server <- function(id, deg_table, geneList) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    filtered_deg <- reactive({
      req(deg_table, geneList)
      deg_table[deg_table$GeneSymbol %in% geneList, ]
    })
    output$deg_table <- renderDT({
      req(filtered_deg())
      datatable(
        filtered_deg(),
        selection = "single",
        options = list(pageLength = 10, scrollX = TRUE)
      )
    })
    
    selected_gene <- reactive({
      sel <- input$deg_table_rows_selected
      if (length(sel) > 0) {
        filtered_deg()[sel, "GeneSymbol"]
      } else {
        NULL
      }
    })
    
    return(selected_gene)
  })
}




