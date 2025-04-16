
#' @title progressPopupUI
#' @description progress popup UI
#' @param id  module id for UI and Server
#' @return progress popup UI
#' @export

progressPopupUI <- function(id) {
  ns <- NS(id)
  
  absolutePanel(
    id = ns("popupPanel"),
    style = "
      position: fixed !important;
      top: 50%;
      left: 50%;
      transform: translate(-50%, -50%);
      
      z-index: 9999;
      background-color: #FFFFFF;
      border: 1px solid #CCC;
      padding: 10px;
      display: none; /* 一開始隱藏 */
      box-shadow: 0 0 5px rgba(0,0,0,0.3);
    ",
    width = "300px",
    
    tags$strong("Just a moment while we load the datasets."),
    br(), br(),
    
    tags$div(
      id = ns("progressBarOuter"),
      class = "progress progress-striped active",
      style = "height: 25px;",
      
      tags$div(
        id = ns("progressBarInner"),
        class = "progress-bar",
        style = "width: 0%; color: black;",
        "0%"
      )
    )
  )
}



#' @title progressPopupServer
#' @description progress popup Server
#' @param id module id for UI and Server
#' @return progress popup Server
#' @export
progressPopupServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    tasksTotal <- reactiveVal(0)
    tasksDone  <- reactiveVal(0)
    
    addPromise <- function(prom, label = NULL) {
      tasksTotal(tasksTotal() + 1)
      prom %...>% (function(res) {
        tasksDone(tasksDone() + 1)
        res
      }) %...!% (function(e) {
        tasksDone(tasksDone() + 1)
        stop(e)
      })
    }
    
    observe({
      total <- tasksTotal()
      done  <- tasksDone()
      
      if (total == 0) {
        runjs(sprintf("$('#%s').hide();", session$ns("popupPanel")))
        return()
      }
      
      runjs(sprintf("$('#%s').show();", session$ns("popupPanel")))
      
      pct <- if (total > 0) round(done / total * 100) else 0
      if (total > 0 && done == 0) {
        pct <- 5
      }
      inner_id <- session$ns("progressBarInner")
      runjs(sprintf("
        $('#%s').css('width', '%d%%');
        $('#%s').text('%d%%');
        if (%d === 0) {
          $('#%s').css('color', 'black');
        } else {
          $('#%s').css('color', 'white');
        }
      ",
      inner_id, pct,
      inner_id, pct,
      pct,
      inner_id,
      inner_id
      ))
      
      if (done >= total) {
        shinyjs::delay(1000, {
          runjs(sprintf("$('#%s').hide();", session$ns("popupPanel")))
          tasksTotal(0)
          tasksDone(0)
        })
      }
    })
    
    list(addPromise = addPromise)
  })
}