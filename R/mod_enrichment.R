#' @title gseaeFCModuleUI
#' @description GSEA module
#' @param id module id for UI and Server
#' @return GSEA module UI
#' @export

gseaFCModuleUI <- function(id) {
  ns <- NS(id)
  tagList(
    # numericInput(ns("pvalCutoff"), "P-value Cutoff:",
    #             value = 0.05, min = 0, max = 1, step = 0.001),
    # actionButton(ns("runGSEA"), "Run GSEA"),
    DTOutput(ns("gseaTable")),
    plotOutput(ns("gseaPlot"))
  )
}


#' @title gseaeFCModuleServer
#' @description GSEA module
#' @param id  module id for UI and Server
#' @param DEG_table DEG table
#' @param direction direction of GSEA (up or down)
#' @return GSEA module Server
#' @export
#'
gseaFCModuleServer <- function(id, DEG_table, direction = c("up", "down"), enrichment_db = "org.Hs.eg.db") {
  direction <- match.arg(direction)

  moduleServer(id, function(input, output, session) {
    result_GSEA_FC <- reactiveVal(NULL)
    local_pvalCutoff <- 0.5
    kegg_organism <- if(identical(enrichment_db, "org.Mm.eg.db")) {
      "mmu"
    } else {
      "hsa"
    }
    if (direction == "up") {
      future_promise(
        {
          start_time <- Sys.time()
          local_deg <- isolate(DEG_table()) # total genes without filtering
          deg_subset <- local_deg
          conv <- bitr(deg_subset$GeneSymbol,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = enrichment_db,
            drop = FALSE
          )
          # deg_subset <- merge(deg_subset, conv, by.x = "GeneSymbol", by.y = "SYMBOL")
          # deg_subset <- deg_subset[!is.na(deg_subset$ENTREZID), ]
          # deg_subset <- deg_subset[!duplicated(deg_subset$ENTREZID), ]
          geneList <- deg_subset %>%
                dplyr::inner_join(conv, by = c("GeneSymbol" = "SYMBOL")) %>%
                dplyr::filter(!is.na(ENTREZID)) %>%                            
                dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
                dplyr::select(ENTREZID, logFC) %>%
                deframe() %>%
                sort(decreasing = TRUE)
          gsea_res <- gseKEGG(
            geneList = geneList,
            organism = kegg_organism,
            scoreType = "pos",
            minGSSize = 1,
            maxGSSize = 500,
            pvalueCutoff = local_pvalCutoff,
            verbose = TRUE
          )
          end_time <- Sys.time()
          list(
            r = gsea_res,
            start_time = start_time,
            end_time = end_time,
            elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
          )
        },
        seed = TRUE
      ) %...>% {
        result_GSEA_FC(.$r)
        cat("GSEA_FC (upregulated) finished. Elapsed:", .$elapsed, "seconds\n")
      }
    } else { # direction == "down"
      future_promise(
        {
          start_time <- Sys.time()
          local_deg <- isolate(DEG_table())
          deg_subset <- local_deg
          conv <- bitr(deg_subset$GeneSymbol,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = enrichment_db,
            drop = FALSE
          )
          #deg_subset <- merge(deg_subset, conv, by.x = "GeneSymbol", by.y = "SYMBOL")
          #deg_subset <- deg_subset[!is.na(deg_subset$ENTREZID), ]
          #deg_subset <- deg_subset[!duplicated(deg_subset$ENTREZID), ]
          geneList <- deg_subset %>%
                dplyr::inner_join(conv, by = c("GeneSymbol" = "SYMBOL")) %>%
                dplyr::filter(!is.na(ENTREZID)) %>%                            
                dplyr::distinct(ENTREZID, .keep_all = TRUE) %>%
                dplyr::select(ENTREZID, logFC) %>%
                deframe() %>%
                sort(decreasing = TRUE)
          gsea_res <- gseKEGG(
            geneList = geneList,
            organism = kegg_organism,
            scoreType = "pos",
            minGSSize = 1,
            maxGSSize = 500,
            pvalueCutoff = local_pvalCutoff,
            verbose = TRUE
          )
          end_time <- Sys.time()
          list(
            r = gsea_res,
            start_time = start_time,
            end_time = end_time,
            elapsed = as.numeric(difftime(end_time, start_time, units = "secs"))
          )
        },
        seed = TRUE
      ) %...>% {
        result_GSEA_FC(.$r)
        cat("GSEA_FC (downregulated) finished. Elapsed:", .$elapsed, "seconds\n")
      }
    }

    output$gseaTable <- renderDT({
      req(result_GSEA_FC())
      as.data.frame(result_GSEA_FC())
    })
    outputOptions(output, "gseaTable", suspendWhenHidden = FALSE)
    output$gseaPlot <- renderPlot({
      req(result_GSEA_FC())
      if (nrow(result_GSEA_FC()@result) > 0) {
        dotplot(result_GSEA_FC(), showCategory = 10)
      } else {
        plot.new()
        text(0.5, 0.5, "No significant pathway found.")
      }
    })
    outputOptions(output, "gseaPlot", suspendWhenHidden = FALSE)

    return(list(result = result_GSEA_FC))
  })
}
