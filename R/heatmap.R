#BiocManager::install("InteractiveComplexHeatmap")
#BiocManager::install("originalHeatmapOutput")

library(airway)
library(DESeq2)
library(shiny)
library(MultiAssayExperiment)
"======================================"
data(airway)
se = airway
dds = DESeqDataSet(se, design = ~ dex)
keep = rowSums(counts(dds)) >= 10
dds = dds[keep, ]
 
dds$dex = relevel(dds$dex, ref = "untrt")

dds = DESeq(dds)
res = results(dds)
res = as.data.frame(res)

"======================================"
# 取得 normalized counts 作為後續分析的表現矩陣
# 取得 normalized counts 作為後續分析的表現矩陣
assay_data <- as.matrix(counts(dds, normalized = TRUE))
# 取得 sample 資訊，並確保 rownames 與 assay_data 的 colnames 一致
sample_info <- as.data.frame(colData(dds))
rownames(sample_info) <- colnames(assay_data)

# 若原本沒有 group 欄位，這邊可直接使用 dex 作為群組資訊
sample_info$group <- sample_info$dex

# 建立 SummarizedExperiment 物件
se_expression_matrix <- SummarizedExperiment(
  assays = list(normCount = assay_data),
  colData = sample_info
)

# 建立 MultiAssayExperiment (MAE) 物件
mae <- MultiAssayExperiment(
  experiments = list(RNAseq = se_expression_matrix),
  colData = sample_info
)

# 定義以 MAE 為資料來源的 heatmap 生成函數
make_heatmap_mae <- function(mae, top_n = 500) {
  # 從 MAE 物件中取得 RNAseq 的 SummarizedExperiment
  se <- experiments(mae)[["RNAseq"]]
  
  # 取得 normalized counts (假設 assay 名稱為 "normCount")
  m <- assay(se, "normCount")
  
  # 以基因變異數篩選前 top_n 變異最大的基因
  variances <- apply(m, 1, var)
  top_genes <- order(variances, decreasing = TRUE)[1:min(top_n, nrow(m))]
  m_sub <- m[top_genes, ]
  
  # 取得 sample 資訊（此處使用 group 欄位）
  sample_info <- as.data.frame(colData(se))
  
  # 建立 heatmap
  ht <- Heatmap(t(scale(t(m_sub))), name = "z-score",
                top_annotation = HeatmapAnnotation(
                  group = sample_info$group  # 這裡如果有其他欄位也可加入
                ),
                show_row_names = FALSE, show_column_names = FALSE,
                column_title = paste0(nrow(m_sub), " top variable genes"),
                show_row_dend = FALSE)
  ht <- draw(ht, merge_legend = TRUE)
  ht
}

# 定義 UI，僅顯示 heatmap 本身
ui <- fluidPage(
  originalHeatmapOutput("ht", height = 800, containment = TRUE)
)

# 定義 Server 邏輯
server <- function(input, output, session) {
  ht <- make_heatmap_mae(mae, top_n = 500)
  
  if (!is.null(ht)) {
    makeInteractiveComplexHeatmap(input, output, session, ht, "ht")
  } else {
    output$ht_heatmap <- renderPlot({
      grid::grid.newpage()
      grid::grid.text("No data available.")
    })
  }
}

# 執行 Shiny App
shinyApp(ui, server)
