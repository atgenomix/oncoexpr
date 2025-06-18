# RNAseqShinyAppSpark

An interactive R Shiny application for **RNA-seq data analysis** built with **Spark** backend support. It offers high-performance DEG analysis, real-time visualization (volcano plots, heatmaps), PCA, and functional enrichment (GO, KEGG, GSEA) in a modular, reactive framework.

This platform is designed to integrate seamlessly with the [**Atgenomix Seqslab platform**](https://docs.atgenomix.com/home.html). Spark connections and delta table queries are executed through the **Seqslab DataHub**, making the app ideal for large-scale biomedical applications in cloud-based environments.

---

## 🔍 Project Overview

**RNAseqShinyAppSpark** allows users to explore and analyze large-scale RNA-seq expression datasets directly from Spark-based delta tables. The app is optimized for both public demo usage and internal data upload scenarios, with full support for asynchronous computation, gene filtering, and biological insight discovery through visual analytics.

---

## ✨ Features

- **Interactive DEG Analysis** – Built-in pipelines for differential expression using `edgeR` and `limma`, with adjustable fold change and p-value thresholds.
- **Dynamic Visualization** – Real-time rendering of volcano plots, violin plots, scatter plots, and heatmaps.
- **Gene Set Enrichment** – GO (BP, CC, MF), KEGG, and GSEA visualization using `clusterProfiler` and `enrichplot`.
- **PCA Explorer** – Perform and visualize Principal Component Analysis with optional clustering.
- **Modular UI** – Shiny modules for plot controls, filtering, and enrichment make extension easy.
- **Asynchronous Execution** – Responsive UI powered by `future_promise` and task progress popups.

> The app manages all Spark and data access processes behind the scenes, so users can focus entirely on analysis and visualization.


---

## ⚙️ Installation

### Requirements

- **R >= 4.1**
- **Apache Spark >= 3.0** with delta table support
- The following R packages are either required or suggested based on the `DESCRIPTION` file:

#### Imports

```r
install.packages(c(
  "shiny", "sparklyr", "promises", "future", "DT", "ggplot2", "shinyjs",
  "shinycssloaders", "InteractiveComplexHeatmap", "clusterProfiler",
  "enrichplot", "org.Hs.eg.db", "org.Mm.eg.db", "MultiAssayExperiment", "SummarizedExperiment",
  "factoextra", "ggiraph", "shinydashboard", "pcaPP", "dplyr", "tidyr", "viridis",
  "reshape2", "stringr", "readr", "readxl", "tibble", "RColorBrewer", "pheatmap",
  "ggrepel", "bslib"
))
```

#### Suggested (Optional)

```r
install.packages(c(
  "airway", "edgeR", "limma", "DESeq2", "ComplexHeatmap"
))
```

---

## 🚀 Usage Guide

### Launch the App

```r
RNAseqShinyAppSpark()
```

The application will automatically connect to the configured Spark cluster and initialize required datasets via Seqslab DataHub.

### Main Workflow

1. Select Spark database (`*_cus_username`).
2. Automatically loads:
   - `normcounts_*` table (normalized expression)
   - `exacttest_*` table (DEG results)
   - `coldata_*` table (sample metadata)
3. Visualize:
   - DEG tables and download CSV
   - Volcano + violin + scatter plots
   - Heatmaps by gene list
   - GO/KEGG enrichment
   - GSEA analysis (up/down-regulated)

---

## 📁 Directory Structure

```bash
.
├── RNAseqShinyApp.R              # Main Shiny app (UI + server)
├── mod_enrichment.R              # GSEA module (GO/KEGG support)
├── mod_volcanoplot.R            # Volcano, scatter, violin plot module
├── mod_sample_selection.R       # Sample filtering and gene selector module
├── mod_progress_popup.R         # Popup progress UI module
├── mod_spark.R                  # Spark connection + DB browser module
├── plot_enrichment.R            # GO/KEGG enrichment plotting
├── plot_volcano.R               # Static and interactive volcano plotting utils
├── plot_heatmap.R               # Heatmap generation from MAE
├── plot_pca.R                   # PCA plot with clustering
├── spark_query.R                # Async Spark query and pattern matching
├── utils.R                      # GTF conversion, expr comparison, helper functions
```

---

## 🤝 Contributing

We welcome contributions via pull requests or issues. To contribute:

1. Fork this repo
2. Create a new branch (`feature/your-feature`)
3. Test your changes
4. Submit a PR

---

## 📄 License

This project is licensed under the **Apache License 2.0**.

Copyright © 2025 Charles Chuang, atgenomix.

See the [LICENSE](./LICENSE) file for details.

---

## 📬 Contact

Please create a GitHub Issue for bug reports, questions, or feature requests.

