# RNAseqShinyAppSpark

RNAseqShinyAppSpark is an R Shiny application designed for RNA sequencing (RNA-seq) data analysis using Spark. The app offers interactive panels for visualizing gene expression profiles, performing differential expression analysis (including volcano plots), conducting gene-enrichment analyses (GO and KEGG), and examining target gene expression.

## Overview

This application connects to a Spark cluster and retrieves data from a user-selected database. It processes the data for various RNA-seq analyses. The user interface is organized into multiple tabs, each dedicated to a specific analysis step, while the server function manages the Spark connection, data querying, and visualization.

## Features

- **Gene Expression Profile:** Browse analysis runs and view normalized count data.
- **Differential Expression Analysis:** Configure thresholds (fold change, p-value, etc.) and generate volcano plots alongside a differential expression gene (DEG) table.
- **Gene-Enriched Analysis:** Automatically produce GO (MF, BP, CC) and KEGG dot plots for up- and downregulated genes.
- **Target Gene Expression:** Input a list of target genes to retrieve and display their expression profiles.

## Installation and Dependencies

Before running the app, ensure you have installed the following R packages:

- **shiny:** For building the interactive web application.
- **sparklyr:** To establish a connection with the Spark cluster.
- **DBI:** For database interactions (executing SQL commands).
- **DT:** For rendering interactive data tables.
- **ggplot2:** For creating custom volcano plots.
- **SummarizedExperiment** and **MultiAssayExperiment:** For handling expression matrices and sample metadata.
- Additional helper functions (e.g., `add_busy_spinner`, `go_enrich_dotplot`, `kegg_enrich_dotplot`, and `target_exprofile`) should also be available in your environment.

## Application Structure

### User Interface (UI)

The UI is built using a `fluidPage` that contains a `navbarPage` with the following tabs:

1. **Gene Expression Profile:**
   - **Sidebar:** Displays analysis runs via a database browser (`dbBrowserUI`) and provides a selection menu for choosing among different analysis runs. This enables users to select the specific batch of data they wish to analyze.
   - **Main Panel:** Shows the normalized count table using an interactive DT data table corresponding to the selected analysis run.

2. **Differential Expression Analysis:**
   - **Sidebar:** Contains multiple slider inputs for setting fold change thresholds, p-value thresholds, point size, transparency, and gene label sizes. Additionally, there is a numeric input to control the number of genes labeled, a checkbox for using adjusted p-values, and a button to run the DEG analysis.
   - **Main Panel:** Includes a tabset panel that displays an interactive volcano plot and a DEG table.

3. **Gene-Enriched Analysis:**
   - **Sidebar:** Features an action button to trigger the enrichment analysis.
   - **Main Panel:** Organized into two sections—one for upregulated genes and one for downregulated genes—with individual tab panels for Molecular Function (MF), Biological Process (BP), Cellular Component (CC), and KEGG pathway enrichment plots.

4. **Target Gene Expression:**
   - **Sidebar:** Provides a text input for entering a comma-separated list of target genes along with a confirmation button.
   - **Main Panel:** Displays the target gene expression table as an interactive DT data table.

### Server Workflow

The server function orchestrates the backend processes:

- **Spark Connection:**  
  On initialization, the app establishes a connection to a Spark cluster using `sparklyr::spark_connect` (configured with the appropriate master URL, connection method, and Spark version). This connection is stored in a reactive value.

- **Database Interaction:**  
  Once connected, the app retrieves database information using `dbBrowserServer`. After the user selects the desired database, it executes SQL queries to load tables with specific prefixes (e.g., `normcountsgene`, `exacttestgene`, `coldata`) into reactive values.

- **Data Preparation:**
  - The normalized count table (`normcount_data`) and DEG results table (`exacttest_data`) are retrieved and processed.
  - A pseudo column data is generated via `generate_colData_random` to create an appropriate sample information table.
  - A `SummarizedExperiment` is created for the RNA-seq expression data and then encapsulated into a `MultiAssayExperiment` (MAE) object for downstream analysis.

- **Differential Expression Analysis:**
  - The app renders a volcano plot based on the DEG table using customizable thresholds specified by the user.
  - A reactive DEG table is provided for review.

- **Gene-Enriched Analysis:**
  - When triggered, the app filters the DEG table using the defined p-value and fold change thresholds.
  - It then separates genes into upregulated and downregulated lists.
  - For each gene list, the app generates GO enrichment (for MF, BP, and CC) and KEGG enrichment dot plots, which are displayed in separate tabs.

- **Target Gene Expression:**
  - The input list of target genes is processed, and their expression profiles are extracted from the normalized count table.
  - A DT table then displays the expression profiles of the selected target genes.

## How to Run the Application

1. **Set Up the Environment:**
   - Ensure that Spark is running and properly configured.
   - Install all required R packages.

2. **Execute the Function:**
   - Source the script containing the `RNAseqShinyAppSpark` function.
   - Run the function in R:
     ```r
     RNAseqShinyAppSpark()
     ```
   - The Shiny app will launch in your default web browser.

3. **Interact with the App:**
   - Navigate through the tabs to perform various analyses.
   - Follow the on-screen instructions (e.g., selecting an analysis run, adjusting thresholds, or entering target gene lists).

## Dependencies Summary

- **R Packages:** shiny, sparklyr, DBI, DT, ggplot2, SummarizedExperiment, MultiAssayExperiment.
- **Spark Connection:** Ensure that Spark is running and properly configured.
- **Custom Functions:** Additional helper functions such as `add_busy_spinner`, `generate_colData_random`, `ggvolcano_custom`, `go_enrich_dotplot`, `kegg_enrich_dotplot`, and `target_exprofile` must be available in your environment.

## License
