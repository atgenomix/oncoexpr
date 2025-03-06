# RNAseqShinyAppSpark

RNAseqShinyAppSpark is an R Shiny application designed for RNA sequencing (RNA-seq) data analysis using Spark. The app provides interactive panels for visualizing gene expression profiles, performing differential expression analysis (including volcano plots), executing gene-enrichment analyses (GO and KEGG), and examining target gene expression.

## Overview

This application connects to a Spark cluster and retrieves data from a user-selected database. It then processes the data for various RNA-seq analyses. The UI is organized into multiple tabs, each addressing a specific analysis step, while the server function handles the Spark connection, data querying, and visualization.

## Features

- **Gene Expression Profile:** Browse analysis runs and display normalized count data from the database.
- **Differential Expression Analysis:** Configure thresholds (fold change, p-value, etc.) and generate volcano plots alongside a differential expression gene (DEG) table.
- **Gene-enriched Analysis:** Automatically generate GO (MF, BP, CC) and KEGG dot plots based on up- and downregulated genes.
- **Target Gene Expression:** Enter a list of target genes to retrieve and display their expression profiles.

## Installation and Dependencies

Before running the app, ensure the following R packages are installed:

- **shiny:** For building the interactive web application.
- **sparklyr:** To establish a connection to the Spark cluster.
- **DBI:** For database interactions (executing SQL commands).
- **DT:** To render interactive data tables.
- **ggplot2:** For generating custom volcano plots.
- **SummarizedExperiment & MultiAssayExperiment:** For handling expression matrices and sample metadata.
- Additional helper functions (e.g., `add_busy_spinner`, `go_enrich_dotplot`, `kegg_enrich_dotplot`, and `target_exprofile`) should be available or defined within your environment.

## Application Structure

### User Interface (UI)

The UI is defined using a `fluidPage` that contains a `navbarPage` with the following tabs:

1. **Gene Expression Profile:**
   - **Sidebar:** Displays analysis runs via a database browser (`dbBrowserUI`) and includes a button (`Confirm DataBase`) to load the selected database.
   - **Main Panel:** Shows the normalized count table using an interactive DT data table.

2. **Differential Expression Analysis:**
   - **Sidebar:** Contains multiple slider inputs for setting fold change thresholds, p-value thresholds, point size, transparency, and gene label sizes. There is also a numeric input to control how many genes are labeled, a checkbox for using adjusted p-values, and a button to run DEG analysis.
   - **Main Panel:** Consists of a tabset panel displaying a volcano plot (interactive) and a DEG table.

3. **Gene-enriched Analysis:**
   - **Sidebar:** Features an action button to trigger the enrichment analysis.
   - **Main Panel:** Organized into two sections (for upregulated and downregulated genes) with individual tab panels for Molecular Function (MF), Biological Process (BP), Cellular Component (CC), and KEGG pathway enrichment plots.

4. **Target Gene Expression:**
   - **Sidebar:** Includes a text input for entering a comma-separated list of target genes and a confirmation button.
   - **Main Panel:** Displays the target gene expression table as an interactive DT data table.

### Server Workflow

The server function orchestrates the backend processes:

- **Spark Connection:**  
  Upon initialization, the app establishes a connection to a Spark cluster using `sparklyr::spark_connect` (configured with a master URL, connection method, and Spark version). This connection is stored in a reactive value.

- **Database Interaction:**  
  Once connected, the app retrieves database information using `dbBrowserServer`. Upon confirming the selected database, it executes SQL queries to load tables with specific prefixes (e.g., `normcountsgene`, `exacttestgene`, `coldata`) into reactive values.

- **Data Preparation:**
  - The normalized count table (`normcount_data`) and DEG results table (`exacttest_data`) are fetched and processed.
  - A pseudo column data is generated via `generate_colData_random` to create a proper sample information table.
  - A `SummarizedExperiment` is created for the RNA-seq expression data and then encapsulated into a `MultiAssayExperiment` (MAE) object for downstream analysis.

- **Differential Expression Analysis:**
  - The app renders a volcano plot based on the DEG table using customizable thresholds provided by the user.
  - A reactive DEG table is displayed for user review.

- **Gene-Enrichment Analysis:**
  - When triggered, the app filters the DEG table by the set p-value and fold change thresholds.
  - It then splits genes into upregulated and downregulated lists.
  - For each gene list, GO enrichment (for MF, BP, and CC) and KEGG enrichment dot plots are generated and rendered in separate tabs.

- **Target Gene Expression:**
  - The entered list of target genes is processed, and the expression profile is extracted from the normalized count table.
  - A DT table displays the expression profile of the target genes.

## How to Run the Application

1. **Set Up the Environment:**
   - Ensure Spark is running and accessible at the specified master address
   - Install all required R packages.

2. **Execute the Function:**
   - Source the script containing the `RNAseqShinyAppSpark` function.
   - Run the function in R:
     ```r
     RNAseqShinyAppSpark()
     ```
   - This will launch the Shiny app in your default web browser.

3. **Interact with the App:**
   - Navigate through the tabs to perform different analyses.
   - Follow on-screen instructions (e.g., confirming the database, adjusting thresholds, or inputting target gene lists).

## Dependencies Summary

- **R Packages:** shiny, sparklyr, DBI, DT, ggplot2, SummarizedExperiment, MultiAssayExperiment.
- **Spark Connection:** Ensure the Spark master is running and accessible.
- **Custom Functions:** Additional helper functions like `add_busy_spinner`, `generate_colData_random`, `ggvolcano_custom`, `go_enrich_dotplot`, `kegg_enrich_dotplot`, and `target_exprofile` must be available in your environment.

## License