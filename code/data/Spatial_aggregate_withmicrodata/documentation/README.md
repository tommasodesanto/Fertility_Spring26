# How to Generate the County-Level Analysis Report

This document provides step-by-step instructions to reproduce the PDF summary report (`county_level_summary_v3.pdf`) from the raw analysis outputs.

## Prerequisites

1.  **R Environment**: A working R installation with the required packages (the script will attempt to install them).
2.  **LaTeX Distribution**: A complete LaTeX distribution (e.g., TeX Live, MiKTeX, MacTeX) is required to compile the `.tex` report into a PDF.

## Step-by-Step Guide

### Step 1: Run the County-Level Analysis Script

First, you must generate all the necessary tables and figures. The primary script handles this.

-   **Script Location**: `Spatial_aggregate_withmicrodata/master_county.R`
-   **Action**: Run this R script from the project's root directory. It will generate all the `.tex` tables and `.png` figures inside the `Spatial_aggregate_withmicrodata/analysis_output_county_restored/` directory.

### Step 2: Navigate to the Documentation Directory

All report compilation is done from within this `documentation` directory.

```bash
cd Spatial_aggregate_withmicrodata/documentation
```

### Step 3: IMPORTANT - Manually Edit the Panel Regressions Table

This is the most critical manual step. The R script that generates the aggregate panel regression table (`agg_panel_regressions.tex`) wraps it in a full `\begin{table}` environment. However, our main report (`county_level_summary_v3.tex`) needs to include this table inside a special `\begin{sidewaystable}` environment to display it in landscape mode.

Nesting a `table` inside another `table` (or `sidewaystable`) causes a `! LaTeX Error: Not in outer par mode`.

**To fix this, you must edit the file:**

`Spatial_aggregate_withmicrodata/analysis_output_county_restored/tables/agg_panel_regressions.tex`

**You need to remove the following lines from the beginning and end of the file:**

1.  **At the beginning of the file, delete these four lines:**
    ```latex
    \begin{table}[htbp]
    \centering
    \caption{Aggregate County Panel Regression Models of Fertility and Housing Costs}
    \label{tab:county_panel_regressions}
    ```

2.  **At the end of the file, delete these two lines:**
    ```latex
    \par\raggedright\footnotesize Notes: All models are weighted by county population. Controlled models include the following covariates: Pct. Graduate+, Log(Pop. Density), Median Age, Pct. Pop. 25-34
    \end{table}
    ```
By doing this, you are left with only the core `\begingroup ... \endgroup` content, which can be safely included in the main report.

### Step 4: Compile the LaTeX Report

Once the table file is corrected, you can compile the final PDF. You need to run `pdflatex` twice to ensure all cross-references (e.g., "Table 1," "Figure 2") are correctly generated.

From within the `Spatial_aggregate_withmicrodata/documentation` directory, run the following commands:

```bash
pdflatex -interaction=nonstopmode county_level_summary_v3.tex
pdflatex -interaction=nonstopmode county_level_summary_v3.tex
```

### Step 5: Verify the Output

After the commands complete successfully, you will find the final report in the current directory:

-   `county_level_summary_v3.pdf`

You can now open this file to see the finished report. 