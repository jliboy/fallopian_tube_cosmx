# Fallopian Tube CosMx Analysis

This repository contains an R Markdown workflow for processing CosMx Spatial Molecular
Imaging (SMI) data from the Fallopian Tube Project in the Howitt Lab at Stanford.
The workflow builds Seurat objects from CosMx output, performs spatial and cellular
quality control, and supports normalization, batch correction, decontamination,
resegmentation, cell-type annotation, DeSpotX processing, and differential expression.

## Repository layout

```text
.
├── main.rmd                 # End-to-end analysis workflow
└── src/
		├── spatial_qc.r          # CosMx loading and spatial/cell QC helpers
		├── spatial_uproc.R       # Normalization and downstream processing helpers
		└── spatial_dproc.r       # Differential expression helpers
```

## Requirements

- R (a current R 4.x installation is recommended)
- RStudio or another environment that can knit R Markdown documents
- Sufficient memory for the full CosMx object; `main.rmd` configures a 30 GB
	`future.globals.maxSize` limit
- CosMx flat files organized into one subdirectory per slide/TMA
- A Python environment with `anndata`, `pandas`, and `scipy` only if the DeSpotX
	conversion steps are used
- The DeSpotX command-line tool and its environment for running DeSpotX itself

The R workflow uses packages including `Seurat`, `SeuratObject`, `Matrix`, `data.table`,
`dplyr`, `ggplot2`, `ggrepel`, `ggthemes`, `patchwork`, `readxl`, `here`, `scPearsonPCA`,
`decontX`, `FastReseg`, `HieraType`, `reticulate`, `anndata`, `lme4`, `smiDE`, `dbscan`,
`pheatmap`, `UpSetR`, and `fs`. Some of these packages may need to be installed from
their upstream or GitHub repositories according to their own installation instructions.

## Configuration

Create a project-local `.Renviron` file and set paths for the data and output locations:

```text
input_path=/path/to/cosmx/flat-files
output_path=/path/to/analysis-output
fov_positions_path=/path/to/fov_positions_file.csv.gz
tma_map_path=/path/to/tma-map.xlsx
RETICULATE_PYTHON=/path/to/despotx/python
```

The checked-in `.Renviron` is ignored by Git and contains example machine-specific
paths. Update it for your environment rather than committing personal or shared data
locations.

The input directory should contain one immediate subdirectory per slide/TMA. Each
slide directory is expected to contain:

- a file matching `*_metadata_file.csv.gz`
- a file matching `*_exprMat_file.csv.gz`
- optionally, a file matching `*-polygons.csv.gz`

The TMA map must contain the sheet `FT- OTHER CA TMA Map`, as referenced in
`main.rmd`. Review the FOV exclusion list and QC thresholds in that file for each
experiment before running the complete workflow.

## Running the workflow

1. Open `fallopian_tube.Rproj` in RStudio.
2. Configure `.Renviron` and restart R so the environment variables are available.
3. Install and load the required R packages.
4. Open `main.rmd` and run the chunks sequentially, or knit the document after
	 confirming all experiment-specific settings.

The workflow performs these main steps:

1. Load and merge CosMx slide data into a Seurat object.
2. Generate FOV-level QC plots and a slide QC table.
3. Compute FOV integrity, signal-to-background ratio, and split-ratio metrics.
4. Flag low-quality cells and optionally remove flagged cells for downstream analysis.
5. Assign TMA cores and sample IDs, then exclude known mislabeled FOVs.
6. Run scPearsonPCA normalization, Leiden clustering, and optional batch correction.
7. Estimate contamination with decontX and optionally run FastReseg.
8. Generate HieraType marker summaries and prepare spatial data for DeSpotX.
9. Run differential expression using the helpers in `src/spatial_dproc.r`.

Most plots and tables are written to `output_path` with the current date in their
filenames. The default QC thresholds and experiment-specific values are defined near
the top of `main.rmd`.

## DeSpotX workflow

When analyzing multiple slides, run the tissue-offsetting and h5ad conversion chunks
in `main.rmd` before invoking DeSpotX. The command is run from a terminal, for example:

```bash
despotx --h5ad /path/to/input.h5ad \
	--cell-type-col cluster \
	--spatial-key spatial \
	--out /path/to/despotx-output/ \
	--device cuda
```

After DeSpotX finishes, set `h5ad_path` to the generated `.h5ad` file and run the
import chunk to add the result back to a Seurat object.

## Notes

- The workflow flags cells during QC before the explicit subsetting step; inspect the
	QC plots and thresholds before removing cells.
- Keep raw CosMx data and generated outputs outside the repository when possible.
- The scripts are analysis helpers rather than an installable R package. Source them
	from `main.rmd` with `here::here()` as shown in the workflow.
