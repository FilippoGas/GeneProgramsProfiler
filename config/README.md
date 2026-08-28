# Configuration

This workflow is configured via `config/config.yaml`. All parameters are validated at startup against the schema defined in `workflow/schemas/config.schema.yaml`.

The test dataset (Natri et al., IPF vs Control) is configured in `.test/config/config.yaml` and used automatically when running CI or dry-runs from the repository root.

## Input data

### Seurat object (`.rds`)

A Seurat object with the following columns in its metadata:

| Metadata column | Description |
|---|---|
| Sample name column | One column identifying the biological sample each cell belongs to (set via `preprocess.annotate_and_save.sample_column`) |
| Condition column | One column with the case/control condition per sample (set via `preprocess.annotate_and_save.condition_column`) |
| Cell-cycle phase column | One column with the cell-cycle phase per cell (set via `preprocess.annotate_and_save.cell_cycle_phase_column`) |
| Cell-type annotation column | One column with the cell-type label per cell (set via `preprocess.annotate_and_save.celltype_annotation_colname`) |

### Cytopus cell-type dictionary (`.json`)

A JSON file mapping the cell-type labels present in your Seurat object to [cytopus](https://github.com/wallet-maker/cytopus) cell-type identifiers. Example:

```json
{
  "Macrophages": "mac",
  "Fibroblasts": "fib",
  "AT1": "at1",
  "AT2": "at2"
}
```

## Global settings

| Key | Type | Description |
|---|---|---|
| `scRNAseq` | string | Path to the input Seurat `.rds` file |
| `celltype_conversion_dictionary` | string | Path to the cytopus cell-type conversion JSON |
| `analysis_name` | string | Name of this analysis. All outputs are written to `results/<analysis_name>/` and logs to `logs/<analysis_name>/` |
| `case_condition` | string | Label for the case/disease condition in the condition column |
| `control_condition` | string | Label for the control condition in the condition column |
| `queues.cpu` | string | HPC queue name for CPU jobs (used by the cluster launcher) |

## Module 1: Preprocessing (`preprocess`)

### `preprocess.annotate_and_save`

Annotates the Seurat object with cytopus cell types, performs cell-cycle scoring, and saves the dataset as `.rds`, 10X Genomics `.mtx` format, and AnnData `.h5ad`.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `rstudio_memory` | integer | yes | Memory (MB) for loading the Seurat object in R |
| `celltype_annotation_colname` | string | yes | Name of the cell-type annotation column in the Seurat metadata |
| `sample_column` | string | yes | Column of sample name in the Seurat metadata |
| `condition_column` | string | yes | Column of condition name in the Seurat metadata |
| `cell_cycle_phase_column` | string | yes | Column of cell-cycle phase in the Seurat metadata |
| `time` | string | no | Job walltime (`HH:MM:SS`), required for HPC schedulers |

## Module 2: Spectra gene program discovery (`spectra`)

Implements [Spectra](https://github.com/dpeerlab/spectra) for identifying cell-type-specific gene programs using expression data and cytopus gene sets.

### `spectra.prepare_cytopus_list`

Downloads cytopus gene sets for the cell types present in the dataset.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `global_celltype` | string | yes | Cell type to use as global cell type in the cytopus list |
| `time` | string | no | Job walltime |

### `spectra.run_spectra`

Runs Spectra to quantify gene program activation in single cells.

| Key | Type | Required | Description |
|---|---|---|---|
| `lambda` | float | yes | Weighs the relative contribution of cytopus list vs expression loss functions (range: 0.0001--0.5) |
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `time` | string | no | Job walltime |

### `spectra.rename_programs`

Labels unlabeled factors via ORA enrichment of marker genes against cytopus gene sets.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `time` | string | no | Job walltime |

### `spectra.spectra_WMW` / `spectra.spectra_LMM`

Differential activation testing of spectra gene programs between conditions. WMW uses Wilcoxon-Mann-Whitney U-test; LMM uses Linear Mixed Models to correct for cell-cycle phase.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `active_cell_thresh` | float | yes | Activation threshold to consider a program active in a cell |
| `time` | string | no | Job walltime |

### `spectra.spectra_WMW_plots` / `spectra.spectra_LMM_plots`

Plots from the differential activation analysis (volcano plots, heatmaps).

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `effect_size_thresh` | float | yes (WMW) | Effect size threshold (rank-biserial correlation) |
| `log2FC_thresh` | float | yes (LMM) | Log2 fold-change threshold |
| `FDR_thresh` | float | yes | False Discovery Rate threshold |
| `time` | string | no | Job walltime |

## Module 3: Differential expression analysis (`DE_analysis`)

### `DE_analysis.run_DE_analysis`

Runs differential expression analysis between the case and control conditions.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `logFC` | float | yes | Log fold-change threshold to consider a gene differentially expressed |
| `FDR` | float | yes | FDR threshold to consider a gene differentially expressed |
| `time` | string | no | Job walltime |

### `DE_analysis.DEA_plots`

Generates diagnostic plots (p-value overlap, correlation plots).

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `time` | string | no | Job walltime |

## Module 4: Functional enrichment (`functional_enrichment`)

### `functional_enrichment.run_gsea`

Runs Gene Set Enrichment Analysis using [fgsea](https://github.com/ctlab/fgsea) on the differential expression results.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `padj_thresh` | float | yes | Adjusted p-value threshold for significance |
| `time` | string | no | Job walltime |

### `functional_enrichment.run_ora`

Runs Over-Representation Analysis on the differential expression results.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `padj_thresh` | float | yes | Adjusted p-value threshold for significance |
| `time` | string | no | Job walltime |

## Module 5: cNMF gene program discovery (`cNMF`)

Implements [consensus NMF](https://github.com/dylkot/cNMF) as an alternative gene program discovery method, with automatic k-selection and consensus clustering.

### `cNMF.cNMF_prepare`

Normalizes the count matrix and prepares the factorization step. Defines the range of k values to evaluate.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `max_nmf_iter` | integer | yes | Maximum NMF optimization iterations per replicate |
| `k_min` | integer | yes | Minimum value of k to try |
| `k_max` | integer | yes | Maximum value of k to try |
| `k_step` | integer | yes | Step size for k sweep |
| `n_iter` | integer | yes | Number of factorization iterations for each k |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_factorize_worker`

Runs a single factorization worker. The `cores` parameter here sets the **number of parallel workers** (each worker runs with `threads: 1`), not CPUs per job.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) per worker |
| `cores` | integer | yes | Number of parallel worker jobs to spawn |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_combine`

Combines factorization results across all k values.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_k_selection_plot`

Generates a plot estimating the trade-off between higher k, stability, and error. Used for diagnostics; the actual k is selected automatically.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `time` | string | no | Job walltime |

### `cNMF.extract_best_k`

Selects the k value with the best stability-error tradeoff from the k-selection statistics.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_consensus`

Generates program usage tables for the selected k. Filters out unstable outlier programs before consensus clustering.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `local_density_threshold` | float | yes | Maximum distance threshold to nearest neighbors for filtering unstable programs |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_rename_programs`

Runs ORA on cNMF program markers to assign biological labels.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_WMW` / `cNMF.cNMF_LMM`

Differential activation testing of cNMF gene programs. Same statistical approaches as the spectra equivalents.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `active_cell_thresh` | float | yes | Activation threshold to consider a program active in a cell |
| `time` | string | no | Job walltime |

### `cNMF.cNMF_WMW_plots` / `cNMF.cNMF_LMM_plots`

Plots from cNMF differential activation analysis.

| Key | Type | Required | Description |
|---|---|---|---|
| `cores` | integer | yes | Number of threads |
| `mem_mb` | integer | yes | Memory (MB) |
| `effect_size_thresh` | float | yes (WMW) | Effect size threshold (rank-biserial correlation) |
| `log2FC_thresh` | float | yes (LMM) | Log2 fold-change threshold |
| `FDR_thresh` | float | yes | False Discovery Rate threshold |
| `time` | string | no | Job walltime |

## Module 6: Collect results (`collect_results`)

### `collect_results.make_comp_table`

Combines results from spectra, cNMF, and functional enrichment into a single comparative table.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `padj_thresh` | float | yes | Adjusted p-value threshold for enrichments to include in the table |
| `time` | string | no | Job walltime |

### `collect_results.comp_table_plots`

Generates UpSet plots showing concordance of detected deregulations across methods.

| Key | Type | Required | Description |
|---|---|---|---|
| `mem_mb` | integer | yes | Memory (MB) |
| `cores` | integer | yes | Number of threads |
| `FDR_thresh` | float | yes | FDR threshold for WMW and LMM results in the plots |
| `effect_size_thresh` | float | yes | Effect size threshold for WMW results in the plots |
| `log2FC_thresh` | float | yes | Log2FC threshold for LMM results in the plots |
| `time` | string | no | Job walltime |

## Resource notes

- `time` fields follow the pattern `HH:MM:SS` and are only needed when running on HPC clusters with job scheduling (PBS/SLURM). They can be left empty or omitted for local execution.
- `mem_mb` values in `.test/config/config.yaml` are tuned for the test dataset; production values may be significantly higher, especially for `cNMF.cNMF_factorize_worker` and `spectra.run_spectra`.
