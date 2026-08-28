# Snakemake workflow: GeneProgramsProfiler

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/FilippoGas/GeneProgramsProfiler/actions/workflows/main.yaml/badge.svg?branch=main)](https://github.com/FilippoGas/GeneProgramsProfiler/actions/workflows/main.yaml?query=branch%3Amain)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![workflow catalog](https://img.shields.io/badge/Snakemake%20workflow%20catalog-darkgreen)](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/FilippoGas/GeneProgramsProfiler)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A Snakemake workflow for characterizing cell type-specific gene program deregulations from scRNA-seq datasets. The workflow takes a Seurat `.rds` object and a [cytopus](https://github.com/wallet-maker/cytopus) cell-type conversion dictionary as inputs, discovers gene programs using two complementary methods ([Spectra](https://github.com/dpeerlab/spectra) and [cNMF](https://github.com/dylkot/cNMF)), tests each program for differential activation between conditions, performs functional enrichment on differentially expressed genes, and produces a comparative cross-method analysis.

## Workflow overview

The pipeline is organized into six modules that run in dependency order:

**1. Preprocessing** — A user-provided cell-type conversion dictionary maps dataset-specific cell-type labels to a standardized [cytopus](https://github.com/wallet-maker/cytopus) ontology. The annotated dataset is saved in three formats (AnnData `.h5ad`, Seurat `.rds`, 10X MTX) to enable cross-tool compatibility with downstream modules.

**2. Gene program discovery — Spectra** — [Spectra](https://github.com/dpeerlab/spectra) ([Peer et al., 2023](https://doi.org/10.1038/s41587-023-01940-3)) decomposes gene expression into interpretable gene programs by combining user-provided gene sets (cytopus cellular identities) with data-driven discovery. It explicitly models cell type and represents input gene sets as a gene–gene knowledge graph, guiding factorization toward biologically meaningful programs while adaptively allocating novel factors for unexplained variation. Programs already labeled by the prior gene sets are retained; remaining unlabeled programs are annotated via ORA enrichment of their marker genes against cytopus gene sets. Programs that remain unlabeled after ORA are discarded as likely technical artifacts.

**3. Gene program discovery — cNMF** — [cNMF](https://github.com/dylkot/cNMF) ([Kotliar et al., 2019](https://doi.org/10.7554/eLife.43803)) performs unsupervised consensus non-negative matrix factorization across a range of K values. Unlike Spectra, cNMF does not use prior gene-set information, making it an independent validation axis: programs detected by both methods are more likely to represent true biology, while method-specific programs may reflect complementary aspects of cellular heterogeneity. The optimal K is selected automatically from the stability–error tradeoff, and programs are labeled using the same cytopus ORA procedure.

**4. Differential expression** — Single-cell and pseudobulk differential expression analyses are performed on the Seurat object, identifying genes that are differentially expressed between case and control conditions. The intersection of both methods is used as a robust DEG list for downstream functional enrichment.

**5. Differential activation — WMW + LMM** — Each program method's cell-level activation scores are tested for case/control differences using two complementary statistical frameworks, both of which account for sample-level variation as a covariate: the Wilcoxon–Mann–Whitney test (non-parametric, with cell scores aggregated to sample-level means) and linear mixed models (`lmer(program ~ diagnosis + phase + (1|sample_name))`, which additionally corrects for cell-cycle phase as a fixed effect and treats sample as a random intercept). Reporting both methods provides transparency and captures different aspects of differential activation.

**6. Functional enrichment — GSEA + ORA** — Gene Set Enrichment Analysis and Over-Representation Analysis test the intersecting DEGs against cytopus gene sets, providing pathway-level context for the detected gene programs.

**7. Comparative analysis** — Results from all methods (Spectra WMW/LMM, cNMF WMW/LMM, GSEA, ORA) are merged into a single comparative table and visualized as UpSet plots showing concordance across methods, highlighting which deregulations are robustly detected.

See the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog/docs/workflows/FilippoGas/GeneProgramsProfiler) for a visual workflow diagram and full parameter table.

## Prerequisites

- **Python ≥ 3.12**
- **Snakemake ≥ 8.0**
- **Conda** with strict channel priorities (`conda-forge` + `bioconda`)
- **`snakemake-executor-plugin-cluster-generic`** (only required for cluster execution — see [Running on a cluster](#running-on-a-cluster))

## Quick start

```bash
# Clone the repository
git clone git@github.com:FilippoGas/GeneProgramsProfiler.git
cd GeneProgramsProfiler

# Edit the config file with your dataset paths and parameters
nano config/config.yaml

# Dry run to verify the DAG
snakemake --snakefile workflow/Snakefile --configfile config/config.yaml --dry-run

# Run the workflow
snakemake --snakefile workflow/Snakefile --sdm conda --cores <N>
```

All output paths derive from `config["analysis_name"]` and are written to `results/<analysis_name>/`. Logs mirror the same structure under `logs/<analysis_name>/`.

## Test dataset

A built-in test dataset ([Natri et al.](https://doi.org/10.1016/j.cell.2023.03.011), IPF vs Control) is committed in `.test/data/`. The Snakefile hardcodes the test config at `.test/config/config.yaml`, so you can run it directly from the repository root:

```bash
snakemake --snakefile workflow/Snakefile --sdm conda --cores 2
```

To run your own analysis, create a new config file in `config/` (one per analysis is recommended) and point the workflow to it:

```bash
cp config/config.yaml config/my_analysis.yaml
nano config/my_analysis.yaml

snakemake --snakefile workflow/Snakefile --configfile config/my_analysis.yaml \
  --sdm conda --cores <N>
```

## Running on a cluster

Every rule carries `mem_mb`, `time`, and `queue` resources that are forwarded to the scheduler. Each launcher script relies on the parameter `--executor cluster-generic`, which requires the **`snakemake-executor-plugin-cluster-generic`** plugin to be installed **in the environment that runs snakemake**, regardless of whether you use PBS Pro or Slurm. Install it with either:

```bash
# conda (recommended, matches the rest of the setup)
conda install -c conda-forge -c bioconda snakemake-executor-plugin-cluster-generic

# or pip
pip install snakemake-executor-plugin-cluster-generic
```

This step is **not needed for local runs** (`snakemake --snakefile workflow/Snakefile --sdm conda --cores <N>`); it is only required when using a cluster launcher.

### PBS Pro

A PBS Pro launcher script is provided at `workflow/scripts/launchers/PBS.sh`. The launcher uses a custom status-check script at `workflow/scripts/cluster_status/PBS_status.py`.

```bash
bash workflow/scripts/launchers/PBS.sh
```
**Remember** to either add `--configfile path/to/your/config.yaml` to the `snakemake` command in the PBS launcher script, or modify `workflow/Snakefile` to point to your config file.

### Slurm

A Slurm launcher script is provided at `workflow/scripts/launchers/slurm.sh`. The `time` and `queue` resources map to Slurm's `--time` and `--partition` options. The launcher uses a custom status-check script at `workflow/scripts/cluster_status/slurm_status.py`, which queries `sacct`/`squeue` for the job state.

```bash
bash workflow/scripts/launchers/slurm.sh
```
**Remember** to either add `--configfile path/to/your/config.yaml` to the `snakemake` command in the Slurm launcher script, or modify `workflow/Snakefile` to point to your config file.

### Other schedulers

The `snakemake-executor-plugin-cluster-generic` plugin is scheduler-agnostic and will typically work with most schedulers out of the box. If your scheduler is not PBS Pro or Slurm, you will need to write your own launcher script (defining the `--cluster-generic-submit-cmd`) and your own cluster status script (defining the `--cluster-generic-status-cmd`, which must print `running`, `success`, or `failed` for a given job id, see `workflow/scripts/cluster_status/PBS_status.py` for reference). Download an additional executor plugin only if your scheduler requires a dedicated one.

## Configuration

See [`config/README.md`](config/README.md) for a full reference of all configuration parameters, including module-specific settings for preprocessing, Spectra, cNMF, differential expression, functional enrichment, and result collection. Config values are validated at startup against `workflow/schemas/config.schema.yaml`.

## Output structure

All outputs are written under `results/<analysis_name>/`, organized by module:

```
results/<analysis_name>/
├── preprocess/           # Annotated dataset (.h5ad, .rds, .mtx)
├── spectra/              # Spectra gene programs + WMW/LMM results + plots
├── DE_analysis/          # DEGs + correlation plots
├── functional_enrichment/# GSEA + ORA results
├── cNMF/                 # cNMF factorization + WMW/LMM results + plots
└── collect_results/      # Comparative table + UpSet concordance plots
```

Logs follow a parallel structure under `logs/<analysis_name>/`. Some rules produce `.done.txt` sentinel files to mark completion of modules with unbounded outputs; delete these to force a re-run.

## Authors

- Filippo Gastaldello
  - Department of Cellular, Computational and Integrative Biology (CIBIO), University of Trento, Trento 38123, Italy
  - Fondazione The Microsoft Research—University of Trento Centre for Computational and Systems Biology (COSBI), Rovereto 38068, Italy
  - [ORCID profile](https://orcid.org/0009-0004-1676-8898)

## Citations

If you use this workflow, please cite:

- **GeneProgramsProfiler**: [https://github.com/FilippoGas/GeneProgramsProfiler](https://github.com/FilippoGas/GeneProgramsProfiler)
- **Spectra**: Kunes, R.Z., Walle, T., Land, M. et al. Supervised discovery of interpretable gene programs from single-cell data. Nat Biotechnol 42, 1084–1095 (2024). https://doi.org/10.1038/s41587-023-01940-3
- **cNMF**: Kotliar et al. (2019) Identifying gene expression programs of cell-type identity and cellular activity with single-cell RNA-Seq eLife 8:e43803 doi.org/10.7554/eLife.43803
- **Cytopus**: [https://github.com/wallet-maker/cytopus](https://github.com/wallet-maker/cytopus) — [Zenodo DOI](https://zenodo.org/badge/latestdoi/389175717)
- **Snakemake**: Köster, J., Mölder, F., et al. Sustainable data analysis with Snakemake. _F1000Research_, 2021. https://doi.org/10.12688/f1000research.29032.2

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
