"""
Script: run_spectra.py
Author: Filippo Gastaldello
Date: 18/07/26
Description:
    Run spectra on the scRNAseq dataset to identify the activated gene programs in the dataset
Expected Snakemake variables:
    snakemake.input.sc_dataset: Path to scRNAseq dataset in .h5ad format 
    snakemake.input.cytopus_list: Path to the cytopus list containing all gene programs associated
                                  to the celltypes present in the dataset
    snakemake.params.lam: Lambda parameter for spectra.est_spectra(). For more details visit https://github.com/dpeerlab/spectra
    snakemake.output.gene_scores: Path for gene score matrix
    snakemake.output.cell_scores: Path for cell score matrix
    snakemake.output.factor_markers: Path for factors markers matrix
"""
import sys

# 1. Setup logging -----------------------------------------------------------------------------------------------
# Redirect all print statement and error tracebacks to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

# 2. Import libraries --------------------------------------------------------------------------------------------
import json
import scipy.sparse as sp
import scanpy as sc
import pandas as pd
import Spectra as spc

# 3. Load dataset and cytopus list -------------------------------------------------------------------------------
# Load h5ad of scRNAseq dataset
print("Loading h5ad dataset ...")
adata = sc.read_h5ad(snakemake.input.sc_dataset)
# Set all genes as variable genes  # FIXME: make configurable via config.yaml
adata.var["highly_variable"] = True 

# RECOVERY CHECK: If adata.X was stripped, restore it from raw or layers
if adata.X is None:
    print("adata.X is missing. Attempting to recover...")
    if adata.raw is not None and adata.raw.X is not None:
        adata.X = adata.raw.X.copy()
    elif "counts" in adata.layers:  # Change "counts" if you used a different layer name
        adata.X = adata.layers["counts"].copy()
    else:
        raise ValueError("Fatal: adata.X is empty and no backup was found in .raw or .layers!")

# Check if matrix is sparse, if so, convert to dense array
if sp.issparse(adata.X):
    adata.X = adata.X.toarray()
print("Done")

# Load cytopus list
print("Loading cytopus list ...")
with open(snakemake.input.cytopus_list, 'r') as f:
    cytopus_list = json.load(f)
print("Done")

# 4. Fit Spectra model -------------------------------------------------------------------------------------------
print("Fitting Spectra ...")
trained_model = spc.est_spectra(adata=adata,
                                gene_set_dictionary=cytopus_list,
                                L=None,
                                use_highly_variable=True,
                                cell_type_key="celltype",
                                use_weights=True,
                                lam=snakemake.params.lam,
                                delta=0.001,
                                kappa=None,
                                rho=0.001,
                                use_cell_types=True,
                                n_top_vals=50,
                                label_factors=True,
                                overlap_threshold=0.2,
                                clean_gs=True,
                                min_gs_num=3,
                                num_epochs=10000) 
print("Done")

# 5. Save spectra's output --------------------------------------------------------------------------------------
print("Saving result matrices ...")
# Save result matrices
gene_scores = pd.DataFrame(adata.uns['SPECTRA_factors'], index = adata.uns['SPECTRA_overlap'].index, 
                           columns = list(adata.var.axes[0]))
markers = pd.DataFrame(adata.uns['SPECTRA_markers'].T,
                       columns = adata.uns['SPECTRA_overlap'].index)
cell_scores = pd.DataFrame(adata.obsm['SPECTRA_cell_scores'], 
                           index = adata.obs_names, 
                           columns = adata.uns['SPECTRA_overlap'].index)
gene_scores.to_csv(snakemake.output.gene_scores)
markers.to_csv(snakemake.output.factor_markers)
cell_scores.to_csv(snakemake.output.cell_scores)
print("Done")