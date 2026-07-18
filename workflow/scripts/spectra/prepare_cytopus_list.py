"""
Script: prepare_cytopus_list.py
Author: Filippo Gastaldello
Date: 17/07/26
Description:
    Connect to the cytopus db with the cytopus python package and download
    the gene sets associated to the celltypes in the conversion dictionary
Expected Snakemake variables:
    snakemake.input[0]: Path to celltype conversion dictionary used to map celltypes in
                                                    the dataset to celltypes in cytopus
    snakemake.output[0]: Destination path for cytopus list
    snakemake.params.global_celltype: The cytopus celltype to use as global celltype
                                     (Refer to cytopus Docs for more info: https://github.com/wallet-maker/cytopus)
    snakemake.log[0]: Path to the script log file
"""

import sys

# 1. Setup logging -----------------------------------------------------------------------------------------------
# Redirect all print statement and error tracebacks to the log file
sys.stderr = sys.stdout = open(snakemake.log[0], 'w')

# 2. Import libraries --------------------------------------------------------------------------------------------
import cytopus as cp
import json

# 3. Load conversion dictionary ----------------------------------------------------------------------------------
print("Loading cell type conversion dictionary ...")
with open(snakemake.input[0]) as f:
    conversion_dict = json.load(f)
print("Done")

# 4. Load Cytopus and get gene sets ------------------------------------------------------------------------------
# Load resource
print("Loading cytopus ...")
G = cp.KnowledgeBase()
print("Done")

# Extract cell types of interest from conversion dictionary
print("Extracting cell types of interest from conversion dictionary ...")
celltypes_of_interest = list(dict.fromkeys(conversion_dict.values()))
global_celltype = list(snakemake.params.global_celltype)
print("Done")

# Get spectra dictionary
print("Getting spectra dictionary ...")
G.get_celltype_processes(celltypes_of_interest, global_celltypes=global_celltype)
processes_dict = G.celltype_process_dict
print("Done")

# Replace process_dict keys with the matching celltype in the dataframe
print("Converting cytopus celltype using celltype conversion dictionary ...")
spectra_dict = {}
for key, value in conversion_dict.items():
    spectra_dict[key] = processes_dict[value]
spectra_dict["global"] = processes_dict["global"]
with open(snakemake.output[0], 'w') as f:
    json.dump(spectra_dict, f)
print("Done")