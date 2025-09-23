import os
import numpy as np
import pandas as pd
import scanpy as sc

file_path_env = 'projects/sparsesampler_reproducibility/data'
OBS_FEATURES = ['celltype']

REFERENCES = [1, 5, 10, 20, 34]
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
SIZES = [50000, 100000, 200000]
REPS = [i for i in range(5)]


directory = "lcmv/benchmark"
PATH = os.path.join(file_path_env, directory)
os.makedirs(PATH, exist_ok=True)

def create_adata(df, obs_features=None):
    """
    Create an AnnData object from a DataFrame.

    Parameters:
    df (pd.DataFrame): Input DataFrame.
    obs_features (list, optional): List of observation features. Defaults to ["Refined_clustering"].

    Returns:
    sc.AnnData: The created AnnData object.
    """
    # Set default observation features if none are provided
    obs_features = ["Refined_clustering"] if obs_features is None else obs_features

    # Create observation DataFrame
    obs = df[obs_features].astype('category')

    # Create variable DataFrame
    var = pd.DataFrame(index=df.drop(columns=obs_features).columns.values)

    # Extract data matrix
    X = df.drop(columns=obs_features).values

    # Create AnnData object
    adata = sc.AnnData(X, obs=obs, var=var)

    return adata


def create_ref_adata_unaligned():
    
    input_file = os.path.join(file_path_env,'lcmv', "2024-02-27_LCMV_all_cells.csv")

    # Read the data from the CSV file
    data = pd.read_csv(input_file, low_memory=False)
    data = data.drop(columns=DROP_FEATURES)

    # Create an AnnData object from the data
    adata = create_adata(data, obs_features=OBS_FEATURES)
    
    # # Apply arcsinh transformation to the data
    # x = np.arcsinh(adata.X/500)
    # adata.X = x

    # Construct the output file path
    output_adata_ref = os.path.join(file_path_env,'lcmv', "adata_lcmv.h5ad")

    # Write the AnnData object to a file
    adata.write(output_adata_ref)

    return adata


### Generate the directories

for reference in REFERENCES:
    for method in METHODS:
        for size in SIZES:
            for rep in REPS:
                path_rep = os.path.join(PATH, f"{reference}/{method}/{size}/{rep}")
                os.makedirs(path_rep, exist_ok=True)



### prepare data


adata = create_ref_adata_unaligned()
# address_adata_ref = os.path.join(file_path_env,'lcmv', "adata_lcmv.h5ad")
# adata = sc.read_h5ad(address_adata_ref)

label_key = 'celltype'

adata.obs[label_key] = adata.obs[label_key].astype('category')
adata.var.index = adata.var.index.astype('object')



### Generate the references

for ref in REFERENCES:
    np.random.seed(164 + ref)
    output_address = os.path.join(PATH, f"{ref}/adata.h5ad")
    
    if ref ==34:
        adata.write(output_address)
    else:
        random_indices = np.random.choice(adata.shape[0], size=ref*1000000, replace=False)
        sampled_adata = adata[random_indices].copy()
        sampled_adata.write(output_address)
