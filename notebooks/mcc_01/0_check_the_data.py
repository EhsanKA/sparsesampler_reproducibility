import os
import numpy as np
import pandas as pd
import scanpy as sc

# Get PROJECT_ROOT from environment or derive from script location
project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
    # Derive project root from script location (script is in notebooks/mcc_01/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))

file_path_env = os.path.join(project_root, 'data')

OBS_FEATURES = ['celltype']

REFERENCES = [5, 10, 20, 25, 30]
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
SIZES = [50000, 100000, 200000, 300000]
REPS = [i for i in range(5)]
label_key = 'celltype'


directory = "mcc_01/benchmark"
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



### Generate the directories


for reference in REFERENCES:
    for method in METHODS:
        for size in SIZES:
            for rep in REPS:
                path_rep = os.path.join(PATH, f"{reference}/{method}/{size}/{rep}")
                os.makedirs(path_rep, exist_ok=True)




### prepare data

address_adata_ref = os.path.join(file_path_env,'mcc_01', "adata.h5ad")
adata = sc.read_h5ad(address_adata_ref)


print(adata.uns['title'])

import scanpy as sc
import scipy.sparse as sp


# 1. Keep only ['author_cell_type', 'cell_type'] in obs
adata.obs = adata.obs[['author_cell_type', 'cell_type']]

# 2. Assign raw.X to X
if adata.raw is not None:
    adata.X = adata.raw.X.copy()
    # Remove raw.X from the AnnData object
    adata.raw = None
else:
    raise ValueError("The 'raw' layer is missing from the AnnData object.")

# 3. Compress X using CSR format (if it's not already sparse)
if not sp.issparse(adata.X):
    adata.X = sp.csr_matrix(adata.X)



sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_genes(adata, min_cells=3)

print(adata.X.min(), adata.X.max() )

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

print(adata.X.min(), adata.X.max() )

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
adata = adata[:, adata.var.highly_variable].copy()


adata.var = adata.var[['gene_short_name']]  # Keep only gene_short_name in var
adata.uns.clear()  # Remove unnecessary unstructured data
adata.obsm.clear()  # Remove all multi-dimensional data



# # 4. Find highly variable genes and subset to keep 2000 genes
# sc.pp.highly_variable_genes(adata, n_top_genes=2000)

# # Subset the adata object to keep only highly variable genes
# adata = adata[:, adata.var.highly_variable].copy()

# # 5. Remove unnecessary information (var and other unwanted layers)
# adata.var = adata.var[['gene_short_name']]  # Keep only gene_short_name in var
# adata.uns.clear()  # Remove unnecessary unstructured data
# adata.obsm.clear()  # Remove all multi-dimensional data

# # Save the new AnnData object (optional)
# adata.write_h5ad('filtered_adata.h5ad')


adata.obs.rename(columns={'cell_type': 'celltype'}, inplace=True)
adata.obs.to_csv('cell_type_mapping.csv', header=True)

label_key = 'celltype'

adata.obs['celltype'] = adata.obs['celltype'].astype('category')
adata.obs['author_cell_type'] = adata.obs['author_cell_type'].astype('category')
adata.var.index = adata.var.index.astype('object')

address_adata = os.path.join(file_path_env,'mcc_01', "adata_clean.h5ad")
adata.write_h5ad(address_adata)


### Generate the references

# adata_clean is the backup of the proccessed adata.
address_adata_ref = os.path.join(file_path_env,'mcc_01', "adata_clean.h5ad")
adata = sc.read_h5ad(address_adata_ref)

print(adata.obs['author_cell_type'].value_counts())

print(adata.obs['celltype'].value_counts())


import numpy as np
np.random.seed(12)

# Identify the index of cells that are Pre-osteoblasts (Sp7+)
osteoblasts_idx = adata.obs[adata.obs['celltype'] == 'osteoblast'].index

# Randomly sample fewer than 3,000 cells from the osteoblasts group
sampled_idx = np.random.choice(osteoblasts_idx, size=3000, replace=False)

# Keep the rest of the cells that are not Pre-osteoblasts (Sp7+)
remaining_idx = adata.obs[adata.obs['celltype'] != 'osteoblast'].index

# Combine the sampled Pre-osteoblasts (Sp7+) with the remaining cells
final_idx = np.concatenate([sampled_idx, remaining_idx])

# Subset the adata with the final indices
adata_sampled = adata[final_idx].copy()


print(adata_sampled.obs['celltype'].value_counts())



for ref in REFERENCES:
    np.random.seed(164 + ref)
    output_address = os.path.join(PATH, f"{ref}/adata.h5ad")
    
    if ref ==30:
        adata_sampled.write(output_address)
    else:
        random_indices = np.random.choice(adata_sampled.shape[0], size=int(ref*100000), replace=False)
        sampled_adata = adata_sampled[random_indices].copy()
        sampled_adata.write(output_address)
    print(ref)