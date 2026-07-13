import os
import numpy as np
import pandas as pd
import scanpy as sc

project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
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



for reference in REFERENCES:
    for method in METHODS:
        for size in SIZES:
            for rep in REPS:
                path_rep = os.path.join(PATH, f"{reference}/{method}/{size}/{rep}")
                os.makedirs(path_rep, exist_ok=True)





address_adata_ref = os.path.join(file_path_env,'mcc_01', "adata.h5ad")
adata = sc.read_h5ad(address_adata_ref)


print(adata.uns['title'])

import scanpy as sc
import scipy.sparse as sp


adata.obs = adata.obs[['author_cell_type', 'cell_type']]

if adata.raw is not None:
    adata.X = adata.raw.X.copy()
    adata.raw = None
else:
    raise ValueError("The 'raw' layer is missing from the AnnData object.")

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


adata.var = adata.var[['gene_short_name']]
adata.uns.clear()
adata.obsm.clear()








adata.obs.rename(columns={'cell_type': 'celltype'}, inplace=True)
adata.obs.to_csv('cell_type_mapping.csv', header=True)

label_key = 'celltype'

adata.obs['celltype'] = adata.obs['celltype'].astype('category')
adata.obs['author_cell_type'] = adata.obs['author_cell_type'].astype('category')
adata.var.index = adata.var.index.astype('object')

address_adata = os.path.join(file_path_env,'mcc_01', "adata_clean.h5ad")
adata.write_h5ad(address_adata)



address_adata_ref = os.path.join(file_path_env,'mcc_01', "adata_clean.h5ad")
adata = sc.read_h5ad(address_adata_ref)

print(adata.obs['author_cell_type'].value_counts())

print(adata.obs['celltype'].value_counts())


import numpy as np
np.random.seed(12)

osteoblasts_idx = adata.obs[adata.obs['celltype'] == 'osteoblast'].index

sampled_idx = np.random.choice(osteoblasts_idx, size=3000, replace=False)

remaining_idx = adata.obs[adata.obs['celltype'] != 'osteoblast'].index

final_idx = np.concatenate([sampled_idx, remaining_idx])

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