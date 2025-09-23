import os
import pandas as pd
import scanpy as sc

file_path_env = 'projects/sparsesampler_reproducibility/data'
REFERENCES = [1, 5, 10, 20, 34]
directory = "lcmv/benchmark"
PATH = os.path.join(file_path_env, directory)

def convert_h5ad_to_csv(ref):
    adata_address = os.path.join(PATH, f"{ref}/adata.h5ad")
    adata = sc.read_h5ad(adata_address)
    
    adata_x_address = os.path.join(PATH,f"{ref}/adata_x.csv")
    obs_address = os.path.join(PATH,f"{ref}/obs.csv")
    print("Reading data...")

    print("Writing data...")

    pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).to_csv(adata_x_address)
    # comma separated csv with header, no rownames
    pd.DataFrame(adata.obs).to_csv(obs_address)



for ref in REFERENCES:
    convert_h5ad_to_csv(ref)

    