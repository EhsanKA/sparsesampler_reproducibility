import os
import pandas as pd
import scanpy as sc
import scipy.sparse

# file_path_env = 'projects/sparsesampler_reproducibility/data'
# Get PROJECT_ROOT from environment or derive from script location
project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
    # Derive project root from script location (script is in notebooks/mcc_05/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))

file_path_env = os.path.join(project_root, 'data')
REFERENCES = [5, 10, 20, 25, 30]
directory = "mcc_05/benchmark"
PATH = os.path.join(file_path_env, directory)

import os
import pandas as pd
import scanpy as sc
import scipy.sparse

def convert_h5ad_to_csv(ref, path=PATH):
    adata_address = os.path.join(path, f"{ref}/adata.h5ad")
    adata_x_address = os.path.join(path, f"{ref}/adata_x.csv")
    obs_address = os.path.join(path, f"{ref}/obs.csv")

    # Check if both output files already exist
    if os.path.exists(adata_x_address) and os.path.exists(obs_address):
        print(f"CSV files already exist for {ref}, skipping conversion.")
        print(adata_x_address)
        print(obs_address)
        return
    
    print("Reading data...")
    adata = sc.read_h5ad(adata_address)

    if isinstance(adata.X, scipy.sparse.csr_matrix):
        adata.X = adata.X.toarray()
        
    print("Writing data...")
    pd.DataFrame(data=adata.X, index=adata.obs_names, columns=adata.var_names).to_csv(adata_x_address)
    pd.DataFrame(adata.obs).to_csv(obs_address)

    print(f"CSV files saved for {ref}.")



for ref in REFERENCES:
    convert_h5ad_to_csv(ref)