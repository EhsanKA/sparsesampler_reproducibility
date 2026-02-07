import os
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sparsesampler.sampling import sample
from scsampler import scsampler
import time

# Get the file path from environment variable
# Get PROJECT_ROOT from environment or derive from script location
project_root = os.environ.get('PROJECT_ROOT')
if project_root is None:
    # Derive project root from script location (script is in notebooks/lcmv/)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(os.path.dirname(script_dir))

file_path_env = os.path.join(project_root, 'data')
OBS_FEATURES = ['prediction','organ','sample_group','label','group','celltype','sample_id']
DROP_FEATURES = ['SSC-B-H','SSC-B-A']

REFERENCES = [1, 5, 10, 20, 34]
METHODS = ['random', 'sps', 'hopper', 'atomic', 'scsampler']
methods = ['random', 'sps', 'hopper', 'scsampler']
SIZES = [50000, 100000, 200000]
REPS = [i for i in range(5)]
label_key = 'celltype'


directory = "lcmv/benchmark"
PATH = os.path.join(file_path_env, directory)



import numpy as np
import pandas as pd

def generate_sps(adata, size, seed = 1234):
    return sample(X=adata.X, size=size, seed=seed, auto_k=False, k=None, feature_index=12)

def generate_scsampler(adata, size, seed = 1234):
    print(f'********* #Start# *********')
    start_time = time.time()
    
    sc.tl.pca(adata, svd_solver='arpack')
    arr = adata.obsm['X_pca']
    arr = arr - arr.min(axis=0)
    arr = arr / np.abs(arr).max(axis=0)
    adata.obsm['X_pca'] = arr
    del arr

    mat = adata.obsm['X_pca']
    res = scsampler(mat, n_obs= size, copy = True, random_split = 16)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

    return res[1].tolist(), elapsed_time

import sys
import io
from contextlib import redirect_stdout

def generate_hopper(adata, size, seed=1234):
    
    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)

    from hopper.treehopper.hoppers import hopper, treehopper, PCATreePartition
    print(f'********* #Start# *********')

    np.random.seed(seed)
    start_time = time.time()

    scaler = StandardScaler()
    data_standardized = scaler.fit_transform(adata.X)
    #data_standardized = scaler.fit_transform(adata_shuffled.X)
    X = data_standardized

    N = size # Number of samples to obtain from the data set.
    with io.StringIO() as buf, redirect_stdout(buf):
        th = treehopper(data_standardized, partition=PCATreePartition, max_partition_size=1000)
        sketch = th.hop(size)

    samples = th.path[:size]

    # X_sketch = X_dimred[sketch_index]
    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time


def generate_random(adata, size, seed=1234):

    print(f'********* #Start# *********')
    np.random.seed(seed)
    start_time = time.time()


    # samples = np.random.randint(0, adata.shape[0], size=size)
    samples = np.random.choice(adata.shape[0], size=size, replace=False)

    elapsed_time = time.time() - start_time
    print(f"Elapsed time: {elapsed_time} seconds")

        
    return samples, elapsed_time




import argparse
import numpy as np
import pickle

# Assuming the sampling functions are defined elsewhere and imported
# from sampling_functions import generate_cubic, generate_random, generate_hopper

def main(ref, method, size, rep, seed):
    # Define a dictionary to store the results
    results = []
    
    address = os.path.join(PATH, f"{ref}/adata.h5ad")
    adata = sc.read_h5ad(address)
    adata.obs[label_key] = adata.obs[label_key].astype('category')
    adata.var.index = adata.var.index.astype('object')
    
    method_dict = {
    "sps": (generate_sps, {"adata": adata, "size": size, "seed": seed}),
    "hopper": (generate_hopper, {"adata": adata, "size": size, "seed": seed}),
    "random": (generate_random, {"adata": adata, "size": size, "seed": seed}),
    "scsampler": (generate_scsampler, {"adata": adata, "size": size, "seed": seed}),
    }
    
    if method in method_dict:
        func, args = method_dict[method]
        results = func(**args)
    else:
        print(f"No function associated with {method}")
        return
        
    output_address = os.path.join(PATH, f"{ref}/{method}/{size}/{rep}/results.pkl")
    
    with open(output_address, 'wb') as handle:
        pickle.dump(results, handle, protocol=pickle.HIGHEST_PROTOCOL)

        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run sampling methods in parallel for a given Reference, Method, Size, Replicate, and Seed.")
    parser.add_argument("--ref", type=int, required=True, help="Reference to process")
    parser.add_argument("--method", type=str, required=True, help="Method to process")
    parser.add_argument("--size", type=int, required=True, help="Size to process")
    parser.add_argument("--rep", type=int, required=True, help="Replicate to process")
    parser.add_argument("--seed", type=int, required=True, help="Seed to process")

    args = parser.parse_args()
    print("###############################")
    print("************ New run **********")
    
    # ref = 1
    # method = 'sps'
    # size = 50000
    # rep = 0
    # seed = 6547
    
    # main(args.batch_id, args.size)
    # main(ref, method, size, rep, seed)
    main(args.ref, args.method, args.size, args.rep, args.seed)