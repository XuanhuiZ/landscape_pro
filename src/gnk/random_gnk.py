"""
Generating GNK data for PTE and thermobodies.
"""

import pickle
from pathlib import Path
from tqdm import tqdm
from utils.landscape_helpers import to_df
from fitness_landscape.models.nk import create_nk_binary_landscape

def make_random_gnk(n_min: int,
                    n_max: int,
                    seed: int,
                    output_dir: Path) -> None:
    """
    Generate random GNK binary fitness landscapes for varying N and K values,
    and save both the landscape objects and associated DataFrames.
    Parameters
    ----------
    n_min : int
    Minimum number of variable sites (N) to start from.
    n_max : int
    Maximum number of variable sites (N) to go up to (inclusive).
    seed : int
    Random seed used for reproducibility.
    output_dir : Path
    Directory where output .csv and .pickle files will be saved.
    """
    # starting from n_min to n_max (+ 1 for Python indexing)
    for N in tqdm(range(n_min, n_max + 1), 
                  desc="Generating random NK for N"):
        k_max = N - 1
        # starting from k_min to k_max (+ 1 for Python indexing)
        for K in tqdm(range(k_max + 1),
                      desc=f"Generating random NK for N={N} for K",
                      leave=False):
            # make binary landsacpe
            nk_binary_landscape = create_nk_binary_landscape(
                N= N,
                K= K,
                seed = seed
            )
            # make dataframe
            df = to_df(nk_binary_landscape,
                        f"nk_k={K}")
            # save dataframe
            df.to_csv(output_dir / f'nk_random_n{N}_k{K}_sites.csv',
                      index=False)
            # save class
            with open(output_dir /f'nk_random_n{N}_k{K}_sites.pickle', 'wb') as f:
                pickle.dump(nk_binary_landscape, f)

