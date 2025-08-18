"""
Generating GNK data for PTE and thermobodies.
"""

import pickle

from pathlib import Path

from utils.landscape_helpers import to_df

from fitness_landscape.models.nk import create_nk_binary_landscape

def make_random_gnk(n_min: int,
                    n_max: int,
                    seed: int,
                    output_dir: Path) -> None:
    """
    python docstring format
    """
    # starting from n_min to n_max
    for N in range(n_min, n_max + 1):
        k_max = N -1
        for K in range(k_max + 1):
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
            df.to_csv(output_dir / f'Random_n{N}_k{K}_sites.csv',
                      index=False)
            # save class
            with open(output_dir /f'Random_n{N}_k{K}_sites.pickle', 'wb') as f:
                pickle.dump(nk_binary_landscape, f)

