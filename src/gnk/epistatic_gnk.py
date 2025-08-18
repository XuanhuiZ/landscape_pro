
import pickle
from pathlib import Path
import json
import pandas as pd
import numpy as np
from fitness_landscape.models.nk import create_nk_binary_landscape
from utils.landscape_helpers import to_df, binary_to_seq



def make_epistatic_gnk(prot_sys: str,
                        seed: int,
                        output_dir: Path) -> None:
    """
    
    """
    # load required information about system
    with open(f"./data/information/{prot_sys}_information.json", "r") as f:
        info_file = json.load(f)
    wt_seq = info_file["wt"]
    n_sites = info_file["n_sites"]
    mut_ls = info_file["mutations"]

    # load distance matrix
    interaction_arr = pd.read_csv(f"./results/epistatic_interaction/{prot_sys}_emperical.csv").to_numpy()  
    #  convert to adjacency matrix
    inter_arr = np.zeros_like(interaction_arr)
    inter_arr[interaction_arr > np.max(interaction_arr)/2] = 1
    np.fill_diagonal(inter_arr, 1)

    # make gnk landscape
    nk_binary_landscape = create_nk_binary_landscape(
        N= n_sites,
        K= n_sites - 1,
        seed = seed,
        adj_mat= inter_arr
    )

    # reconstruct sequences
    seq_ls = binary_to_seq(nk_binary_landscape.sequences,
                           mut_ls,
                           wt_seq)
    # make df
    df = to_df(nk_binary_landscape,
                        f"nk_k={n_sites - 1}")
    df["full_seq"]=seq_ls
    
    # save dataframe
    df.to_csv(output_dir / f'Epistatic_n{n_sites}_sites.csv',
                      index=False)
    # save class
    with open(output_dir /f'Epistatic_n{n_sites}_sites.pickle', 'wb') as f:
                pickle.dump(nk_binary_landscape, f)









