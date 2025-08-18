"""
Generating GNK data for PTE and thermobodies.
"""
import fitness_landscape

import pandas as pd
import pickle

from tqdm import tqdm

# # load pte data
# pte_df = pd.read_csv("./data/raw_data/PTE_2NH_ActvSite.csv")
# wt_seq = pte_df.sequence.tolist()[0]

# ordered_sites = [254, 233, 272, 271, 313, 306, 308, 172,203, 174, 258, 130]

# variable_site_ls = []
# for k in range(3, 13):
#     variable_site_ls.append(
#         [site - 1 for site in ordered_sites[:k]]
#     )

# wt_seq_modified = wt_seq
# for site in ordered_sites:
#     wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]

# prepare thermobody data
wt_seq = "KEVKAAELAAAKEAAKAELKALNLSEGQKDFYIKKINDAKTVEGVKALLEEALKLNDAKKEI"
wt_marked = "KEVKAAELAAAKXXAXXELXXLNLSXXQXXFYIXKIXDAKTVEGVKALLEEALKLNDAKKEI"

ordered_sites = [13, 14, 16, 17, 20, 21, 26, 27, 29, 30, 34, 37]
for i, site in enumerate(wt_marked):
    if site == "X":
        ordered_sites.append(i)

variable_site_ls = []
for k in range(3, 13):
    variable_site_ls.append(
        [site - 1 for site in ordered_sites[:k]]
    )
wt_seq_modified = wt_seq
for site in ordered_sites:
    wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]


for varsite in tqdm(variable_site_ls,
                    desc="Number of sites varied"):

    n_sites = len(varsite)

    for k_val in range(0, n_sites):

        # get landscape
        gnk_landscape = fitness_landscape.models.create_gnk_landscape(
            N=n_sites,
            K=k_val,
            alphabet=['A', 'C'], # where A == WT and C = mutation
            seed=0,
            base_sequence=list(wt_seq_modified),
            variable_sites = varsite
        )

        # save class
        with open(f'data/nk_data/Thermo/Thermo_random_n{n_sites}_k{k_val}.pickle', 'wb') as f:
            pickle.dump(gnk_landscape, f)

        # process results
        gnk_seqs = gnk_landscape.sequences
        proseq = []

        for seq_arr in gnk_seqs:
            joint_seq = "".join(seq_arr.sequence)
            proseq.append(joint_seq)

        gnk_fitness = gnk_landscape.fitness_layers[f"nk_k={k_val}"]._replicates
        profitness = []

        for fitness in gnk_fitness:
            fitness_val = fitness[0].item()
            profitness.append(fitness_val)

        # make dataframe
        data = {'seq': proseq,
                'fitness': profitness}

        # save dataframe
        df = pd.DataFrame(data)
        df.to_csv(f'data/nk_data/Thermo/Thermo_random_n{n_sites}_k{k_val}.sites.csv', index=False)