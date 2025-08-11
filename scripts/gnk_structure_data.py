""""
Generating GNK based on episatatic interactions in the PTE data.
"""
import fitness_landscape

import json
import numpy as np
import pandas as pd
import pickle

prot_sys = "thermo"

# load required information about system
with open(f"./data/information/{prot_sys}_information.json", "r") as f:
    info_file = json.load(f)
wt_seq = info_file["wt"]
ordered_sites = info_file["ordered_sites"]
n_sites = info_file["n_sites"]


# load distance matrix
dist_arr = pd.read_csv(f"./results/structure_interaction/{prot_sys}_dist.csv").to_numpy()

# convert to adjacency matrix
inter_arr = np.zeros_like(dist_arr)
inter_arr[dist_arr < 4.5] = 1
np.fill_diagonal(inter_arr, 1)

# make gnk landscape
variable_sites = [site - 1 for site in ordered_sites]
variable_sites = variable_sites[:n_sites]

# convert to binary form (A == WT, C == mutant)
wt_seq_modified = wt_seq
for site in variable_sites:
    wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]

gnk_landscape = fitness_landscape.models.create_gnk_landscape(
    N=len(variable_sites),
    K=len(variable_sites) - 1,
    alphabet=['A', 'C'],
    seed=0,
    base_sequence=list(wt_seq_modified),
    variable_sites = variable_sites,
    adj_mat=inter_arr
)

 # save class
with open(f'data/datasets/nk_data/{prot_sys}/{prot_sys}_structure_n{len(variable_sites)}.pickle', 'wb') as f:
    pickle.dump(gnk_landscape, f)

# process results
gnk_seqs = gnk_landscape.sequences
pro_seq = []
for seq_arr in gnk_seqs:
    joint_seq = "".join(seq_arr.sequence)
    pro_seq.append(joint_seq)

gnk_fitness = gnk_landscape.fitness_layers[f"nk_k={len(variable_sites) - 1}"]._replicates
pro_fitness = []

for fitness in gnk_fitness:
    fitness_val = fitness[0].item()
    pro_fitness.append(fitness_val)

# make dataframe
data = {'seq': pro_seq,
        'fitness': pro_fitness}

# save dataframe
df = pd.DataFrame(data)
df.to_csv(f'data/datasets/nk_data/{prot_sys}/{prot_sys}_structure_n{len(variable_sites)}.csv', index=False)