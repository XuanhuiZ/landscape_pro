""""
Generating GNK based on episatatic interactions in the PTE data.
"""
import fitness_landscape

import numpy as np
import pandas as pd
import pickle

# load pte data
pte_df = pd.read_csv("./data/raw_data/PTE_2NH_ActvSite.csv")
wt_seq = pte_df.sequence.tolist()[0]

# load epistatic interactions
pte_ep_arr = pd.read_csv("./results/epistatic_interaction/PTE_emperical.csv").to_numpy()

# convert to adjacency matrix
inter_arr = np.zeros_like(pte_ep_arr)
inter_arr[pte_ep_arr > 0.00367] = 1
np.fill_diagonal(inter_arr, 1)

# make gnk landscape
variable_sites = [233,254,271,272,306,313]
variable_sites = [site - 1 for site in variable_sites]

wt_seq_modified = wt_seq
for site in variable_sites:
    wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]

gnk_landscape = fitness_landscape.models.create_gnk_landscape(
    N=6,
    K=3,
    alphabet=['A', 'C'],
    seed=0,
    base_sequence=list(wt_seq_modified),
    variable_sites = variable_sites,
    adj_mat=inter_arr
)
print()
 # save class
with open(f'data/nk_data/PTE/PTE_epistasis_n6.pickle', 'wb') as f:
            pickle.dump(gnk_landscape, f)

             # process results
gnk_seqs = gnk_landscape.sequences
proseq = []
for seq_arr in gnk_seqs:
            joint_seq = "".join(seq_arr.sequence)
            proseq.append(joint_seq)

gnk_fitness = gnk_landscape.fitness_layers[f"nk_k=3"]._replicates
profitness = []

for fitness in gnk_fitness:
            fitness_val = fitness[0].item()
            profitness.append(fitness_val)

        # make dataframe
data = {'seq': proseq,
                'fitness': profitness}

        # save dataframe
df = pd.DataFrame(data)
df.to_csv(f'data/nk_data/PTE/PTE_epistasis_n6_sites.csv', index=False)