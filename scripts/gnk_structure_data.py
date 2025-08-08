""""
Generating GNK based on episatatic interactions in the PTE data.
"""
import fitness_landscape

import numpy as np
import pandas as pd
import pickle

# # load pte data
pte_df = pd.read_csv("./data/raw_data/PTE_2NH_ActvSite.csv")
wt_seq = pte_df.sequence.tolist()[0]

# # load epistatic interactions
# dist_arr = pd.read_csv("./results/structure_interaction/PTE_structure.csv").to_numpy()

# load thermobody data

# wt_seq = "KEVKAAELAAAKEAAKAELKALNLSEGQKDFYIKKINDAKTVEGVKALLEEALKLNDAKKEI"

# load distance matrix
dist_arr = pd.read_csv("/Users/zhouxuanhui/Desktop/landscape_pro/results/structure_interaction/PTE_dis.csv").to_numpy()[:, 1:]

# convert to adjacency matrix
inter_arr = np.zeros_like(dist_arr)
inter_arr[dist_arr < 4.5] = 1
np.fill_diagonal(inter_arr, 1)

# make gnk landscape
variable_sites = [233,254,271,272,306,313]
variable_sites = [site - 1 for site in variable_sites]
#variable_sites = [13, 14, 16, 17, 20, 21, 26, 27, 29, 30, 34, 37]



wt_seq_modified = wt_seq
for site in variable_sites:
    wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]

gnk_landscape = fitness_landscape.models.create_gnk_landscape(
    N=len(variable_sites),
    K=6,
    alphabet=['A', 'C'],
    seed=0,
    base_sequence=list(wt_seq_modified),
    variable_sites = variable_sites,
    adj_mat=inter_arr
)
print()
 # save class
with open(f'data/nk_data/PTE/PTE_structure_n6.pickle', 'wb') as f:
            pickle.dump(gnk_landscape, f)

             # process results
gnk_seqs = gnk_landscape.sequences
proseq = []
for seq_arr in gnk_seqs:
            joint_seq = "".join(seq_arr.sequence)
            proseq.append(joint_seq)

gnk_fitness = gnk_landscape.fitness_layers[f"nk_k=6"]._replicates
profitness = []

for fitness in gnk_fitness:
            fitness_val = fitness[0].item()
            profitness.append(fitness_val)

        # make dataframe
data = {'seq': proseq,
                'fitness': profitness}

        # save dataframe
df = pd.DataFrame(data)
df.to_csv(f'data/nk_data/PTE/PTE_structure_n6_sites.csv', index=False)