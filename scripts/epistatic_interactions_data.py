""""
Generating interaction matrices for PTE and PTE control.
"""
import numpy as np
import pandas as pd
import pickle

from fitness_landscape.core import BinarySequence, FitnessLandscape
from fitness_landscape.core.fitness import NumericFitness
from fitness_landscape.core.graph import create_hamming_graph
from fitness_landscape.analysis.epistasis import get_epistasis_matrix

# load raw data
pte_norm_df    = pd.read_csv("./data/raw_data/PTE_norm.csv", dtype={"sequence": str})
pte_control_df = pd.read_csv("./data/raw_data/PTE_smooth_con.csv", dtype={"sequence": str})

# convert control binary to str
name_seq_list = pte_control_df['binary_sequence'].astype(str)
str_len = name_seq_list.str.len().max()
padded_seqs = name_seq_list.str.zfill(str_len)
pte_control_df["binary_sequence"] = padded_seqs


# add fake noise based on error in real data
name_fit_list = pte_control_df['fitness_mean'].tolist()

## get mean and stdev from emperical
pte_norm_mean = pte_norm_df['fitness_stdev'].mean()
pte_norm_stdev = pte_norm_df['fitness_stdev'].std()
## make noise of length 64 based on mean and stdev
noise = np.random.normal(loc = pte_norm_mean, scale = pte_norm_stdev, 
                         size = len(name_fit_list))
## add noise to control data
noisy_fitness = [[name_fit_list[i] + noise[i].item()] for i in range(len(name_fit_list))]


# process sequences 
seq_arr = np.array([list(seq) for seq in padded_seqs]).astype(int)
base_np_seq = [BinarySequence(sequence=seq) for seq in seq_arr]

# process fitness data
base_fit = NumericFitness("control",
                          noisy_fitness)

# make and save fitness landscape
hamming_graph = create_hamming_graph(base_np_seq, name_fit_list)
fit_land = FitnessLandscape (base_np_seq,
                             hamming_graph,
                             fitness_layers= {"noisy": base_fit},
                              )
with open('./data/nk_data/PTE/PTE_control.pickle', 'wb') as f:
    pickle.dump(fit_land, f)


# get epistatic interactions for emperical
epi_matrix = get_epistasis_matrix(landscape=fit_land)
df = pd.DataFrame(epi_matrix)
df.to_csv('./results/epistatic_interaction/PTE_control.csv', index=False)
print()