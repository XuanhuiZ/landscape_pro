"""
Generate NK landscape data.
"""

import numpy as np
import os
import pandas as pd
from pathlib import Path
import pickle
from fitness_landscape import (create_hamming_graph, 
                               BinarySequence, 
                               FitnessLandscape, 
                               NumericFitness)
from gnk.epistatic_gnk import make_epistatic_gnk
from gnk.random_gnk import make_random_gnk
from gnk.structural_gnk import make_structural_gnk
from utils.landscape_helpers import pad_binary_ints


# make synthetic
print("Making random NK landscapes...")
os.makedirs("./data/datasets/nk_data_new/synthetic_random", exist_ok=True)
make_random_gnk(n_min=3,
                n_max=12,
                seed=0,
                output_dir=Path("./data/datasets/nk_data_new/synthetic_random"))

# make emperical
print("Making emperical landscapes...")
os.makedirs("./data/datasets/emperical", exist_ok=True)
df = pd.read_csv("./data/raw_data/PTE_norm.csv")
sequences = [BinarySequence(seq) for seq in
             pad_binary_ints(df["binary_sequence"].tolist())]
fitness = {"arylesterase_activity": NumericFitness(
    name="arylesterase_activity",
    values=[[item] for item in df["fitness_mean"].tolist()]
)}
graph = create_hamming_graph(sequences)
landscape = FitnessLandscape(sequences=sequences, 
                             fitness_layers=fitness, 
                             graph=graph)
with open("./data/datasets/emperical/pte_emperical.pickle", "wb") as f:
    pickle.dump(landscape, f)

# make epistatic
print("Making epistatic NK landscapes...")
os.makedirs("./data/datasets/nk_data_new/pte_epistatic", exist_ok=True)
make_epistatic_gnk(prot_sys="pte",
                   seed=0,
                   output_dir=Path("./data/datasets/nk_data_new/pte_epistatic"))


# make structural
print("Making structural NK landscapes...")
os.makedirs("./data/datasets/nk_data_new/pte_structural", exist_ok=True)
make_structural_gnk(prot_sys="pte",
                    seed=0,
                    distance_threshold=4.5,
                    output_dir = Path("./data/datasets/nk_data_new/pte_structural"))
os.makedirs("./data/datasets/nk_data_new/thermo_structural", exist_ok=True)
make_structural_gnk(prot_sys="thermo",
                    seed=0,
                    distance_threshold=4.5,
                    output_dir=Path("./data/datasets/nk_data_new/thermo_structural"))