"""calculating the ruggedness"""
import json
import os
import pandas as pd
import pickle

from argparse import ArgumentParser
from tqdm import tqdm

from fitness_landscape.analysis.dirichlet_energy import calculate_ruggedness_dirichlet_energy
from fitness_landscape.analysis.graph import calculate_ruggedness_local_optima
from fitness_landscape.analysis.random_walk import calculate_ruggedness_autocorrelation_analytical



# landscape_dir = "./data/nk_data/PTE/"
# n_sites = 6
# output_name = "PTE_ruggedness"

# get user supplied args
parser = ArgumentParser()
parser.add_argument("--landscape_dir", type = str, help="Path to folder containing pickled landscapes.")
parser.add_argument("--n_sites", type= int, help= "Default number of variable sites.")
parser.add_argument("--output_name", type= str, help= "Output save name.")

args = parser.parse_args()

landscape_dir = args.landscape_dir
n_sites = args.n_sites
output_name = args.output_name


# data storage
result_dict = {
    "neighbourhood_scheme": [],
    "N": [],
    "K": [],
    "dirichlet_energy": [],
    "local_optima_energy": [],
    "autocorrelation_energy": [],
}

full_autocorrelation_res = {}

landscape_ls = os.listdir(landscape_dir)

# iterate through all landscapes
for item in tqdm(landscape_ls):

    # check file extension
    if not item.endswith(".pickle"):
        continue
    
    item = item[:-7]

    # extract landscape properties
    item_elements = item.split("_")

    neighbourhood_scheme = item_elements[1]

    if neighbourhood_scheme == "random":
        n = int(item_elements[2][1:])
        k = int(item_elements[3][1:])

    else:
        n = n_sites
        k = None

    # load premade landscape
    with open(landscape_dir + item + ".pickle", 'rb') as f:
        loaded_data = pickle.load(f)

    #dirichlet energy
    dirichlet_energy_res = calculate_ruggedness_dirichlet_energy(loaded_data)["total_dirichlet_energy"].item()

    #local
    local_energy_res = calculate_ruggedness_local_optima(loaded_data)["local_optima_count"]


    #auto energy data
    auto_energy_data = calculate_ruggedness_autocorrelation_analytical(loaded_data)["autocorrelation"]
    auto_1 = auto_energy_data[1].item()

    # add to results
    result_dict["neighbourhood_scheme"].append(neighbourhood_scheme)
    result_dict["N"].append(n)
    result_dict["K"].append(k)
    result_dict["dirichlet_energy"].append(dirichlet_energy_res)
    result_dict["local_optima_energy"].append(local_energy_res)
    result_dict["autocorrelation_energy"].append(auto_1)

    # add all autocorr res
    res_key = item.split(".")[0]
    full_autocorrelation_res[res_key] = [val.item() for val in auto_energy_data]


# make csv
df = pd.DataFrame(result_dict)
df.to_csv(f'./data/results/{output_name}.csv', index=False)

# save json
with open(f"./data/results/{output_name}_autocorrelation.json", "w") as f:
    json.dump(full_autocorrelation_res, f)