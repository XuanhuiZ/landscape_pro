"""order graph"""
import os
import pandas as pd
import pickle

from tqdm import tqdm

from fitness_landscape.analysis.epistasis import calculate_epistasis_walsh
from fitness_landscape.core.sequence import BinarySequence

## system parameters
landscape_dir = "./data/nk_data/PTE/"
n_sites = 6

pte_df = pd.read_csv("./data/raw_data/PTE_2NH_ActvSite.csv")
wt_seq = pte_df.sequence.tolist()[0]
wt_seq_modified = wt_seq
ordered_sites = [254, 233, 272, 271, 313, 306, 308, 172,203, 174, 258, 130]
for site in ordered_sites:
    wt_seq_modified = wt_seq_modified[:site -1] + "A" + wt_seq_modified[site:]


result_dict = {
    "neighbourhood_scheme": [],
    "N": [],
    "K": [],
}
for order in range(1, 13):
    result_dict[f"{order}_order_variance_explained"] = []


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
            
    # convert sequences in landscape into base sequences
    if not isinstance(loaded_data.sequences[0], BinarySequence):
        relevant_sites = ordered_sites[:n]
        binary_seq_ls = []
        for seq in loaded_data.sequences:
            # get binary sequence given WT
            binary_seq = []
            seq_basenp = seq.sequence
            for site in relevant_sites:
                site_indicator = seq_basenp[site - 1]
                if site_indicator == "A":
                    binary_seq.append(0)
                elif site_indicator == "C":
                    binary_seq.append(1)
                else:
                    raise ValueError(f"Unexpected AA: {site_indicator} at site: {site}")
            # make binary sequence class
            binary_seq_object = BinarySequence(binary_seq)
            binary_seq_ls.append(binary_seq_object)

        ## overwrite landscapes sequences with binary
        loaded_data.sequences = binary_seq_ls

    # get variance explained by orders
    order_res = calculate_epistasis_walsh(loaded_data, 
                                        order=n)
    
    # update results
    for order in range(1, 13):
        if order in order_res["variance_explained"]:
            result_dict[f"{order}_order_variance_explained"].append(order_res["variance_explained"][order].item())
        else:
            result_dict[f"{order}_order_variance_explained"].append(None)
        
    result_dict["neighbourhood_scheme"].append(neighbourhood_scheme)
    result_dict["N"].append(n)
    result_dict["K"].append(k)

# make csv
df = pd.DataFrame(result_dict)
df.to_csv(f'./data/results/PTE_order.csv', index=False)
