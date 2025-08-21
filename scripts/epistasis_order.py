"""order graph"""
import os
import pandas as pd
import pandas.io.formats.format as _fmt
pd.io.formats.string = _fmt
from argparse import ArgumentParser
import pickle
from tqdm import tqdm
from fitness_landscape.analysis.epistasis import calculate_epistasis_walsh

# get user supplied args
parser = ArgumentParser()
parser.add_argument("--landscape_dir", type = str, help="Path to folder containing pickled landscapes.")
parser.add_argument("--n_sites", type= int, help= "Default number of variable sites.")
parser.add_argument("--prot_sys", type= str, help= "Name of protein system.")

args = parser.parse_args()

landscape_dir = args.landscape_dir
n_sites = args.n_sites
prot_sys = args.prot_sys

# Store results
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
        
    # get variance explained by orders
    order_res = calculate_epistasis_walsh(loaded_data, 
                                          order=n)
    
    # update results
    for order in range(1, 13):
        if order in order_res["variance_explained"]:
            result_dict[f"{order}_order_variance_explained"].append(order_res["variance_explained"][order])
        else:
            result_dict[f"{order}_order_variance_explained"].append(None)
        
    result_dict["neighbourhood_scheme"].append(neighbourhood_scheme)
    result_dict["N"].append(n)
    result_dict["K"].append(k)

# make csv
df = pd.DataFrame(result_dict)
df.to_csv(f'./results/epistasis_orders/{prot_sys}_epistatic_orders.csv', index=False)