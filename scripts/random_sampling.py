
from argparse import ArgumentParser
from tqdm import tqdm
import os
import pandas as pd

# get user supplied args
parser = ArgumentParser()
parser.add_argument("--landscape_dir", type = str, help="Path to folder containing pickled landscapes.")
parser.add_argument("--running_cycle", type= int, help= "number of iterations")
args = parser.parse_args()
# load random data CSV

args = parser.parse_args()
landscape_dir = args.landscape_dir



random_ls = os.listdir(landscape_dir)

for item in tqdm(random_ls):

    for batsize in [1, 25, 50, 100, 150, 250, 350, 500]:

        # check file extension
        if not item.endswith(".csv"):
            continue

        df = pd.read_csv(landscape_dir + item)
       


        #data storage
        data_storage = {
            "cycle": [],
            "sequence": [],
            "fitness": [],
            "max_fitness": [],
        }

        for cycle in tqdm(range(args.running_cycle)):
            #check length of the dataframe is big enough for batch size
            df_len = len(df)
            if df_len < batsize:
                break

            # perform random search
            sample_row = df.sample(n=batsize, random_state= 10)
            drop_inx = sample_row.index
            df = df.drop(drop_inx)

            # extract data/results
            sample_seq = sample_row["seq"].tolist()
            sample_fitness = sample_row["fitness"].tolist()
            cycle_ls = [cycle + 1]*batsize
            max_fit_list = [max(sample_fitness)]*batsize

            # add to data storage
            data_storage["cycle"] += cycle_ls
            data_storage["sequence"] += sample_seq
            data_storage["fitness"] += sample_fitness
            data_storage["max_fitness"] += max_fit_list

        # save csv
        df = pd.DataFrame(data_storage)
        df.to_csv(f"./results/random_searching/random_search_{item}{batsize}_batchsize.csv")



        







