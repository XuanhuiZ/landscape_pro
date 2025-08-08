import pandas as pd



df = pd.read_csv("./data/raw_data/PTE_2NH_ActvSite.csv")

seqs = df.sequence.tolist()

first_seq = seqs[0]

mutations = [233, 254, 271, 272, 306, 313]

for mut in mutations:
    print(first_seq[mut - 1])
    