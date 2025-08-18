from fitness_landscape import FitnessLandscape
from fitness_landscape import BinarySequence
import pandas as pd

from typing import List

def to_df(landscape: FitnessLandscape,
          fitness_layer: str) -> pd.DataFrame:
    """
    
    """
    # process results
    gnk_seqs = landscape.sequences
    proseq = []

    for seq_arr in gnk_seqs:
        joint_seq = "".join(list(seq_arr.sequence))
        proseq.append(joint_seq)

    gnk_fitness = landscape.fitness_layers[fitness_layer]._replicates
    profitness = []

    for fitness in gnk_fitness:
        fitness_val = fitness[0].item()
        profitness.append(fitness_val)

    # make dataframe
    data = {'seq': proseq,
            'fitness': profitness}

    # save dataframe
    df = pd.DataFrame(data) 
    return df
    

"""binary seq to full seq"""
def binary_to_seq(binary_seqs: List[BinarySequence],
                  mutation_seq: List[str],
                  WT_seq: str,
                  ) -> List[str]:
    """
    
    """
    # check seq len and mutations correspond
    assert len(binary_seqs[0]) == len(mutation), \
        "Sequence length does not correspond to number of mutations provided"

    # extract mutation info
    mutation_tuples = []
    for mutation in mutation_seq:
        aa_from = str(mutation[0])
        site = int(mutation[1:-1])
        aa_to = str(mutation[-1])
        mutation_tuples.append((aa_from, site, aa_to))

    # make sequences
    full_seq= []
    for seq in binary_seqs:
        new_seq = WT_seq
        for idx, mutation in enumerate(mutation_tuples):
            # change wt seq with mutation
            if seq[idx] == 1:
                change_site = mutation[1]
                new_aa = mutation[2]
                new_seq = new_seq[:change_site] + new_aa + new_seq[change_site + 1:]
        full_seq.append(new_seq)

        return full_seq