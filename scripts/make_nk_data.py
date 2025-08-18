

from gnk.epistatic_gnk import make_epistatic_gnk
from gnk.random_gnk import make_random_gnk
from gnk.structural_gnk import make_structural_gnk

# make synthetic
print("Making random NKs")
make_random_gnk(n_min=3,
                n_max=12,
                seed=0,
                output_dir="./data/datasets/nk_data_new/synthetic_random")

# make epistatic
print("Making epistatic")
make_epistatic_gnk(prot_sys = "pte",
                        seed=0,
                        output_dir="/Users/zhouxuanhui/Desktop/ML_ASC_code/landscape_pro/data/datasets/nk_data_new/pte_epistatic")

# make structural
print("Making structural")
make_structural_gnk(prot_sys = "pte",
                        seed = 0,
                        distance_thershold = 4.5,
                        output_dir = "/Users/zhouxuanhui/Desktop/ML_ASC_code/landscape_pro/data/datasets/nk_data_new/pte_structural")

make_structural_gnk(prot_sys = "thermo",
                        seed = 0,
                        distance_thershold = 4.5,
                        output_dir = "/Users/zhouxuanhui/Desktop/ML_ASC_code/landscape_pro/data/datasets/nk_data_new/thermo_structural")