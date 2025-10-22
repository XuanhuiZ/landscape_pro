# landscape_pro

**landscape_pro** is a Python toolkit for generating, analyzing, and benchmarking protein fitness landscapes using generalized NK (GNK) models. It supports synthetic, empirical-informed, and structure-informed landscape construction, and provides quantitative tools to assess ruggedness and active learning performance.

## Overview

This repository enables the construction and analysis of both idealized and biologically informed protein fitness landscapes using GNK models. The toolkit supports:

- Generation of GNK fitness landscapes with tunable interaction complexity
- Construction of empirical-informed and structure-informed landscapes based on interaction networks
- Ruggedness quantification using Dirichlet energy, autocorrelation, and local optima count
- Epistasis decomposition via Walsh–Hadamard transforms
- Evaluation of active learning strategies (e.g., ALDE) against random sampling baselines

## Methodological Background

### 1. Synthetic GNK Landscapes
Epistatic interactions are randomly assigned across `N` sites with interaction degree `K ∈ [0, N−1]`. This allows systematic exploration of ruggedness and search difficulty under controlled settings.

### 2. Empirical-Informed GNK Landscapes
To approximate realistic epistatic structures, pairwise interaction matrices are derived from background-averaged fitness effects or related metrics. These matrices are thresholded to retain the strongest interaction signals and binarized to form adjacency matrices, which define the GNK neighbourhoods for each site. This approach preserves major non-additive contributions while discarding weaker, potentially noisy interactions.

### 3. Structure-Informed GNK Landscapes
Spatial interaction networks are constructed using residue-residue distances (e.g., from PDB structures). Sites within a defined spatial cutoff (e.g., 4.5 Å) are considered interacting. This generates adjacency matrices that reflect local structural constraints and allow GNK landscapes to incorporate spatial proximity among mutable positions.

### 4. Epistasis Decomposition
The Walsh–Hadamard transform is applied to decompose each landscape into components corresponding to different epistasis orders. The proportion of total fitness variance explained by each order is computed, allowing for direct comparisons of interaction complexity across landscapes.

### 5. Ruggedness Metrics
Ruggedness is quantified using multiple complementary measures:
- **Dirichlet energy**: local fitness gradient magnitude
- **Autocorrelation**: decay in fitness similarity over mutational steps
- **Number of local optima**: adaptive peaks in the landscape

All metrics are computed using the `Landscapy` package and stored for comparative analysis.

##  Script Descriptions

| Script | Description |
|--------|-------------|
| `make_nk_data.py` | Generates synthetic GNK landscapes for given `N` and `K` |
| `epistatic_interactions_data.py` | Constructs empirical-informed adjacency matrices from pairwise interaction data |
| `pairwise_distance.py` | Calculates residue-residue distance matrices from structural data, used in structure-informed GNK construction |
| `epistasis_order.py` | Computes Walsh–Hadamard decomposition and quantifies variance by epistasis order |
| `ruggedness_calculation.py` | Computes ruggedness metrics (Dirichlet energy, autocorrelation, local optima) for GNK landscapes |
| `random_sampling.py` | Evaluates baseline performance of random sampling on different landscapes |
| `ALDE_expected_ratio_calculation.py` | Calculates ALDE vs random sampling performance ratios under varied conditions |

## Quick Start

1. **Clone the repository**

```bash
git clone https://github.com/XuanhuiZ/landscape_pro.git
cd landscape_pro
