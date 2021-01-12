# Plot neutralization curves and pseudovirus titers

The input neutralizaiton data are in [./data/](data).
These are plate reader data processed to give the fraction infectivity at each concentration, taken from Amin's BloomLab computational notebooks.
Each file is named by the date of the assay.

The input titer data are in [./titer_data/](titer_data). This file contains the normalized RLU/µL of pseudovirus for each of the constructs.

The Python Jupyter notebook [analyze_neuts.ipynb](analyze_neuts.ipynb) plots these neut data.
