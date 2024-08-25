# Quantifying age specific household contacts for infectious disease modelling 

This repository provides data, code and outputs to accompany the article "Quantifying age specific household contacts for infectious disease modelling" by Sullivan et al.

A preprint of this article is available [here](https://arxiv.org/abs/2404.04300).
Results in the first version were generated using the commit tagged `v1.0`.
Results in the revised version were generated using the commit tagged `v2.1`.

All analyses were run in Matlab R2022b.


# Structure of this repository

The repository contains the following folders:
- `data` - the raw data that is used as model inputs and associated metdata (see [data](#data)).
- `functions` - Matlab functions that are called by the top-level script `main.m`.
- `results` - model outputs and associated metadata (see [results](#results)).


# How to use this repository

Open Matlab and run the top-level script `main.m` to reproduce the results that are in the article. This will read in raw data from the [data](#data) directory and save model outputs in the [results](#results) directory.

By default, `main.m` is set up to run a single realisation of the agent-based model for each combination of infection rate parameters (nReps = 1). You can change this by changing the value of nReps (set nReps = 100 to reproduce the results in the paper).  

You can specify the values of the household and non-household infection rate parameters to run in the vectors aRateHouse and aRateNonhouse.




# Data

The `data` folder contains the following files:
- pop_size_by_age_2018census_usually_resident.csv - 2018 census usually resident population estimates in 5-year age bands.
- Table_1_rounded.csv - the Stats NZ household composition data. Each row of this CSV file corresponds to one household type. The first 8 columns show the number of people in that household type in each age gorup, the column `Count` contains the number of households of that type, and the column `Count_people` contains the total number of people in households of that type. See Table 1a in the article for an example.
- JOB-12355 metadata.xslx - Stats NZ metadata and conditions of supply.
- nz_contacts_xxxx.xlsx (where xxxx is either home, school, work or other) - synthetic contact matrices projected onto the New Zealand population by [Prem et al. (2017)](https://doi.org/10.1371/journal.pcbi.1005697).


# Results

The `results` folder contains the following files:
- household_matrix_10.csv - contact matrix derived from the New Zealand household composition data in 10-year age bands.
- popsize.csv - population size derived from the New Zealand household composition data (after imputation).
- household_matrix_metadata.xlsx - metadata for the household contact matrix files.
- synpop.zip - compressed folder containing synpop.csv - the synthetic population derived from the New Zealand household composition data in 10-year age bands.
- synpop_metadata.xlsx - metadata for the synthetic population file.

It also contains the figure files that are written by main.m and appear in the article.



# Acknowledgements

This work includes customised Stats NZ data which are licensed by Stats NZ for re-use under the Creative Commons Attribution 4.0 International licence.
All data and other material produced by Stats NZ constitutes Crown copyright administered by Stats NZ. 
The authors are grateful to two anonymous reviewers for helpful comments on a previous version of the manuscript.



