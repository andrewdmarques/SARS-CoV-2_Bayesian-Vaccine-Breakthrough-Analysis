# SARS-CoV-2_Bayesian-Vaccine-Breakthrough-Analysis

This project includes R and Stan code to analyse point mutations in SARS-CoV-2 virus. The overall goal of the project is to understand the evolution of the virus and how different mutations appear and propagate over time.

## Project Structure

- `2021-09-23_mutations.xlsx`: This file contains the raw mutation data.
- `2021-09-25_genomeMetaData.xlsx`: This file contains metadata about the genomes that were analyzed.
- `2021-10-14_vbt-mutation-enrichment_v5.32.R`: This is the main analysis script.
- `functions.R`: This script contains various helper functions that are used by the other scripts.
- `model_mutations.rds`: This is a saved Stan model.
- `model_mutations.stan`: This is the Stan model script used for Bayesian inference.

## Instructions

1. Ensure you have the necessary R packages installed. You can install them using the following command:


install.packages(c("tidyverse", "reshape2", "ggplot2", "rstan"))

2. Run the `2021-10-14_vbt-mutation-enrichment_v5.32.R` script for the main analysis.


# Running the analysis.
source('2021-10-14_vbt-mutation-enrichment_v5.32.R')


The results of the analysis will be saved in the 'results' directory.

## Additional Information

Please make sure that you have all these files and packages installed before running the scripts. The scripts are set up to run with relative paths from the main project directory.

