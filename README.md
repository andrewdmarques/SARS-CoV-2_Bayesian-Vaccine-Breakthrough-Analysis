# SARS-CoV-2 Bayesian Vaccine Breakthrough Analysis

View our peer-review publication [SARS-CoV-2 Variants Associated with Vaccine Breakthrough in the Delaware Valley through Summer 2021](https://pubmed.ncbi.nlm.nih.gov/35130727/).

The project includes R and Stan code to analyse point mutations in SARS-CoV-2 virus. The overall goal of the project is to understand the evolution of the virus and how different mutations appear and propagate over time.



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

### Data Block

```stan
data {
  int<lower=1> nTime;
  int counts[nTime];
  int vaccine[nTime];
  int drop[nTime];
  int nCounts[nTime];
  int nVaccine[nTime];
  int nDrop[nTime];
}
```
This block includes the variables that represent the data that are input into the model. The variables include counts (total counts over time), vaccine (vaccine rates over time), and drop (drop rates over time) with their respective numbers (nCounts, nVaccine, nDrop). 'nTime' is the total number of time periods.

### Parameters Block

```stan
parameters{
  real meanWeek1;
  vector[nTime-1] changes;
  real<lower=0> changeSigma;
  real vaccineChange;
  real dropChange;
}
```
This block includes the parameters that the model will estimate. It includes the mean value for the first week (`meanWeek1`), a vector of changes (`changes`), the standard deviation of the changes (`changeSigma`), the rate of change for the vaccine (`vaccineChange`), and the rate of change for the drop (`dropChange`).

### Transformed Parameters Block

```stan
transformed parameters{
  vector[nTime] means;
  means[1]=meanWeek1;
  for(ii in 2:nTime) means[ii]=means[ii-1]+changes[ii-1]*changeSigma;
}
```
This block includes the transformed parameters. It contains a vector `means` which holds the means for each time point, calculated based on the previous mean and the corresponding change times the change standard deviation (`changeSigma`).

### Model 

```stan
model {
  meanWeek1 ~ normal(0,10);
  changes[1] ~ normal(0,1);
  for(ii in 2:(nTime-1)) changes[ii]~normal(changes[ii-1],1);
  counts ~ binomial_logit(nCounts,means);
  vaccine ~ binomial_logit(nVaccine,means+vaccineChange);
  drop ~ binomial_logit(nDrop,means+dropChange);
  changeSigma~gamma(1,2);
  vaccineChange~double_exponential(0,1);
  dropChange~double_exponential(0,1);
}
```
This block includes the model specifications, including prior distributions for the parameters and the likelihood functions for the observations.


First, the prior distributions for the parameters are defined as follows:

```stan
meanWeek1 ~ Normal(0, 10)
changes[1] ~ Normal(0, 1)
changes[ii] ~ Normal(changes[ii-1], 1) for ii in 2:(nTime-1)
changeSigma ~ Gamma(1, 2)
vaccineChange ~ Laplace(0, 1)
dropChange ~ Laplace(0, 1)
```

Second, the transformed parameters are computed as:

```stan
means[1] = meanWeek1
means[ii] = means[ii-1] + changes[ii-1]*changeSigma for ii in 2:nTime
```

Finally, the likelihood functions for the data are defined as:

```stan
counts ~ Binomial_Logit(nCounts, means)
vaccine ~ Binomial_Logit(nVaccine, means + vaccineChange)
drop ~ Binomial_Logit(nDrop, means + dropChange)
```

In this model, `Binomial_Logit(n, x)` represents a binomial distribution with `n` trials and log-odds of success `x`. `Normal(mu, sigma)` represents a normal distribution with mean `mu` and standard deviation `sigma`. `Gamma(alpha, beta)` represents a gamma distribution with shape parameter `alpha` and rate parameter `beta`. `Laplace(mu, b)` represents a Laplace (or double exponential) distribution with location `mu` and scale `b`.

This Bayesian model thus uses the prior distributions to inform the parameters' values, which are then used to compute the likelihood of observing the given data. This information is combined in a posterior distribution that provides updated estimates of the parameters given the data.

## Additional Information

Please make sure that you have all these files and packages installed before running the scripts. The scripts are set up to run with relative paths from the main project directory.

