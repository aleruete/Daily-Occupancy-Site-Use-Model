Exploiting Opportunistic Observations with Seasonal Site Use Models (Daily Occupancy Models)
=======================================================================================================================
Supplementary data and scripts

This repository accompanies the text: Ruete et al. 2017. Exploiting Opportunistic Observations to Estimate Changes in Seasonal Site Use: An Example with Wetland Birds. JX X:XX-XX

Here you can download all data and scripts needed to replicate the analyses done in the article. 
To run it you will need to install R <http://www.r-project.org/> and JAGS <http://mcmc-jags.sourceforge.net/>.

### NOTE: 
To run these scripts you will need at least 16GB of RAM, and the code is written so that models are fitted in 4 independent chains, therefore you will need at least 4 cores (you can change that if you want, though).
Be warned, model fitting takes about 2 hours per species (CPU i7 @3.5Ghz).

## Content
Original Observation data: Uppland90180.csv
To collate the data: Jags Compile data Optimized 20160125.R (returns OccDataUppland.RData)
The JAGS model: Write Jags Occupancy Model Saturation 20150815.R
To fit the models: Jags Occupancy Model Parallel Saturation 20150815.R 
+ other files to plot the results

#/Results
a folder to store the results. Example of results for the first species

#/Simulations
All scripts required to created the simulated data and fit the model to it and by products. 
Running order:
1. Scenarios Simulate Occ Data.R
2. Scenarios Fit Models.R
3. Scenarios Check Model Fit.R (or Scenarios Check Model Fit 3x3.R for a shorter version)



### Authors and Contributors
Scripts and Repository created by Alejandro Ruete in Apr 2017.
DOI: dx.doi.org/
### Licence GNU v.3