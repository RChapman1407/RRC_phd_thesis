This folder contains MATLAB scripts to produce the heatmap of AMOC collapse probabilities in Chapter 5

HadGEM3MM_run_parameters_multiples_decadal_all.m: Contains parameter file to load in with noise and forcing scenarios for each scenario in the heatmap
HadGEM3MM_MC_tipping_probability_auto.m: Main script to calculate probability of collapse for different noise/ forcing scenarios using Monte Carlo method. Each scenario run for 1000 years for 1000 noise realisations.
Heatmap_tipping_probability.m: Script to plot heatmap fig from output of above files.

tipping_probabilities_decadal_all.txt: Output of HadGEM3MM_MC_tipping_probability_auto.m used to plot the heatmap.
