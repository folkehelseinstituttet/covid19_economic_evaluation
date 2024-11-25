# Optimal intervention strategies for a new pandemic: The case of SARS-CoV-2

This repository contains all the code and data to reproduce the results in the paper "Optimal intervention strategies for a new pandemic:  The case of SARS-CoV-2". We implement an infectious disease model that describes the transmission of a new COVID-19 variant in Norway in a pipeline that based on the model output can estimate the total cost of epidemic including the cost of interventions. Using this model we perform a numerical optimisation to find the government intervention strategy that minimises the total cost. 
## Structure of the code

The transmission model itself is implemented in the (metapop)[https://github.com/Gulfa/metapop] package. 
- run_analysys_sbatch_array.R - Running the main analysis as an sbatch array
- run_functions.R - Includes most of the code for setting up and running the model, calculating costs
- optimise_beta.R - Code for optimising beta(t)
- QALY.R - Code for estimating QALY loss, production loss and loss of exeeding hospitalisation capacity
- paramerters.r - Includes key parameters
- plot_utils.R - Plotting utilty functions
- plot_paper.R - Plotting the figures in the paper based on the optimal scenarios
- plot_individual_runs.R - Plotting figures
- install_reqs.R - Install the required R-packages
- parameter_files/ - Folder for key transmission parameters
