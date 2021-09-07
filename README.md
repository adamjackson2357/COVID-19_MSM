## Introduction

This documentation details options and features implemented in QOMiC and in the companion simulation tool QOMiC_simul.

This is an adaptation of QOMiC: Quasi-observed Markov Chain; An MCMC sampler for smoking-induced lung cancer [1]

The QOMiC tool has been adapted into a multi-state discrete-time model for COVID-19. The report associated with this model can be found in the documentation.

## Methodology

Metropolis Hastings (all parameters are samples from a multivariate proposal)

- QOMiC will estimate the parameters of the model using a Metropolis-Hastings algorithm, assuming a uniform prior distribution and a random walk multivariate proposal distribution. 

Adaptive procedure to automatically calibrate the variances of the proposal distibution [2]

- In order to ensure a smooth estimation of the tail of the posterior distribution, the variance of the proposal distribution must be tuned such that the acceptance rate lies between 20% and 60%

Simulation

- Of individual trajectories between each of the four states based on the joint posterior distribution of the parameters
- In practice, it samples sets of parameters from the joint posterior distribution, computes the transition probabilities for each participant at each calendar year, and samples transitions from this set of probabilities.

## Installation

1. QOMiC is a C++ program that uses packages from the GNU Scientific Library (GSL), available at:

    ftp://ftp.gnu.org/gnu/gsl/

    For Linux users, the installation is straightforward and documented in the downloaded archive.

    QOMiC has been developed and tested using the g++ compiler on several Linux platforms. Compilation using alternative compiler may require modifications to the Makefile.

2. Clone the project from [git_lab_link]

3. Create the msm environment:

    module load anconda3/personal
    conda env create -f configs/msm.yml

4. Using the qomic_demo as an example, compile the MCMC scripts and simulation scripts using:

    make -C qomic_simul/main
    
    make -C qomic_simul/QOMiC_simul
    
5. Run the code by navigating to the scripts directory, and running ./run_mcmc_demo.sh, taking care to change any filepaths

## Contents

### documentation

Contains the associated report.

### scripts

Contains bash scripts to run the MCMC and simulations for each qomic.

### configs

Contains:

- input data for the MCMC
- yaml file to create an anaconda environment

### qomic

The code has yet to be generalised, so has been split using the sequential adjustment approach detailed in the COVID-19 MSM paper:

- qomic_demo: model with only demographic parameters
- qomic_social: demographic + social
- qomic_health: demographic + social + health
- qomic_medical: demographic + social + health + medical
- qomic_time: demographic + social + health + medical + time

Each qomic folder contains all C++ codes:

- main: This directory contains the QOMiC source code and runs the MCMC. The main program is QOMiC.cc. In most files, a variable called DEBUG is defined and set to 0 by default. If the user changes it to 1 and recompiles the program, step-by-step details of the performed operations will be printed on the standard output device (or log file). This option is useful to understand how QOMiC works, but usually generates very large log files.
- QOMiC_simul: This directory contains the simulation tool source code. The main program is QOMIc_simul.cc.
- Routines: This directory contains the main routines to extract information from the input files, to calculate the likelihood, to perform the MCMC estimation procedure, and to output the results.
- Classes: This directory contains the classes (here only matrices) used in the code.

### visualisations

R codes to create visualisations for MCMC and simulation diagnostics.

## Running QOMiC

### Run configuration

The bash scripts are currently configured for use in an HPC environment.

The bash scripts run the:

- MCMC
- diagnostic scripts
- simulation

Steps: 

1. Change file paths and input names
2. Set run parameters
    - name: run name should correspond to the name of the folder in the configs directory
    - n_iter: number of iterations for the MCMC
    - burn_in: MCMC burn-in iterations
    - adaptive: set as 1 for adaptive iterations
    - iter_adaptive: the adaptive memory iterations
    - update_frequency: the update frequency for the adaptive
    - n_params: number of parameters to include in the model
    - n_trans: number of transitions in the model
    - simul_iter: number of iterations for the simulation
    - nb_iter: n_iter - burn_in
    - type: model name corresponds to the variables included
    - seed: set random seed

### Input files

Matrix structures are defined as C++ objects in the Classes directory: the first two elements read from the matrix file are the number of rows and columns, respectively.

**Date file (-dates file_name option)**
This input file contains, for each participant, the day that a participant entered a given state, such that each state is a column. Days from the start are formatted as an integer (with 0 being the start of the time space). If the participant never entered a state, then set as -1. The four columns correspond to:

1.  O, Outpatient COVID-19 positive
2.  I, Inpatient COVID-19 positive
3.  R, Recovered from COVID-19
4.  D, Death from COVID-19;

**Status file (-status file_name option)**
A binarised version of the dates file

**Covars file (-covars file_name option)**
Time independent covariates for each participant

**Hospital file (-hospital file_name option)**
qomic_time only. Contains the UK-wide number of COVID-19 hospitalisations at time (t-1).

**Trans file (-trans file_name option)**
A transition matrix indicating which transitions occurred. Rows = from, columns = to.

**theta file (-theta file_name option)**
Vector of starting parameters

**sigma file (-sigma file_name option)**
Starting covariance matrix

### Output files

Results from the **MCMC** are automatically outputed to results/mcmc:

- History file: compiling the retained values for each parameter at each iteration, as well as the value of the log-likelihood estimated at each parameter update (i.e., at each step of the Metropolis-Hastings algorithm);
- Posterior file: giving, for several values of the burn-in, the estimates of the posterior mean of each parameter, as well as the log-likelihood obtained fixing all parameters to these posterior means.
- Log file: information about the run (including estimates of the acceptance rates)
- Figures: diagnostic visualisations:
    - Trace plots
    - Posterior distribution
    - Scatter plots for values of each parameter combination
    - Correlation plots of all parameters
    
Results from the **QOMiC_simul** are automatically outputed to results/probabilities:

- P_: Average probability matrices (individual-day) for each transition
- Transitions: The average number of transitions for each transition
- Conditional: The average number of days spent in a state before transitioning to a new state
- log file
    

## References

[1] Chadeau-Hyam,  M;Tubert-Bitter, P.   Dynamics  of  the  risk  of  smoking-induced  lung  cancer:a  compartmental  hidden  Markov  model  for  longitudinal  analysis.Epidemiology,  25(1):28–34,2014.

[2] Heikki  Haario,  Eero  Saksman,  and  Johanna  Tamminen.    Adaptive  proposal  distribution  forrandom walk Metropolis algorithm.Computational Statistics, 14(3):375–395, 1999.