# Selecting pharmacies for COVID-19 testing
This repository contains code for a facility location optimization problem and helper functions to plot and prepare input data. Although the code is general, its intent was to evaluate the access to COVID-19 tests if they were made available through pharmacies in the US. We present this study in the paper *Selecting pharmacies for COVID-19 testing to ensure access*, which is openly available as a preprint on [medRxiv](https://www.medrxiv.org/content/10.1101/2020.09.17.20185090v1). Please see the paper for further information about the methods and assumptions.

Although the pharmacy data used in the study is not publicly available, we provide alternative sources in the *Data* section below. The repository contains a small fictional example with three states, state ONE, TWO, and THREE, and three pharmacies, Pharma A, B, and C. Independent pharmacies are labeled as *Other*. Simply run ```run_national.jl``` and ```run_states.jl``` to solve the fictional example. The ```results``` folder contains results from the solved fictional example. For results on realistic and larger case studies, please see *Selecting pharmacies for COVID-19 testing to ensure access* on [medRxiv](https://www.medrxiv.org/content/10.1101/2020.09.17.20185090v1) and particularly its online supplement.

To run your own study, simply replicate the format of the ```areas.csv``` and ```pharmacies.csv``` files with your own data and change the ```STATES``` and ```INCLUDED_STATES``` variables in ```run_national.jl``` and ```run_states.jl``` to contain the states in your dataset and run ```run_national.jl``` and ```run_states.jl```.

## Repository content
The repository contains the following files and folders:
* ```data``` contains ```areas.csv``` and ```pharmacies.csv``` of a fictional example.
* ```model``` contains the facility location problem in ```model.jl```, helper fuctions in ```utils.jl```, and functions to plot results in ```plotting.jl```.
* ```results``` contains results as ```csv``` files and plots for the fictional example.
* ```run_national.jl``` runs code for the national case, i.e. selects pharamcies to maximize access in all states.
* ```run_states.jl``` runs code that selects pharmacies to maximize access for individual states.

## Requirements to run code
The code is in the Julia programming language, using the following open source packages (versions in parentheses and Julia v1.4.2 were used to run ```run_national.jl``` and ```run_states.jl``` to get the results in ```results```):
* [JuMP](https://jump.dev/JuMP.jl/stable/) (v0.21.3) for mathematical optimization.
* [JSON](https://github.com/JuliaIO/JSON.jl) (v0.21.0) for parsing and printing ```JSON``` files.
* [CSV](https://juliadata.github.io/CSV.jl/stable/) (v0.7.5) for utilities for working with ```CSV``` files.
* [DataFrames](http://juliadata.github.io/DataFrames.jl/stable/) (v0.21.3) for tabular data manipulation.
* [DataStructures](https://juliacollections.github.io/DataStructures.jl/latest/) (v0.17.19) for the ```PriorityQueue``` data structure.
* [Plots](http://docs.juliaplots.org/latest/) (v1.5.7) to plot results.
* A solver. The code currently uses the open source [GLPK](http://www.gnu.org/software/glpk/) solver (v0.13.0). Other solvers are possible as long as they are compatible with JuMP. The user must remember to change GLPK to the solver of their choice in ```model.jl```.

## Data
The code is general, and users can run their own code by using their own ```areas.csv``` and ```pharmacies.csv``` files. Note that altough we use pharmacies in our study, one can estimate access and travel distances to other facilties than pharmacies.

```areas.csv``` requires the following columns:
* ZCTA: ZIP Code Tabulation Areas (ZCTAs), i.e. generalized areal representations of ZIP codes.
* POP: Population at ZCTA.
* LANDAREA: The land area of a ZCTA.
* LAT: Latitude of ZCTA centroid.
* LON: Longitude of ZCTA centroid.
* STATE: State that ZCTA belong to.

```pharmacies.csv``` requires the following columns:
* STATE: State that pharmacy belong to.
* ZCTA: The ZCTA where the pharmacy is located.
* CHAIN: The chain the pharmacy belong to. Note that all pharmacies not belonging to a chain are labeled by *Other*.

Several openly accessible data sources are available. [The US Census](https://www.census.gov/) provides population and land area at ZCTA detail for the US. For pharmacies, [Rx Open](http://rxopen.org/) makes operating statuses of healthcare facilities available for the whole US during the COVID-19 pandemic. Although the Rx Open list is quite comprehensive, we note that they only include facilities enrolled in Rx Open.
