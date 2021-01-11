This is the code developed for our paper "[Public policy and economic dynamics of COVID-19 spread: a mathematical modeling study](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0244174)":

Uri Goldsztejn, David Schwartzman, Arye Nehorai\
Washington University in St. Louis, 2020

*If you find this code useful in your research, please consider citing.*

## Content
* Overview
* Files
* Further materials
* Citation format
<!--* Contact-->

## Overview

We have developed a predictive model for COVID-19 that considers, for the first time, its intercoupled effect on both economic and health (deaths and hospitalizations) outcomes for different quarantine policies.

We analyzed three different scenarios and concluded that the best policy to balance health and economic outcomes after an initial lockdown, would be to keep seniors in place and release non-seniors only gradually and after the first pandemic wave is over.

Furthermore, we performed a sensitivity analysis to study the effect of reducing the contagiousness of asymptomatic infected individuals. Our results suggest that public health measures that restrict disease spread outside of quarantine (e.g., social distancing, masks, avoiding crowded places) are important for policymakers to invest in.

## Files

* *code/baseline.m* - The code to simulate the baseline scenario described in our mansucript.

* *code/sudden_release.m* - The code to simulate the sudden release scenario described in our mansucript.

* *code/gradual_release.m* - The code to simulate the gradual release scenario described in our mansucript.

* *code/sensitivity_analysis.m* - The code to analyze the sensitivity of the parameters that regulate public policy.

* *code/diff_system.m* - The system of differential equations described in our manuscript.

* *code/diff_system_sensitivity.m* - The system of differential equations used for the sensitivity analysis.

* *code/initial_conditions_simulation.m* - The code used to seed the pandemic and determine the initial conditions.


## Further materials
A seminar about our work can be found [here](https://www.youtube.com/watch?v=a1qZjUVoe_E&t=1s)

A summary and a news article about our work can be found [here](https://www.ese.wustl.edu/~nehorai/research/Covid-19/Goldsztejn_Schwartzman_Nehorai_MedRxiv_2020.html)

## Citation format

Prefered reference format:\
Goldsztejn U, Schwartzman D, Nehorai A (2020) Public policy and economic dynamics of COVID-19 spread: A mathematical modeling study. PLoS ONE 15(12): e0244174. https://doi.org/10.1371/journal.pone.0244174
