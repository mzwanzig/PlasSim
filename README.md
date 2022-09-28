# PlasSim
Plasmid population dynamics simulator - a population-level model approach considering the amelioration of plasmid costs by chromosomal and plasmid mutations

The model ensemble comprises two sets of ordinary differential equations (ODE), which are structurally different in conjugation and adaptation:

## PLASMID ADAPTATION model 'PlasSim.PA'
* adaptations are generated through conjugation AND fission (mutations occur in relation to the frequency of plasmid replications)
* adaptations are inherited through conjugation AND fission (the adaptation is transfered with the plasmid)

## CHROMOSOMAL ADAPTATION model 'PlasSim.CA'
-> major difference to PlasSim.PA: transition and mutation of adaptation through cell fission only
* adaptations are generated through fission (mutations occur in relation to the frequency of chromosomal replications)
* adaptations are inherited through fission (the chromosmal adaptation is inherited to the daughter cells, but not to infected cells (conjugation))

A comparison of both models (modes of amelioration) is presented in Zwanzig et al. 2019 (see references below). This work demonstrated that plasmid-located compensatory evolution is better at enhancing plasmid persistence, even when its effects are smaller than those provided by chromosomal compensation.

## Implementation
The model ensemble is implemented in R, the free software environment for statistical computing and graphics, which is freely available from https://www.r-project.org/. It was created under R Version 3.4.2 and executed using aligned package 'deSolve' and 'rootSolve' for solving and analyzing the
steady state of ordinary differential equations. It was last tested for r version 4.2.1. It is recommended to run the code using the development plattform RStudio which is freely available under https://www.rstudio.com/

## Shiny app
Downloading and running the file 'launchPlasSim.R' with R allows you to run and explore the models. The file will execute 'ui.R' (defining the user interface) and 'server.R' (defining the calculations, including the model equations) stored in the subdirectory PlasSim_v1.0. To run 'launchPlasSim.R', an installation of the packages deSolve and shiny are required. The last is used to create a user interface with which model parameters can be easily changed with sliders and the resulting population dynamics be automatically viewed and downloaded as a file so that even non-experts in R are able to use the models with self-defined experimental conditions.

All cell densities and related processes in this model can be either defined in relative terms or by absolute values. We suggest to explore the model behavior as given with relative values. If you wish to use cells per ml, you have to make three adjustments: define a maximal attainable cell density (carrying capacity) in cells/ml, an initial density in cells/ml and a value for the conjugation rate that was for instance estimated by the endpoint-method (e.g. 10^-11)

## Reference
A full description of the model and its major results is given in this publication:
Zwanzig, M., Harrison, E., Brockhurst, M. A., Hall, J. P. J., Berendonk, T. U. & Berger, U. (2019): Mobile compensatory mutations promote plasmid survial. mSystems 4(1), https://doi.org/10.1128/mSystems.00186-18

Please indicate any use of this model that contributes to a publication with a reference to this article. If you plan to use or modify this model for your own research, please contact martin.zwanzig@tu-dresden.de to avoid working on the same or very similar research questions in parallel.
