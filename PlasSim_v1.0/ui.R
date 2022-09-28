###########################################################################################
### Simulation and Shiny Application of PlasSim (Plasmid population dynamics Simulator) ###
#
# © Martin Zwanzig, 2018
#
# model compartment abbreviations:
# F ... plasmid-free bacteria
# P ... non-adapted plasmid-bearing bacteria
# A ... adapted plasmid-bearing bacteria (either by chromosomal or plasmid mutations)
###########################################################################################

library(shiny)
library(deSolve)
library(RColorBrewer)

########################################################################
####### DEFINE UI FOR PARAMETERIZATION AND LAYOUT OF THE RESULTS #######
########################################################################

fluidPage(
  
  #  Application title
  headerPanel("PlasSim - Plasmid population dynamics Simulator"),
  
  wellPanel(
    fluidRow(
      tags$body("This web application is part of a publication:
                Zwanzig, M., Harrison, E., Brockhurst, M. A., Hall, J. P. J., Berendonk, T. U. & Berger, U. (2019): Mobile compensatory mutations promote plasmid survial. mSystems 4(1), https://doi.org/10.1128/mSystems.00186-18
                All cell densities and related processes in this model can be either defined in relative terms or by absolute values. We suggest to explore the model behavior as given with relative values. If you wish to use cells per ml, you have to make three adjustments: define a maximal attainable cell density (carrying capacity) in cells/ml, an initial density in cells/ml and a value for the conjugation rate that was for instance estimated by the endpoint-method (e.g. 10^-11)"
      ))),
  
  
  # Sidebar with sliders that enable the user to define parameter values
  fluidRow(
    
    column(4, wellPanel(
      
      #h3("Scales"),
      
      sliderInput("simTime", "simulation time [hours]:", 
                  min=1, max=10000, value=4000, step=1),
      
      numericInput("cc", "maximal attainable cell density (either 1 for relative considerations or in cells/ml)",
                   value=1),
      
      p(textOutput("cc_scinot")),
      
      numericInput("id", "initial cell density (as above either relative or in cells/ml)",
                   value = 1),
      
      p(textOutput("id_scinot")),
      
      sliderInput("iP", "initial proportion of (non-adapted) plasmid-bearing bacteria:", 
                  min=0, max=1, value=.5, step=0.01)
    )),
    
    wellPanel(
      h3("Parameterization"),
      fluidRow(
        column(4,
               sliderInput("gr", "maximum growth rate:", 
                           min=0, max=3, value=1, step=.01), 
               
               sliderInput("pb", "plasmid costs:", 
                           min=0, max=1, value=.2, step=.01), 
               
               sliderInput("br", "burden reduction by compensatory mutation:", 
                           min=0, max=1, value=.9, step=.01),
               
               sliderInput("mr", "dilution rate:", 
                           min=0, max=1, value=.1, step=.01)
        ), 
        
        column(4,
               numericInput("cr", "conjugation rate:",
                            min = 0, value = 0.02, step=.01),
               
               p(textOutput("cr_scinot")),
               
               numericInput("s", "segregation rate:", 
                            min=0, max=0.5, value=1e-03, step=1e-03), 
               
               p(textOutput("s_scinot")),
               
               numericInput("ar", "mutation rate:",
                            min = 0, max=1, value=1e-06, step=1e-06), 
               
               p(textOutput("ar_scinot")),
               
               numericInput("aa", "antibiotic action:", 
                            min=0, max=1, value=0, step=1e-03),
               
               p(textOutput("aa_scinot"))
        )
      )
    )
  ),
  
  wellPanel(
    fluidRow(
      column(3,
             checkboxInput("check_reldens", label = "Plot relative densities (cells / carrying capacity)", value = TRUE)
      ),
      column(3,
             checkboxInput("check_ylogscale", label = "Plot y-axis on log10 scale", value = TRUE)
      ),
      column(3,
             checkboxInput("check_fszoomin", label = "Zoom in fitness space plot", value = FALSE)
      ),
      column(3,
             downloadButton("downloadData", "Download simulation results") # Download the simulation results
      )
    )),
  
  # Show the simulation results of the parameterization entered
  mainPanel(
    plotOutput("graph", width = "1600px")),
  
  
  # Show copyright
  wellPanel(
    fluidRow(
      column(12,
             tags$a("Created by: Martin Zwanzig (martin.zwanzig@tu-dresden.de)", 
                    href="https://www.tu-dresden.de/forst/sysan/die-professur/staff/academic/martin-zwanzig")
      )))
)
