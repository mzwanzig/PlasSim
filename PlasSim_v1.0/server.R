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

##########################################################################################################
####### DEFINE SERVER LOGIC TO PERFROM SIMULATIONS USING THE ORDINARY DIFFERENTIAL EQUATIONS ABOVE #######
##########################################################################################################

function(input, output) {
  
  simulation <- reactive({
    # Model Parameters:
    # Initialization with parameter values of input sliders
    
    simTime <- input$simTime  # Simulation time in hours
    cc <- input$cc            # carrying capacity (defaults to 1, if unkown)
    id <- input$id            # initial total cell density
    iP <- input$iP            # initial proportion of plasmid-bearers
    gr <- input$gr            # maximum growth rate
    pb <- input$pb            # plasmid burden
    br <- input$br            # burden reduction (for adapted cells)
    cr <- input$cr            # conjugation rate
    s <- input$s              # segregation rate
    ar <- input$ar            # adaptation rate; according to the given exponent for the desired adaptation rate
    mr <- input$mr            # mortality rate
    aa <- input$aa            # antibiotic action (addition to mortility rate)
    
    reldens <- input$check_reldens     # show results as relative densities? [boolean variable]
    ylogscale <- input$check_ylogscale # show results on log10 scaled y-axis? [boolean variable]
    fszoomin <- input$check_fszoomin   # zoom in plasmid fitness space plot? [boolean variable]
    
    default <- FALSE
    if(default){
      simTime <- 500; cc <- 1; id <- 1; iP <- .5; gr <- 1; pb <- .2
      br <- .9; ec <- 1; s <- 0.5^10; ar <- 1e-03; mr <- .1; aa <- 0
      reldens <- TRUE; ylogscale <- TRUE
    }
    
    # merge single parameters to the parameter set used for simulation
    parms.def <- c(mr = mr, gr = gr, pb = pb, cr = cr, br = br, s = s, ar = ar, aa = aa, cc = cc)
    parms.na <- parms.def; parms.na["ar"] <- 0 # the adaptation rate is set to 0 for no adaptation mode
    
    # defining initial state
    iP <- iP * as.numeric(id)
    iF <- (as.numeric(id) - iP)
    iA <- 0
    y0 <- c(F = iF, P = iP, A = iA)
    
    # define simulation period and temporal discretization
    times <- seq(from = 0, to = simTime, by = 1) # numerical integration in 1h steps until simTime is reached
    
    # SET OF ORDINARY DIFFERENTIAL EQUATIONS  
    # ---------------------------------------------------------
    # PLASMID DYNAMICS CONSIDERING ADAPTATION AND BACTERICIDAL ANTIBIOTIC ACTION
    
    # The follwing two models PlasSim.PA and PlasSim.CA are structurally different in conjugation and adaptation:
    
    # PLASMID ADAPTATION:
    # * adaptations are generated through conjugation AND fission (mutations occur in relation to the frequency of plasmid replications)
    # * adaptations are inherited through conjugation AND fission (the adaptation is transfered with the plasmid)
    PlasSim.PA <- function(t,y,p) {
      with(as.list(c(y, p)),{
        f <- 1 - ((F + P + A) / cc)
        #        GROWTH                SEGREGATION                                CONJUGATION       ADAPTATION                        MORTALITY/WASHOUT
        dF  <- f*(gr*F                 + s*(gr*(1-pb))*P + s*(gr*(1-pb*(1-br)))*A - cr*F*P - cr*F*A                                 ) - (mr+aa)*F
        dP  <- f*((gr*(1-pb))*P        - s*(gr*(1-pb))*P                          + cr*F*P          - ar*cr*F*P  - ar*(gr*(1-pb))*P ) - mr*P
        dA  <- f*((gr*(1-pb*(1-br)))*A - s*(gr*(1-pb*(1-br)))*A	                  + cr*F*A          + ar*cr*F*P  + ar*(gr*(1-pb))*P ) - mr*A
        list(c(dF, dP, dA),
             "PA prop." = (P+A) / (F+P+A)
        )
      })
    }
    
    # CHROMOSOMAL ADAPTATION: transition and mutation of adaptation through fission
    # * adaptations are generated through fission (mutations occur in relation to the frequency of chromosomal replications)
    # * adaptations are inherited through fission (the chromosmal adaptation is inherited to the daughter cells, but not to infected cells (conjugation))
    PlasSim.CA <- function(t,y,p) {
      with(as.list(c(y, p)),{
        f <- 1 - ((F + P + A) / cc)
        #        GROWTH                SEGREGATION                                CONJUGATION       ADAPTATION           MORTALITY/WASHOUT
        dF  <- f*(gr*F                 + s*(gr*(1-pb))*P + s*(gr*(1-pb*(1-br)))*A - cr*F*P - cr*F*A                    ) - (mr+aa)*F
        dP  <- f*((gr*(1-pb))*P        - s*(gr*(1-pb))*P                          + cr*F*P + cr*F*A - ar*(gr*(1-pb))*P ) - mr*P
        dA  <- f*((gr*(1-pb*(1-br)))*A - s*(gr*(1-pb*(1-br)))*A	                                    + ar*(gr*(1-pb))*P ) - mr*A
        
        list(c(dF, dP, dA),
             "PA prop." = (P+A) / (F+P+A)
        )
      })
    }
    
    ###################
    # RUN SIMULATIONS #
    ###################
    
    out.NA <- ode(y0, times, PlasSim.PA, parms.na)  # simulation without adaptation (using PA or CA makes no difference here, since ar = 0)
    out.CA <- ode(y0, times, PlasSim.CA, parms.def) # simulation considering chromosomal compensatory adaptation
    out.PA <- ode(y0, times, PlasSim.PA, parms.def) # simulation considering plasmid compensatory adaptation
    
    SimResults <- data.frame(NA.time = out.NA[,1], NA.F = out.NA[,2], NA.P = out.NA[,3], NA.A = out.NA[,4],
                             CA.time = out.CA[,1], CA.F = out.CA[,2], CA.P = out.CA[,3], CA.A = out.CA[,4],
                             PA.time = out.PA[,1], PA.F = out.PA[,2], PA.P = out.PA[,3], PA.A = out.PA[,4])
    
    return(SimResults)
    
  })
  
  output$graph <- renderPlot({
    
    simres <- simulation()
    out.NA <- simres[,1:4]
    out.CA <- simres[,5:8]
    out.PA <- simres[,9:12]
    
    simTime <- input$simTime  # Simulation time in hours
    cc <- input$cc            # carrying capacity (defaults to 1, if unkown)
    id <- input$id            # initial total cell density
    iP <- input$iP            # initial proportion of plasmid-bearers
    gr <- input$gr            # maximum growth rate
    pb <- input$pb            # plasmid burden
    br <- input$br            # burden reduction (for adapted cells)
    cr <- input$cr            # conjugation rate
    s <- input$s              # segregation rate
    ar <- input$ar            # adaptation rate; according to the given exponent for the desired adaptation rate
    mr <- input$mr            # mortality rate
    aa <- input$aa            # antibiotic action (addition to mortility rate)
    
    reldens <- input$check_reldens     # show results as relative densities? [boolean variable]
    ylogscale <- input$check_ylogscale # show results on log10 scaled y-axis? [boolean variable]
    fszoomin <- input$check_fszoomin   # zoom in plasmid fitness space plot? [boolean variable]
    
    FPA.col <- brewer.pal(5, "Set1")[c(2,1,5)]
    a.col <- .9 # alpha factor for FPA.col transparency
    cex.ef <- 1.2 #character size expansion factor
    hcol <- grey(.1)
    ca.col <- "purple"
    pa.col <- "plum"
    
    pit <- function(x, title.name){
      
      X <- x[,2:4]; adaptive.ylab <- "density [cells / ml]"; adaptive.ylim <- c(0, cc); adaptive.yaxt <- NULL
      
      if(reldens){
        X <- X / cc
        adaptive.ylab <- "density [cells / max. cell density]"
        adaptive.ylim <- c(0, 1)
        if(ylogscale){adaptive.ylim <- c(-9, 0)}
        adaptive.yaxt <- "n"
      }
      if(ylogscale){
        X <- log10(X)
        adaptive.ylab <- paste(adaptive.ylab,"(log10 scale)")
        if(reldens == FALSE){adaptive.ylim <- c(0, log10(cc))}
      }
      
      matplot(X, type = "l", lty = 1, lwd = 3, cex.lab = cex.ef, cex.axis = cex.ef,
              ylim = adaptive.ylim, yaxt = adaptive.yaxt,
              col = scales::alpha(FPA.col, a.col), las = 1,
              main = title.name, xlab = "time [hours]", ylab = adaptive.ylab)
      if(reldens){
        if(ylogscale){
          axis(2, at = seq(-9,0,3), labels = c(expression(10^-9),expression(10^-6),expression(10^-3),expression(1)), las = 2)
        }else{
          axis(2, at = seq(0,1,0.2), seq(0,1,0.2), las = 2)
        }
      }
    }
    
    fitnessspaceplot <- function(){
      
      pht <- (cr*cc) / gr
      pvt_P <- ((gr * (1 - pb)) - (gr * (1 - pb) * s)) / gr
      pvt_A <- ((gr * (1 - pb * (1 - br))) - (gr * (1 - pb * (1 - br)) * s)) / gr
      
      pht.v <- c(0, pht, 0, pht)
      pvt.v <- c(1, pvt_P, pvt_A, pvt_A)
      
      if(fszoomin){
        my.xlim <- c(min(pht.v) - 0.1, max(pht.v) + 0.1)
        my.ylim <- c(min(pvt.v) - 0.1, max(pvt.v) + 0.1)
        ex <- 0.05 * (max(c(abs(diff(range(pvt.v))),abs(diff(range(pht.v))))) / 1.2)
      }else{
        my.xlim <- c(-0.1,1.1)
        my.ylim <- c(-0.1,1.1)
        ex <- 0.04
      }
      
      plot(cbind(pht.v, pvt.v), cex.lab = cex.ef, cex.axis = cex.ef,
           pch = 19, cex = 1.5, col = scales::alpha(FPA.col[c(1:3,3)], a.col), pty = "s",
           xlab = "HORIZONTAL TRANSMISSION FITNESS",
           ylab = "VERTICAL TRANSMISSION FITNESS",
           xlim = my.xlim, ylim = my.ylim, xaxs = "i", yaxs = "i", las = 1,
           main = "PLASMID FITNESS SPACE", bty = "n")
      lines(c(0,1),c(1,0), col = grey(.5, .9), lwd = 0.5, lty = 1)  #power-isoline
      lines(c(0,1,1,0,0), c(1,1,0,0,1), col = grey(.9, .9), lwd = 0.5, lty = 1)
      box()
      text(x = c(0 - ex, pht + ex, 0 - ex, pht + ex), y = c(1 + ex, pvt_P - ex, pvt_A - ex, pvt_A + ex),
           labels = expression(F,P,A[C],A[P]), col = scales::alpha(FPA.col[c(1:3,3)], a.col), cex = 2)
    }
    
    ode_term_plot <- function(){
      # GROWTH TERM
      barplot(c(gr, gr * (1 - pb), gr * (1 - pb * (1- br))),
              main = "GROWTH", names.arg = expression(F,P,A),
              col = FPA.col, border = NA, las = 1, cex.lab = cex.ef, cex.axis = cex.ef)
      abline(h = 0, col = hcol)
      # SEGREGATION TERM
      barplot(c(s*(gr*(1-pb))+s*(gr*(1-pb*(1-br))), -s*(gr*(1-pb)), -s*(gr*(1-pb*(1-br)))),
              main = "SEGREGATION", names.arg = expression(F,P,A),
              col = FPA.col[c(3,2,3)], border = NA, las = 1, cex.lab = cex.ef, cex.axis = cex.ef)
      barplot(c(s*(gr*(1-pb)), 0, 0),
              add = TRUE, col = c(FPA.col[2], NA, NA), border = NA, yaxt="n", xaxt="n")
      abline(h = 0, col = hcol)
      # CONJUGATION TERM
      barplot(c(-cr-cr, +cr+cr, +cr),
              main = "CONJUGATION", names.arg = expression(F,P,A),
              col = c(FPA.col[3], ca.col, pa.col), border = NA, las = 1, cex.lab = cex.ef, cex.axis = cex.ef)
      barplot(c(-cr, +cr, 0),
              add = TRUE, col = c(FPA.col[c(2,2)], NA), border = NA, yaxt="n", xaxt="n")
      barplot(c(-cr-cr, +cr+cr, +cr), add = TRUE, border = NA, 
              col = c(FPA.col[c(1,1,1)]), density = 25, angle = 45, yaxt="n", xaxt="n")
      abline(h = 0, col = hcol)
      # ADAPTATION TERM
      barplot(c(0, -ar*cr-ar*(gr*(1-pb)), +ar*cr+ar*(gr*(1-pb))),
              main = "MUTATION", names.arg = expression(F,P,A),
              col = c(1,pa.col,pa.col), border = NA, las = 1, cex.lab = cex.ef, cex.axis = cex.ef)
      barplot(c(0, -ar*(gr*(1-pb)), +ar*(gr*(1-pb))),
              add = TRUE, col = FPA.col[c(2,2,2)], border = NA, yaxt="n", xaxt="n")
      abline(h = 0, col = hcol)
      # MORTALITY TERM
      barplot(c(-aa-mr, -mr, -mr),
              main = "MORTALITY", names.arg = expression(F,P,A),
              col = FPA.col, border = NA, las = 1, cex.lab = cex.ef, cex.axis = cex.ef)
      abline(h = 0, col = hcol)
    }
    
    layout(matrix(c(1:9),1,9), width = c(2,2,2,1.5,rep(0.5,5)))
    par(pty = "m", oma = c(6.5, 1, 3, 1))
    pit(out.NA, "NO MUTATION")
    legend("bottomleft", legend = expression(F,P), bty = "n", col = scales::alpha(FPA.col[1:2], a.col), lty = 1, lwd = 3)
    pit(out.CA, "CHROMOSOMAL MUTATION")
    legend("bottomleft", legend = expression(F,P,A[C]), bty = "n", col = scales::alpha(FPA.col, a.col), lty = 1, lwd = 3)
    pit(out.PA, "PLASMID MUTATION")
    legend("bottomleft", legend = expression(F,P,A[P]), bty = "n", col = scales::alpha(FPA.col, a.col), lty = 1, lwd = 3)
    par(pty = "s"); fitnessspaceplot()
    par(pty = "m"); ode_term_plot()
    
    #define legend function
    add_legend <- function(...) {
      opar <- par(fig=c(0, 1, 0, 1), oma=c(0, 0, 0, 0), 
                  mar=c(0, 0, 0, 0), new=TRUE)
      on.exit(par(opar))
      plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n')
      legend(...)
    }
    
    add_legend("bottomright", legend = c("F", "P", "A", "(A) only for chromosomal mutation", "(A) only for plasmid mutation", "second order reaction involving F cells"),
               title = "Process rates (y-axis) affecting compartments F, P, A (x-axis); colors indicate the causative compartments:",
               fill = c(FPA.col, ca.col, pa.col, FPA.col[1]), bty = "n", title.adj = 0, cex = cex.ef,
               border = NA, density = c(NA,NA,NA,NA,NA,25), angle = c(NA,NA,NA,NA,NA,25))
    
    add_legend("bottomleft", legend = expression(paste(F," plasmid-free cells"),
                                                 paste(P," non-adapted plasmid-bearing cells"),
                                                 paste(A[C]," plasmid-bearing cells with chromosomal mutation"),
                                                 paste(A[P]," plasmid-bearing cells with plasmid mutation")),
               col = scales::alpha(FPA.col[c(1:3,3)], a.col), bty = "n", title.adj = 0, cex = cex.ef, lwd = 3, lty = 1)
    
    mtext("SINGLE REACTION RATES", side = 3, outer = TRUE, adj = 1, cex = 2)
    mtext("POPULATION DYNAMICS", side = 3, outer = TRUE, adj = 0, cex = 2)
  })
  
  output$cc_scinot <- renderText({
    formatC(input$cc, format = "e", digits = 1)
  })
  
  output$id_scinot <- renderText({
    formatC(input$id, format = "e", digits = 1)
  })
  
  output$cr_scinot <- renderText({
    formatC(input$cr, format = "e", digits = 1)
  })
  
  output$s_scinot <- renderText({
    formatC(input$s, format = "e", digits = 1)
  })
  
  output$ar_scinot <- renderText({
    formatC(input$ar, format = "e", digits = 1)
  })
  
  output$aa_scinot <- renderText({
    formatC(input$aa, format = "e", digits = 1)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste('PlasSim_data-', Sys.Date(), '.csv', sep='')
    },
    content = function(file) {
      write.csv(simulation(), file)
    }
  )
}