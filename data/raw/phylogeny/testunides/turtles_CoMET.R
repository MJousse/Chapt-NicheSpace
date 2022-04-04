library(TESS)
setwd("~/Desktop/TESS/")
turtles <- read.tree("data/tess_input_tree.tre")

numExpectedRateChanges <-1
samplingFraction <- 0.85

tess.analysis(turtles,
              empiricalHyperPriors = TRUE,
              samplingProbability = samplingFraction,
              numExpectedRateChanges = numExpectedRateChanges,
              estimateNumberMassExtinctions = FALSE,
              MAX_ITERATIONS = 2000000,
              MIN_ESS = 500,
              priorOnly = FALSE,
              dir = "comet_auto_stop_3")


output <- tess.process.output("comet_auto_stop_3",
                              numExpectedRateChanges = numExpectedRateChanges
                              )

layout.mat <- matrix(1:6,nrow=3,ncol=2,byrow=TRUE)
layout(layout.mat)
tess.plot.output(output,
                 fig.types = c("speciation rates",
                               "speciation shift times",
                               "extinction rates",
                               "extinction shift times",
                               "net-diversification rates",
                               "relative-extinction rates"),
                las=2)
