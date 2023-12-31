## ----load_libraries-----------------------------------------------------------
library(covidStateSird)
library(foreach)


# the last day of data to use
endDate <- "2020-11-23"
minCase <- 100

set.seed(525600)

covidDir <- "../"

doParallel::registerDoParallel(cores=5)

outputPath <- file.path(covidDir, "/Output", Sys.Date())
dir.create(outputPath)
dir.create(file.path(outputPath, "/Tables/"))
dir.create(file.path(outputPath, "/Plots/"))
dir.create(file.path(outputPath, "/Data/"))

stateInterventions <- read.csv(paste0(covidDir, "Data/StateInterventionDates.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)

stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"

stateCovidData <- read.csv(stateDataCSV)

covariates <- read.csv(paste0(covidDir, "Data/COVID-19 Location Covariates - analyticStates.csv"))
covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] <-
  covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] / 100000

states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", 
            "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", 
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", 
            "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", 
            "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")
 
velocLogCases <- velocLogDeaths <- data.frame()
loc <- 0
for(i in 1:length(states)) {
  loc <- loc + 1

  velocLoc <- velocitiesState(stateCovidData, states[i], stateInterventions, minCases = minCase, endDate = endDate)
  population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == states[i]]
  
  velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases,  loc, row.names = NULL))
}

velocLogCasesList <- as.list(velocLogCases)
velocLogCasesList$N <- nrow(velocLogCases)
velocLogCasesList$nLoc <- length(unique(velocLogCasesList$loc))

velocLogCasesList$y[velocLogCasesList$y <= 0] <- NA

params <- c("mu_a", "mu_b", "tau", "a", "b", "g", "d", "mu_g", "mu_d", "mu", "alpha", "beta", "mu_alpha", "mu_beta")
 
start <- Sys.time()

velocModel <- R2jags::jags(data = velocLogCasesList, inits = NULL, parameters.to.save = params,
  model.file = riasJagsModel, n.chains = 3, n.iter = 100, n.burnin = 10, n.thin = 10, DIC = F)

timeElapsed <- (Sys.time() - start)

plotVelocityFit(velocModel$BUGSoutput$mean,
                stateCovidData, states, stateInterventions,
                fileName = paste0(outputPath, "/Plots/velocityModelFit.pdf"))
