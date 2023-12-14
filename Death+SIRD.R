# SIRD Analysis
library(covidStateSird)
library(foreach)
library(doParallel)


# set final day of training data -----------------------------------------------
endDate <- "2021-03-09"
minCase <- 100

# toggles "demo" mode, when TRUE exectues an abbreviated run of the model ------
DEMO <- TRUE

if(DEMO) {
  n.chains = 1
  n.iter = 200
  n.burnin = 10
  n.thin = 20
  ntree = 5
} else{
  n.chains = 3
  n.iter = 200000
  n.burnin = 10000
  n.thin = 500
  ntree = 500
}

# set RNG seed -----------------------------------------------------------------
set.seed(525600)

doParallel::registerDoParallel(cores=5)

covidDir <- "/Users/zavier/Desktop/Master/time\ series/Time_Series_Project"
outputPath <- file.path(covidDir, "Output")
dir.create(outputPath)
outputPath <- file.path(covidDir, "Output", Sys.Date())
dir.create(outputPath)
dir.create(file.path(outputPath, "Tables"))
dir.create(file.path(outputPath, "Plots"))
dir.create(file.path(outputPath, "Data"))

stateInterventions <- read.csv(paste0(covidDir,
                                      "/Data/StateInterventionDates.csv"),
                               stringsAsFactors = FALSE,
                               header = TRUE)

stateDataCSV <- "https://raw.githubusercontent.com/COVID19Tracking/covid-public-api/master/v1/states/daily.csv"

stateCovidData <- read.csv(stateDataCSV)

# covariates for death model ---------------------------------------------------
covariates <- read.csv(paste0(covidDir, 
                              "/Data/COVID-19_Location_Covariates_AnalyticStates.csv"))
covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] <-
  covariates[, c("heartDisease", "lungCancer", "diabetes", "copd")] / 100000

# Exclude stroke, because it is missing for 5 states 
covariates <- covariates[, -which(names(covariates) == "stroke")]

states <- c("AK", "AL", "AR", "AZ", "CA", "CO", "CT", "DE", "FL", "GA", 
            "HI", "IA", "ID", "IL", "IN", "KS", "KY", "LA", "MA", "MD", 
            "ME", "MI", "MN", "MO", "MS", "MT", "NC", "ND", "NE", "NH", 
            "NJ", "NM", "NV", "NY", "OH", "OK", "OR", "PA", "RI", "SC", 
            "SD", "TN", "TX", "UT", "VA", "VT", "WA", "WI", "WV", "WY")

velocLogCases <- velocLogDeaths <- data.frame()
loc <- 0
for(i in 1:length(states)) {
  loc <- loc + 1
  
  velocLoc <- velocitiesState(stateCovidData, states[i], stateInterventions, 
                              minCases = minCase, endDate = endDate)
  population <- stateInterventions$statePopulation[stateInterventions$stateAbbreviation == states[i]]
  
  velocLogCases <- rbind(velocLogCases, cbind(velocLoc$cases, loc, 
                                              row.names = NULL))
}

velocLogCases$y[velocLogCases$y <= 0] <- NA
velocLogCases <- velocLogCases[complete.cases(velocLogCases),]


locs <- unique(velocLogCases $loc)
firstObs <- rep(NA, length(locs))
for(i in 1:length(locs)) {
  firstObs [i] <- min(which(velocLogCases$loc == locs[i]))
}

velocLogCasesList = list(y = velocLogCases$y, 
                         N = table(velocLogCases$loc), 
                         nLoc = length(locs), 
                         firstObs = firstObs)

params = c("tau", "u", "mu", "phi")

start <- Sys.time()

# define ar1jagsModel
#' @export
ar1JagsModel <- function() {
  mu_mu ~ dnorm(0, .1)
  mu_phi ~ dunif(0, 1)
  mu_tau ~ dgamma(0.001, 0.001)
  sig2_phi <- .05
  sig2_tau <- 1
  for(j in 1:nLoc) {
    logmu[j] ~ dnorm(mu_mu, .01)
    mu[j] <- -exp(logmu[j])
    # truncate the beta to avoid numerical instability in the slice sampler near 0 and 1
    phi[j] ~ dbeta(((1 - mu_phi)/(sig2_phi) - (1/mu_phi)) * mu_phi^2,
                   (((1 - mu_phi)/(sig2_phi) - (1/mu_phi)) * mu_phi^2) * ((1/mu_phi) - 1));T(.001,.999)
    tau[j] ~ dgamma(mu_tau^2 / sig2_tau, mu_tau / sig2_tau)
    sd[j] <- 1 / sqrt(tau[j])
  }
  
  for(j in 1:nLoc) {
    
    u[firstObs[j]] <- log(y[firstObs[j]])
    
    for(i in (firstObs[j]+1):(firstObs[j] + N[j] - 1)) {
      u[i] <- mu[j] + phi[j] * log(y[i - 1])
      y[i] ~ dlnorm(u[i], tau[j])
    }
  }
}


velocModel <- R2jags::jags(data = velocLogCasesList, 
                           inits = NULL, 
                           parameters.to.save = params,
                           model.file = ar1JagsModel, 
                           n.chains = n.chains, 
                           n.iter = n.iter, 
                           n.burnin = n.burnin, 
                           n.thin = n.thin, 
                           DIC = F)

timeElapsed <- (Sys.time() - start)

# plotVelocityFit(velocModel$BUGSoutput$mean,
#                stateCovidData, states, stateInterventions,
#                fileName = paste0(outputPath, "/Plots/velocityModelFit.pdf"))


posteriorSamples <- velocModel$BUGSoutput$sims.list

save(posteriorSamples,
     file = paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))


# train death model ------------------------------------------------------------
randomForestDeathModel <- deathForest(stateCovidData, 
                                      states, 
                                      covariates, 
                                      21, 
                                      fileOut = paste0(outputPath, 
                                                       "/randomForestDeathModel.Rdata"),
                                      ntree = ntree)

load(paste0(outputPath, "/CasePosteriorSamples", endDate, ".Rdata"))

# test on CA
# stateSird("CA", covariates, stateInterventions, stateCovidData, 
# randomForestDeathModel, posteriorSamples, rfError = T, plots = T, 
# endPlotDay = "2021-03-01")

# run state SIRD models --------------------------------------------------------
foreach(i = 1:length(states)) %dopar% {
  stateSird(states[i],
            states,
            covariates,
            stateInterventions,
            velocLogCases,
            stateCovidData,
            randomForestDeathModel,
            posteriorSamples,
            rfError = TRUE,
            plots = TRUE,
            n.t = 120,
            lagDays = 21,
            minCases = 100,
            endDay = endDate,
            endPlotDay = "2021-05-01")
}