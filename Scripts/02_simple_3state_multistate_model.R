library(tidyverse)
library(lubridate)

rm(list=ls())

bc <- readRDS("BadgersSG2024.rds")
#capture_history <- read_csv("CHexample.csv")
#capture_history <- capture_history[-1]

## Arrange the data by ID and date
badger_data <- bc %>%
  arrange(tattoo, date)

badger_data <- badger_data %>%
  group_by(tattoo, date) %>%
  distinct(tattoo, date, .keep_all = TRUE) %>% # Keep only the first entry for each individual on any given date
  ungroup() %>%
  filter(pm != "Yes") %>% # Filter out rows where pm is "Yes"
  filter(!is.na(sett)) %>%
  group_by(tattoo) %>%
  filter(n() > 1) %>% # Keep individuals that were caught more than once
  ungroup() %>%
  select(-c(9:19)) %>% # Select columns to keep (removing columns 9 to 19)
  arrange(tattoo, occ) # Arrange by tattoo and occ

## Create a new column to indicate if the badger was caught at the same sett as the previous occasion
badger_data <- badger_data %>%
  group_by(tattoo) %>%
  mutate(prev_sett = lag(sett),
         state = if_else(sett == prev_sett, 1, 2, missing = 1)) %>%
  ungroup() %>%
  select(-prev_sett)

## Have a look
summary(as.factor(badger_data$state))
length(unique(as.factor(badger_data$tattoo)))

## Construct CH
## Get the unique list of individuals and occasions
badger_data <- arrange(badger_data, occ)
individuals <- unique(badger_data$tattoo)
occasions <- as.numeric(1:177)

## Initialize the capture history matrix with zeros (not captured)
capture_history <- matrix(0, nrow = length(individuals), ncol = length(occasions))
rownames(capture_history) <- individuals
colnames(capture_history) <- occasions

## Fill the capture history matrix with the state values
for (i in 1:nrow(badger_data)) {
  ind <- badger_data$tattoo[i]
  occ <- badger_data$occ[i]
  state <- badger_data$state[i]
  capture_history[ind, occ] <- state
}

## Find first capture event
get.first <- function(x) min(which(x != 0))
f <- apply(capture_history, 1, get.first)

#capture_history <- capture_history[1:10,]
#f <- f[1:10]
capture_history <- as.matrix(capture_history)

## Build the model
code  <- nimbleCode({
  #phiA : survival at same site
  #phiB : survival at new site
  #psiAB : movememnt probability to a new site
  #pA : recapture probability
  
  # States (z)
  # 1 = Alive at same site
  # 2 = Alive at new site
  # 3 = Dead
  # Observations (y)
  # 1 = Not seen
  # 2 = Seen at same site
  # 3 = Seen at new site
  
  # Priors and constraints
    phiA ~ dunif(0, 1)
    phiB ~ dunif(0, 1)
    psiAB ~ dunif(0, 1)
    psiBA ~ dunif(0, 1)
    pA ~ dunif(0, 1)
    pB ~ dunif(0, 1)
    
    # Probabilities of state z(t + 1) given z(t)
    # gamma[from state, to state]
      gamma[1, 1] <- phiA * (1 - psiAB)
      gamma[1, 2] <- phiA * psiAB
      gamma[1, 3] <- 1 - phiA
      gamma[2, 1] <- phiB * psiBA
      gamma[2, 2] <- phiB * (1 - psiBA)
      gamma[2, 3] <- 1 - phiB
      gamma[3, 1] <- 0
      gamma[3, 2] <- 0
      gamma[3, 3] <- 1
      
    # omega[z(t), y(t)] = omega[State, Observation]
      omega[1, 1] <- 1 - pA
      omega[1, 2] <- pA
      omega[1, 3] <- 0
      omega[2, 1] <- 1 - pB
      omega[2, 2] <- 0
      omega[2, 3] <- pB
      omega[3, 1] <- 1
      omega[3, 2] <- 0
      omega[3, 3] <- 0
  
  # Likelihood
  for (i in 1:nind){
    z[i, f[i]] <- y[i, f[i]] - 1
    for(t in (f[i] + 1):n.occasions){
      # State process: z(t) given z(t - 1)
      z[i, t] ~ dcat(gamma[z[i, t-1], 1:3])
      # Observation process: y(t) given z(t)
      y[i, t] ~ dcat(omega[z[i, t], 1:3]) 
    }
  }
})

consts <- list(f = f, n.occasions = ncol(capture_history), nind = nrow(capture_history))

data <- list(y = capture_history + 1)

zinits <- capture_history
zinits[zinits == 0] <- sample(c(1, 2), sum(zinits==0), replace = TRUE)

inits <- list(phiA = runif(1, 0, 1),
              phiB = runif(1, 0, 1),
              pA = runif(1, 0, 1),
              pB = runif(1, 0, 1),
              psiAB = runif(1, 0, 1),
              psiBA = runif(1, 0, 1),
              z = zinits)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)
config <- configureMCMC(model)
rMCMC <- buildMCMC(config)
cMCMC <- compileNimble(rMCMC, project = model)

# Run the MCMC
system.time(run <- runMCMC(cMCMC, 
                           niter = 100000, 
                           nburnin = 24000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

run$summary
plot(run$samples)
