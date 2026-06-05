library(tidyverse)
library(nimble)

rm(list=ls())

bc <- readRDS("BadgersSG2024.rds")
load("SpatialCMR/spatialInfo.RData")
setts <- settGrid$Sett
levels(as.factor(bc$sett))
bc$Sett <- iconv(bc$sett, from = "latin1", to = "UTF-8", sub = "")
bc$Sett <- gsub(" ", "", bc$Sett)
settGrid <- select(settGrid, -geometry)
studyArea <- select(studyArea, - geometry)
bc <- bc %>%
  filter(Sett %in% setts)

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
  arrange(tattoo, occ) %>% # Arrange by tattoo and occ
  mutate(Sett = as.numeric(as.factor(sett))) %>%
  filter(Sett <= 17, .preserve = TRUE)

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
  state <- badger_data$Sett[i]
  capture_history[ind, occ] <- state
}

## Find first capture event
get.first <- function(x) min(which(x != 0))
f <- apply(capture_history, 1, get.first)

code  <- nimbleCode({
  #----------------------------------------------
  # Parameters:
  # phi[k]: survival probability at site k (k from 1 to 17)
  # psi[k, l]: movement probability from site k to site l (k, l from 1 to 17)
  # p[k]: recapture probability at site k (k from 1 to 17)
  #----------------------------------------------
  # States (S):
  # 1 alive at site 1
  # 2 alive at site 2
  # ...
  # 17 alive at site 17
  # 18 dead
  #
  # Observations (O):
  # 1 seen at site 1
  # 2 seen at site 2
  # ...
  # 17 seen at site 17
  # 18 not seen
  #----------------------------------------------
  # Priors and constraints
  for (k in 1:17) {
    phi[k] ~ dunif(0, 1) #survival
    p[k] ~ dunif(0, 1) #recapture
  }
  
  # Transition probabilities
  for (k in 1:17) { # for each site
    for (l in 1:16) {
        lpsi[k, l] ~ dnorm(0, 0.001)
      }
    }
  
  # Define the transition probabilities
  for (k in 1:17) {
    for (l in 1:16) {
        psi[k, l] <- exp(lpsi[k, l]) / (1 + sum(exp(lpsi[k, 1:16])))
    }
    # last transition probability
    psi[k, 17] <- 1 - sum(psi[k, 1:16])
  }
  
  # Define state-transition and observation matrices
  for (i in 1:nind) { 
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions - 1)) {
      
      for (k in 1:17) {
        for (l in 1:17) {
          ps[k, i, t, l] <- phi[k] * psi[k, l]
        }
        ps[k, i, t, 18] <- 1 - phi[k]  # transition to dead state
      }
      
      ps[18, i, t, 18] <- 1  # once dead, stay dead
      
      for (k in 1:17) {
        ps[18, i, t, k] <- 0  # dead cannot transition to alive
      }
      
      # Define probabilities of O(t) given S(t)
      for (k in 1:17) {
        for (l in 1:17) {
          po[k, i, t, l] <- 0  # initialize all to 0
        }
        po[k, i, t, k] <- p[k]  # set the diagonal to p[k]
        po[k, i, t, 18] <- 1 - p[k]  # not seen probability
      }
      po[18, i, t, 18] <- 1  # if dead, not seen
      
      for (k in 1:17) {
        po[18, i, t, k] <- 0  # if dead, not seen at any site
      }
    } # t
  } # i
  
  # Likelihood
  for (i in 1:nind) {
    # Define latent state at first capture
    z[i, f[i]] <- Y[i, f[i]]
    for (t in (f[i] + 1):n.occasions) {
      # State process: draw S(t) given S(t-1)
      z[i, t] ~ dcat(ps[z[i, t-1], i, t-1, 1:18])
      # Observation process: draw O(t) given S(t)
      y[i, t] ~ dcat(po[z[i, t], i, t-1, 1:18])
    } # t
  } # i
})

consts <- list(f = f, n.occasions = ncol(capture_history), nind = nrow(capture_history))

data <- list(y = capture_history)

zinits <- capture_history
zinits[zinits == 0] <- sample(c(1:17), sum(zinits==0), replace = TRUE)

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

      