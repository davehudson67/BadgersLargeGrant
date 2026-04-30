library(sf)
library(nimble)
library(lubridate)
library(tidyverse)

bc <- readRDS("../BadgersSG2024.rds")
load("spatialInfo.RData")
SGs <- settGrid$Sett
levels(as.factor(bc$socg))
bc$SG <- iconv(bc$socg, from = "latin1", to = "UTF-8", sub = "")
bc$SG <- gsub(" ", "", bc$SG)

settGrid <- settGrid %>%
  st_drop_geometry() %>%
  rename(SG = Sett)

studyArea <- studyArea %>%
  st_drop_geometry()

bc <- bc %>%
  filter(SG %in% SGs) %>%
  left_join(settGrid, by = "SG") %>%
  mutate(date = ymd(date)) %>%
  filter(year(date) > 1981) %>%
  arrange(date) %>%
  mutate(primary = factor(year(date), levels = c("1982", "1983", "1984", "1985", "1986", "1987", 
                                                 "1988", "1989", "1990", "1991", "1992", "1993", 
                                                 "1994", "1995", "1996", "1997", "1998", "1999", 
                                                 "2000", "2001", "2002", "2003", "2004", "2005", 
                                                 "2006", "2007", "2008", "2009", "2010", "2011", 
                                                 "2012", "2013", "2014", "2015", "2016", "2017", 
                                                 "2018", "2019", "2020"))) %>%
  dplyr::select(tattoo, date, sett, pm, socg, SG, sex, age, primary, trap_season, row_index, col_index) %>%
  ungroup()

bc$primary <- as.numeric(bc$primary)

bc <- bc %>%
  group_by(tattoo) %>%
  mutate(minPrimary = min(primary)) %>%
  mutate(maxPrimary = max(primary)) %>%
  group_by(tattoo, primary) %>%
  mutate(lastSecondary = max(trap_season)) %>%
  ungroup() %>%
  arrange(tattoo)

n.prim <- max(bc$primary)
n.sec <- max(bc$trap_season)
dt <- rep(1, times = n.prim - 1) #yearly primary sessions
nind <- length(unique(bc$tattoo))

# Get unique individual IDs and time occasions
ids <- unique(bc$tattoo)

# K[i]: Last primary session for individual i (allows censoring)
K <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(maxPrimary)

first <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(minPrimary)

# J[i, k]: Last secondary session for individual i in primary session k 
J <- matrix(NA, nind, n.prim)
rownames(J)<- ids

for(i in 1:nrow(bc)){
  ind <- bc$tattoo[i]
  J[ind, bc$primary[i]] <- bc$lastSecondary[i]
}

# H[i,j,k]: Trap index for capture of individual i in secondary session j
#             within primary session k. 0 = not captured, other values are r+1
# Create a named vector for mapping sett names to numbers
sett_map <- setNames(seq_along(SGs), SGs)
# Replace sett names with corresponding numbers in the dataframe
bc$SGnum <- sett_map[bc$SG]

H <- array(NA, dim = c(nind, n.sec, n.prim))#, dimnames = list(ids, NULL, NULL))




bc$tat <- as.numeric(as.factor(bc$tattoo))





# Populate H capture history array
for (i in 1:nrow(bc)) {
  p <- bc$primary[i]
  s <- bc$trap_season[i]
  ind <- bc$tat[i]
  H[ind, s, p] <- bc$SGnum[i]
}

# R: Number of sett locations
R <- length(unique(bc$SG))

# X[r,]:    Location of trap r = 1..R (coordinate)
setts <- settGrid$SG
X <- settGrid %>%
  dplyr::select(row_index, col_index)

# Ones[i,j,k]: An array of all ones used in the "Bernoulli-trick"
Ones <- array(1, dim = c(nind, n.sec, n.prim))

# Priors for first centre of activity:
# Assuming that first centre of activity is always within +/- maxd
# from the mean capture location in both x and y direction.
maxd = 2*4
mean.x = apply(H, c(1, 3), function(i) mean(X[i - 1, 1], na.rm=T))
mean.y = apply(H, c(1, 3), function(i) mean(X[i - 1, 2], na.rm=T))
first.mean.x = apply(mean.x, 1, function(i) i[min(which(!is.na(i)))])
first.mean.y = apply(mean.y, 1, function(i) i[min(which(!is.na(i)))])
xlow = first.mean.x - maxd
xupp = first.mean.x + maxd
ylow = first.mean.y - maxd
yupp = first.mean.y + maxd
# xlow[i]: Lower bound in uniform prior for first x-location of individual i
# xupp[i]: Upper bound in uniform prior for first x-location of individual i
# ylow[i]: Lower bound in uniform prior for first y-location of individual i
# yupp[i]: Upper bound in uniform prior for first y-location of individual i

# N[1]: Number of individuals that are only seen in the last primary session
#         or the session they are censored (these individuals must be sorted
#         first in the data)
# N[2]: Total number of individuals
N = c(sum((K-first)==0), length(first))

# For the model, individuals must be sorted such that the individuals
# that are only seen in the last primary session or the session they are
# censored comes first:
ord = order(K-first)
K = K[ord]
J = J[ord,]
first = first[ord]
xlow = xlow[ord]


###################################################################
### JAGS MODEL SPECIFICATION - Spatial robust design model with ###
### (zero-inflated) exponential dispersal distribution          ###
###################################################################
#
# DATA:

# - n.prim:   Number of primary sessions
# - dt[k]:    Length of interval k
# - first[i]: First primary session of individual i
# - K[i]:     Last primary session for individual i (allows censoring)
# - J[i,k]:   Last secondary session for individual i in primary session k 
#             (allows censoring)
# - R:        Number of traps
# - X[r,]:    Location of trap r = 1..R (coordinate)
# - H[i,j,k]: Trap index for capture of individual i in secondary session j
#             within primary session k. 1 = not captured, other values are r+1
# - Ones[i,j,k]: An array of all ones used in the "Bernoulli-trick" (see 
#                OpenBUGS user manual) in the observation likelihood. This saves
#                computation time as it is not necessary to compute the complete 
#                capture probabilities for every individual*trap*session (it is
#                sufficient to compute the g-variable for every
#                ind*trap*primary)
# - xlow[i]: Lower bound in uniform prior for first x-location of individual i
# - xupp[i]: Upper bound in uniform prior for first x-location of individual i
# - ylow[i]: Lower bound in uniform prior for first y-location of individual i
# - yupp[i]: Upper bound in uniform prior for first y-location of individual i
#
# PARAMETERS:
# - kappa, sigma and lambda: Space use and recapture probability parameters
#                            (eq. 5 in paper)
# - beta: additive effects on log(lambda):
#           - beta[1]: effect of tod==2
#           - beta[2]: effect of gr==2
# - Phi[gr[i]]:   Survival for group gr[i] per time-unit 
# - phi[gr[i],k]: Survival probability for individual i belonging to group (sex) 
#                 gr[i] over the k'th interval given that the individual is
#                 alive at the beginning of the interval.
# - dmean[gr[i]]: Mean dispersal distance (exponential distribution) given
#                 dispersal for group gr[i]
# - Psi[gr[i]]:   Probability of not emigrating for group gr[i] per time-unit
#                 (only for zero-inflated model, commented out)
# - psi[gr[i],k]: Probability of dispersal (1-zero-inflation probability) for
#                 group gr[i] in interval k (only for zero-inflated model, 
#                 commented out)
# - d.offset[gr[i]]: Minimum dispersal distance given dispersal for group gr[i]
#                    (only for zero-inflated model, commented out)                 
# STATE VARIABLES:
# - z[i,k]: Alive state of individual i in primary session k (1=alive, 0=dead)
# - S[i,,k]: Centre of activity of individual i in primary session k (x- and
#            y-coordinate)
# - u[i,k]: Dispersal state for individual i in interval k (1=dispersal,
#           0=no dispersal)(only for zero-inflated model, commented out)

# Writing the model specification to a text file:
code <- nimbleCode({
  ## PRIORS AND CONSTRAINTS: ##
  #  Space use and recapture probability parameters:
  kappa ~ dunif(0,50)
  sigma ~ dunif(0.1,20)
  lambda <- exp(log.lambda0)
  PL ~ dunif(0.01, 0.99)
  log.lambda0 <- log(-log(1-PL)) 
  
  # Survival parameters:
  Phi ~ dunif(0,1)
    for(k in 1:(n.prim - 1)){
      phi[k] <- pow(Phi, dt[k])
    }

  # Dispersal parameters:
  dmean ~ dunif(0,100)
  dlambda <- 1/dmean
  # For zero-inflated model:  
  #  d.offset[sex] ~ dunif(5,100)
  #  Psi0[sex] ~ dunif(0,1)
  #  for(k in 1:(n.prim-1)){
  #    psi[sex,k] <- 1 - pow(Psi0[sex], dt[k])
  #  }
  
  ## MODEL: ##
  
  # Loop over individuals that are only seen in the last primary session or the
  # session they are censored
  for(i in 1:N[1]){
    z[i, first[i]] ~ dbern(1)
    S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
    S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
    g[i, first[i], 1] <- 0
    for(r in 1:R){ # trap
        D[i, r, first[i]] <- sqrt(pow(S[i, 1, first[i]] - X[r, 1], 2) + pow(S[i, 2, first[i]] - X[r, 2], 2))
        g[i, first[i], r + 1] <- exp(-pow(D[i, r, first[i]] / sigma)) # Trap exposure
    }
    G[i, first[i]] <- sum(g[i, first[i], 1:R]) # Total trap exposure
    for(j in 1:J[i, first[i]]){
        P[i, j, first[i]] <- 1 - exp(-lambda * G[i, first[i]]) # Probability of being captured
      	PI[i, first[i], j] <- step(H[i, j, first[i]] - 2) * (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) * P[i, j, first[i]] +
      	  (1 - step(H[i, j, first[i]] - 2)) * (1 - P[i, j, first[i]])
      	Ones[i, j, first[i]] ~ dbern(PI[i, first[i], j])
    }
  }
  
  # Loop over all other individuals
  for(i in (N[1] + 1):N[2]){
    z[i, first[i]] ~ dbern(1)
    S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
    S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
    # First primary session:
    g[i, first[i], 1] <- 0
    for(r in 1:R){ # trap
        D[i, r, first[i]] <- sqrt(pow(S[i, 1, first[i]] - X[r, 1], 2) + pow(S[i, 2, first[i]] - X[r, 2], 2))
        g[i, first[i], r + 1] <- exp(-pow(D[i, r, first[i]] / sigma)) # Trap exposure
    }
    G[i, first[i]] <- sum(g[i, first[i], 1:R]) # Total trap exposure
    for(j in 1:J[i, first[i]]){
        P[i, j, first[i]] <- 1 - exp(-lambda * G[i, first[i]]) # Probability of being captured
      	PI[i, first[i], j] <- step(H[i, j, first[i]] - 2) * (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) * P[i, j, first[i]] + 
      	  (1 - step(H[i, j, first[i]] - 2)) * (1 - P[i, j, first[i]])
      	Ones[i, j, first[i]] ~ dbern(PI[i, first[i], j])
    }
    
    ## Later primary sessions
    for(k in (first[i] + 1):K[i]){ # primary session
      theta[i, k - 1] ~ dunif(-3.141593,3.141593) # Prior for dispersal direction 
      z[i, k] ~ dbern(Palive[i, k - 1])
      Palive[i, k - 1] <- z[i, k - 1] * phi[k - 1] # Pr(alive in primary session k) gr[i] = sex
      d[i, k - 1] ~ dexp(dlambda)
      # For zero-inflated model, replace line above with:
      #u[i,k-1] ~ dbern(psi[gr[i],k-1])
      #dd[i,k-1] ~ dexp(dlambda[gr[i]])
      #d[i,k-1] <- u[i,k-1]*(d.offset[gr[i]] + dd[i,k-1])
      S[i, 1, k] <- S[i, 1, k - 1] + d[i, k - 1] * cos(theta[i, k - 1])
      S[i, 2, k] <- S[i, 2, k - 1] + d[i, k - 1] * sin(theta[i, k - 1])
      g[i, k, 1] <- 0
      for(r in 1:R){ # trap
        D[i, r, k] <- sqrt(pow(S[i, 1, k] - X[r, 1], 2) + pow(S[i, 2, k] - X[r, 2], 2))  # Squared distance to trap
        g[i, k, r + 1] <- exp(-pow(D[i, r, k] / sigma)) # Trap exposure
      }
      G[i, k] <- sum(g[i, k, 1:R]) # Total trap exposure
      for(j in 1:J[i, k]){
        P[i, j, k] <- (1 - exp(-lambda * G[i, k])) * z[i, k] # Probability of being captured
      	PI[i, k, j] <- step(H[i, j, k] - 2) * (g[i, k, H[i, j, k]] / (G[i, k] + 0.000000001)) * P[i, j, k] + (1 - step(H[i, j, k] - 2)) * (1 - P[i, j, k])
      	Ones[i, j, k] ~ dbern(PI[i, k, j])
      }
    }
  }
  })


consts <- list(R = nrow(X),
               # Prior parameters for initial centre of activity:
               xlow = xlow,
               xupp = xupp,
               ylow = ylow,
               yupp = yupp,
               N = N,
               K = K,
               J = J,
               first = first,
               X = as.matrix(X),
               n.prim = n.prim,
               dt = dt)
               
H[is.na(H)] <- 1

data = list(H = H,
            Ones = Ones)
  
model <- nimbleModel(code, constants = consts, data = data)
model$initializeInfo()

cModel <- compileNimble(model)
config <- configureMCMC(model)
rMCMC <- buildMCMC(config)
cMCMC <- compileNimble(rMCMC)

# Run the MCMC
system.time(run <- runMCMC(cMCMC,
                           niter = 10000, 
                           nburnin = 2400, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE))
