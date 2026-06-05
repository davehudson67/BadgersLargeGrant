library(sf)
library(nimble)
library(lubridate)
library(tidyverse)
library(mcmcplots)

# Load and preprocess data
bc <- readRDS("../BadgersSG2024.rds") %>%
  droplevels()
load("spatialInfo.RData")

bc$brock[is.na(bc$brock)] <- 0
bc$GAMMA[is.na(bc$GAMMA)] <- 0
bc$statpak[is.na(bc$statpak)] <- 0

# Create a binary "test_positive" column
bc <- bc %>%
  rowwise() %>%
  mutate(test_positive = ifelse(any(c_across(starts_with("Cult_")) > 0, na.rm = TRUE) | 
                                  brock > 0 | statpak > 0 | GAMMA  > 0, 1, 0)) %>%
  ungroup()

# Categorize individuals as "never_positive" or "cub_positive"
bc <- bc %>%
  group_by(tattoo) %>%
  mutate(
    first_year = min(year_fc, na.rm = TRUE),  # Identify the first year of life for each individual
    cub_positive = ifelse(any(test_positive == 1 & year_fc == first_year, na.rm = TRUE), 1, 0),
    never_positive = ifelse(all(test_positive == 0, na.rm = TRUE), 1, 0)
  ) %>%
  ungroup()

# Step 3: Summarize the results for each individual
individual_summary <- bc %>%
  group_by(tattoo) %>%
  summarize(
    cub_positive = max(cub_positive),  # If cub_positive is 1 in any row, it will be 1 for the individual
    never_positive = max(never_positive)  # If never_positive is 1 in all rows, it will be 1 for the individual
  )

# Convert and clean data
bc$SG <- iconv(bc$socg, from = "latin1", to = "UTF-8", sub = "")
bc$SG <- gsub(" ", "", bc$SG)

settGrid <- settGrid %>%
  st_drop_geometry() %>%
  rename(SG = Sett)

# Remove unnecessary data
rm(studyArea)

SGs <- settGrid$SG

# Filter and join datasets, mutate necessary columns
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

# Remove pm only individuals
bc <- bc %>%
  group_by(tattoo) %>%
  filter(!all(pm == "Yes")) %>%
  ungroup() %>%
  arrange(tattoo)

# Further filter to keep only the "pm == Yes" entry if it occurs in the same year/primary session as another observation
#bc <- bc %>%
#  group_by(tattoo, primary) %>%
#  filter(!(pm == "Yes" & n() > 1 & !duplicated(pm == "Yes"))) %>%
#  ungroup()

# repeat this step to remove individuals who now only have one entry and its a pm
#bc <- bc %>%
#  group_by(tattoo) %>%
#  filter(!all(pm == "Yes")) %>%
#  ungroup() %>%
#  arrange(tattoo)

# Randomly select 50 individuals
#set.seed(42)  # for reproducibility
#sampled_ids <- sample(unique(bc$tattoo), 100)
#bc <- bc %>%
#  filter(tattoo %in% sampled_ids)

# Convert primary to numeric
bc$primary <- as.numeric(bc$primary) + 1

# Set parameters
n.prim <- max(bc$primary)
n.sec <- max(bc$trap_season)
dt <- rep(1, times = n.prim - 1) # yearly primary sessions
nind <- length(unique(bc$tattoo))

# Calculate individual-specific parameters
bc <- bc %>%
  group_by(tattoo) %>%
  mutate(minPrimary = min(primary),
         maxPrimary = max(primary)) %>%
  group_by(tattoo, primary) %>%
  mutate(lastSecondary = max(trap_season)) %>%
  ungroup() %>%
  group_by(tattoo) %>%
  mutate(firstSecondary = ifelse(primary == minPrimary, min(trap_season), NA)) %>%
  fill(firstSecondary, .direction = "downup") %>%  # Propagate first Secondary value to all primary occasions
  mutate(lastSecondary = ifelse(primary == maxPrimary, lastSecondary, NA)) %>%
  fill(lastSecondary, .direction = "downup") %>%  # Propagate lastSecondary to all previous primary sessions
  mutate(# Adjust death.occasion if pm == "Yes"
    death.occasion = ifelse(pm == "Yes", primary, n.prim + 1),
    death.occasion = min(death.occasion, na.rm = TRUE),
    death.secondary = ifelse(pm == "Yes", trap_season, 0)) %>%
  mutate(lastSecondary = ifelse(primary == maxPrimary & pm == "Yes", lastSecondary - 1, lastSecondary)) %>%
  mutate(maxPrimary = ifelse(pm == "Yes" & lastSecondary == 0, maxPrimary - 1, maxPrimary)) %>%
  mutate(lastSecondary = ifelse(lastSecondary == 0 & pm == "Yes", 4, lastSecondary)) %>%
  arrange(tattoo) %>%
  mutate(maxPrimary = min(maxPrimary)) %>%
  ungroup()

# Make sure individuals are alive and dead in same primary session
bc <- bc %>%
  group_by(tattoo) %>%
  mutate(maxPrimary = ifelse(maxPrimary == death.occasion, maxPrimary - 1, maxPrimary))

death.primary <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.occasion)
death.primary[death.primary > n.prim] <- NA

death.secondary <- bc %>%
  group_by(tattoo) %>%
  mutate(death.secondary = max(death.secondary)) %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.secondary)
death.secondary[death.secondary == 0] <- NA

# Get unique SGs and filter settGrid
SGs <- levels(as.factor(bc$SG))

settGrid <- settGrid %>%
  filter(SG %in% SGs)

# Get unique individual IDs
ids <- unique(bc$tattoo)

# K: Last primary session for individual
#K <- bc %>%
#  distinct(tattoo, .keep_all = TRUE) %>%
#  pull(maxPrimary)

# First primary session for each individual
first <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(minPrimary)

# J[i, k]: possible sampling occasions from first to last capture 
J <- matrix(n.sec, nind, n.prim)

# Death occasion
death.occasion <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(death.occasion)

# K matrix
K <- rep(n.prim, nind)
K[death.occasion < n.prim] <- death.occasion[death.occasion < n.prim]

# Adjust K and J based on death information
for (i in 1:length(death.primary)) {
  if (!is.na(death.primary[i])) {
    K[i] <- death.primary[i] - 1
    if (!is.na(death.secondary[i])) {
      J[i, death.primary[i]] <- death.secondary[i]
    }
  }
}

# Fill in J matrix
#for (i in 1:nind) {
#  individual_data <- bc %>% 
#    arrange(desc(primary)) %>%
#    filter(tattoo == unique(bc$tattoo)[i]) %>%
#    slice(1)  # Ensures only one row per individual
#  
#  minPrimary <- individual_data$minPrimary
#  maxPrimary <- individual_data$maxPrimary
#  lastSecondary <- individual_data$lastSecondary
#  firstSecondary <- individual_data$firstSecondary
#  
#  # Set J for the first primary session
#  J[i, minPrimary] <- 4 - firstSecondary + 1
#  
#  # Set J to the maximum number of secondary sessions for all sessions between first and last
#  if (maxPrimary > minPrimary) {
#    J[i, (minPrimary + 1):(maxPrimary - 1)] <- n.sec
#  }
#  
#  # Set J for the last primary session to reflect all possible secondary sessions up to lastSecondary
#  J[i, maxPrimary] <- lastSecondary
#}

# H[i,j,k]: Trap index for capture of individual i in secondary session j within primary session k. 1 = not captured, other values are r+1
# Create a named vector for mapping sett names to numbers
sett_map <- setNames(seq_along(SGs), SGs)
# Replace sett names with corresponding numbers in the dataframe
bc$SGnum <- sett_map[bc$SG]

H <- array(1, dim = c(nind, n.sec, n.prim), dimnames = list(ids, NULL, NULL))

# Populate H capture history array
for (i in 1:nrow(bc)) {
  p <- bc$primary[i]
  s <- bc$trap_season[i]
  ind <- bc$tattoo[i]
  H[ind, s, p] <- bc$SGnum[i] + 1
}

# R: Number of sett locations
R <- length(unique(bc$SGnum))

# X[r,]:    Location of trap r = 1..R (coordinate)
setts <- settGrid$SG

settGrid$SG <- factor(settGrid$SG, levels = SGs)
# Reorder the dataframe based on this factor
settGrid <- settGrid[order(settGrid$SG), ]
rownames(settGrid) <- NULL
X <- settGrid %>%
  dplyr::select(col_index, row_index)

# Ones[i,j,k]: An array of all ones used in the "Bernoulli-trick"
Ones <- array(1, dim = c(nind, n.sec, n.prim))

# Last primary session for individual
last.prim <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  pull(maxPrimary)

# z_data: Create a matrix of NA values
z_data <- matrix(NA, nrow = nind, ncol = n.prim)

death.occNA <- death.occasion
death.occNA[death.occNA == 41] <- NA

# Fill the matrix with 1s based on first_primary and last_primary
for (i in 1:nind) {
  z_data[i, first[i]:last.prim[i]] <- 1
  
  # If a known death is recorded, set all entries after the death to 0
  if (!is.na(death.occNA[i])) {
    z_data[i, death.occNA[i]:n.prim] <- 0
  }
}

which(first == death.occasion)
first[first == death.occasion] <- (first - 1)[first == death.occasion]

# Sex 
sex <- bc %>%
  distinct(tattoo, .keep_all = TRUE) %>%
  mutate(sex = as.numeric(sex)) %>%
  pull(sex)

# For the model, individuals must be sorted such that the individuals
# that are only seen in the last primary session or the session they are
# censored comes first:
ord = order(K - first)
K = K[ord]
J = J[ord,]
first = first[ord]
H <- H[ord, , ]
death.occasion <- death.occasion[ord]
z_data <- z_data[ord,]
sex <- sex[ord]

# Priors for first centre of activity:
# Assuming that first centre of activity is always within +/- maxd
# from the mean capture location in both x and y direction.
maxd = 2*8
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
N = c(sum((K - first) == 0), length(first))

code <- nimbleCode({
  ## PRIORS AND CONSTRAINTS: ##
  
  # Space use and recapture probability parameters:
  for (sx in 1:2){
    kappa[sx] ~ dunif(0, 10)
    sigma[sx] ~ dunif(0, 30)
    lambda[sx] <- exp(log.lambda0[sx])
    log.lambda0[sx] <- log(-log(1 - PL[sx])) 
    PL[sx] ~ dunif(0.01, 0.99)
    
    # Survival parameters:
    for(k in 1:(n.prim - 1)){
      phi[sx, k] ~ dunif(0, 1)
    }
    
    # Dispersal parameters:
    dmean[sx] ~ dunif(0, 100)
    dlambda[sx] <- 1/dmean[sx]
  }
  
  ## MODEL: ##
  
  # Loop over individuals that are only seen in the last primary session or the
  # session they are censored
  for(i in 1:N[1]){
    z[i, first[i]] ~ dbern(1)
    S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
    S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
    g[i, first[i], 1] <- 0 # for non-capture
    for(r in 1:R){ # trap
      D[i, r, first[i]] <- sqrt(pow(S[i, 1, first[i]] - X[r, 1], 2) + pow(S[i, 2, first[i]] - X[r, 2], 2))
      g[i, first[i], r + 1] <- exp(-pow(D[i, r, first[i]] / sigma[sex[i]], kappa[sex[i]])) # Trap exposure
    }
    G[i, first[i]] <- sum(g[i, first[i], 1:(R + 1)]) # Total trap exposure
    for(j in 1:J[i, first[i]]){
      P[i, j, first[i]] <- 1 - exp(-lambda[sex[i]] * G[i, first[i]]) # Probability of being captured
      captureProb[i, first[i], j] <- step(H[i, j, first[i]] - 2) * (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) * P[i, j, first[i]] + (1 - step(H[i, j, first[i]] - 2)) * (1 - P[i, j, first[i]])
      Ones[i, j, first[i]] ~ dbern(captureProb[i, first[i], j])
    }
  }
  
  # Loop over all other individuals
  for(i in (N[1] + 1):N[2]){
    z[i, first[i]] ~ dbern(1)
    S[i, 1, first[i]] ~ dunif(xlow[i], xupp[i]) # Prior for the first x coordinate
    S[i, 2, first[i]] ~ dunif(ylow[i], yupp[i]) # Prior for the first y coordinate
    # First primary session:
    g[i, first[i], 1] <- 0 # for non-capture
    for(r in 1:R){ # trap
      D[i, r, first[i]] <- sqrt(pow(S[i, 1, first[i]] - X[r,1],2) + pow(S[i, 2, first[i]] - X[r, 2], 2))
      g[i, first[i], r + 1] <- exp(-pow(D[i, r, first[i]] / sigma[sex[i]], kappa[sex[i]])) # Trap exposure
    }
    G[i, first[i]] <- sum(g[i, first[i], 1:(R + 1)]) # Total trap exposure
    for(j in 1:J[i, first[i]]){
      P[i, j, first[i]] <- 1 - exp(-lambda[sex[i]] * G[i, first[i]]) # Probability of being captured
      captureProb[i, first[i], j] <- step(H[i, j, first[i]] - 2) * (g[i, first[i], H[i, j, first[i]]] / (G[i, first[i]] + 0.000000001)) * P[i, j, first[i]] + (1 - step(H[i, j, first[i]] - 2)) * (1 - P[i, j, first[i]])
      Ones[i, j, first[i]] ~ dbern(captureProb[i, first[i], j])
    }
    ## Later primary sessions
    for(k in (first[i] + 1):K[i]){ # primary session
      
      # Known death occasion constraint
      z[i, k] ~ dbern(Palive[i, k - 1] * step(death.occasion[i] - k)) 
      #z[i, k] ~ dbern(Palive[i, k - 1])
      
      theta[i, k - 1] ~ dunif(-3.141593, 3.141593) # Prior for dispersal direction 
      Palive[i, k - 1] <- z[i, k - 1] * phi[sex[i], k - 1] # Pr(alive in primary session k) 
      
      d[i, k - 1] ~ dexp(dlambda[sex[i]])
      S[i, 1, k] <- S[i, 1, k - 1] + d[i, k - 1] * cos(theta[i, k - 1])
      S[i, 2, k] <- S[i, 2, k - 1] + d[i, k - 1] * sin(theta[i, k - 1])
      g[i, k, 1] <- 0
      for(r in 1:R){ # trap
        D[i, r, k] <- sqrt(pow(S[i, 1, k] - X[r, 1], 2) + pow(S[i, 2, k] - X[r, 2], 2))  # Squared distance to trap
        g[i, k, r + 1] <- exp(-pow(D[i, r, k] / sigma[sex[i]], kappa[sex[i]])) # Trap exposure
      }
      G[i, k] <- sum(g[i, k, 1:(R + 1)]) # Total trap exposure
      for(j in 1:J[i, k]){
        P[i, j, k] <- (1 - exp(-lambda[sex[i]] * G[i, k])) * z[i, k] # Probability of being captured
        captureProb[i, k, j] <- step(H[i, j, k] - 2) * (g[i, k, H[i, j, k]] / (G[i, k] + 0.000000001)) * P[i, j, k] + (1 - step(H[i, j, k] - 2)) * (1 - P[i, j, k])
        Ones[i, j, k] ~ dbern(captureProb[i, k, j])
      }
    }
  }
}
)

##############################################
### FUNCTION FOR GENERATING INITIAL VALUES ###
##############################################
#dt = rep(1, times = n.prim - 1)  # time between primary sessions
n = length(first)
#gr = rep(1, length(first))          

# make all 1:
#tod = matrix(1, ncol=n.sec, nrow=n.prim, byrow=T)

source("InitsForEGBadgersSex.R")

# Random initial centres of activities
trap.dist <- 4
X.max = Y.max = 5
bord = 2*trap.dist # How far outside the trapping grid can first centre og activity be located?
S.first = cbind(runif(n, -bord, X.max+bord), runif(n, -bord, Y.max+bord))

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
               H = H,
               death.occasion = death.occasion,
               sex = sex)

data = list(Ones = array(1, dim(H)), z = z_data)

#inits = list(Inits(H, X = X, K = K))
#inits$S <- Sinit_array
#inits$z <- z_init

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

#model$initializeInfo()

cModel <- compileNimble(model)

config <- configureMCMC(model, monitors = c("phi", "kappa", "sigma", "PL", "dmean"))

rMCMC <- buildMCMC(config)
cMCMC <- compileNimble(rMCMC)

# Run the MCMC
system.time(run <- runMCMC(cMCMC,
                           niter = 70000, 
                           nburnin = 19000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE))

run$summary
saveRDS(run, "EG_Badgers_SexOut.rds")
run <- readRDS("EG_Badgers_SexOut.rds")


mcmcplot(run$samples, parms = c("phi", "kappa", "sigma", "PL", "dmean"))
plot(run$samples)

