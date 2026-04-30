library(nimble)
library(tidyverse)

# Parameters:
# rel: number of newly released individuals at each occasion
# nyears: number of capture occasions
# s: survival probability (we assume here that it is constant over time)
# p: recapture probability (we assume here that it is constant over time)
# disp.varX: variance of the dispersal distribution in x-direction
# disp.varY: variance of the dispersal distribution in y-direction
# grid.x.study.area: x-coordinate of the grid cells
# grid.y.study.area: y-coordinate of the grid cells
# grid.size: side length of the quadratic grid cells

# Output parameters:
# Y.in: capture-recapture data sampled within the study area
# Y.comp: all capture-recapture data (i.e. not restricted to the study area)
# G.in: X and Y-coordinates of the encountered individuals within the study area
# G.comp: X and Y-coordinates of all encountered individuals (i.e. not restricted to the study area)
# f: vector with the release occasion for each individual

bc <- readRDS("../BadgersSG2024.rds")
load("spatialInfo.RData")
setts <- settGrid$Sett
levels(as.factor(bc$sett))
bc$Sett <- iconv(bc$sett, from = "latin1", to = "UTF-8", sub = "")
bc$Sett <- gsub(" ", "", bc$Sett)

settGrid <- select(settGrid, -geometry)
studyArea <- select(studyArea, - geometry)

bc <- bc %>%
  filter(Sett %in% setts)

bc <- bc %>%
  left_join(settGrid, by = "Sett")

# Get unique individual IDs and time occasions
ids <- unique(bc$tattoo)
occs <- seq(1:(max(bc$occ) + 2))

# Initialize capture history array
capture_history <- matrix(0, nrow = length(ids), ncol = length(occs))
rownames(capture_history) <- ids
colnames(capture_history) <- occs

# Populate capture history array
for (i in 1:nrow(bc)) {
  ind <- as.character(bc$tattoo[i])
  time <- as.character(bc$occ[i])
  capture_history[ind, time] <- 1
}

# Initialize location arrays for x and y coordinates
loc_x <- array(NA, dim = c(length(ids), length(occs)))
loc_y <- array(NA, dim = c(length(ids), length(occs)))
rownames(loc_x) <- ids
colnames(loc_x) <- occs
rownames(loc_y) <- ids
colnames(loc_y) <- occs

# Populate location arrays
for (i in 1:nrow(bc)) {
  ind <- as.character(bc$tattoo[i])
  time <- as.character(bc$occ[i])
  loc_x[ind, time] <- bc$col_index[i]
  loc_y[ind, time] <- bc$row_index[i]
}

grid_size <- 1

# Create the first capture vector
first_capture <- apply(capture_history, 1, function(x) which(x == 1)[1])

# Combine location arrays into a single 3D array G
G <- array(NA, dim = c(length(ids), length(occs), 2))
G[,,1] <- loc_x
G[,,2] <- loc_y

# sCJS-N with irregularly shaped study area 
code <- nimbleCode({
  # Priors and constraints
  for (i in 1:nind){
    for (t in f[i]:(n.occasions-1)){
      phi[i,t] <- mean.phi
      p[i,t] <- mean.p
    } #t
  } #i
  
  mean.phi ~ dunif(0,1)
  mean.p ~ dunif(0,1)
  
  for (i in 1:2){
    tau[i] <- pow(sigma[i], -2)
    sigma[i] ~ dunif(0, 50)
  }
  
  # Likelihood 
  for (i in 1:nind){
    # Define latent state at first capture
    z[i, f[i]] <- 1
    
    for (t in (f[i] + 1):n.occasions){
      # State processes
      # Survival
      z[i, t] ~ dbern(phi[i, t - 1] * z[i, t - 1])
      
      # Dispersal
      G[i, t, 1] ~ dnorm(G[i, t - 1, 1], tau[1])
      G[i, t, 2] ~ dnorm(G[i, t - 1, 2], tau[2])
      
      # Observation process
      # Test whether the actual location is in- or outside the state-space
      for (g in 1:ngrids){
        inside[i, t, g] <- step(G[i, t, 1] - grid.x.study.area[g]) * step(grid.x.study.area[g] + grid.size - G[i, t, 1]) *
          step(G[i, t, 2] - grid.y.study.area[g]) * step(grid.y.study.area[g] + grid.size - G[i, t, 2])
      }
      r[i, t] <- sum(inside[i, t, 1:ngrids])
      Y[i, t] ~ dbern(p[i, t - 1] * z[i, t] * r[i, t])
    } #t
  } #i
})

####################################################
# Run the model

# Package data
G[is.na(G)] <- 0

capture_history <- capture_history[1:100,]
G <- G[1:100,,]
first_capture <- first_capture[1:100]

data <- list(Y = capture_history,  G = G)
consts <- list(nind = nrow(capture_history), f = first_capture, n.occasions = ncol(capture_history), 
               grid.x.study.area = studyArea$col_index, grid.y.study.area = studyArea$row_index, grid.size = grid_size, ngrids = nrow(studyArea))

# Function that is required for creating initial values (from Kéry & Schaub 2012)
known.state.cjs <- function(ch){
  state <- ch
  for (i in 1:dim(ch)[1]){
    n1 <- min(which(ch[i,]==1))
    n2 <- max(which(ch[i,]==1))
    state[i,n1:n2] <- 1
    state[i,n1] <- NA
  }
  state[state==0] <- NA
  return(state)
}

# 4. Initial values
inits <- list(mean.phi = runif(1,0,1), mean.p = runif(1,0,1), sigma = runif(2,0.01,2), z = known.state.cjs(capture_history))

# 5. Parameters monitored
parms <- c("mean.phi", "mean.p", "sigma")

# 6. MCMC settings
ni <- 10000
nt <- 1
nb <- 5000
nc <- 2

model <- nimbleModel(code, constants = consts, data = data, inits = inits)
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
run$summary
plot(run$samples)

