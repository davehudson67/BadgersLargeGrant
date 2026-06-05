library(nimble)
library(tidyverse)
library(MCMCpack)

bc <- readRDS("../BadgersSG2024.rds")
bc$socg <- as.factor(bc$socg)
levels(bc$socg)

bc <- bc %>%
  group_by(tattoo) %>%
  filter(pm != "Yes") %>%
  filter(!is.na(socg)) %>%
  mutate(count = n()) %>%
  filter(count > 2) %>%
  ungroup() %>%
  droplevels() %>%
  mutate(state = as.numeric(socg))

nSG <- max(bc$state)

# Get unique individual IDs and time occasions
ids <- unique(bc$tattoo)
occs <- seq(1:(max(bc$occ) + 2))

# Initialize capture history array
capture_history <- matrix(0, nrow = length(ids), ncol = length(occs))
rownames(capture_history) <- ids
colnames(capture_history) <- occs

## Fill the capture history matrix with the state values
for (i in 1:nrow(bc)) {
  ind <- bc$tattoo[i]
  occ <- bc$occ[i]
  state <- bc$state[i]
  capture_history[ind, occ] <- state
}

## Find first capture event
get.first <- function(x) min(which(x != 0))
f <- apply(capture_history, 1, get.first)

## Define multistate model in NIMBLE
code <- nimbleCode({
  
  # Priors and constraints
  for (k in 1:nSG) {
    phi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  
  # Transition probabilities
  for (j in 1:(nSG - 1)) {
    for (k in 1:14) {
      lpsi[j, k] ~ dnorm(0, 0.001)
    }
  }
  
  # Loop to calculate psi values
  for (i in 1:(nSG - 1)) {
    for (k in 1:14) {
      psi[i, k] <- exp(lpsi[i, k]) / (1 + sum(exp(lpsi[1:(nSG - 1), k])))
    }
  }
  
  # Calculate the last state as 1 minus the sum of others
  for (k in 1:14) {
    psi[nSG, k] <- 1 - sum(psi[1:(nSG - 1), k])
  }
  
  for (i in 1:nind) {
    for (t in f[i]:(n.occasions - 1)) {
      for (j in 1:nSG) {
        for (k in 1:nSG) {
          ps[j, i, t, k] <- phi[j] * psi[j, k]
        }
        ps[nSG + 1, i, t, j] <- 0 # Set the last row except the last column to 0
        ps[j, i, t, nSG + 1] <- 1 - phi[j] # Set the last column for each row
      }
      ps[nSG + 1, i, t, nSG + 1] <- 1 # Set the bottom-right cell to 1
      
      
      # Define probabilities of O(t) given S(t)
      for (j in 1:43) {
        po[j, i, t, j] <- p[j]
        po[44, i, t, j] <- 0
        po[j, i, t, 44] <- 1 - p[j]
      }
      
      for (j in 2:43) {
        po[1, i, t, j] <- 0
      }
      
      po[2, i, t, 1] <- 0
      for (j in 3:43) {
        po[2, i, t, j] <- 0
      }
      
      for (j in 1:2) {
        po[3, i, t, j] <- 0
      }
      for (j in 4:43) {
        po[3, i, t, j] <- 0
      }
      
      for (j in 1:3) {
        po[4, i, t, j] <- 0
      }
      for (j in 5:43) {
        po[4, i, t, j] <- 0
      }
      
      for (j in 1:4) {
        po[5, i, t, j] <- 0
      }
      for (j in 6:43) {
        po[5, i, t, j] <- 0
      }
      
      for (j in 1:5) {
        po[6, i, t, j] <- 0
      }
      for (j in 7:43) {
        po[6, i, t, j] <- 0
      }
      
      for (j in 1:6) {
        po[7, i, t, j] <- 0
      }
      for (j in 8:43) {
        po[7, i, t, j] <- 0
      }
      
      for (j in 1:7) {
        po[8, i, t, j] <- 0
      }
      for (j in 9:43) {
        po[8, i, t, j] <- 0
      }
      
      for (j in 1:8) {
        po[9, i, t, j] <- 0
      }
      for (j in 10:43) {
        po[9, i, t, j] <- 0
      }
      
      for (j in 1:9) {
        po[10, i, t, j] <- 0
      }
      for (j in 11:43) {
        po[10, i, t, j] <- 0
      }
      
      for (j in 1:10) {
        po[11, i, t, j] <- 0
      }
      for (j in 12:43) {
        po[11, i, t, j] <- 0
      }
      
      for (j in 1:11) {
        po[12, i, t, j] <- 0
      }
      for (j in 13:43) {
        po[12, i, t, j] <- 0
      }
      
      for (j in 1:12) {
        po[13, i, t, j] <- 0
      }
      for (j in 14:43) {
        po[13, i, t, j] <- 0
      }
      
      for (j in 1:13) {
        po[14, i, t, j] <- 0
      }
      for (j in 15:43) {
        po[14, i, t, j] <- 0
      }
      
      for (j in 1:14) {
        po[15, i, t, j] <- 0
      }
      for (j in 16:43) {
        po[15, i, t, j] <- 0
      }
      
      for (j in 1:15) {
        po[16, i, t, j] <- 0
      }
      for (j in 17:43) {
        po[16, i, t, j] <- 0
      }
      
      for (j in 1:16) {
        po[17, i, t, j] <- 0
      }
      for (j in 18:43) {
        po[17, i, t, j] <- 0
      }
      
      for (j in 1:17) {
        po[18, i, t, j] <- 0
      }
      for (j in 19:43) {
        po[18, i, t, j] <- 0
      }
      
      for (j in 1:18) {
        po[19, i, t, j] <- 0
      }
      for (j in 20:43) {
        po[19, i, t, j] <- 0
      }
      
      for (j in 1:19) {
        po[20, i, t, j] <- 0
      }
      for (j in 21:43) {
        po[20, i, t, j] <- 0
      }
      
      for (j in 1:20) {
        po[21, i, t, j] <- 0
      }
      for (j in 22:43) {
        po[21, i, t, j] <- 0
      }
      
      for (j in 1:21) {
        po[22, i, t, j] <- 0
      }
      for (j in 23:43) {
        po[22, i, t, j] <- 0
      }
      
      for (j in 1:22) {
        po[23, i, t, j] <- 0
      }
      for (j in 24:43) {
        po[23, i, t, j] <- 0
      }
      
      for (j in 1:23) {
        po[24, i, t, j] <- 0
      }
      for (j in 25:43) {
        po[24, i, t, j] <- 0
      }
      
      for (j in 1:24) {
        po[25, i, t, j] <- 0
      }
      for (j in 26:43) {
        po[25, i, t, j] <- 0
      }
      
      for (j in 1:25) {
        po[26, i, t, j] <- 0
      }
      for (j in 27:43) {
        po[26, i, t, j] <- 0
      }
      
      for (j in 1:26) {
        po[27, i, t, j] <- 0
      }
      for (j in 28:43) {
        po[27, i, t, j] <- 0
      }
      
      for (j in 1:27) {
        po[28, i, t, j] <- 0
      }
      for (j in 29:43) {
        po[28, i, t, j] <- 0
      }
      
      for (j in 1:28) {
        po[29, i, t, j] <- 0
      }
      for (j in 30:43) {
        po[29, i, t, j] <- 0
      }
      
      for (j in 1:29) {
        po[30, i, t, j] <- 0
      }
      for (j in 31:43) {
        po[30, i, t, j] <- 0
      }
      
      for (j in 1:30) {
        po[31, i, t, j] <- 0
      }
      for (j in 32:43) {
        po[31, i, t, j] <- 0
      }
      
      for (j in 1:31) {
        po[32, i, t, j] <- 0
      }
      for (j in 33:43) {
        po[32, i, t, j] <- 0
      }
      
      for (j in 1:32) {
        po[33, i, t, j] <- 0
      }
      for (j in 34:43) {
        po[33, i, t, j] <- 0
      }
      
      for (j in 1:33) {
        po[34, i, t, j] <- 0
      }
      for (j in 35:43) {
        po[34, i, t, j] <- 0
      }
      
      for (j in 1:34) {
        po[35, i, t, j] <- 0
      }
      for (j in 36:43) {
        po[35, i, t, j] <- 0
      }
      
      for (j in 1:35) {
        po[36, i, t, j] <- 0
      }
      for (j in 37:43) {
        po[36, i, t, j] <- 0
      }
      
      for (j in 1:36) {
        po[37, i, t, j] <- 0
      }
      for (j in 38:43) {
        po[37, i, t, j] <- 0
      }
      
      for (j in 1:37) {
        po[38, i, t, j] <- 0
      }
      for (j in 39:43) {
        po[38, i, t, j] <- 0
      }
      
      for (j in 1:38) {
        po[39, i, t, j] <- 0
      }
      for (j in 40:43) {
        po[39, i, t, j] <- 0
      }
      
      for (j in 1:39) {
        po[40, i, t, j] <- 0
      }
      for (j in 41:43) {
        po[40, i, t, j] <- 0
      }
      
      for (j in 1:40) {
        po[41, i, t, j] <- 0
      }
      for (j in 42:43){
        po[41, i, t, j] <- 0
      }
      
      for (j in 1:41) {
        po[42, i, t, j] <- 0
      }
      po[42, i, t, 43] <- 0
      
    } #t
  } #i
  
  # Likelihood
  for (i in 1:nind) {
    # Define latent state at first capture
    z[i, f[i]] <- y[i, f[i]]
    for (t in (f[i] + 1):(n.occasions)) {
      # State process: draw S(t) given S(t-1)
      z[i, t] ~ dcat(ps[z[i, t-1], i, t-1, 1:44])
      # Observation process: draw O(t) given S(t)
      y[i, t] ~ dcat(po[z[i, t], i, t-1, 1:44])
    } # t
  } # i
})

consts <- list(f = f, n.occasions = ncol(capture_history), nind = nrow(capture_history), nSG = nSG)

zinits <- capture_history
zinits[zinits == 0] <- sample(c(1:43), sum(zinits==0), replace = TRUE)
capture_history[capture_history == 0] <- 44

data <- list(y = capture_history)

lpsi_init <- array(NA, dim = c(nSG - 1, 14))

# Assign initial values
for (j in 1:(nSG - 1)) {
  for (k in 1:14) {
    lpsi_init[j, k] <- rnorm(1, 0, 0.001) # Initial values can be set from a normal distribution
  }
}

inits <- list(phi = runif(43, 0, 1),
              p = runif(43, 0, 1),
              lpsi = lpsi_init,
              z = zinits)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)

cModel <- compileNimble(model)
config <- configureMCMC(model)
rMCMC <- buildMCMC(config)
cMCMC <- compileNimble(rMCMC, project = model)

# Run the MCMC
system.time(run <- runMCMC(cMCMC, 
                           niter = 10000, 
                           nburnin = 2400, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
plot(run$samples)
