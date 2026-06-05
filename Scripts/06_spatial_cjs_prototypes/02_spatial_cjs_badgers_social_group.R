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
  filter(as.numeric(as.factor(socg)) <= 14) %>%
  mutate(count = n()) %>%
  filter(count > 2) %>%
  ungroup() %>%
  droplevels() %>%
  mutate(state = as.numeric(socg))

# Select a random sample of 100 unique tattoos
sampled_tattoos <- bc %>%
  distinct(tattoo) %>%
  sample_n(100)

# Filter the original dataset to include only the sampled tattoos
bc <- bc %>%
  filter(tattoo %in% sampled_tattoos$tattoo)

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

# Initialize an empty list to store the vectors
#po_list <- vector("list", 14)

# Loop to fill the list with the required sequences
#for (j in 1:14) {
#  po_list[[j]] <- setdiff(1:14, j)
#}

# Create an empty matrix to store the result
#po_mat <- matrix(NA, nrow = 14, ncol = 13)

# Fill the matrix with the values from the list
#for (j in 1:14) {
#  po_mat[j, ] <- po_list[[j]]
#}

code <- nimbleCode({
  
  # Priors and constraints
  for (k in 1:nSG) {
    phi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  
  # Transition probabilities
  for (j in 1:(nSG - 1)) { # the last site transition will be calculated as 1 minus the sum of others
    lpsiA[j] ~ dnorm(0, 0.001)
    lpsiB[j] ~ dnorm(0, 0.001)
    lpsiC[j] ~ dnorm(0, 0.001)
    lpsiD[j] ~ dnorm(0, 0.001)
    lpsiE[j] ~ dnorm(0, 0.001)
    lpsiF[j] ~ dnorm(0, 0.001)
    lpsiG[j] ~ dnorm(0, 0.001)
    lpsiH[j] ~ dnorm(0, 0.001)
    lpsiI[j] ~ dnorm(0, 0.001)
    lpsiJ[j] ~ dnorm(0, 0.001)
    lpsiK[j] ~ dnorm(0, 0.001)
    lpsiL[j] ~ dnorm(0, 0.001)
    lpsiM[j] ~ dnorm(0, 0.001)
    lpsiN[j] ~ dnorm(0, 0.001)
  }
  
  for (i in 1:(nSG - 1)){
    psiA[i] <- exp(lpsiA[i]) / (1 + sum(exp(lpsiA[1:13])))
    psiB[i] <- exp(lpsiB[i]) / (1 + sum(exp(lpsiB[1:13])))
    psiC[i] <- exp(lpsiC[i]) / (1 + sum(exp(lpsiC[1:13])))
    psiD[i] <- exp(lpsiD[i]) / (1 + sum(exp(lpsiD[1:13])))
    psiE[i] <- exp(lpsiE[i]) / (1 + sum(exp(lpsiE[1:13])))
    psiF[i] <- exp(lpsiF[i]) / (1 + sum(exp(lpsiF[1:13])))
    psiG[i] <- exp(lpsiG[i]) / (1 + sum(exp(lpsiG[1:13])))
    psiH[i] <- exp(lpsiH[i]) / (1 + sum(exp(lpsiH[1:13])))
    psiI[i] <- exp(lpsiI[i]) / (1 + sum(exp(lpsiI[1:13])))
    psiJ[i] <- exp(lpsiJ[i]) / (1 + sum(exp(lpsiJ[1:13])))
    psiK[i] <- exp(lpsiK[i]) / (1 + sum(exp(lpsiK[1:13])))
    psiL[i] <- exp(lpsiL[i]) / (1 + sum(exp(lpsiL[1:13])))
    psiM[i] <- exp(lpsiM[i]) / (1 + sum(exp(lpsiM[1:13])))
    psiN[i] <- exp(lpsiN[i]) / (1 + sum(exp(lpsiN[1:13])))
  }
  
  psiA[14] <- 1 - sum(psiA[1:13])
  psiB[14] <- 1 - sum(psiB[1:13])
  psiC[14] <- 1 - sum(psiC[1:13])
  psiD[14] <- 1 - sum(psiD[1:13])
  psiE[14] <- 1 - sum(psiE[1:13])
  psiF[14] <- 1 - sum(psiF[1:13])
  psiG[14] <- 1 - sum(psiG[1:13])
  psiH[14] <- 1 - sum(psiH[1:13])
  psiI[14] <- 1 - sum(psiI[1:13])
  psiJ[14] <- 1 - sum(psiJ[1:13])
  psiK[14] <- 1 - sum(psiK[1:13])
  psiL[14] <- 1 - sum(psiL[1:13])
  psiM[14] <- 1 - sum(psiM[1:13])
  psiN[14] <- 1 - sum(psiN[1:13])
  
  
  for (i in 1:nind) { 
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions - 1)) {
      for (j in 1:14) {
        ps[1, i, t, j] <- phi[1] * psiA[j]
        ps[2, i, t, j] <- phi[2] * psiB[j]
        ps[3, i, t, j] <- phi[3] * psiC[j]
        ps[4, i, t, j] <- phi[4] * psiD[j]
        ps[5, i, t, j] <- phi[5] * psiE[j]
        ps[6, i, t, j] <- phi[6] * psiF[j]
        ps[7, i, t, j] <- phi[7] * psiG[j]
        ps[8, i, t, j] <- phi[8] * psiH[j]
        ps[9, i, t, j] <- phi[9] * psiI[j]
        ps[10, i, t, j] <- phi[10] * psiJ[j]
        ps[11, i, t, j] <- phi[11] * psiK[j]
        ps[12, i, t, j] <- phi[12] * psiL[j]
        ps[13, i, t, j] <- phi[13] * psiM[j]
        ps[14, i, t, j] <- phi[14] * psiN[j]
        ps[15, i, t, j] <- 0
        ps[j, i, t, 15] <- 1 - phi[j]
      }
      ps[15, i, t, 15] <- 1
      # Define probabilities of O(t) given S(t)
      
      for(j in 1:14){
        po[j, i, t, j] <- p[j]
        po[15, i, t, j] <- 0
        po[j, i, t, 15] <- 1 - p[j]
      }
      
      for(j in 2:14){
        po[1, i, t, j] <- 0
      }
      
      po[2, i, t, 1] <- 0
      for(j in 3:14){
        po[2, i, t, j] <- 0
      }
      
      for(j in 1:2){
        po[3, i, t, j] <- 0
      }
      for(j in 4:14){
        po[3, i, t, j] <- 0
      }
      
      for(j in 1:3){
        po[4, i, t, j] <- 0
      }
      for(j in 5:14){
        po[4, i, t, j] <- 0
      }
      
      for(j in 1:4){
        po[5, i, t, j] <- 0
      }
      for(j in 6:14){
        po[5, i, t, j] <- 0
      }
      
      for(j in 1:5){
        po[6, i, t, j] <- 0
      }
      for(j in 7:14){
        po[6, i, t, j] <- 0
      }
      
      for(j in 1:6){
        po[7, i, t, j] <- 0
      }
      for(j in 8:14){
        po[7, i, t, j] <- 0
      }
      
      for(j in 1:7){
        po[8, i, t, j] <- 0
      }
      for(j in 9:14){
        po[8, i, t, j] <- 0
      }
      
      for(j in 1:8){
        po[9, i, t, j] <- 0
      }
      for(j in 10:14){
        po[9, i, t, j] <- 0
      }
      
      for(j in 1:9){
        po[10, i, t, j] <- 0
      }
      for(j in 11:14){
        po[10, i, t, j] <- 0
      }
      
      for(j in 1:10){
        po[11, i, t, j] <- 0
      }
      for(j in 12:14){
        po[11, i, t, j] <- 0
      }
      
      for(j in 1:11){
        po[12, i, t, j] <- 0
      }
      for(j in 13:14){
        po[12, i, t, j] <- 0
      }
      
      for(j in 1:12){
        po[13, i, t, j] <- 0
      }
      po[13, i, t, 14] <- 0
    } #t
  } #i
  
  # Likelihood
  for (i in 1:nind) {
    # Define latent state at first capture
    z[i, f[i]] <- y[i, f[i]]
    for (t in (f[i] + 1):(n.occasions)) {
      # State process: draw S(t) given S(t-1)
      z[i, t] ~ dcat(ps[z[i, t-1], i, t-1, 1:15])
      # Observation process: draw O(t) given S(t)
      y[i, t] ~ dcat(po[z[i, t], i, t-1, 1:15])
    } # t
  } # i
})

consts <- list(f = f, n.occasions = ncol(capture_history), nind = nrow(capture_history), nSG = 14)

zinits <- capture_history
zinits[zinits == 0] <- sample(c(1:14), sum(zinits==0), replace = TRUE)
capture_history[capture_history == 0] <- 15

data <- list(y = capture_history)

inits <- list(phi = runif(14, 0, 1),
              p = runif(14, 0, 1),
              lpsiA = rnorm(13, 0, 1),
              lpsiB = rnorm(13, 0, 1),
              lpsiC = rnorm(13, 0, 1),
              lpsiD = rnorm(13, 0, 1),
              lpsiE = rnorm(13, 0, 1),
              lpsiF = rnorm(13, 0, 1),
              lpsiG = rnorm(13, 0, 1),
              lpsiH = rnorm(13, 0, 1),
              lpsiI = rnorm(13, 0, 1),
              lpsiJ = rnorm(13, 0, 1),
              lpsiK = rnorm(13, 0, 1),
              lpsiL = rnorm(13, 0, 1),
              lpsiM = rnorm(13, 0, 1),
              lpsiN = rnorm(13, 0, 1),
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
