# Set parameters
num_individuals <- 10
num_occasions <- 20
num_sites <- 3

# Create an empty matrix to hold the CMR data
cmr_data <- matrix(0, nrow = num_individuals, ncol = num_occasions)

# Simulate data
set.seed(123)  # Set seed for reproducibility
for (i in 1:num_individuals) {
  for (j in 1:num_occasions) {
    cmr_data[i, j] <- sample(0:num_sites, 1)
  }
}

## Find first capture event
get.first <- function(x) min(which(x != 0))
f <- apply(cmr_data, 1, get.first)

code <- nimbleCode({
  
  # Priors and constraints
  for (k in 1:num_sites) {
    phi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  # Transitions
  for (i in 1:(num_sites - 1)) {
    lpsiA[i] ~ dnorm(0, 0.001)
    lpsiB[i] ~ dnorm(0, 0.001)
    lpsiC[i] ~ dnorm(0, 0.001)
  }
  # Constrain the transitions
  for (i in 1:(num_sites - 1)) {
    psiA[i] <- exp(lpsiA[i]) / (1 + sum(exp(lpsiA[1:2])))
    psiB[i] <- exp(lpsiB[i]) / (1 + sum(exp(lpsiB[1:2])))
    psiC[i] <- exp(lpsiC[i]) / (1 + sum(exp(lpsiC[1:2])))
  }
  
  psiA[3] <- 1 - sum(psiA[1:2])
  psiB[3] <- 1 - sum(psiB[1:2])
  psiC[3] <- 1 - sum(psiC[1:2])
  
  # Define state transitions and obs matrices
  for (i in 1:nind){
    for(t in f[i]:(n.occasions - 1)){
      ps[1, i, t, 1] <- phi[1] * psiA[1]
      ps[1, i, t, 2] <- phi[1] * psiA[2]
      ps[1, i, t, 3] <- phi[1] * psiA[3]
      ps[1, i, t, 4] <- 1 - phi[1]
      ps[2, i, t, 1] <- phi[2] * psiB[1]
      ps[2, i, t, 2] <- phi[2] * psiB[2]
      ps[2, i, t, 3] <- phi[2] * psiB[3]
      ps[2, i, t, 4] <- 1 - phi[2]
      ps[3, i, t, 1] <- phi[3] * psiC[1]
      ps[3, i, t, 2] <- phi[3] * psiC[2]
      ps[3, i, t, 3] <- phi[3] * psiC[3]
      ps[3, i, t, 4] <- 1 - phi[3]
      ps[4, i, t, 1] <- 0
      ps[4, i, t, 2] <- 0
      ps[4, i, t, 3] <- 0
      ps[4, i, t, 4] <- 1
      
      po[1, i, t, 1] <- p[1]
      po[1, i, t, 2] <- 0
      po[1, i, t, 3] <- 0
      po[1, i, t, 4] <- 1 - p[1]
      po[2, i, t, 1] <- 0
      po[2, i, t, 2] <- p[2]
      po[2, i, t, 3] <- 0
      po[2, i, t, 4] <- 1 - p[2]
      po[3, i, t, 1] <- 0
      po[3, i, t, 2] <- 0
      po[3, i, t, 3] <- p[3]
      po[3, i, t, 4] <- 1 - p[3]
      po[4, i, t, 1] <- 0
      po[4, i, t, 2] <- 0
      po[4, i, t, 3] <- 0
      po[4, i, t, 4] <- 1
      }
  }
  
  # Likelihood
  for (i in 1:nind){
    z[i, f[i]] <- y[i, f[i]]
    for (t in (f[i] + 1):n.occasions){
      z[i, t] ~ dcat(ps[z[i, t-1], i, t-1, 1:4])
      y[i, t] ~ dcat(po[z[i, t], i, t-1, 1:4])
    }
  }
})

consts <- list(f = f, n.occasions = num_occasions, nind = num_individuals, num_sites = num_sites)

zinits <- cmr_data
zinits[zinits == 0] <- sample(c(1:3), sum(zinits==0), replace = TRUE)
cmr_data[cmr_data == 0] <- 4

data <- list(y = cmr_data)

inits <- list(phi = runif(3, 0, 1),
              p = runif(3, 0, 1),
              lpsiA = rnorm(2, 0, 1),
              lpsiB = rnorm(2, 0, 1),
              lpsiC = rnorm(2, 0, 1),
              z = zinits)

model <- nimbleModel(code, constants = consts, data = data, inits = inits)
model$initializeInfo()

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
