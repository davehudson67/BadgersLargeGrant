install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)

# Load required packages
library(sf)
library(dplyr)
library(tidyverse)
library(nimble)
library(coda)
library(spdep)
library(patchwork)
library(GGally)
library(ggmcmc)
library(INLA)
library(kableExtra)

# Import the Scottish lip cancer data
ScotLip <- read_sf("Data/scotlip.shp")

# Calculate the standardized incidence ratios (SIRs)
ScotLip$crudeSIR <- ScotLip$CANCER / ScotLip$CEXP

# Plot the crude SIRs
p1 <- ggplot() + geom_boxplot(data = ScotLip, aes(y = crudeSIR)) + theme_bw() + theme(text = element_text(size = 12))
p2 <- ggplot() + geom_sf(data = ScotLip, aes(fill = crudeSIR), col = "NA") + scale_fill_viridis_c() + theme_bw() + theme(axis.text = element_blank(), text = element_text(size = 12))
(p1 | p2)

# Model specification for Unstructured random effect using NIMBLE
UnstrCode <- nimbleCode({
  for (i in 1:N){
    O[i] ~ dpois(mu[i])
    log(mu[i]) <- log(E[i]) + alpha + theta[i]
    theta[i] ~ dnorm(0, tau = tau.theta)
    SIR[i] <- exp(alpha + theta[i])
    resSIR[i] <- exp(theta[i])
    e[i] <- (O[i]-mu[i])/sqrt(mu[i])
  }
  alpha ~ dflat()
  overallSIR <- exp(alpha)
  tau.theta ~ dgamma(1, 0.01)
  sigma2.theta <- 1/tau.theta
})

# Data and constants for NIMBLE
N <- dim(ScotLip)[1]
ScotLipdata <- list(O = ScotLip$CANCER)
ScotLipConsts <- list(N = N, E = ScotLip$CEXP)

# Initial values
inits <- list(
  list(alpha = 0.01, tau.theta = 10, theta = rep(0.01, times = N)),
  list(alpha = 0.5, tau.theta = 1, theta = rep(-0.01, times = N))
)

# Monitored parameters
params <- c("sigma2.theta", "overallSIR", "resSIR", "SIR", "e", "mu", "alpha", "theta")

# MCMC settings
ni <- 50000
nt <- 10
nb <- 10000
nc <- 2

# Run the MCMC simulations
UnstrCodesamples <- nimbleMCMC(code = UnstrCode, 
                               data = ScotLipdata, 
                               constants = ScotLipConsts, 
                               inits = inits, 
                               monitors = params, 
                               niter = ni, 
                               nburnin = nb, 
                               thin = nt, 
                               nchains = nc, 
                               setSeed = 9, 
                               progressBar = FALSE, 
                               samplesAsCodaMCMC = TRUE, 
                               summary = TRUE, 
                               WAIC = TRUE)

# Check convergence
GR.diag <- gelman.diag(UnstrCodesamples$samples, multivariate = FALSE)
all(GR.diag$psrf[,"Point est."] < 1.1)

# Extract results
ScotLip$unstr_SIR <- UnstrCodesamples$summary$all.chains[paste0("SIR[", 1:N, "]"), "Median"]

# Plot smoothed SIRs
p1 <- ggplot() + geom_sf(data = ScotLip, aes(fill = unstr_SIR), col = NA) + theme_bw() + scale_fill_viridis_c(limits = c(0,5)) + theme(axis.text = element_blank(), text = element_text(size = 12))
p1

# Posterior probability that the spatial SIR per area is larger than 1
postProb <- sapply(paste0("SIR[", 1:N, "]"), function(X) mean(UnstrCodesamples$samples$chain1[,X] > 1))
ScotLip$unstr_postProb <- postProb

p1 <- ggplot() + geom_sf(data = ScotLip, aes(fill = postProb), col = NA) + theme_bw() + scale_fill_viridis_c(limits = c(0,1)) + theme(axis.text = element_blank(), text = element_text(size = 12))
p2 <- ggplot() + geom_boxplot(data = ScotLip, aes(y = postProb)) + scale_fill_viridis_d(alpha = .5) + theme_bw() + xlim(c(-1,1)) + theme(axis.text.x = element_blank(), text = element_text(size = 12))
(p1 | p2) + plot_annotation(title = 'Posterior probability that Pr(SIR>1)')
