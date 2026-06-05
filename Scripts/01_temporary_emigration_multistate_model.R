library(nimble)
library(tidyverse)
library(MCMCpack)

set.seed(11)
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
  mutate(state = 1)

# Select a random sample of 100 unique tattoos
sampled_tattoos <- bc %>%
  distinct(tattoo) %>%
  sample_n(10)

# Filter the original dataset to include only the sampled tattoos
bc <- bc %>%
  filter(tattoo %in% sampled_tattoos$tattoo)

# Get unique individual IDs and time occasions
ids <- unique(bc$tattoo)
occs <- seq(1:(max(bc$occ)))

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
get.last <- function(x) max(which(x != 0))
l <- apply(capture_history, 1, get.last)

capture_history[capture_history == 0] <- 2
nind <- nrow(capture_history)

## SG changes
# Function to count SG changes for each individual
sg_changes <- function(sg) {
  return(sum(sg[-1] != sg[-length(sg)]))
}

# Summarize the DataFrame to count SG changes for each individual
result <- bc %>%
  group_by(tattoo) %>%
  summarize(Captures = n(),
            SG_Changes = sg_changes(socg))

# View the result
print(result)
ggplot(result, aes(x = Captures, y = SG_Changes)) +
  geom_point(position = "jitter", alpha = 0.2)

hist(result$SG_Changes)

## Proportion that dont change SG in capture history
length(result$SG_Changes[result$SG_Changes == 0])/length(result$tattoo)

# Temporary Emigration Model P282 in BPA book
TEcode <- nimbleCode({
  # Priors and constraints
  for (t in 1:(n.occasions-1)){
    phi[t] <- mean.phi
    psiIO[t] <- mean.psiIO
    psiOI[t] <- mean.psiOI
    p[t] <- mean.p
  }
  mean.phi ~ dunif(0, 1)       # Prior for mean survival
  mean.psiIO ~ dunif(0, 1)     # Prior for mean temp. emigration
  mean.psiOI ~ dunif(0, 1)     # Prior for mean temp. immigration
  mean.p ~ dunif(0, 1)         # Prior for mean recapture
  
  # Define state-transition and observation matrices 	
  for (i in 1:nind){
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions-1)){
      ps[1,i,t,1] <- phi[t] * (1-psiIO[t])
      ps[1,i,t,2] <- phi[t] * psiIO[t]
      ps[1,i,t,3] <- 1-phi[t]
      ps[2,i,t,1] <- phi[t] * psiOI[t]
      ps[2,i,t,2] <- phi[t] * (1-psiOI[t])
      ps[2,i,t,3] <- 1-phi[t]
      ps[3,i,t,1] <- 0
      ps[3,i,t,2] <- 0
      ps[3,i,t,3] <- 1
      
      # Define probabilities of O(t) given S(t)
      po[1,i,t,1] <- p[t]
      po[1,i,t,2] <- 1-p[t]
      po[2,i,t,1] <- 0
      po[2,i,t,2] <- 1
      po[3,i,t,1] <- 0
      po[3,i,t,2] <- 1
    } #t
  } #i
  
  # Likelihood
  for (i in 1:nind){
    # Define latent state at first capture
    z[i,f[i]] <- y[i,f[i]]
    for (t in (f[i]+1):(n.occasions)){
      # State process: draw S(t) given S(t-1)
      z[i,t] ~ dcat(ps[z[i,t-1], i, t-1, 1:3])
      # Observation process: draw O(t) given S(t)
      y[i,t] ~ dcat(po[z[i,t], i, t-1, 1:2])
    } #t
  } #i
})

## set up other components of model
consts <- list(nind = nind, f = f, n.occasions = ncol(capture_history))

data <- list(y = capture_history)

# Function to create initial values for unknown z
z_init <- capture_history
for (i in 1:nrow(capture_history)){
  z_init[i, (f[i] : l[i])] <- 1
}

inits <- list(mean.phi = runif(1, 0, 1),
              mean.psiIO = runif(1, 0, 1),
              mean.psiOI = runif(1, 0, 1),
              mean.p = runif(1, 0, 1),
              z = z_init
              )

## define the model, data, inits and constants
model <- nimbleModel(code = TEcode, constants = consts, data = data, inits = inits, calculate = FALSE)
model$initializeInfo()
## compile the model
cmodel <- compileNimble(model)

## try with adaptive slice sampler
config <- configureMCMC(cmodel, monitors = c("mean.phi", "mean.psiIO", "mean.psiOI", "mean.p"), thin = 1)
#config$removeSamplers(c("a", "b"))
#config$addSampler(target = c("a", "b"), type = 'AF_slice')

## check monitors and samplers
config$printMonitors()
#config$printSamplers(c("a", "b", "mean.p"))

## build the model
built <- buildMCMC(config)
cbuilt <- compileNimble(built)

## run the model
system.time(run <- runMCMC(cbuilt, 
                           niter = 500000, 
                           nburnin = 19000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))

saveRDS(run, "TEModelRun.rds")

run$summary
plot(run$samples)
