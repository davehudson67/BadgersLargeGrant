# Function to simulate a data for the spatial Cormack-Jolly-Seber model
# - Dispersal just with a normal distribution
# - The study area is irregularly shaped

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

simul.sCJSN.grid <- function(rel = 100, nyears = 7, s = 0.7, p = 0.8, disp.varX = 0.5, disp.varY = 0.5, grid.x.study.area, grid.y.study.area, grid.size = 2){
  n <- rel * (nyears-1)   # total number of individuals
  Z <- matrix(NA, ncol = nyears, nrow = n)
  G <- G.in <- G.comp <- array(NA, dim=c(n, nyears, 2))
  f <- rep(1:(nyears-1), rep(rel, (nyears-1)))  # vector indicating when individuals were captured for the first time
  
  # Simulating survival process
  for (i in 1:n){
    Z[i,f[i]] <- 1
    for (t in (f[i]+1):nyears){
      Z[i,t] <- rbinom(1, 1, s) * Z[i,t-1]
    } # t
  } # i
  
  # Simulating dispersal process
  # determine site of first capture
  for (i in 1:n){
    # Select randomly a grid
    grid <- sample(1:length(grid.x.study.area), 1)
    # Generate a random location within the selected grid         
    G[i,f[i],1] <- runif(1, grid.x.study.area[grid], grid.x.study.area[grid]+grid.size)
    G[i,f[i],2] <- runif(1, grid.y.study.area[grid], grid.y.study.area[grid]+grid.size)
  } # i
  # dispersal (assuming a normal distribution in both directions)
  for (i in 1:n){
    for (t in (f[i]+1):nyears){   
      G[i,t,1] <- G[i,t-1,1] + rnorm(1, 0, sqrt(disp.varX))
      G[i,t,2] <- G[i,t-1,2] + rnorm(1, 0, sqrt(disp.varY))
    } # t
  } # i
  # only keep coordinates of those that are alive
  G[,,1] <- G[,,1] * Z
  G[,,2] <- G[,,2] * Z
  G[G==0] <- NA
  
  # Simulating sampling 
  # Capture histories
  Y <- Z
  for (i in 1:n){
    for (t in (f[i]+1):nyears){
      Y[i,t] <- rbinom(1, 1, p) * Z[i,t]
    } # t
  } # i
  
  # Exclude locations that are outside the boundary (study area)
  inside <- Z
  for (i in 1:n){
    for (t in f[i]:nyears){
      inside[i,t] <- check.position(G[i,t,1], G[i,t,2], grid.x.study.area, grid.y.study.area, grid.size)
    } # t
  } # i
  
  # Capture-recapture data within study area
  Y.in <- Y * inside
  Y.in[is.na(Y.in)] <- 0
  G.in[,,1] <- G[,,1]*inside*Y
  G.in[,,2] <- G[,,2]*inside*Y
  G.in[G.in==0] <- NA
  
  # Capture-recapture data including all captures (not restrcted to study area)
  Y.comp <- Y
  G.comp[,,1] <- G[,,1]*Y
  G.comp[,,2] <- G[,,2]*Y
  G.comp[G.comp==0] <- NA
  
  # Output
  return(list(Y.in = Y.in, Y.comp = Y.comp, G.in = G.in, G.comp = G.comp, f = f))
}



# Function to check whether a position is within the study area
check.position <- function(X, Y, grid.x.study.area, grid.y.study.area, grid.size){
  ngrids <- length(grid.x.study.area)
  inside <- numeric()
  for (i in 1:ngrids){
    inside[i] <- (X > grid.x.study.area[i]) & (X < (grid.x.study.area[i] + grid.size)) & (Y > grid.y.study.area[i]) & (Y < (grid.y.study.area[i] + grid.size))   
  }
  ins <- sum(inside)
  return(ins)
}

############################################################################

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

##########################################


# Specify analyzing models

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
  
     for (t in (f[i]+1):n.occasions){
        # State processes
        # Survival
        z[i,t] ~ dbern(phi[i,t-1] * z[i,t-1])
        
        # Dispersal
        G[i,t,1] ~ dnorm(G[i,t-1,1], tau[1])
        G[i,t,2] ~ dnorm(G[i,t-1,2], tau[2])
  
        # Observation process
        # Test whether the actual location is in- or outside the state-space
        for (g in 1:ngrids){
           inside[i,t,g] <- step(G[i,t,1]-grid.x.study.area[g]) * step(grid.x.study.area[g]+grid.size-G[i,t,1]) *
                        step(G[i,t,2]-grid.y.study.area[g]) * step(grid.y.study.area[g]+grid.size-G[i,t,2])
           }
        r[i,t] <- sum(inside[i,t,1:ngrids])
        Y[i,t] ~ dbern(p[i,t-1] * z[i,t] * r[i,t])
        } #t
     } #i
  })


####################################################
# Run the model on a simulated data set

# 1. Create study area 
x.space <- 15        # Number of grids (in x direction)
y.space <- 15        # Number of grids (in y direction)

# Grid the study area
grid.size <- 2
gx <- seq(0, x.space, by = grid.size)
gy <- seq(0, y.space, by = grid.size)

grid.matrix <- matrix(0, ncol=length(gx), nrow=length(gy))

# Define the cells to belong to the study area
grid.matrix[2,4:6] <- 1
grid.matrix[3,3:6] <- 1
grid.matrix[4,3:7] <- 1
grid.matrix[5,1:7] <- 1
grid.matrix[6,4:7] <- 1
grid.matrix[7,5:6] <- 1

# List with grids and their min coordinates - the max is determined by grid.size
grid.x <- rep(gx, length(gx))
grid.y <- rep(gy, rep(length(gy), length(gy)))

grid.x.study.area <- grid.x[which(grid.matrix==1)]
grid.y.study.area <- grid.y[which(grid.matrix==1)]

simul.sCJSN.grid(rel = 100, nyears = 7, s = 0.7, p = 0.8, disp.varX = 0.5, disp.varY = 0.5, grid.x.study.area = grid.x.study.area, grid.y.study.area = grid.y.study.area)

# 2. Create data
data <- list(Y = dat$Y.in,  G = dat$G.in)
consts <- list(nind = nrow(dat$Y.in), f = dat$f, n.occasions = ncol(dat$Y.in), grid.x.study.area = grid.x.study.area, grid.y.study.area = grid.y.study.area, grid.size = grid.size, ngrids = length(grid.x.study.area))

# 4. Initial values
inits <- list(mean.phi = runif(1,0,1), mean.p = runif(1,0,1), sigma = runif(2,0.01,2), z = known.state.cjs(dat$Y.in))

# 5. Parameters monitored
parameters <- c("mean.phi", "mean.p", "sigma")

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
cMCMC <- compileNimble(rMCMC, project = model)

# Run the MCMC
system.time(run <- runMCMC(cMCMC,
                           parms
                           niter = 100000, 
                           nburnin = 24000, 
                           nchains = 2, 
                           progressBar = TRUE, 
                           summary = TRUE, 
                           samplesAsCodaMCMC = TRUE, 
                           thin = 1))
run$summary
plot(run$samples)

