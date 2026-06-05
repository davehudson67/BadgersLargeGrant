# Load your actual data
# Assuming 'bc' is your data frame, and 'X' is your trap location matrix
# 'bc' should have columns like 'tattoo' (individual ID), 'primary', 'trap_season', 'x', 'y' for individual locations

S_init <- array(NA, dim = c(nind, 2, n.prim))

# Populate S with the actual x, y coordinates from your data
for (i in 1:nind) {
  # Subset data for this individual
  individual_data <- bc[bc$tattoo == unique(bc$tattoo)[i], ]
  
  for (k in 1:n.prim) {
    # Check if there is an entry for this primary session
    matching_rows <- individual_data[individual_data$primary == k, ]
    
    if (nrow(matching_rows) == 1) {
      # Single entry: Assign coordinates
      S_init[i, 1, k] <- matching_rows$col_index
      S_init[i, 2, k] <- matching_rows$row_index
      
    } else if (nrow(matching_rows) > 1) {
      # Multiple entries: Take the first entry (or modify logic if needed)
      S_init[i, 1, k] <- matching_rows$col_index[1]
      S_init[i, 2, k] <- matching_rows$row_index[1]
      
    } else {
      # No entry for this session, carry forward the last known position
      if (k > 1) {
        S_init[i, 1, k] <- S_init[i, 1, k-1]
        S_init[i, 2, k] <- S_init[i, 2, k-1]
      } else {
        S_init[i, 1, k] <- NA
        S_init[i, 2, k] <- NA
      }
    }
  }
}

S_init <- S_init[ord,,]

# Initialize the distance matrix D: [nind, R, n.prim]
D_init <- array(1, dim = c(nind, R, n.prim))

# Calculate Euclidean distance for each individual, trap, and session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    for (r in 1:R) {
      D_init[i, r, k] <- sqrt((S_init[i, 1, k] - X[r, 1])^2 + (S_init[i, 2, k] - X[r, 2])^2)
    }
  }
}

# Initialize g: [nind, n.prim, R+1]
g_init <- array(1, dim = c(nind, n.prim, R + 1))

kappa = 2.0
sigma = 5.0

# Calculate trap exposure for each individual, trap, and session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    for (r in 1:R) {
      g_init[i, k, r + 1] <- exp(-((D_init[i, r, k] / sigma)^kappa))
    }
  }
}

# Initialize G: [nind, n.prim]
G_init <- matrix(1, nrow = nind, ncol = n.prim)

# Sum trap exposures for each individual and primary session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    G_init[i, k] <- sum(g_init[i, k, 2:(R + 1)])  # Summing from r = 2 to R+1 to skip the non-capture state
  }
}

P_init <- array(1, dim = c(nind, max(J), n.prim))

PL = 0.7
log.lambda0 <- log(-log(1 - PL))
lambda <- exp(log.lambda0)

for (i in 1:nind) {
  for (k in 1:n.prim) {
    for (j in 1:J[i, k]) {
      P_init[i, j, k] <- 1 - exp(-lambda * G_init[i, k])
    }
  }
}

# Extract relevant data for individual i and primary session k
inits <- list(
  PL = 0.7,
  #log.lambda0 <- log(-log(1 - PL)),
  #lambda <- exp(log.lambda0),
  kappa = 1.0,
  sigma = 10.0,
  phi = runif(n.prim - 1, 0.5, 0.8),
  dmean = 10 + runif(1,-3,3),
  S = S_init,
  D = D_init,
  g = g_init,
  G = G_init,
  P = P_init)

