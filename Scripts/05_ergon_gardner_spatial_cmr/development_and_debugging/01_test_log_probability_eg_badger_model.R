# Load your actual data
# Assuming 'bc' is your data frame, and 'X' is your trap location matrix
# 'bc' should have columns like 'tattoo' (individual ID), 'primary', 'trap_season', 'x', 'y' for individual locations

# Example setup for an individual
i <- 1340  # Specific individual
k <- 21    # Specific primary session
j <- 3     # Specific secondary session
r <- 12    # Specific trap

# Extract relevant data for individual i and primary session k
PL <- 0.7  # Replace with actual data or calculated values
log.lambda0 <- log(-log(1 - PL)) 
lambda <- exp(log.lambda0)

kappa <- 2.0  # Replace with actual parameter value if known
sigma <- 5.0  # Replace with actual parameter value if known

# Define number of individuals, traps, and primary sessions
nind <- length(unique(bc$tattoo))  # Number of individuals from your data
n.prim <- max(bc$primary)  # Number of primary sessions
R <- nrow(X)  # Number of traps (assuming X has trap coordinates)

# Initialize S based on actual data
S <- array(NA, dim = c(nind, 2, n.prim))

# Populate S with the actual x, y coordinates from your data
for (i in 1:nind) {
  # Subset data for this individual
  individual_data <- bc[bc$tattoo == unique(bc$tattoo)[i], ]
  
  for (k in 1:n.prim) {
    # Check if there is an entry for this primary session
    matching_rows <- individual_data[individual_data$primary == k, ]
    
    if (nrow(matching_rows) == 1) {
      # Single entry: Assign coordinates
      S[i, 1, k] <- matching_rows$col_index
      S[i, 2, k] <- matching_rows$row_index
      
    } else if (nrow(matching_rows) > 1) {
      # Multiple entries: Take the first entry (or modify logic if needed)
      S[i, 1, k] <- matching_rows$col_index[1]
      S[i, 2, k] <- matching_rows$row_index[1]
      
    } else {
      # No entry for this session, carry forward the last known position
      if (k > 1) {
        S[i, 1, k] <- S[i, 1, k-1]
        S[i, 2, k] <- S[i, 2, k-1]
      } else {
        S[i, 1, k] <- NA
        S[i, 2, k] <- NA
      }
    }
    
    # Debugging output
    #cat("i =", i, "k =", k, "S[i, 1, k] =", S[i, 1, k], "S[i, 2, k] =", S[i, 2, k], "\n")
  }
}

# Initialize the distance matrix D: [nind, R, n.prim]
D <- array(NA, dim = c(nind, R, n.prim))

# Calculate Euclidean distance for each individual, trap, and session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    for (r in 1:R) {
      D[i, r, k] <- sqrt((S[i, 1, k] - X[r, 1])^2 + (S[i, 2, k] - X[r, 2])^2)
    }
  }
}

# Initialize g: [nind, n.prim, R+1]
g <- array(NA, dim = c(nind, n.prim, R + 1))

# Calculate trap exposure for each individual, trap, and session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    for (r in 1:R) {
      g[i, k, r + 1] <- exp(-((D[i, r, k] / sigma)^kappa))
    }
  }
}

# Example: Print the results for the specific individual, primary session, and trap
cat("Distance D[i, r, k] =", D[i, r, k], "\n")
cat("Trap exposure g[i, k, r + 1] =", g[i, k, r + 1], "\n")

# Initialize G: [nind, n.prim]
G <- matrix(NA, nrow = nind, ncol = n.prim)

# Sum trap exposures for each individual and primary session
for (i in 1:nind) {
  for (k in 1:n.prim) {
    G[i, k] <- sum(g[i, k, 2:(R + 1)])  # Summing from r = 2 to R+1 to skip the non-capture state
  }
}

G_i_k <- G[i, k]  # Total trap exposure for individual i in primary session k
z_i_k <- 1#z[i, k]  # Survival status of individual i in primary session k
H_i_j_k <- H[i, j, k]  # Capture history entry for individual i, j, k
g_i_k_H <- g[i, k, H_i_j_k]  # Trap exposure at specific trap for individual i, k, H
Ones_i_j_k <- 1 #Ones[i, j, k]  # Observed capture event

# Manually calculate the capture probability
P_i_j_k <- (1 - exp(-lambda * G_i_k)) * z_i_k
captureProb_i_k_j <- step(H_i_j_k - 2) * (g_i_k_H / (G_i_k + 0.000000001)) * P_i_j_k + (1 - step(H_i_j_k - 2)) * (1 - P_i_j_k)

# Print out the values for inspection
cat("Debugging info for individual", i, "in primary session", k, "and secondary session", j, "\n")
cat("  P[i, j, k] = ", P_i_j_k, "\n")
cat("  g[i, k, H[i, j, k]] = ", g_i_k_H, "\n")
cat("  G[i, k] = ", G_i_k, "\n")
cat("  captureProb[i, k, j] = ", captureProb_i_k_j, "\n")
cat("  Ones[i, j, k] = ", Ones_i_j_k, "\n")

# Check if captureProb[i, k, j] is very small or zero
if (captureProb_i_k_j == 0 && Ones_i_j_k == 1) {
  cat("Warning: captureProb is zero while Ones is 1, leading to -Inf log probability.\n")
}


# Ensure all variables are defined correctly
H_i_j_k <- H[i, j, k]  # Capture history entry
g_i_k_H <- g[i, k, H_i_j_k]  # Trap exposure at specific trap
G_i_k <- G[i, k]  # Total trap exposure for individual i in session k
P_i_j_k <- P[i, j, k]  # Probability of being captured

# Check if these are atomic vectors or not
cat("H_i_j_k =", H_i_j_k, "g_i_k_H =", g_i_k_H, "G_i_k =", G_i_k, "P_i_j_k =", P_i_j_k, "\n")

# Ensure H_i_j_k is not an atomic vector and is correctly evaluated
if (length(H_i_j_k) != 1) {
  stop("H_i_j_k is not a scalar; it should be a single value.")
}

# Calculate capture probability
captureProb_i_k_j <- ifelse(H_i_j_k >= 2, 
                            (g_i_k_H / (G_i_k + 1e-9)) * P_i_j_k, 
                            1 - P_i_j_k)

# Print the result to check
cat("captureProb_i_k_j =", captureProb_i_k_j, "\n")


step(1 - 40)
