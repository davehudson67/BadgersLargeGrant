nimblemodel({
  #----------------------------------------------
  # Parameters:
  # phi[k]: survival probability at site k (k from 1 to 17)
  # psi[k, l]: movement probability from site k to site l (k, l from 1 to 17)
  # p[k]: recapture probability at site k (k from 1 to 17)
  #----------------------------------------------
  # States (S):
  # 1 alive at site 1
  # 2 alive at site 2
  # ...
  # 17 alive at site 17
  # 18 dead
  # Observations (O):
  # 1 seen at site 1
  # 2 seen at site 2
  # ...
  # 17 seen at site 17
  # 18 not seen
  #----------------------------------------------
  
  # Priors and constraints
  for (k in 1:17) {
    phi[k] ~ dunif(0, 1)
    p[k] ~ dunif(0, 1)
  }
  
  # Transition probabilities
  for (k in 1:16) { # the last site transition will be calculated as 1 minus the sum of others
    lpsiA[k] ~ dnorm(0, 0.001)
    lpsiB[k] ~ dnorm(0, 0.001)
    lpsiC[k] ~ dnorm(0, 0.001)
    lpsiD[k] ~ dnorm(0, 0.001)
    lpsiE[k] ~ dnorm(0, 0.001)
    lpsiF[k] ~ dnorm(0, 0.001)
    lpsiG[k] ~ dnorm(0, 0.001)
    lpsiH[k] ~ dnorm(0, 0.001)
    lpsiI[k] ~ dnorm(0, 0.001)
    lpsiJ[k] ~ dnorm(0, 0.001)
    lpsiK[k] ~ dnorm(0, 0.001)
    lpsiL[k] ~ dnorm(0, 0.001)
    lpsiM[k] ~ dnorm(0, 0.001)
    lpsiN[k] ~ dnorm(0, 0.001)
    lpsiO[k] ~ dnorm(0, 0.001)
    lpsiP[k] ~ dnorm(0, 0.001)
    lpsiQ[k] ~ dnorm(0, 0.001)
  }
  
  for (i in 1:16){
    psiA[i] <- exp(lpsiA[i]) / (1 + exp(lpsiA[1]) + exp(lpsiA[2]))
    psiB[i] <- exp(lpsiB[i]) / (1 + exp(lpsiB[1]) + exp(lpsiB[2]))
    psiC[i] <- exp(lpsiC[i]) / (1 + exp(lpsiC[1]) + exp(lpsiC[2]))
    psiD[i] <- exp(lpsiD[i]) / (1 + exp(lpsiD[1]) + exp(lpsiD[2]))
    psiE[i] <- exp(lpsiE[i]) / (1 + exp(lpsiE[1]) + exp(lpsiE[2]))
    psiF[i] <- exp(lpsiF[i]) / (1 + exp(lpsiF[1]) + exp(lpsiF[2]))
    psiG[i] <- exp(lpsiG[i]) / (1 + exp(lpsiG[1]) + exp(lpsiG[2]))
    psiH[i] <- exp(lpsiH[i]) / (1 + exp(lpsiH[1]) + exp(lpsiH[2]))
    psiI[i] <- exp(lpsiI[i]) / (1 + exp(lpsiI[1]) + exp(lpsiI[2]))
    psiJ[i] <- exp(lpsiJ[i]) / (1 + exp(lpsiJ[1]) + exp(lpsiJ[2]))
    psiK[i] <- exp(lpsiK[i]) / (1 + exp(lpsiK[1]) + exp(lpsiK[2]))
    psiL[i] <- exp(lpsiL[i]) / (1 + exp(lpsiL[1]) + exp(lpsiL[2]))
    psiM[i] <- exp(lpsiM[i]) / (1 + exp(lpsiM[1]) + exp(lpsiM[2]))
    psiN[i] <- exp(lpsiN[i]) / (1 + exp(lpsiN[1]) + exp(lpsiN[2]))
    psiO[i] <- exp(lpsiO[i]) / (1 + exp(lpsiO[1]) + exp(lpsiO[2]))
    psiP[i] <- exp(lpsiP[i]) / (1 + exp(lpsiP[1]) + exp(lpsiP[2]))
    psiQ[i] <- exp(lpsiQ[i]) / (1 + exp(lpsiQ[1]) + exp(lpsiQ[2]))
  }
  
  psiA[17] <- 1 - psiA[1] - psiA[2] - psiA[3] - psiA[4] - psiA[5] - psiA[6] - psiA[7] - psiA[8] - psiA[9] - psiA[10] - psiA[11] - psiA[12] - psiA[13] - psiA[14] - psiA[15] - psiA[16]
  psiB[17] <- 1 - psiB[1] - psiB[2] - psiB[3] - psiB[4] - psiB[5] - psiB[6] - psiB[7] - psiB[8] - psiB[9] - psiB[10] - psiB[11] - psiB[12] - psiB[13] - psiB[14] - psiB[15] - psiB[16]
  psiC[17] <- 1 - psiC[1] - psiC[2] - psiC[3] - psiC[4] - psiC[5] - psiC[6] - psiC[7] - psiC[8] - psiC[9] - psiC[10] - psiC[11] - psiC[12] - psiC[13] - psiC[14] - psiC[15] - psiC[16]
  psiD[17] <- 1 - psiD[1] - psiD[2] - psiD[3] - psiD[4] - psiD[5] - psiD[6] - psiD[7] - psiD[8] - psiD[9] - psiD[10] - psiD[11] - psiD[12] - psiD[13] - psiD[14] - psiD[15] - psiD[16]
  psiE[17] <- 1 - psiE[1] - psiE[2] - psiE[3] - psiE[4] - psiE[5] - psiE[6] - psiE[7] - psiE[8] - psiE[9] - psiE[10] - psiE[11] - psiE[12] - psiE[13] - psiE[14] - psiE[15] - psiE[16]
  psiF[17] <- 1 - psiF[1] - psiF[2] - psiF[3] - psiF[4] - psiF[5] - psiF[6] - psiF[7] - psiF[8] - psiF[9] - psiF[10] - psiF[11] - psiF[12] - psiF[13] - psiF[14] - psiF[15] - psiF[16]
  psiG[17] <- 1 - psiG[1] - psiG[2] - psiG[3] - psiG[4] - psiG[5] - psiG[6] - psiG[7] - psiG[8] - psiG[9] - psiG[10] - psiG[11] - psiG[12] - psiG[13] - psiG[14] - psiG[15] - psiG[16]
  psiH[17] <- 1 - psiH[1] - psiH[2] - psiH[3] - psiH[4] - psiH[5] - psiH[6] - psiH[7] - psiH[8] - psiH[9] - psiH[10] - psiH[11] - psiH[12] - psiH[13] - psiH[14] - psiH[15] - psiH[16]
  psiI[17] <- 1 - psiI[1] - psiI[2] - psiI[3] - psiI[4] - psiI[5] - psiI[6] - psiI[7] - psiI[8] - psiI[9] - psiI[10] - psiI[11] - psiI[12] - psiI[13] - psiI[14] - psiI[15] - psiI[16]
  psiJ[17] <- 1 - psiJ[1] - psiJ[2] - psiJ[3] - psiJ[4] - psiJ[5] - psiJ[6] - psiJ[7] - psiJ[8] - psiJ[9] - psiJ[10] - psiJ[11] - psiJ[12] - psiJ[13] - psiJ[14] - psiJ[15] - psiJ[16]
  psiK[17] <- 1 - psiK[1] - psiK[2] - psiK[3] - psiK[4] - psiK[5] - psiK[6] - psiK[7] - psiK[8] - psiK[9] - psiK[10] - psiK[11] - psiK[12] - psiK[13] - psiK[14] - psiK[15] - psiK[16]
  psiL[17] <- 1 - psiL[1] - psiL[2] - psiL[3] - psiL[4] - psiL[5] - psiL[6] - psiL[7] - psiL[8] - psiL[9] - psiL[10] - psiL[11] - psiL[12] - psiL[13] - psiL[14] - psiL[15] - psiL[16]
  psiM[17] <- 1 - psiM[1] - psiM[2] - psiM[3] - psiM[4] - psiM[5] - psiM[6] - psiM[7] - psiM[8] - psiM[9] - psiM[10] - psiM[11] - psiM[12] - psiM[13] - psiM[14] - psiM[15] - psiM[16]
  psiN[17] <- 1 - psiN[1] - psiN[2] - psiN[3] - psiN[4] - psiN[5] - psiN[6] - psiN[7] - psiN[8] - psiN[9] - psiN[10] - psiN[11] - psiN[12] - psiN[13] - psiN[14] - psiN[15] - psiN[16]
  psiO[17] <- 1 - psiO[1] - psiO[2] - psiO[3] - psiO[4] - psiO[5] - psiO[6] - psiO[7] - psiO[8] - psiO[9] - psiO[10] - psiO[11] - psiO[12] - psiO[13] - psiO[14] - psiO[15] - psiO[16]
  psiP[17] <- 1 - psiP[1] - psiP[2] - psiP[3] - psiP[4] - psiP[5] - psiP[6] - psiP[7] - psiP[8] - psiP[9] - psiP[10] - psiP[11] - psiP[12] - psiP[13] - psiP[14] - psiP[15] - psiP[16]
  psiQ[17] <- 1 - psiQ[1] - psiQ[2] - psiQ[3] - psiQ[4] - psiQ[5] - psiQ[6] - psiQ[7] - psiQ[8] - psiQ[9] - psiQ[10] - psiQ[11] - psiQ[12] - psiQ[13] - psiQ[14] - psiQ[15] - psiQ[16]
  
  for (i in 1:nind) { 
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions - 1)) {
      ps[1, i, t, 1] <- phiA * psiA[1]
      ps[1, i, t, 2] <- phiA * psiA[2]
      ps[1, i, t, 3] <- phiA * psiA[3]
      ps[1, i, t, 4] <- phiA * psiA[4]
      ps[1, i, t, 5] <- phiA * psiA[5]
      ps[1, i, t, 6] <- phiA * psiA[6]
      ps[1, i, t, 7] <- phiA * psiA[7]
      ps[1, i, t, 8] <- phiA * psiA[8]
      ps[1, i, t, 9] <- phiA * psiA[9]
      ps[1, i, t, 10] <- phiA * psiA[10]
      ps[1, i, t, 11] <- phiA * psiA[11]
      ps[1, i, t, 12] <- phiA * psiA[12]
      ps[1, i, t, 13] <- phiA * psiA[13]
      ps[1, i, t, 14] <- phiA * psiA[14]
      ps[1, i, t, 15] <- phiA * psiA[15]
      ps[1, i, t, 16] <- phiA * psiA[16]
      ps[1, i, t, 17] <- phiA * psiA[17]
      ps[1, i, t, 18] <- 1 - phiA
      
      } # t
  } # i



    # Probability of staying at the same site
    psi[k, k] <- 1 - sum(psi[k, -k])
  }
  
  # Define state-transition and observation matrices
  for (i in 1:nind) {
    # Define probabilities of state S(t+1) given S(t)
    for (t in f[i]:(n.occasions - 1)) {
      for (k in 1:17) {
        for (l in 1:17) {
          ps[k, i, t, l] <- phi[k] * psi[k, l]
        }
        ps[k, i, t, 18] <- 1 - phi[k]  # transition to dead state
      }
      ps[18, i, t, 18] <- 1  # once dead, stay dead
      
      # Define probabilities of O(t) given S(t)
      for (k in 1:17) {
        po[k, i, t, k] <- p[k]
        for (l in 1:17) {
          if (k != l) {
            po[k, i, t, l] <- 0
          }
        }
        po[k, i, t, 18] <- 1 - p[k]
      }
      po[18, i, t, 18] <- 1  # if dead, not seen
      for (k in 1:17) {
        po[18, i, t, k] <- 0  # if dead, not seen at any site
      }
    } # t
  } # i
  
  # Likelihood
  for (i in 1:nind) {
    # Define latent state at first capture
    z[i, f[i]] <- Y[i, f[i]]
    for (t in (f[i] + 1):n.occasions) {
      # State process: draw S(t) given S(t-1)
      z[i, t] ~ dcat(ps[z[i, t-1], i, t-1,])
      # Observation process: draw O(t) given S(t)
      y[i, t] ~ dcat(po[z[i, t], i, t-1,])
    } # t
  } # i
})
