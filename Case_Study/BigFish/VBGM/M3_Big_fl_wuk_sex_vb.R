start.time <- Sys.time()

library(nimble)
library(MCMCvis)
library(tidyverse)
library(dplyr)
library(abind)

#' Read in the data. 
#------------------------------------------------------------------------#
cap_data <- read.csv("Bigfish_capture_full_wuk_sex.csv")
dd_fl <- read.csv("Bigfish_forklength_full.csv")

y <- cap_data %>%
  select(Fall_2009:Fall_2020) %>%
  #select(year_2014:year_2020) %>%
  as.matrix()
#head(y)

dd_fl <- dd_fl %>%
  select(Fall_2009: Fall_2020) %>%
  as.matrix()
dd_fl[dd_fl == 0] <- NA
#head(dd_fl)

#Scaling the forklength
#-----------------------------------------------------------------------#

# Define a function to scale each row separately
logscale_row <- function(row) {
  # Filter out NA values
  row_without_nas <- row[!is.na(row)]
  # Scale the row
  logscaled_row <- log(row_without_nas)
  # Replace NA positions with NA in scaled row
  logscaled_row_with_nas <- rep(NA, length(row))
  logscaled_row_with_nas[!is.na(row)] <- logscaled_row
  return(logscaled_row_with_nas)
}

# Apply the scaling function to each row of the forklength
d_fl <- t(apply(dd_fl, 1, logscale_row))

#d_fl <- dd_fl

#head(d_fl)

#' Get occasion of first capture.
#------------------------------------------------------------------------#
get.first <- function(x) min(which(x != 0))
first <- apply(y, 1, get.first)
#first

# Nimble Code
#------------------------------------------------------------------------#

code_m <- nimbleCode({
  # Priors
  
  alpha[1] ~ dnorm(0,0.1)
  alpha[2] ~ dnorm(0,0.1)
  alpha[3] ~ dnorm(0,0.1)
  
  beta[1] ~ dnorm(0,0.1)
  beta[2] ~ dnorm(0,0.1)
  beta[3] ~ dnorm(0,0.1)
  
  p.male ~ dbeta(1,1)
  
  a_sigma <- 2.1
  b_sigma <- 1.1
  
  sigma ~ dinvgamma(a_sigma, b_sigma)
  
  ## priors for VB growth parameters
  
  Linf_f ~ T(dnorm(6.67, 0.1),6.54, 6.8)
  Linf_m ~ T(dnorm(6.67, 0.1),6.54, 6.8)
  
  # Growth Coefficient
  for(t in 1:nyears){
    
    VBk_m[t] ~ dgamma(shape = 1.6, scale = 2.5)
    VBk_f[t] ~ dgamma(shape = 0.65, scale = 4)
    
  }
  
  for (i in 1:N) {
    
    # sex of unknown individuals
    male[i] ~ dbern(p.male)
    
    # asymptotic size - male or female
    a[i] <- Linf_f*(1-male[i]) + (male[i]*Linf_m) 
    
    # growth coefficient - male or female
    
    for (t in (first[i]+1):nyears) {
      
      # growth coefficient depending on reproduction
      
      k[i,t] <-  VBk_f[t]*(1-male[i]) + (male[i]*VBk_m[t])
    }
    
    for (t in (first[i]+1):nyears) {
      
      fl[i,t] ~ T(dnorm(increment[i,t],sigma), fl[i,t-1],a[i])
      increment[i,t] <- fl[i,t-1] + (a[i] - fl[i,t-1])*(1-exp(-k[i,t]))  # VB Model
    }
    
    # Multistate-Modelling
    
    # Initial state
    z[i,first[i]] <- y[i,first[i]] - 1
    
    for (t in (first[i]+1):nyears) {
      # State transition
      z[i, t] ~ dcat(px[z[i, t - 1], 1:2, i, t])
      
      # Observation process
      y[i, t] ~ dcat(po[z[i, t], 1:2, i, t])
    }
  }
  
  # probabilities of state z(t+1) given z(t)
  for (i in 1:N) {
    
    for (t in first[i]:nyears) {
      
      logit(phi[i, t]) <- alpha[1] + (alpha[2]*fl[i, t]*male[i]) +
        (alpha[3]*fl[i,t]*(1-male[i]))
      
      logit(p[i, t]) <- beta[1] + (beta[2]*fl[i,t]*male[i]) +
        (beta[3]*fl[i,t]*(1-male[i]))
      
      # Latent state process probabilities
      px[1, 1, i, t] <- phi[i,t]
      px[1, 2, i, t] <- 1 - phi[i,t]
      
      px[2, 1, i, t] <- 0
      px[2, 2, i, t] <- 1
      
      # Observation process probabilities
      po[1, 1, i, t] <- 1 - p[i,t]
      po[1, 2, i, t] <- p[i,t]
      
      po[2, 1, i, t] <- 1
      po[2, 2, i, t] <- 0
      
    }
  }
})

#Input Data
#-----------------------------------------------------------------------#
sex.st<-as.vector(cap_data$Sex)
my.data = list(y=y + 1, fl = d_fl, male = sex.st)

#Constants in a list. 
#-----------------------------------------------------------------------#
N = nrow(y)
nyears = ncol(y)

my.constants <- list(N = N, nyears = nyears, first = first)

# Generate inits for the latent states
#-----------------------------------------------------------------------#

# Initializing fl

fl_inits_temp <- dd_fl

replace_nas_increasing_order <- function(row) {
  na_indices <- which(is.na(row))
  
  for (i in na_indices) {
    if (i == 1) {
      row[i] <- row[i + 1] + 1
    } else if (i == length(row)) {
      row[i] <- row[i - 1] + 1
    } else {
      # Adjust the range to handle consecutive NAs
      start_range <- ifelse(!is.finite(row[i - 1]), 0, row[i - 1] + 1)
      end_range <- ifelse(!is.finite(row[i + 1]), start_range + 1, row[i + 1] - 1)
      
      if (!is.finite(start_range) || !is.finite(end_range)) {
        row[i] <- row[i - 1] + 1
      } else {
        row[i] <- sample(seq(start_range, end_range), 1)
      }
    }
  }
  
  return(row)
}

fl_inits_1 <- apply(fl_inits_temp, 1, replace_nas_increasing_order)

fl_inits_2 <- t(fl_inits_1)

#head(fl_inits_1)
#head(fl_inits_2)

# Function to replace values in Matrix B with NAs based on NAs in Matrix A
replace_values_with_nas <- function(matrix_a, matrix_b) {
  # Get the indices where Matrix A has NAs
  na_indices <- which(!is.na(matrix_a), arr.ind = TRUE)
  
  # Replace values in Matrix B with NAs at the corresponding positions
  matrix_b[na_indices] <- NA
  
  return(matrix_b)
}

d_fl_inits <- replace_values_with_nas(fl_inits_temp, fl_inits_2)

d_fl_inits <- as.matrix(d_fl_inits)

for (i in 1:N){
  if (first[i] == 1) next
  if (first[i] > 1) d_fl_inits[i,1:(first[i]-1)] <- NA
}

#head(d_fl_inits)

# Apply the scaling function to each row of the initial forklength
fl_inits <- t(apply(d_fl_inits, 1, logscale_row))
#fl_inits <- d_fl_inits
#head(fl_inits)
#tail(fl_inits)
#------------------------------------------------------------------------#
# Initializing alive state

x.init <- matrix(1, N, nyears)  # Set to 1 or another appropriate value for initial state
for (i in 1:N) {
  if (first[i] > 1) x.init[i, 1:(first[i]-1)] <- NA
}

z_inits <- x.init
#head(z_inits)
#------------------------------------------------------------------------#

increment_inits <- matrix(runif(N * nyears), nrow = N, ncol = nyears)

VBk_m_inits <- runif(nyears)
VBk_f_inits <- runif(nyears)

#------------------------------------------------------------------------#

#H_m_inits <- matrix(runif(N * nyears, 5.5, 6.75), nrow = N, ncol = nyears)


#------------------------------------------------------------------------#

s_inits <- sex.st #ifelse(!is.na(firth_sex)==NA, 1)

s_inits[is.na(sex.st)] <- sample(c(0,1),sum(is.na(sex.st)), replace = TRUE)
s_inits[!is.na(sex.st)] <- NA 

#------------------------ Final initialization --------------------------#

initial.values <- list(alpha = runif(3, 0, 1), beta = runif(3,0,1),
                       Linf_f = 6.67, Linf_m = 6.67, 
                       VBk_m = VBk_m_inits, VBk_f = VBk_f_inits,
                       increment = increment_inits,sigma = runif(1, 0, 10), 
                       a_sigma = 2.1, b_sigma = 1.1, male = s_inits, p.male=0.3,
                       fl = fl_inits, z = z_inits)


#' Parameters to monitor. 
#------------------------------------------------------------------------#
parameters.to.save <- c("alpha", "beta", "Linf_m", "Linf_f", 
                        "VBk_f", "VBk_m", "sigma","p.male", "fl")

#' MCMC settings.
#------------------------------------------------------------------------#
n.iter <- 600000
n.burnin <- 75000
n.chains <- 3

#--- CHECKED NIMBLE MODEL ---#

mcmc.multistate <- nimbleMCMC(code = code_m,
                              constants = my.constants,
                              data = my.data,
                              inits = initial.values,
                              monitors = parameters.to.save,
                              niter = n.iter,
                              nburnin = n.burnin,
                              nchains = n.chains,
                              summary = TRUE, WAIC = TRUE)


waic <- mcmc.multistate$WAIC
waic
#Model Summary
samples<- mcmc.multistate$samples

combined_samples <- do.call(rbind, samples)

write.csv(combined_samples, file = "Big_m3_samples_vb.csv", row.names = FALSE)

pdf(file = "Big_mp_m3_vb.pdf")
MCMCplot(samples, HPD = T)
dev.off()

s <- MCMCsummary(samples, round = 5, params = c("alpha", "beta", "Linf_m", "Linf_f", 
                                                "VBk_f", "VBk_m", "sigma","p.male"))

MCMCtrace(samples,pdf = T,open_pdf = F,filename = "Big_m3_vb", ind = TRUE,
          Rhat = FALSE, n.eff = FALSE)

write.csv(s, file = "Big_m3_sum_vb.csv")

# Check column names related to fl
fl_columns <- grep("^fl\\[", colnames(samples[[1]]), value = TRUE)

# If you want to combine across chains
fl_all_chains <- abind::abind(samples[[1]][, fl_columns],
                              samples[[2]][, fl_columns],
                              samples[[3]][, fl_columns],
                              along = 1)

# Posterior mean of fl for each [i, t]
fl_means <- apply(fl_all_chains, 2, mean)

# Convert to matrix form: assume you know N and nyears
fl_matrix <- matrix(fl_means, nrow = N, ncol = nyears, byrow = FALSE)

f <- exp(fl_matrix)

write.csv(f, file = "Big_m3_length_VB.csv")

end.time <- Sys.time()

time.taken <- end.time - start.time
time.taken
