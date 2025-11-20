start.time <- Sys.time()

library(nimble)
library(MCMCvis)
library(tidyverse)
library(dplyr)
library(abind)

#' Read in the data. 
#------------------------------------------------------------------------#
cap_data <- read.csv("JoeCreek_capture_full_wuk_sex.csv")
dd_fl <- read.csv("JoeCreek_forklength_full.csv")

y <- cap_data %>%
  select(Fall_2014:Fall_2020) %>%
  as.matrix()
#head(y)

dd_fl <- dd_fl %>%
  select(Fall_2014: Fall_2020) %>%
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
  
  a_sigma <- 2.1 # For Inv Gamma variance 10
  b_sigma <- 1.1
  
  # For Gamma variance 10, the rate of growth is higher of males then females,
  # Exp value for males is taken higher than female.
  mu_m_o <- 1.6
  mu_f_o <- 0.65
  
  sigma_mum_o <- 2.5
  sigma_muf_o <- 4
  
  sigma ~ dinvgamma(a_sigma, b_sigma)
  
  for (t in 1:nyears){
    
    k_m[t] ~ dgamma(shape = mu_m_o, scale = sigma_mum_o)
    k_f[t] ~ dgamma(shape = mu_f_o, scale = sigma_muf_o)
    
  }
  
  Linf_f ~ T(dnorm(6.76, 0.1), 6.73, 6.8)
  Linf_m ~ T(dnorm(6.76, 0.1), 6.73, 6.8)
  
  for (i in 1:N) {
    
    male[i] ~ dbern(p.male)
    
    # asymptotic size - male or female
    a[i] <- Linf_f*(1-male[i]) + (male[i]*Linf_m) 
    
    for (t in (first[i]+1):nyears) {
      
      k[i,t] <- k_m[t]*male[i] + k_f[t]*(1-male[i])
      
      fl[i, t] ~ T(dnorm(fl[i, t-1] + k[i,t], sd = sigma), fl[i,t-1], a[i])
      
    }
    
    # Modelling
    
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

#head(d_fl_inits)

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

k_m_inits <- runif(nyears)
k_f_inits <- runif(nyears)

#head(mu_inits)


s_inits <- sex.st #ifelse(!is.na(firth_sex)==NA, 1)

s_inits[is.na(sex.st)] <- sample(c(0,1),sum(is.na(sex.st)), replace = TRUE)
s_inits[!is.na(sex.st)] <- NA 

#------------------------ Final initialization --------------------------#

initial.values <- list(alpha = runif(3, 0, 1), beta = runif(3,0,1),
                       Linf_f = 6.76, Linf_m = 6.76,
                       mu_m_o = 1.6, sigma_mum_o = 2.5,
                       mu_f_o = 0.65, sigma_muf_o = 4,
                       sigma = 1, a_sigma = 2.1, b_sigma = 1.1, 
                       male = s_inits,p.male = 0.2,
                       fl = fl_inits, z = z_inits, k_m = k_m_inits, 
                       k_f = k_f_inits)


#' Parameters to monitor. 
#------------------------------------------------------------------------#
parameters.to.save <- c("alpha", "beta", "sigma", "Linf_f",
                        "Linf_m","k_m","k_f","p.male","fl")

#' MCMC settings.
#------------------------------------------------------------------------#
n.iter <- 850000
n.burnin <- 100000
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

# Combine MCMC samples from all chains
combined_samples <- do.call(rbind, samples)

# Remove columns related to 'fl'
combined_samples_nofl <- combined_samples[, !grepl("^fl\\[", colnames(combined_samples))]

# Write CSV without fl
write.csv(combined_samples_nofl, file = "Joe_m3_samples_S.csv", row.names = FALSE)

# Filter samples for MCMCplot and MCMCtrace (exclude fl)
samples_nofl <- lapply(samples, function(chain) chain[, !grepl("^fl\\[", colnames(chain))])

# Plot posterior distributions (no fl)
pdf(file = "Joe_mp_m3_S.pdf")
MCMCplot(samples_nofl, HPD = TRUE)
dev.off()

# MCMC summary (no fl)
s <- MCMCsummary(samples_nofl, round = 5, params =  c("alpha", "beta", "sigma", "Linf_f",
                                                      "Linf_m","k_m","k_f","p.male"))
write.csv(s, file = "Joe_m3_sum_S.csv")

# Trace plot (no fl)
MCMCtrace(samples_nofl, pdf = TRUE, open_pdf = FALSE, filename = "Joe_m3_S", ind = TRUE,
          Rhat = FALSE, n.eff = FALSE)

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

write.csv(f, file = "Joe_m3_length_s.csv")

end.time <- Sys.time()

time.taken <- end.time - start.time
time.taken
