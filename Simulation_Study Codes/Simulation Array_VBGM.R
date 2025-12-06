library(nimble)
library(MCMCvis)
library(tidyverse)
library(dplyr)

# Define the model code

sim_function <- function(params){
  
  iter = params$iter
  
  code_sim_simple <- nimbleCode({
    
    alpha[1] <- params$alpha0
    alpha[2] <- params$alpha1
    
    beta[1] <- params$beta0
    beta[2] <- params$beta1
    
    a_sigma <- 2.1
    b_sigma <- 1.1
    
    mu_o <- 1.1
    sigma_mu_o <- 3
    
    sigma ~ dinvgamma(a_sigma, b_sigma)
    
    a ~  T(dnorm(6.76, 0.1), 6.73, 6.8)
    
    for (i in 1:N) {
      
      fl[i, first[i]] ~ dunif(log(300), log(860))
    }
    
    for (t in 1:nyears){
      
      k[t] ~ dgamma(shape = mu_o, scale = sigma_mu_o)
    }
    
    for (i in 1:N) {
      
      for (t in (first[i]+1):nyears) {
        
        fl[i, t] ~ T(dnorm(fl[i, t-1] + (a - fl[i,t-1])*(1-exp(-k[t])), sd = sigma),
                     fl[i,t-1], a)
      }
      
      y[i, first[i]] <- 2
      z[i, first[i]] <- y[i, first[i]] - 1
      
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
        
        logit(phi[i, t]) <- alpha[1] + (alpha[2]*fl[i, t])
        logit(p[i, t]) <- beta[1] + (beta[2]*fl[i,t])
        
        # Latent state process probabilities
        px[1, 1, i, t] <- phi[i,t]
        px[1, 2, i, t] <- 1 - phi[i,t]
        
        px[2, 1, i, t] <- 0
        px[2, 2, i, t] <- 1
        
        # Observation process probabilities
        po[1, 1, i, t] <- 1 - p[i,t]   # alive but not cap
        po[1, 2, i, t] <- p[i,t]       # alive and cap
        
        po[2, 1, i, t] <- 1            # dead and not cap
        po[2, 2, i, t] <- 0            # dead but cap
        
      }
    }
  })
  
  N = 5000        # Number of individuals
  nyears = 5  # Number of sampling years
  
  first <- rep(NA, N)
  for (j in 1:N) {
    first[j] <- sample(1:nyears,1)
  }
  # Initialize fl values outside the model to ensure different starting values
  initial_fl <- log(runif(N, 300, 860))
  inits <- list(alpha = runif(2, 0, 1), beta = runif(2, 0, 1),
                mu_o = 1.1, sigma_mu_o = 3, sigma = 1,
                a_sigma = 2.1, b_sigma = 1.1, a = 6.76,
                fl = matrix(0, nrow = N, ncol = nyears))
  
  
  for (i in 1:N) {
    inits$fl[i, first[i]] <- initial_fl[i]
  }
  
  sim_model <- nimbleModel(code_sim_simple, 
                           constants = list(N = N, nyears = nyears, first = first),
                           inits = inits, calculate = TRUE)
  
  # Get nodes for simulation
  sim_nodes <- sim_model$getDependencies(c("alpha", "beta", "mu_o",
                                           "sigma_mu_o", "sigma",
                                           "a_sigma", "b_sigma", "a"), self = FALSE,
                                         downstream = TRUE)
  
  # Simulate the data
  sim_model$simulate(sim_nodes)
  
  cap_sim2 <- sim_model$y      # Saving the simulated capture matrix
  fl_sim3 <- sim_model$fl
  
  cap_sim2[cap_sim2==1] <- 0
  cap_sim2[cap_sim2==2] <- 1
  
  cap_sim2[is.na(cap_sim2)] <- 0  # NA's to zero's (not captured)
  fl_sim3[is.na(fl_sim3)] <- 0
  
  fl_sim2 <- fl_sim3*cap_sim2
  
  cap_sim2 = data.frame(cap_sim2)
  
  fl_sim2 = data.frame(fl_sim2)
  
  # Removing all the rows containing only zeros
  cap_sim1<- cap_sim2[!(apply(cap_sim2, 1, function(y) all(y == 0))),]
  fl_sim1<- fl_sim2[!(apply(fl_sim2, 1, function(y) all(y == 0))),]
  
  m_sample <- 1000
  
  # Randomly select 3 unique row indices from the capture matrix
  selected_indices <- sort(sample(nrow(cap_sim1), m_sample))
  
  # Subset the capture matrix with the selected indices
  cap_sim <- cap_sim1[selected_indices, ]
  
  # Subset the forklength matrix with the same indices
  fl_sim <- fl_sim1[selected_indices, ]
  
  for (t in 1:nyears){
    names(cap_sim)[t]<- t 
    names(fl_sim)[t]<- t
  }
  
  capture_mat <- cap_sim  # Final capture matrix
  fl_mat <- round(fl_sim, 3)
  n_ind = nrow(capture_mat)
  
  ####--------####
  
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
  
  # Function to replace values in Matrix B with NAs based on NAs in Matrix A
  replace_values_with_nas <- function(matrix_a, matrix_b) {
    # Get the indices where Matrix A has NAs
    na_indices <- which(!is.na(matrix_a), arr.ind = TRUE)
    
    # Replace values in Matrix B with NAs at the corresponding positions
    matrix_b[na_indices] <- NA
    
    return(matrix_b)
  }
  
  get.first <- function(x) min(which(x != 0))
  
  # -- fl data preparation --- #
  
  # Exponentiate and round the values in fl_mat[[s]]
  a_fl <- round(exp(fl_mat))
  
  # Replace any occurrences of 1 with 0
  a_fl[a_fl == 1] <- 0
  
  ddd_fl <- a_fl
  ddd_fl[ddd_fl == 0] <- NA
  
  # Apply the scaling function to each row of the forklength
  dd_fl <- t(apply(ddd_fl, 1, logscale_row))
  
  
  code_run_simple <- nimbleCode({
    # Priors
    
    alpha[1] ~ dnorm(0,0.1)
    alpha[2] ~ dnorm(0,0.1)
    
    beta[1] ~ dnorm(0,0.1)
    beta[2] ~ dnorm(0,0.1)
    
    a_sigma <- 2.1
    b_sigma <- 1.1 
    
    mu_o <- 1.1 
    sigma_mu_o <- 3
    
    a ~  T(dnorm(6.76, 0.1), 6.73, 6.8)
    
    sigma ~ dinvgamma(a_sigma, b_sigma)
    
    for (t in 1:nyears){
      
      k[t] ~ dgamma(shape = mu_o, scale = sigma_mu_o)
    }
    
    for (i in 1:N) {
      
      for (t in (first_r[i]+1):nyears) {
        
        fl[i, t] ~ T(dnorm(fl[i, t-1] + (a - fl[i,t-1])*(1-exp(-k[t])), sd = sigma),
                     fl[i,t-1], a)
      }
      
      # Modelling
      
      # Initial state
      z[i,first_r[i]] <- y[i,first_r[i]] - 1
      
      for (t in (first_r[i]+1):nyears) {
        # State transition
        z[i, t] ~ dcat(px[z[i, t - 1], 1:2, i, t])
        
        # Observation process
        y[i, t] ~ dcat(po[z[i, t], 1:2, i, t])
      }
    }
    
    # probabilities of state z(t+1) given z(t)
    for (i in 1:N) {
      
      for (t in first_r[i]:nyears) {
        
        logit(phi[i, t]) <- alpha[1] + (alpha[2]*fl[i, t])
        logit(p[i, t]) <- beta[1] + (beta[2]*fl[i,t])
        
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
  
  y <- capture_mat
  fl <- dd_fl
  
  #Input Data
  #-----------------------------------------------------------------------#
  my.data = list(y=y + 1, fl = fl)
  
  #Constants in a list. 
  #-----------------------------------------------------------------------#
  N = nrow(y)
  nyears = ncol(y)
  
  first_r <- apply(y, 1, get.first)
  
  my.constants <- list(N = N, nyears = nyears, first_r = first_r)
  
  fl_inits_temp <- ddd_fl
  fl_inits_1 <- apply(fl_inits_temp, 1, replace_nas_increasing_order)
  
  fl_inits_2 <- t(fl_inits_1)
  
  d_fl_inits <- replace_values_with_nas(fl_inits_temp, fl_inits_2)
  
  d_fl_inits <- as.matrix(d_fl_inits)
  
  for (i in 1:N){
    if (first_r[i] == 1) next
    if (first_r[i] > 1) d_fl_inits[i,1:(first_r[i]-1)] <- NA
  }
  
  # Apply the scaling function to each row of the initial forklength
  fl_inits <- t(apply(d_fl_inits, 1, logscale_row))
  #------------------------------------------------------------------------#
  # Initializing alive state
  
  x.init <- matrix(1, N, nyears)  # Set to 1 or another appropriate value for initial state
  for (i in 1:N) {
    if (first_r[i] > 1) x.init[i, 1:(first_r[i]-1)] <- NA
  }
  
  z_inits <- x.init
  #head(z_inits)
  #------------------------------------------------------------------------#
  
  k_inits <- runif(nyears)
  #head(mu_inits)
  
  #------------------------ Final initialization --------------------------#
  
  initial.values <- list(alpha = runif(2, 0, 1), beta = runif(2,0,1),
                         mu_o = 1.1, sigma_mu_o = 3,a = 6.76,
                         sigma = 1, a_sigma = 2.1, b_sigma = 1.1,
                         fl = fl_inits, z = z_inits, k = k_inits)
  
  ## ---------------------------------------------------------------------------------------------------------
  parameters.to.save <- c("alpha","beta")
  
  #' MCMC settings.
  ## ---------------------------------------------------------------------------------------------------------
  n.iter <- 300000
  n.burnin <- 70000
  n.chains <- 2
  
  mcmc.multistate <- nimbleMCMC(code = code_run_simple, 
                                constants = my.constants,
                                data = my.data,              
                                inits = initial.values,
                                monitors = parameters.to.save,
                                niter = n.iter,
                                nburnin = n.burnin, 
                                nchains = n.chains,
                                summary = TRUE)
  
  samples_sim <- mcmc.multistate$samples
  
  s_alpha0 <- (MCMCsummary(samples_sim, params = "alpha", round=5)$mean)[1]
  s_alpha1 <- (MCMCsummary(samples_sim, params = "alpha", round=5)$mean)[2]
  s_beta0 <- (MCMCsummary(samples_sim, params = "beta", round=5)$mean)[1]
  s_beta1 <- (MCMCsummary(samples_sim, params = "beta", round=5)$mean)[2]
  
  d_b = data.frame(s_alpha0, s_alpha1, s_beta0, s_beta1)
  
  write.csv(d_b, paste0("alpha0", s_alpha0, "_alpha1",s_alpha1,
                        "_beta0",s_beta0,"_beta1",s_beta1,
                        "_Sim_",iter,".csv"), row.names = TRUE)
  
  return(d_b)
  
}

args <- commandArgs(trailingOnly = TRUE)
iter <- as.numeric(args[1])
params= list(alpha0=6.225, alpha1= -1,
             beta0=17, beta1= -2.5,
             iter=iter)
sim_function(params)
