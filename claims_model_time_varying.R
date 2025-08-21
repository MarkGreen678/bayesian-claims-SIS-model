# ----------------------------------------
# Bayesian Inference: Time-Varying SIS Model
# ----------------------------------------

## Step 1: Simulate Time-Varying Claim Data

simulate_claim_SIS_tv <- function(N, beta_vec, gamma_vec, IO, T) {
  S <- numeric(T)
  I <- numeric(T)
  new_claims <- numeric(T)
  resolved_claims <- numeric(T)
  
  S[1] <- N - IO
  I[1] <- IO
  
  for (t in 1:(T - 1)) {
    p_claim <- 1 - exp(-beta_vec[t] * I[t] / N)
    p_resolve <- 1 - exp(-gamma_vec[t])
    
    new_I <- rbinom(1, S[t], p_claim)
    new_S <- rbinom(1, I[t], p_resolve)
    
    new_claims[t + 1] <- new_I
    resolved_claims[t + 1] <- new_S
    
    S[t + 1] <- S[t] - new_I + new_S
    I[t + 1] <- I[t] + new_I - new_S
  }
  
  return(data.frame(
    time = 0:(T - 1), S = S, I = I,
    new_claims = new_claims,
    resolved_claims = resolved_claims
  ))
}

## Step 2: Log-likelihood for Time-Varying Params

loglik_claim_SIS_tv <- function(beta_vec, gamma_vec, data) {
  N <- data$S[1] + data$I[1]
  loglik <- 0
  
  for (t in 1:(nrow(data) - 1)) {
    St <- data$S[t]
    It <- data$I[t]
    
    new_I <- data$new_claims[t + 1]
    new_S <- data$resolved_claims[t + 1]
    
    p_claim <- 1 - exp(-beta_vec[t] * It / N)
    p_resolve <- 1 - exp(-gamma_vec[t])
    
    # Bound probabilities to avoid log(0)
    p_claim <- min(max(p_claim, 1e-10), 1 - 1e-10)
    p_resolve <- min(max(p_resolve, 1e-10), 1 - 1e-10)
    
    loglik <- loglik +
      dbinom(new_I, St, p_claim, log = TRUE) +
      dbinom(new_S, It, p_resolve, log = TRUE)
  }
  
  return(loglik)
}

## Step 3: Random Walk Prior

log_prior_rw <- function(beta_vec, gamma_vec, sigma_beta, sigma_gamma, drift) {
  beta_diffs <- diff(beta_vec) - drift  # subtract expected increase
  gamma_diffs <- diff(gamma_vec)
  
  logp_beta <- sum(dnorm(beta_diffs, mean = 0, sd = sigma_beta, log = TRUE))
  logp_gamma <- sum(dnorm(gamma_diffs, mean = 0, sd = sigma_gamma, log = TRUE))
  
  return(logp_beta + logp_gamma)
}

## Step 4: MH Sampler for Time-Varying Params 

MH_sampler_tv <- function(data, n_iter, beta_init, gamma_init, lambda_beta, lambda_gamma, proposal_sd, drift) {
  T <- nrow(data)
  beta <- rep(beta_init, T - 1)
  gamma <- rep(gamma_init, T - 1)
  samples_beta <- matrix(NA, n_iter, T - 1)
  samples_gamma <- matrix(NA, n_iter, T - 1)
  accepted <- 0
  total <- 0
  
  for (i in 1:n_iter) {
    # Update each beta_t individually
    for (t in 1:(T - 1)) {
      beta_prop <- beta
      beta_prop[t] <- beta[t] + rnorm(1, 0, proposal_sd)
      
      if (beta_prop[t] > 0) {
        loglik_curr <- loglik_claim_SIS_tv(beta, gamma, data)
        loglik_prop <- loglik_claim_SIS_tv(beta_prop, gamma, data)
        
        logprior_curr <- log_prior_rw(beta, gamma, lambda_beta, lambda_gamma, drift)
        logprior_prop <- log_prior_rw(beta_prop, gamma, lambda_beta, lambda_gamma, drift)
        
        log_alpha <- (loglik_prop + logprior_prop) - (loglik_curr + logprior_curr)
        
        if (log(runif(1)) < log_alpha) {
          beta <- beta_prop
          accepted <- accepted + 1
        }
        total <- total + 1
      }
    }
    
    # Update each gamma_t individually
    for (t in 1:(T - 1)) {
      gamma_prop <- gamma
      gamma_prop[t] <- gamma[t] + rnorm(1, 0, proposal_sd)
      
      if (gamma_prop[t] > 0) {
        loglik_curr <- loglik_claim_SIS_tv(beta, gamma, data)
        loglik_prop <- loglik_claim_SIS_tv(beta, gamma_prop, data)
        
        logprior_curr <- log_prior_rw(beta, gamma, lambda_beta, lambda_gamma, drift)
        logprior_prop <- log_prior_rw(beta, gamma_prop, lambda_beta, lambda_gamma, drift)
        
        log_alpha <- (loglik_prop + logprior_prop) - (loglik_curr + logprior_curr)
        
        if (log(runif(1)) < log_alpha) {
          gamma <- gamma_prop
          accepted <- accepted + 1
        }
        total <- total + 1
      }
    }
    
    samples_beta[i, ] <- beta
    samples_gamma[i, ] <- gamma
  }
  
  cat("Acceptance rate:", round(accepted / total, 3), "\n")
  return(list(beta = samples_beta, gamma = samples_gamma))
}

## Step 5: Run Simulation and Inference 

set.seed(42)
T <- 30 # Set time variable

# Initiate true beta, gamma vectors
true_beta_vec <- seq(0.1, 0.8, length.out = T - 1)
true_gamma_vec <- rep(0.1, T - 1)

# Calculate true drift for beta vector
beta_drift <- (max(true_beta_vec) - min(true_beta_vec)) / (T - 2)

# Simulate time-varying data
data <- simulate_claim_SIS_tv(N = 500, beta_vec = true_beta_vec, gamma_vec = true_gamma_vec, IO = 20, T = T)

# Run time-varying MH sampler
samples <- MH_sampler_tv(
  data = data,
  n_iter = 1000,
  beta_init = 0.1,
  gamma_init = 0.2,
  lambda_beta = 0.01,
  lambda_gamma = 0.1,
  proposal_sd = 0.03,
  drift = beta_drift
)

## Step 7: Posterior Summary

# Extract posterior matrices
beta_mat <- samples$beta  # 5000 x 99
gamma_mat <- samples$gamma  # 5000 x 99

# Compute summary statistics (mean and 95% CI) at each time step
beta_mean <- apply(beta_mat, 2, mean)
beta_lower <- apply(beta_mat, 2, quantile, probs = 0.025)
beta_upper <- apply(beta_mat, 2, quantile, probs = 0.975)

gamma_mean <- apply(gamma_mat, 2, mean)
gamma_lower <- apply(gamma_mat, 2, quantile, probs = 0.025)
gamma_upper <- apply(gamma_mat, 2, quantile, probs = 0.975)

# Plot of Posterior Mean and 95% CI for beta(t)
plot(beta_mean, type = "l", ylim = range(c(beta_lower, beta_upper)), 
     main = expression(paste("Posterior Mean and 95% CI of ", beta(t))), 
     ylab = expression(beta(t)), xlab = "Time")
lines(beta_lower, col = "gray", lty = 2)
lines(beta_upper, col = "gray", lty = 2)
lines(true_beta_vec, col = "blue", lwd = 2, lty = 3)
legend("topleft", 
       legend = c("Posterior Mean", "95% CI", "True β(t)"),
       col = c("black", "gray", "blue"), 
       lty = c(1, 2, 3), 
       lwd = c(1, 1, 2))


# Plot of posterior Mean and 95% CI for for gamma(t)
plot(gamma_mean, type = "l", ylim = range(c(gamma_lower, gamma_upper)), 
     main = expression(paste("Posterior Mean and 95% CI of ", gamma(t))), 
     ylab = expression(gamma(t)), xlab = "Time")
lines(gamma_lower, col = "gray", lty = 2)
lines(gamma_upper, col = "gray", lty = 2)
lines(true_gamma_vec, col = "blue", lwd = 2, lty = 3)
legend("topright", 
       legend = c("Posterior Mean", "95% CI", "True γ(t)"),
       col = c("black", "gray", "blue"), 
       lty = c(1, 2, 3), 
       lwd = c(1, 1, 2))

## Step 8: Predictive Checks

# Plot of Observed vs Expected new claims over time 
expected_claims <- simulate_claim_SIS_tv(
  N = 500,
  beta_vec = beta_mean,
  gamma_vec = gamma_mean,
  IO = 20,
  T = T
)

plot(data$time, data$new_claims, type = "l",
     col = "darkred", lwd = 2,
     xlab = "Time", ylab = "Number of New Claims",
     main = "Observed vs Expected New Claims (Time-Varying SIS Model)")

lines(expected_claims$time, expected_claims$new_claims, 
      col = "blue", lwd = 2, lty = 2)

legend("topleft", legend = c("Observed", "Expected (Posterior Mean)"),
       col = c("darkred", "blue"), lty = c(1, 2), lwd = 2)
