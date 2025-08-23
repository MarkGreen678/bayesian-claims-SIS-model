# ----------------------------------------
# Bayesian Inference: Static SIS Model
# ----------------------------------------

## Step 1: Simulate Policyholder Claim Data

simulate_claim_SIS <- function(N, beta, gamma, I0, T) {
  S <- numeric(T)
  I <- numeric(T)
  new_claims <- numeric(T)
  resolved_claims <- numeric(T)
  
  S[1] <- N - I0
  I[1] <- I0
  
  for (t in 1:(T - 1)) {
    p_claim <- 1 - exp(-beta * I[t] / N)
    p_resolve <- 1 - exp(-gamma)
    
    new_I <- rbinom(1, S[t], p_claim)
    new_S <- rbinom(1, I[t], p_resolve)
    
    new_claims[t + 1] <- new_I
    resolved_claims[t + 1] <- new_S
    
    S[t + 1] <- S[t] - new_I + new_S
    I[t + 1] <- I[t] + new_I - new_S
  }
  
  return(data.frame(
    time = 0:(T - 1), S = S, I = I,
    new_claims = new_claims, resolved_claims = resolved_claims
  ))
}

## Step 2: Log-likelihood for Claims

loglik_claim_SIS <- function(params, data) {
  beta <- params[1]
  gamma <- params[2]
  if (beta <= 0 || gamma <= 0) return(-Inf)
  
  N <- data$S[1] + data$I[1]
  loglik <- 0
  
  for (t in 1:(nrow(data) - 1)) {
    St <- data$S[t]
    It <- data$I[t]
    new_I <- data$new_claims[t + 1]
    new_S <- data$resolved_claims[t + 1]
    
    p_claim <- min(max(1 - exp(-beta * It / N), 1e-10), 1 - 1e-10)
    p_resolve <- min(max(1 - exp(-gamma), 1e-10), 1 - 1e-10)
    
    loglik <- loglik +
      dbinom(new_I, St, p_claim, log = TRUE) +
      dbinom(new_S, It, p_resolve, log = TRUE)
  }
  
  return(loglik)
}

## Step 3: Log-prior (Exponential)

log_prior <- function(beta, gamma, lambda_beta = 1, lambda_gamma = 1) {
  dexp(beta, rate = lambda_beta, log = TRUE) +
    dexp(gamma, rate = lambda_gamma, log = TRUE)
}

## Step 4: Metropolis-Hastings Inference

MH_sampler <- function(data, n_iter, beta_init, gamma_init, lambda_beta, lambda_gamma, proposal_sd) {
  beta <- beta_init
  gamma <- gamma_init
  samples <- matrix(NA, n_iter, 2)
  accepted <- 0
  
  for (i in 1:n_iter) {
    beta_prop <- rnorm(1, beta, proposal_sd)
    gamma_prop <- rnorm(1, gamma, proposal_sd)
    
    if (beta_prop > 0 && gamma_prop > 0) {
      loglik_curr <- loglik_claim_SIS(c(beta, gamma), data)
      loglik_prop <- loglik_claim_SIS(c(beta_prop, gamma_prop), data)
      
      logprior_curr <- log_prior(beta, gamma, lambda_beta, lambda_gamma)
      logprior_prop <- log_prior(beta_prop, gamma_prop, lambda_beta, lambda_gamma)
      
      log_alpha <- (loglik_prop + logprior_prop) - (loglik_curr + logprior_curr)
      
      if (log(runif(1)) < log_alpha) {
        beta <- beta_prop
        gamma <- gamma_prop
        accepted <- accepted + 1
      }
    }
    
    samples[i, ] <- c(beta, gamma)
  }
  
  cat("Acceptance rate:", accepted / n_iter, "\n")
  samples <- as.data.frame(samples)
  colnames(samples) <- c("beta", "gamma")
  return(samples)
}

## Step 5: Simulate & Run

set.seed(42)
true_beta <- 0.25
true_gamma <- 0.1

data <- simulate_claim_SIS(N = 500, beta = true_beta, gamma = true_gamma, I0 = 20, T = 100)

samples <- MH_sampler(
  data = data,
  n_iter = 5000,
  beta_init = 0.2,
  gamma_init = 0.2,
  lambda_beta = 1,
  lambda_gamma = 1,
  proposal_sd = 0.02
)

## Step 6: Analyze Posterior

# Histogram of Posterior of Beta
hist(samples$beta, col = "lightblue", breaks = 50,
     main = "Posterior of Beta", xlab = "beta")
abline(v = true_beta, col = "red", lwd = 2)
legend("topleft", legend = c("True Beta"),
       col = c("red"), lwd = 2)

# Histogram of Posterior of Gamma
hist(samples$gamma, col = "lightgreen", breaks = 50, main = "Posterior of Gamma", xlab = "gamma")
abline(v = true_gamma, col = "red", lwd = 2)
legend("topright", legend = c("True Gamma"),
       col = c("red"), lwd = 2)

# Trace plot of Beta
plot(samples$beta, type = "l", main = "Trace of Beta", xlab = "Iteration" )

# Trace plot of Gamma
plot(samples$gamma, type = "l", main = "Trace of Gamma", xlab = "Iteration")

# Autocorrelation plot of Beta
acf(samples$beta, main = "Autocorrelation of Beta Samples")

# Autocorrelation plot of Gamma
acf(samples$gamma, main = "Autocorrelation of Gamma Samples")

# Posterior Estimates of Parameters
mean(samples$beta)
quantile(samples$beta, probs = c(0.025, 0.975))  # 95% credible interval
# ie. “There is a 95% probability that the true claim rate β lies between ... and ...”

## Step 7: Forecast Future Claims

# Sample one posterior draw:
burn_in <- 500
beta_draw <- samples$beta[sample((burn_in+1):nrow(samples), 1)]
gamma_draw <- samples$gamma[sample((burn_in+1):nrow(samples), 1)]

# Simulate claim path from model
simulated <- simulate_claim_SIS(N = 1000, beta = beta_draw, gamma = gamma_draw, 5, 150)

# Plot claims over time
plot(simulated$time, simulated$new_claims, type = "l", main = "Simulated Claims Path", xlab = "time",
     ylab = "new claims")

# Plot number of Active and Susceptible Policyholders over Time
matplot(data$time, cbind(data$S, data$I), type = "l", lty = 1, 
        col = c("blue", "red"),
        xlab = "Time", ylab = "Number of Policyholders",
        main = "SIS Claim Dynamics: Susceptible vs Active Claims")

legend("topright", 
       legend = c("Susceptible (No Claim)", "Active Claims (I)"),
       col = c("blue", "red"), lty = 1)


## Step 8: Multiple simulations

n_sims <- 100

# Initiate matrix for future claims for all sims
future_claims <- matrix(0, nrow = 150, ncol = n_sims) 

# Run simulations
for (i in 1:n_sims) {
  beta_i <- samples$beta[sample(1:nrow(samples), 1)]
  gamma_i <- samples$gamma[sample(1:nrow(samples), 1)]
  
  sim <- simulate_claim_SIS(N = 1000, beta = beta_i, gamma = gamma_i, I0 = 5, T = 150)
  future_claims[, i] <- sim$new_claims
}

# Plot 5 simulated paths
matplot(future_claims[, 1:5], type = "l", lty = 1, col = rainbow(5), 
        main = "Posterior Predictive Claim Paths",  xlab = "time")
