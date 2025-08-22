# Bayesian Claims SIS Model

This project implements a Bayesian framework for modeling and forecasting insurance claim dynamics using both:
- A static SIS model with fixed transmission (β) and resolution (γ) rates
- A time-varying SIS model, where β(t) and γ(t) evolve over time via random walk priors

Metropolis-Hastings sampling is used for parameter inference, and posterior predictive simulations are conducted to assess model fit and forecast future claims.

Note: This code has been adapted from SIS modeling examples presented during a short course at the LMS Summer School 2025.

## Overview

The model tracks the dynamics of policyholders:
- **Susceptible (S):** Not currently claiming
- **Infected (I):** Currently claiming

Two core parameters:
- `β(t)`: Claim initiation rate (can vary over time)
- `γ(t)`: Claim resolution rate (can vary over time)

This project implements both:
- A **static** version with constant `β`, `γ`
- A **time-varying** version with vector-valued `β(t)`, `γ(t)`

---

## Contents

| r code and README               | Description                                                        |
|--------------------|--------------------------------------------------------------------|
| `claims_model_static.R`    | Static SIS model with fixed β, γ and MH inference                  |
| `claims_model_time_varying.R` | Time-varying β(t), γ(t) inference via component-wise MH sampling |
| `README.md`         | Project overview (this file)                                      |

| Figures                                                                                   | Description                                           |
| -------------------------------------------------------------------------------------- | ----------------------------------------------------- |
| `Histogram_of_Posterior_of_Beta_(Static_SIS_model).png`                                | Posterior distribution of fixed β                     |
| `Histogram_of_Posterior_of_Gamma_(Static_SIS_model).png`                               | Posterior distribution of fixed γ                     |
| `Trace_of_Beta_(Static_SIS_model).png`                                                 | MCMC trace plot for fixed β                           |
| `Trace_of_Gamma_(Static_SIS_model).png`                                                | MCMC trace plot for fixed γ                           |
| `Plot_of_Simulated_claims_path_(Static_SIS_model).png`                                 | Simulated claims path over time (Static SIS model)    |
| `Plot_of_Posterior_Predictive_Claim_Paths_for_five_sims_(Static_SIS_model).png`        | Posterior predictive claim paths (static model, five simulations)       |
| `Plot_number_of_Active_and_Susceptible_Policyholders_over_Time_(Static_SIS_model).png` | SIS model dynamics (S and I)                          |
| `Plot_of_Posterior_Mean_and_95pc_CI_for_beta(t).png`                                   | Posterior mean and 95% CI for β(t)                    |
| `Plot_of_Posterior_Mean_and_95pc_CI_for_gamma(t).png`                                  | Posterior mean and 95% CI for γ(t)                    |
| `Plot_Observed_vs_Expected_New_Claims_(Time-Varying_SIS_Model).png`                    | Observed vs predicted new claims (time-varying model) |

---

## Features

- Bayesian inference using custom Metropolis-Hastings samplers
- Exponential priors (static) and random walk priors (time-varying)
- Posterior summaries with credible intervals
- Posterior predictive simulations
- Visualizations: parameter traces, posteriors, autocorrelation plots and claim paths.

---

## Example Outputs

 - Posterior distributions of static parameters β, γ (histograms, trace plots)
 - Time-varying posterior paths of β(t), γ(t) with 95% credible intervals
 - Simulated SIS dynamics showing susceptible and active policyholders over time
 - Observed vs expected new claims under the time-varying SIS model
 - Posterior predictive simulations: multiple forward-simulated claim trajectories
 - Simulated claims path using posterior draws (static SIS model)
 - Autocorrelation plots for both static and time-varying models

---

## How to Run

1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/bayesian-claims-SIS-model.git
   cd bayesian-claims-SIS-model

---

## Summary of Findings

- **Posterior inference successfully recovers true parameters**  
  For both static and time-varying models, the Metropolis-Hastings algorithm accurately estimates the transmission (`β`) and resolution (`γ`) rates from simulated claim data.

- **Time-varying model captures dynamic trends**  
  The model tracks evolving claim rates over time by allowing `β(t)` and `γ(t)` to vary, improving realism and forecasting accuracy in non-stationary environments.

- **Posterior credible intervals are well-calibrated**  
  The 95% credible intervals for `β(t)` and `γ(t)` consistently include the true simulated paths, validating the uncertainty quantification.

- **Predictive performance aligns with observations**  
  Simulated claim trajectories based on posterior estimates closely match the observed number of new claims, both in shape and scale.

- **Component-wise MH sampling is effective**  
  Updating parameters one at a time allows flexible inference for high-dimensional latent paths (e.g., full trajectories of `β(t)`), with reasonable acceptance rates.

- **Drift priors improve inference precision** \
  Including a prior drift term helps regularize the posterior, leading to narrower credible intervals and improved tracking of underlying parameter trends over time.

- **Posterior samples exhibit low autocorrelation** \
  The MCMC chains for both static and time-varying parameters demonstrate rapidly decaying autocorrelation, indicating efficient exploration of the posterior space and yielding high effective sample sizes.


