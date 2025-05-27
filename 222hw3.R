rm(list = ls())
setwd("C:/Users/lix23/Desktop/")
dat = y = unlist(read.table("C:/Users/lix23/Desktop/hwk3-data.txt"))
get_info <- function(theta){
  n_star <- length(unique(theta))
  theta_star <- as.numeric(names(table(theta)))
  nj <- table(theta)
  return(list("nstar" = n_star, "theta_star" = theta_star,"nj" = nj))
}
# Hyper parameter choice
amu = mean(dat)
aphi = atau = 2
bphi = bmu = btau = var(dat)/3
aalpha = 3
balpha = 1
burn = 0
thin = 1
tot_sample = 500
n = length(dat)


DPM_Sampling <- function(y = y, amu = amu, bmu = bmu, aphi = aphi, bphi = bphi, aalpha = aalpha, balpha = balpha,
                         atau = atau, btau = btau, burn = burn, thin = thin, tot_sample = tot_sample){
  # Initialize
  n <- length(y)
  theta_all <- matrix(NA,nrow = n, ncol = tot_sample)
  mu_all = phi_all = tau_all = alpha_all = rep(NA, tot_sample)
  curr_idx = save_idx = 0
  curr_theta <- rep(0,n)
  curr_phi = curr_tau = 1
  curr_alpha = 10
  curr_mu = mean(y)
  
  while(TRUE){
    curr_idx <- curr_idx + 1
    print(curr_idx)
    # Sample theta
    for (i in 1:n) {
      # Calculate current q0
      theta_info <- get_info(curr_theta[-i])
      q0 <- 1/sqrt(2*pi*(curr_phi+curr_tau))*exp(-0.5*(y[i]-curr_mu)^2/(curr_phi + curr_tau))
      denom <- curr_alpha * q0 + sum(theta_info$nj*dnorm(y[i],mean = theta_info$theta_star, sd = sqrt(curr_phi)))
      normal_weights <- curr_alpha*q0/denom
      dirac_weights <- theta_info$nj*dnorm(y[i], mean = theta_info$theta_star, sd = sqrt(curr_phi))/denom
      selection <- rmultinom(1,1,prob = c(normal_weights,dirac_weights))
      if(selection[1]==1){ curr_theta[i] <- rnorm(1, 
                                                  mean = (curr_tau*y[i]+ curr_phi*curr_mu)/(curr_phi+curr_tau), 
                                                  sd = 1/sqrt(1/curr_tau+1/curr_phi) )
      }else{
        curr_theta[i] <- theta_info$theta_star[(which(selection==1) - 1)]
      }
    }
    
    # Sample phi
    curr_phi <- 1/rgamma(1, shape = aphi+n/2, rate = bphi+sum((y-curr_theta)^2)/2)
    
    # Sample mu
    theta_info <- get_info(curr_theta)
    curr_mu <- rnorm(1, 
                     mean = (curr_tau*amu+theta_info$nstar*mean(theta_info$theta_star)*bmu)/(curr_tau+theta_info$nstar*bmu),
                     sd = 1/sqrt(1/bmu+theta_info$nstar/curr_tau) )
    # Sample tau
    curr_tau <- 1/rgamma(1, shape = atau + theta_info$nstar/2, rate = btau + sum((theta_info$theta_star-curr_mu)^2)/2)
    
    # Sample alpha
    eta <- rbeta(1, curr_alpha+1, n)
    eps <- (aalpha + theta_info$nstar-1)/(n*(balpha-log(eta))+aalpha+theta_info$nstar-1)
    l <- runif(1)
    if(l <= eps){
      curr_alpha = rgamma(1, shape = aalpha+theta_info$nstar, rate = balpha-log(eta))
    }else{
      curr_alpha = rgamma(1, shape = aalpha+theta_info$nstar-1, rate = balpha-log(eta))
    }
    
    if(curr_idx >= burn & (curr_idx - burn)%%thin == 0 ){
      save_idx = save_idx + 1
      theta_all[,save_idx] <- curr_theta
      mu_all[save_idx] <- curr_mu
      tau_all[save_idx] <- curr_tau
      phi_all[save_idx] <- curr_phi
      alpha_all[save_idx] <- curr_alpha
      if(save_idx == tot_sample){break}
    }
  }
  return(list("theta" = theta_all, "mu" = mu_all, "tau" = tau_all, "phi" = phi_all, "alpha" = alpha_all))
}

S1 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 3, balpha = 1,
             atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S2 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 20, bphi = var(dat)/3*19, aalpha = 3, balpha = 1,
             atau = 20, btau = var(dat)/3*19, burn = 100, thin = 1, tot_sample = 500)

S3 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 200, bphi = var(dat)/3*199, aalpha = 3, balpha = 1,
             atau = 200, btau = var(dat)/3*199, burn = 100, thin = 1, tot_sample = 500)
save(S1,S2,S3, file = "226hw3.Rdata")
saveRDS(list(S1,S2,S3), "225hw3.Rda")

prev_res <- readRDS("225hw3.Rda")

library(ggplot2)
library(dplyr)

phi_df <- bind_rows(
  data.frame(phi = prev_res[[1]]$phi, group = "a=2"),
  data.frame(phi = prev_res[[2]]$phi, group = "a=20"),
  data.frame(phi = prev_res[[3]]$phi, group = "a=200")
)

p1 <- 
  ggplot(phi_df, aes(x = phi, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Density Comparison of phi",
       x = "phi", y = "Density",
       color = "Group", fill = "Group") +
  theme(legend.position = "bottom")


tau_df <- bind_rows(
  data.frame(tau = prev_res[[1]]$tau, group = "a=2"),
  data.frame(tau = prev_res[[2]]$tau, group = "a=20"),
  data.frame(tau = prev_res[[3]]$tau, group = "a=200")
)

p2 <- 
  ggplot(tau_df, aes(x = tau, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Density Comparison of tau",
       x = "tau", y = "Density",
       color = "Group", fill = "Group") +
  theme(legend.position = "bottom")

mu_df <- bind_rows(
  data.frame(mu = prev_res[[1]]$mu, group = "a=2"),
  data.frame(mu = prev_res[[2]]$mu, group = "a=20"),
  data.frame(mu = prev_res[[3]]$mu, group = "a=200")
)

p3 <- 
  ggplot(mu_df, aes(x = mu, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Density Comparison of mu",
       x = "mu", y = "Density",
       color = "Group", fill = "Group") +
  theme(legend.position = "bottom")

cowplot::plot_grid(p1,p2,p3, nrow = 1)



S4 <- DPM_Sampling(y = y, amu = -5, bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 3, balpha = 1,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S5 <- DPM_Sampling(y = y, amu = 0, bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 3, balpha = 1,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S6 <- DPM_Sampling(y = y, amu = 5, bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 3, balpha = 1,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)


mu_df <- bind_rows(
  data.frame(mu = S4$mu, group = "a=-5"),
  data.frame(mu = S5$mu, group = "a=0"),
  data.frame(mu = S6$mu, group = "a=5")
)

p3 <- 
  ggplot(mu_df, aes(x = mu, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  theme_minimal() +
  labs(title = "Density Comparison of mu",
       x = "mu", y = "Density",
       color = "Group", fill = "Group") +
  theme(legend.position = "bottom")

full_results <- list(S4,S5,S6)

# Predictive Density
pred_samples <- matrix(NA, nrow = 3, ncol = tot_sample)
for (i in 1:tot_sample) {
  for (j in 1:3) {
    curr_info <- get_info(full_results[[j]]$theta[,i])
    theta_star <- curr_info$theta_star
    weight_vector <- c(full_results[[j]]$alpha[i]/(full_results[[j]]$alpha[i]+n),curr_info$nj/(full_results[[j]]$alpha[i]+n))
    idx <- rmultinom(1, 1, weight_vector)
    if(idx[1]==1){
      print("A new table")
      pred_samples[j,i] <- rnorm(1, mean = rnorm(1, mean = full_results[[j]]$mu[i], sd = sqrt(full_results[[j]]$tau[i])),
                                 sd = sqrt(full_results[[j]]$phi[i]))
    }else{
      print("Existing table")
      pred_samples[j,i] <- rnorm(1, mean = theta_star[which(idx==1) - 1], sd = sqrt(full_results[[j]]$phi[i]))
  }
  }
}

group_labels <- c("mu = -5", "mu = 0", "mu = 5")

df_long <- as.data.frame(pred_samples) %>%
  mutate(group = group_labels) %>%
  pivot_longer(
    cols = -group,
    names_to = "sample_id",
    values_to = "value"
  )

x_grid <- seq(min(df_long$value, na.rm = TRUE) - 1,
              max(df_long$value, na.rm = TRUE) + 1,
              length.out = 1000)
true_density <- data.frame(
  x = x_grid,
  y = 0.2 * dnorm(x_grid, mean = -5, sd = 1) +
    0.5 * dnorm(x_grid, mean = 0, sd = 1) +
    0.3 * dnorm(x_grid, mean = 3.5, sd = 1)
)

ggplot(df_long, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  geom_line(data = true_density, aes(x = x, y = y),
            color = "red", size = 1.2, inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Posterior Predictive Sensitivity Analysis",
       x = "Value", y = "Density",
       color = NULL, fill = NULL) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.text = element_text(size = 11)
  )

S7 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 2, balpha = 15,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S8 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 2, balpha = 4,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S9 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 2, balpha = 0.9,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)

S10 <- DPM_Sampling(y = y, amu = mean(y), bmu = var(dat)/3, aphi = 2, bphi = var(dat)/3, aalpha = 2, balpha = 0.1,
                   atau = 2, btau = var(dat)/3, burn = 100, thin = 1, tot_sample = 500)


library(ggplot2)
library(dplyr)

alpha_df <- bind_rows(
  data.frame(alpha = S7$alpha, label = "aalpha=2, balpha=15"),
  data.frame(alpha = S8$alpha, label = "aalpha=2, balpha=4"),
  data.frame(alpha = S9$alpha, label = "aalpha=2, balpha=0.9"),
  data.frame(alpha = S10$alpha, label = "aalpha=2, balpha=0.1")
)

alpha_df$label <- factor(alpha_df$label,
                         levels = c("aalpha=2, balpha=15",
                                    "aalpha=2, balpha=4",
                                    "aalpha=2, balpha=0.9",
                                    "aalpha=2, balpha=0.1"))
pa <- 
ggplot(alpha_df, aes(x = alpha, fill = label)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 100, color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Posterior Alpha under Different Priors",
       x = expression(alpha), y = "Count",
       fill = "Prior") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.text = element_text(size = 11)
  )

unique_counts7 <- apply(S7$theta, 2, function(col) length(unique(col)))
unique_counts8 <- apply(S8$theta, 2, function(col) length(unique(col)))
unique_counts9 <- apply(S9$theta, 2, function(col) length(unique(col)))
unique_counts10 <- apply(S10$theta, 2, function(col) length(unique(col)))

df_counts <- bind_rows(
  data.frame(counts = unique_counts7, label = "aalpha=2, balpha=15"),
  data.frame(counts = unique_counts8, label = "aalpha=2, balpha=4"),
  data.frame(counts = unique_counts9, label = "aalpha=2, balpha=0.9"),
  data.frame(counts = unique_counts10, label = "aalpha=2, balpha=0.1")
)

df_counts$label <- factor(df_counts$label,
                          levels = c("aalpha=2, balpha=15",
                                     "aalpha=2, balpha=4",
                                     "aalpha=2, balpha=0.9",
                                     "aalpha=2, balpha=0.1"))
pns <- 
ggplot(df_counts, aes(x = counts, fill = label)) +
  geom_histogram(alpha = 0.4, position = "identity", bins = 30, color = "black") +
  theme_minimal() +
  labs(title = "Histogram of Unique Theta Counts per Sample",
       x = "Number of Unique Theta Values", y = "Count",
       fill = "Prior Setting") +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.text = element_text(size = 11)
  )
cowplot::plot_grid(pa,pns, ncol = 1)



full_results = list(S7,S8,S9,S10)
pred_samples <- matrix(NA, nrow = 4, ncol = tot_sample)
for (i in 1:tot_sample) {
  for (j in 1:4) {
    curr_info <- get_info(full_results[[j]]$theta[,i])
    theta_star <- curr_info$theta_star
    weight_vector <- c(full_results[[j]]$alpha[i]/(full_results[[j]]$alpha[i]+n),curr_info$nj/(full_results[[j]]$alpha[i]+n))
    idx <- rmultinom(1, 1, weight_vector)
    if(idx[1]==1){
      print("A new table")
      pred_samples[j,i] <- rnorm(1, mean = rnorm(1, mean = full_results[[j]]$mu[i], sd = sqrt(full_results[[j]]$tau[i])),
                                 sd = sqrt(full_results[[j]]$phi[i]))
    }else{
      print("Existing table")
      pred_samples[j,i] <- rnorm(1, mean = theta_star[which(idx==1) - 1], sd = sqrt(full_results[[j]]$phi[i]))
    }
  }
}

group_labels <- c("aalpha=2, balpha=15",
                  "aalpha=2, balpha=4",
                  "aalpha=2, balpha=0.9",
                  "aalpha=2, balpha=0.1")

library(lubridate)
library(dplyr)
library(tidyr)
df_long <- as.data.frame(pred_samples) %>%
  mutate(group = group_labels) %>%
  pivot_longer(
    cols = -group,
    names_to = "sample_id",
    values_to = "value"
  )

x_grid <- seq(min(df_long$value, na.rm = TRUE) - 1,
              max(df_long$value, na.rm = TRUE) + 1,
              length.out = 1000)
true_density <- data.frame(
  x = x_grid,
  y = 0.2 * dnorm(x_grid, mean = -5, sd = 1) +
    0.5 * dnorm(x_grid, mean = 0, sd = 1) +
    0.3 * dnorm(x_grid, mean = 3.5, sd = 1)
)

ggplot(df_long, aes(x = value, color = group, fill = group)) +
  geom_density(alpha = 0.3) +
  geom_line(data = true_density, aes(x = x, y = y),
            color = "red", size = 1.2, inherit.aes = FALSE) +
  theme_minimal() +
  labs(title = "Posterior Predictive Sensitivity Analysis",
       x = "Value", y = "Density",
       color = NULL, fill = NULL) +
  theme(
    legend.position = "bottom",
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.text = element_text(size = 11)
  )


# Posterior Predictive Distribution for Theta
pp_thetas <- rep(NA, tot_sample)
for (i in 1:tot_sample) {
  theta_info = get_info(S9$theta[,i])
  curr_alpha = S9$alpha[i]
  curr_mu = S9$mu[i]
  curr_tau = S9$tau[i]
  idx = rmultinom(1, 1, c(curr_alpha, theta_info$nj))
  if(idx[1]==1){pp_thetas[i] <- rnorm(1, mean = curr_mu, sd = sqrt(curr_tau)) }else{
    pp_thetas[i] <- theta_info$theta_star[which(idx == 1) - 1]
  }
}

df_pp <- data.frame(pp_thetas = pp_thetas)

ggplot(df_pp, aes(x = pp_thetas)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "gray80", color = "black", alpha = 0.6) +
  geom_density(color = "blue", fill = "skyblue", alpha = 0.4, linewidth = 1) +
  labs(title = "Posterior Predictive Density with Histogram",
       x = expression(theta),
       y = "Density") +
  theme_minimal()

sample_DP <- function(alpha, n, mu, tau,  theta, K = 100) {
  # Stick-breaking weights
  V <- rbeta(K, 1, alpha)
  pi <- numeric(K)
  pi[1] <- V[1]
  for (k in 2:K) {
    pi[k] <- V[k] * prod(1 - V[1:(k-1)])
  }
  
  # Atom locations from base distribution
  idx <- rmultinom(K,1,c(alpha/(alpha+n), rep(1/(alpha+n), n)))
  theta_out <- rep(NA, K)
  for (k in 1:K) {
    if(idx[1,k]==1){theta_out[k] <- (rnorm(1, mean = mu, sd = sqrt(tau)))}else{
      theta_out[k] <- (theta[which(idx[,k]==1)-1])
  }
  }
  return(list(weights = pi, atoms = theta_out))
}

out <- list()
for (i in 1:tot_sample) {
  print(i)
  out[[i]] <- sample_DP(S10$alpha[i], n = n, mu = S10$mu[i], tau = S10$tau[i], theta = S10$theta[,i],
                        K = 1000)
}



plot_data <- list()

for (i in 1:tot_sample) {
  dfi <- data.frame(
    theta = out[[i]]$atoms,
    weight = out[[i]]$weights,
    draw = paste0("Draw ", i)
  )
  plot_data[[i]] <- dfi
}

df_all <- bind_rows(plot_data)

df_all_sorted <- df_all %>% 
  arrange(draw, theta) %>%
  group_by(draw) %>%
  mutate(cum_weight = cumsum(weight))

ggplot(df_all_sorted, aes(x = theta, y = cum_weight, group = draw)) +
  geom_step(linewidth = 0.7, direction = "hv", alpha = 0.2) +
  labs(title = "Posterior Draws from DP (Step CDF)",
       x = expression(theta), y = "CDF") +
  theme_minimal() + 
  theme(legend.position = "")

grid_y <- seq(from = -10, to = 10, by = 0.01)
all_density <- matrix(NA, nrow = tot_sample, ncol = length(grid_y))
for (i in 1:tot_sample) {
  print(i)
  for(j in 1:length(grid_y)){
    all_density[i,j] <- sum(out[[i]]$weights*dnorm(grid_y[j], mean = out[[i]]$atoms, sd = sqrt(S10$phi[i])))
  }
}


mean_density <- apply(all_density, 2, mean)
lower_density <- apply(all_density, 2, quantile, probs = 0.025)
upper_density <- apply(all_density, 2, quantile, probs = 0.975)


true_density <- 0.2 * dnorm(grid_y, mean = -5, sd = 1) +
  0.5 * dnorm(grid_y, mean = 0, sd = 1) +
  0.3 * dnorm(grid_y, mean = 3.5, sd = 1)


summary_df <- data.frame(
  y = grid_y,
  mean = mean_density,
  lower = lower_density,
  upper = upper_density,
  truth = true_density
)


ggplot(summary_df, aes(x = y)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.4) +
  geom_line(aes(y = mean), color = "blue", linewidth = 1) +
  geom_line(aes(y = truth), color = "red", linewidth = 1, linetype = "dashed") +
  labs(
    title = "Posterior Density with 95% Interval vs. True Density",
    x = expression(y),
    y = "Density"
  ) +
  theme_minimal()



