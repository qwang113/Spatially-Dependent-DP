rm(list = ls())
library(MASS)
library(ggplot2)
library(sp)
library(mvtnorm)
sim_size <- 50
long_grid <- lat_grid <- seq(from = -2, to = 2, length.out = sim_size)
grid <- expand.grid(long = long_grid, lat = lat_grid)
coords <- as.matrix(grid)
n <- nrow(coords)

A <- matrix(c(1, 0, 0, 0.2), nrow = 2)
A_inv <- solve(A)

dists <- array(NA, dim = c(n, n, 2))
for (k in 1:2) {
  dists[,,k] <- outer(coords[,k], coords[,k], "-")
}
mahal_sq <- A_inv[1,1]*dists[,,1]^2 +
  (A_inv[1,2] + A_inv[2,1])*dists[,,1]*dists[,,2] +
  A_inv[2,2]*dists[,,2]^2
sigma2 <- 1
Sigma <- sigma2 * exp(-mahal_sq)


set.seed(123)
nobs <- 50
z <- mvrnorm(n = nobs, mu = rep(0, n), Sigma = Sigma)

# 
ggplot(grid, aes(x = long, y = lat, fill = z[1,])) +
  geom_raster() +
  scale_fill_viridis_c() +
  coord_fixed() +
  labs(title = "Anisotropic Gaussian Process Simulation")


get_info <- function(mat) {
  unique_rows <- unique(mat)  # Unique rows directly
  
  counts <- apply(unique_rows, 1, function(ur) {
    sum(apply(mat, 1, function(r) all(r == ur)))
  })
  
  unique_columns <- t(unique_rows)  # Each column is a unique row (transposed)
  n_unique <- ncol(unique_columns)  # Total number of unique rows
  
  return(list(
    unique_rows = unique_columns,
    counts = counts,
    n_unique = n_unique
  ))
}

train_grids <- sample(1:ncol(z), 100)


# Hyperparameters
anu <- 1
bnu <- 1
abeta <- mean(z)
bbeta <- var(as.vector(z))/3
asig = atau = 2
bsig = btau = var(as.vector(z))/3
bphi = 10
# Initialization

curr_theta <- matrix(0, nrow = nrow(z), ncol = length(train_grids))
curr_beta <- 0
curr_nu <- 1
curr_phi <- 5 
curr_tau <- 1
curr_sig <- 1
curr_idx <- 0
save_idx <- 0

burn = 0
thin = 1
tot_sample = 1000

all_theta <- array(NA, dim = c(tot_sample, dim(curr_theta)))
all_beta = all_nu = all_phi = all_tau = all_sig <- rep(NA, tot_sample)

paired_dist <- spDists(as.matrix(grid[train_grids,]))
train_dat <- z[,train_grids]
while (TRUE) {
  curr_idx = curr_idx + 1
  print(curr_idx)
  # Sample theta
  for (i in 1:nrow(z)) {
    LAM <- solve(diag( 1/curr_tau, ncol = length(train_grids) , nrow = length(train_grids)) +
                   1/curr_sig * solve(exp(-curr_phi*paired_dist)))
    q0 <- curr_nu * sqrt( det(LAM) ) *
      exp(-0.5/curr_tau* 
            t(train_dat[i,] - curr_beta)%*%
            (diag(1, nrow = length(train_grids), ncol = length(train_grids)) - 1/curr_tau*LAM)%*%
            (train_dat[i,] - curr_beta)
      )*(2*pi*curr_tau*curr_sig)^(-length(train_grids)/2)/sqrt(det(exp(-curr_phi*paired_dist)))
    
    theta_info <- get_info(curr_theta[-i,])
    qj_vec <- rep(NA, theta_info$n_unique)
    for (j in 1:theta_info$n_unique) {
      qj_vec[j] <- theta_info$counts[j]*dmvnorm(train_dat[i,], mean = curr_beta + theta_info$unique_rows[,j],
                                                sigma = diag(curr_tau, nrow = length(train_grids), ncol = length(train_grids)))
    }
    idx <- rmultinom(1,1,prob = c(q0,qj_vec))
    if(idx[1]==1){curr_theta[i,] = rmvnorm(1, mean = LAM%*%matrix(train_dat[i,]-curr_beta, ncol = 1)/curr_tau,
                                           sigma = LAM )}else{
                                             curr_theta[i,] <- theta_info$unique_rows[,which(idx==1)-1]}
  }
  
  # Sample beta
  xt <- rep(1, length(train_grids))
  sig_beta <- 1/(1/bbeta+1/curr_tau*nrow(train_dat)*sum(t(xt)%*%(xt)) )
  mu_beta <- sig_beta*(1/bbeta*abeta+1/curr_tau*sum(train_dat-curr_theta))
  curr_beta <- rnorm(1, mean = mu_beta, sd = sqrt(sig_beta))
  
  # Sample Tau
  curr_tau = 1/rgamma(1, shape = atau+length(train_dat)/2, rate = btau+1/2*sum((train_dat-curr_theta-curr_beta)^2) )
  
  # Sample nu
  temp_eta <- rbeta(1, curr_nu+1, nrow(train_dat))
  theta_info <- get_info(curr_theta)
  nstar <- theta_info$n_unique
  p <- (anu + nstar - 1)/(nrow(train_dat)*(bnu-log(temp_eta))+anu+nstar-1)
  l <- runif(1)
  if(l<p){curr_nu <- rgamma(1, shape = anu + nstar, rate = bnu-log(temp_eta))}else{
    curr_nu <- rgamma(1, shape = anu + nstar - 1, rate = bnu - log(temp_eta))
  }
  
  # Sample sigma
  shape_pos <- asig + 1/2*theta_info$n_unique*ncol(train_dat)
  rate_pos <- bsig + 1/2*sum(diag(t(theta_info$unique_rows)%*%solve(exp(-curr_phi*paired_dist))%*%theta_info$unique_rows))
  curr_sig <- 1/rgamma(1, shape_pos, rate_pos)
  
  # Sample phi
  curr_llh <- -theta_info$n_unique/2*log(det(exp(-curr_phi*paired_dist)))-
    sum(diag(t(theta_info$unique_rows)%*%solve(exp(-curr_phi*paired_dist))%*%theta_info$unique_rows))/(2*curr_sig)
  
  logit_curr <- log(curr_phi/bphi)-log(1-curr_phi/bphi)
  new_logit <- rnorm(1, logit_curr, sd = 0.1)
  new_phi <- exp(new_logit)/(1+exp(new_logit)) * bphi
  
  new_llh <- -theta_info$n_unique/2*log(det(exp(-new_phi*paired_dist)))-
    sum(diag(t(theta_info$unique_rows)%*%solve(exp(-new_phi*paired_dist))%*%theta_info$unique_rows))/(2*curr_sig)
  ptrans <- exp(new_llh - curr_llh)
  l = runif(1)
  if(l<ptrans){curr_phi = new_phi}
  
  if(curr_idx >= burn & curr_idx %% thin ==0){
    save_idx = save_idx + 1
    all_theta[save_idx,,] <- curr_theta
    all_beta[save_idx] <- curr_beta
    all_nu[save_idx] <- curr_nu
    all_sig[save_idx] <- curr_sig
    all_tau[save_idx] <- curr_tau
    all_phi[save_idx] <- curr_phi
    
    if(save_idx == tot_sample){break}
  }
  
}
library(ggplot2)

plot_pred <- function(idx){
  
  # Calculate the density of the first row of train_dat
  real_density_vals <- density(train_dat[idx, ])$y
  real_density_xvals <- density(train_dat[idx, ])$x
  
  # Create a data frame for the real density curve
  real_density_df <- data.frame(
    value = real_density_vals,
    type = "True",
    x = real_density_xvals
  )
  
  # Create a data frame to store estimated density curves
  all_theta_long <- data.frame(
    value = c(),
    type = c(),
    x = c()
  )
  # all_theta <- all_theta[500:1000,,]
  # Calculate the density for each row of all_theta and add it to the data frame
  for (i in 100:nrow(all_theta)) {
    density_vals <- density(all_theta[i, idx, ])$y
    x_vals <- density(all_theta[i, idx, ])$x
    
    all_theta_long <- rbind(all_theta_long, data.frame(value = density_vals, type = paste("Estimate ", i), x = x_vals))
  }
  
  # Combine the two data frames
  combined_df <- rbind(real_density_df, all_theta_long)
  
  ggplot() +
    geom_line(aes(x = all_theta_long$x, y = all_theta_long$value, group = all_theta_long$type), color = "blue", alpha = 0.2) +
    geom_line(aes(x = real_density_df$x, y = real_density_df$value), color = "red", linewidth = 0.8) +
    theme_minimal() +
    labs(x = paste("Theta",train_grids[idx]), y = "Density") +
    theme(legend.title = element_blank())  # 可选：去除图例标题
}
