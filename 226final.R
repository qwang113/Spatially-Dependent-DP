rm(list = ls())
library(sp)
library(mvtnorm)
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
# Prior Simulation
b_phi = 0.1
a_sig = a_tau = b_sig = b_tau = 1
a_nu = b_nu = 1
coords <- as.matrix(expand.grid(1:50,1:50))
L <- 1000
nsim <- 25

prior_sim <- function(a_sig=1, b_sig = 1, a_tau = 1, b_tau = 1, a_nu = 1, b_nu = 1, b_phi = 0.1,
                      coords, L = 1000, nsim = 9){
  full_surface <- array(NA, dim = c(nsim,nrow(coords),3))
  for (sim in 1:nsim) {
    print(paste("Simulation",sim,"out of",nsim))
    phi <- runif(1,0,b_phi)
    nu <- rgamma(1, a_nu, b_nu)
    sigma2 <- 1/rgamma(1, a_sig, b_sig)
    dists <- spDists(coords)
    H <- sigma2*exp(-phi*dists)
    
    V <- rbeta(L, 1, nu)
    omegas <- rep(NA, L)
    omegas[1] <- V[1]
    for (i in 2:L) {
      omegas[i] <- V[i]*prod(1-V[1:(i-1)])
    }
    theta_0s <- rmvnorm(L, mean = rep(0,nrow(coords)), sigma = H)
    sim_surface <- theta_0s[rmultinom(1,1,omegas),]
    full_surface[sim,,] <- cbind(coords, sim_surface)
  }
  plot_df <- NULL
  for (sim in 1:nsim) {
    temp <- as.data.frame(full_surface[sim,,])
    colnames(temp) <- c("Lon", "Lat", "Value")
    temp$Sim <- sim
    plot_df <- bind_rows(plot_df, temp) 
  }
  plots <- lapply(1:nsim, function(sim) {
    ggplot(filter(plot_df, Sim == sim), aes(x = Lon, y = Lat, fill = Value)) +
      geom_tile() +
      scale_fill_viridis_c() +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none") 
    
  })
  
  return(wrap_plots(plots, ncol = sqrt(nsim)))
}




p1 <- prior_sim(coords = coords, b_phi = 0.01, nsim = 25, a_nu = 1)
p2 <- prior_sim(coords = coords, b_phi = 1, nsim = 25, a_nu = 1)
p3 <- prior_sim(coords = coords, b_phi = 100, nsim = 25, a_nu = 1)

p4 <- prior_sim(coords = coords, b_phi = 1, nsim = 25, a_nu = 0.1)
p5 <- prior_sim(coords = coords, b_phi = 1, nsim = 25, a_nu = 1)
p6 <- prior_sim(coords = coords, b_phi = 1, nsim = 25, a_nu = 100)

