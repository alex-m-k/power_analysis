# This function simulates the values derived in closed form solution
library(fishmethods)

simulate <- function(cluster_autocorr){
  
  # k, 24 clusters
  k <- 24
  
  # m, number of associates in each cluster
  m_assoc <- 122
  n <- sum(m_assoc)
  
  # Assign treatment to cluster
  treatment <- if((k %% 2) == 0){
    sample(rep(0:1, k/2))
  } else {
    c(sample(rep(0:1, k/2)), 
      sample(0:1, 1))  
  }  
  
  # Put associates in clusters
  fc <- rep(1:k, m_assoc)
  
  # sigma, overall SD, coming from Kossek et al.
  sigma <- 4.809
  
  # ICC, setting equal to 0.02
  rho <- 0.02
  
  # Generate cluster variance
  rho_c <- cluster_autocorr
  # Set time invariant variance
  sigma_c <- rnorm(k, mean = 0, sd = sigma)
  # Set variance in period 1
  sigma_1 <- rnorm(k, mean = 0, sd = sigma)
  # Set variance in period 2
  sigma_2 <- rnorm(k, mean = 0, sd = sigma)
  # Set overall variance
  sigma_c1 <- (sqrt(rho_c) * sigma_c) + (sqrt(1-rho_c) * sigma_1)
  sigma_c2 <- (sqrt(rho_c) * sigma_c) + (sqrt(1-rho_c) * sigma_2)
  
  # Generate associate variance, setting rho_p equal to 0.55
  rho_p <- 0.55
  # Set time invariant variance
  sigma_p <- rnorm(n, mean = 0, sd = sigma)
  # Set variance in period 1
  sigma_1 <- rnorm(n, mean = 0, sd = sigma)
  # Set variance in period 2
  sigma_2 <- rnorm(n, mean = 0, sd = sigma)
  # Set overall variance
  sigma_p1 <- (sqrt(rho_p) * sigma_p) + (sqrt(1-rho_p) * sigma_1)
  sigma_p2 <- (sqrt(rho_p) * sigma_p) + (sqrt(1-rho_p) * sigma_2)
  
  # Calculate rho?
  var(sigma_c) / (var(sigma_c1) + var(sigma_c) + var(sigma_p1) + var(sigma_p))
  var(sigma_c) / (var(sigma_c2) + var(sigma_c) + var(sigma_p2) + var(sigma_p))
  
  sqrt(.55 * (1-rho) * sigma*sigma)
  sqrt((1-.55)* (1-rho) * sigma*sigma)
  
  # r value
  r <- ((m * rho) / (1 + ((m-1)*rho))*rho_c) + ((1- rho) / (1 + ((m-1)*rho))*rho_p)
  print(paste0("r equals ", r))
  
  # Generate outcome values
  # Baseline, with mean equal to 9.25 from Kossek et al.
  y0 <- 9.5 + (sqrt(rho) * sigma_c1[fc]) + (sqrt(1-rho) * sigma_p1)
  # Post-treatment, adding a 0.07 decrease for common time trend 
  # This was the decrease among controls in Kossek et al. 2019
  y1 <- 9.5 + (sqrt(rho) * sigma_c2[fc]) + (sqrt(1-rho) * sigma_p2)
  
  # Add treatment effect, setting value to -0.527
  y1 <- ifelse(treatment[fc] == 1,
               y1 + 0.527,
               y1)
  
  # Create table
  dat_dd <- data.frame(y = c(y0, y1), 
                       fc = rep(fc, 2),
                       assoc = rep(1:n, 2),
                       period = rep(0:1, c(n, n)),
                       d = rep(treatment[fc], 2))
  
  clus.rho(dat_dd$y, dat_dd$fc)
  clus.rho(dat_dd$y, dat_dd$assoc)
  clus.rho(dat_cov$y0, dat_cov$fc)
  clus.rho(dat_cov$y1, dat_cov$fc)
  
  summary_aov = summary(aov(y ~ as.factor(fc),data=dat_dd))
  summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2])
  summary_aov = summary(aov(y ~ as.factor(assoc),data=dat_dd))
  summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2])
  summary_aov = summary(aov(y0 ~ as.factor(fc),data=dat_cov))
  summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2])
  summary_aov = summary(aov(y1 ~ as.factor(fc),data=dat_cov))
  summary_aov[[1]][1,2]/sum(summary_aov[[1]][,2])
  
  

  # Verify that ICCs were correctly set
  
  y0_mean <- mean(y0)
  
  # Generate cluster component
  # Period 1
  c_mu0 <- data.frame(fc, y0) %>%
    group_by(fc) %>%
    mutate(group_mean = mean(y0),
           mean_i = y0-group_mean,
           sum_square = sum(mean_i^2))
  
  y0_between_var <- var(unique(c_mu0$group_mean))
  
  y0_within_var <- sum(c_mu0$sum_square) / (length(y0) - length(unique(fc)))
    
  y0_between_var / (y0_between_var + y0_within_var)
  
  # Period 2
  c_mu1 <- data.frame(fc, y1) %>%
    group_by(fc) %>%
    mutate(group_mean = mean(y1),
           mean_i = y1-group_mean,
           sum_square = sum(mean_i^2))
  
  y1_between_var <- var(unique(c_mu1$group_mean))
  
  y1_within_var <- sum(c_mu1$sum_square) / (length(y1) - length(unique(fc)))
  
  y1_between_var / (y1_between_var + y1_within_var)
  
  # rho_c estimate
  # Not sure if this is working
  cor(c_mu1$group_mean, c_mu0$group_mean)
  
  # Generate associate component
  p0 <- y0 - c_mu0$group_mean[fc]
  p1 <- y1 - c_mu1$group_mean[fc]
  
  var(p0)
  var(p1)
  
  
  aov(y0~as.factor(fc)+as.factor(assoc), dat_cov)
  aov(y1~as.factor(fc)+as.factor(assoc), dat_cov)
  
  
  # rho_p estimate
  # Not sure if this is working
  cor(p1, p0)
  
  dat_cov <- data.frame(y0 = y0,
                        y1 = y1,
                        fc = fc,
                        assoc = 1:n,
                        d = treatment[fc])
  
  return(list("dd" = dat_dd, "cov" = dat_cov))
}

power_test <- function(sims, cluster_autocorr){
  n.sims <- sims
  signif <- cbind("dd" = rep(NA, n.sims),
                  "cov" = rep(NA, n.sims))
  for(s in 1:n.sims){
    gen <- simulate(cluster_autocorr)
    
    # I don't think the errors should be clustered
    felm.dd <- felm(y~d*period | 0 | 0 | 0, gen$dd)
    # beta.hat <- coef(felm.dd)["d:period"]
    # beta.se <- felm.dd$cse["d:period"]
    # signif[s] <- ifelse((beta.hat - 2*beta.se) > 0,
    #                     1,
    #                     0)
    signif[s, "dd"] <- ifelse(felm.dd$pval["d:period"] < 0.05,
                              1,
                              0)
    felm.cov <- felm(y1~d+y0 | 0 | 0 | 0, gen$cov)
    signif[s, "cov"] <- ifelse(felm.cov$pval["d"] < 0.05,
                               1,
                               0)
  }
  power.dd <- mean(signif[, "dd"])
  power.cov <- mean(signif[, "cov"])
  return(c(power.dd, power.cov))
}

simulations <- 100

results  <- mapply(power_test, 
                   sim = simulations, 
                   cluster_autocorr = c(0.48))


library(swCRTdesign)
design <- swDsn(clusters = c(24, 24),
                tx.effect.frac = c(0, 0.5))

swPwr(design, 
      distn = "gaussian", 
      n = 2928,
      m0 = 9,
      m1 = 9.6,
      sigma = 4.8,
      tau = 4.8,
      eta = 4.8,
      rho = 0.02,
      icc = 0.3)
