library(tidyverse)
fit_bg_comp <-function(dat_obs){

  ##
  sl_libs <- c('SL.glm')

  dat_obs_ <- dat_obs
  colnames(dat_obs_) <- c("Y","A","W")
  dat_obs <- dat_obs_
  Y <- dat_obs$Y
  W_A <- dat_obs %>% select(-Y) # remove the outcome to make a matrix of predictors (A, W1, W2, W3, W4) for SuperLearner


  Q <- stan_glm(Y ~ A + W, data = dat_obs,
           family = binomial(link = "logit"))

  W_A1 <- W_A %>% mutate(A = 1)  # data set where everyone received treatment
  W_A0 <- W_A %>% mutate(A = 0) # data set where no one received treatment

  D_Q_A <- ((rstanarm::posterior_epred(Q,dat_obs,type = "response"))) # obtain predictions for everyone using the treatment they actually received
  D_Q_1 <- ((posterior_epred(Q,newdata = W_A1,type = "response")))# predict on that everyone-exposed data set
  D_Q_0 <- (posterior_epred(Q, newdata = W_A0,type="response"))

  ve <- rep(0,1000)
  for (mcmc in 1:1000){
    Q_A_update <- D_Q_A[mcmc,] # updated expected outcome given treatment actually received
    Q_1_update <- D_Q_1[mcmc,] # updated expected outcome for everyone receiving treatment
    Q_0_update <- D_Q_0[mcmc,]# updated expected outcome for everyone not receiving treatment
    ve[mcmc] <- 1- mean(Q_1_update)/mean(Q_0_update) # mean diff in updated expected outcome estimates

  }
  return (mean(ve))
}
