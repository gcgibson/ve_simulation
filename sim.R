library(ggplot2)
library(cowplot)
library(rjags)
library(R2jags)

estimate_ve <- function(transmission_pairs_df){
  model_code <- "
model
{
  # Likelihood
  for (i in 1:T) {
    y[i] ~ dbern(p[i])
    logit(p[i]) <- alpha + beta_1*vac_bin_contact[i]  + beta_2*vac_bin_case[i]
  }

  for (i in 1:T) {
      E_Y_A_1[i] <- p_y_a_1[i]
      logit(p_y_a_1[i]) <- alpha + beta_1 + beta_2*vac_bin_case[i]
  }
  for (i in 1:T) {
      E_Y_A_0[i] <- p_y_a_0[i]
      logit(p_y_a_0[i]) <- alpha  + beta_2*vac_bin_case[i]
  }
  # Priors
  alpha ~ dnorm(0.0, .01)
  beta_1 ~ dnorm(0.0, .01)
  beta_2 ~ dnorm(0.0, .01)

}
"

  # Set up the data
  model_data <- list(T = nrow(transmission_pairs_df), y = transmission_pairs_df$inf,
                     vac_bin_contact= transmission_pairs_df$vax,
                     vac_bin_case= transmission_pairs_df$id_vax)
  model_parameters <- c("E_Y_A_0", "E_Y_A_1")

  model_run <- jags(
    data = model_data,
    parameters.to.save = model_parameters,
    model.file = textConnection(model_code),
    n.chains = 4,
    n.iter = 1000,
    n.burnin = 200,
    n.thin = 2
  )
  ve_hat <- (1-mean(model_run$BUGSoutput$sims.list$E_Y_A_1)/mean(model_run$BUGSoutput$sims.list$E_Y_A_0))
  return (ve_hat)

}





inner_loop_sim <- function(rho,p_v,p_vv,p_uv,p_vu,p_uu,obs_p){
  multiplier <- 2
  inf_status <- rep(0,1000)
  vac_status <- rep(0,1000)
  g <- sample_fitness_pl(no.of.nodes = 1000,no.of.edges = 200,exponent.out = 100)
  mat <- as_adjacency_matrix(g,sparse = F)

  ### randomly choose 10 index cases
  index_case_ids <- sample(1000,100)
  contact_ids <- setdiff(1:1000,index_case_ids)
  #### update status of indiiduals to infected if index case
  inf_status[index_case_ids ] <- 1.0

  initial_inf_status <- inf_status

  #### randomly assign index case vaccination status
  index_case_status <- rbinom(n=length(index_case_ids),size=1,p=p_v)
  while (sum(index_case_status)  == length(index_case_status)){
    index_case_status <- rbinom(n=length(index_case_ids),size=1,p=p_v)
  }
  vac_status[index_case_ids ] <- rbinom(n=length(index_case_ids),size=1,p=p_v)

  #### randomly assign contact vax status
  #vac_status[contact_ids] <- rbinom(n = length(contact_ids),size = 1,p=.75)


  transmission_pairs_df <- matrix(NA,ncol=3)

  #### propagate a single time step
  for (i in index_case_ids){
    ### iterate through contacts
    for (nbors in which(mat[i,]==1)){

      ### introduce correlation
      if (vac_status[i] == 1){
        vac_status[nbors]  = rbinom(1,1,prob=rho)
      } else{
        vac_status[nbors]  = rbinom(1,1,prob=1-rho)
      }

      if (vac_status[nbors] == 1 & vac_status[i] == 1){
        inf_status[nbors] <-  rbinom(1,1,p_vv*multiplier)
        if (inf_status[nbors] == 1){
          transmission_pairs_df <- rbind(transmission_pairs_df,c(1,1,1))
        } else{
          transmission_pairs_df <- rbind(transmission_pairs_df,c(0,1,1))
        }
      } else if(vac_status[nbors] == 1 & vac_status[i] == 0){
        inf_status[nbors] <-  rbinom(1,1,p_uv*multiplier)
        if (inf_status[nbors] == 1){
          transmission_pairs_df <- rbind(transmission_pairs_df,c(1,1,0))
        } else{
          transmission_pairs_df <- rbind(transmission_pairs_df,c(0,1,0))
        }
      } else if(vac_status[nbors] == 0 & vac_status[i] == 1){
        inf_status[nbors] <-  rbinom(1,1,p_vu*multiplier)
        if (inf_status[nbors] == 1){
          transmission_pairs_df <- rbind(transmission_pairs_df,c(1,0,1))
        } else{
          transmission_pairs_df <- rbind(transmission_pairs_df,c(0,0,1))
        }
      } else if(vac_status[nbors] == 0 & vac_status[i] == 0){
        inf_status[nbors] <-  rbinom(1,1,p_uu*multiplier)


        if (inf_status[nbors] == 1){
          if (rbinom(1,1,obs_p) > 0 ){

            transmission_pairs_df <- rbind(transmission_pairs_df,c(1,0,0))
          }
        } else{
          transmission_pairs_df <- rbind(transmission_pairs_df,c(0,0,0))
        }



      }
    }
  }

  new_infections = inf_status[contact_ids] - initial_inf_status[contact_ids]

  new_infections_vac = new_infections[which(vac_status[contact_ids] == 1)]
  new_infections_unvac = new_infections[which(vac_status[contact_ids] == 0)]

  num <- sum(new_infections_vac)/length(new_infections_vac)
  denom <- sum(new_infections_unvac)/length(new_infections_unvac)

  ve_hat_unadj  <- 1- num/denom

  transmission_pairs_df <- transmission_pairs_df[2:nrow(transmission_pairs_df),]
  transmission_pairs_df <- data.frame(inf = transmission_pairs_df[,1],vax=transmission_pairs_df[,2],
                                      id_vax=transmission_pairs_df[,3])
  fit2 <- glm(inf~ vax ,family = "poisson",data =transmission_pairs_df )
  fit <- glm(inf~ vax + id_vax ,family = "poisson",data =transmission_pairs_df )
#1-exp(fit$coefficients)[2]#
  ret_list <- list()
  ret_list[[1]] <-1-exp(fit2$coefficients)[2]
  ret_list[[2]] <-estimate_ve(transmission_pairs_df)

  return (ret_list)
}

run_sim_of_ve <- function(rho,obs_p){

  ### marginal probability of vaccination of index case
  p_v = .75
  p_u = 1-p_v

  ### conditional probability of infection given contact with
  ## vaccinated and unvaccinated
  p_vv = .28
  p_uv = .29

  ### conditional probability of infection given
  ## vaccinated and unvaccinated
  p_vu = .3
  p_uu = .4

  #### true VE as of EQ 3
  true_ve = 1- (p_vv*p_v + p_uv*p_u)/(p_vu*p_v + p_uu*p_u)


  #### simulate random graph over contacts and index cases
  library(igraph)
  ve_unadj_hats <- c()
  ve_adj_hats <- c()

  for (mcmc in 1:100){
    ret_list <- inner_loop_sim(rho,p_v,p_vv,p_uv,p_vu,p_uu,obs_p)
    ve_unadj_hats <- c(ve_unadj_hats,ret_list[[1]])
    ve_adj_hats <- c(ve_adj_hats,ret_list[[2]])

  }

  ve_unadj <- median(ve_unadj_hats[!is.infinite(ve_unadj_hats) & !is.nan(ve_unadj_hats)],na.rm=T)
  ve_unadj_var <- var(ve_unadj_hats[!is.infinite(ve_unadj_hats) & !is.nan(ve_unadj_hats)],na.rm=T)

  ve_adj <- median(ve_adj_hats[!is.infinite(ve_adj_hats) & !is.nan(ve_adj_hats)],na.rm = T)
  ve_adj_var <- var(ve_adj_hats[!is.infinite(ve_adj_hats) & !is.nan(ve_adj_hats)],na.rm = T)


  ret_list <- list()
  ret_list[[1]] <- true_ve
  ret_list[[2]] <- ve_unadj
  ret_list[[3]] <- ve_adj
  ret_list[[4]] <- ve_unadj_var
  ret_list[[5]] <- ve_adj_var

  return (ret_list)
}


### set observation probability of unvaccinated contact to test
## positive
obs_p <-1

## define sequence of homophily to test
rho_to_test  <- seq(.1,.9,by=.1)


true_ve <- c()
ve_hat_adj <- c()
ve_hat_adj_var <- c()

ve_hat_unadj <- c()
ve_hat_unadj_var <- c()

for (rho in rho_to_test){
  sim_res <- run_sim_of_ve(rho,obs_p)
  true_ve <- c(true_ve, sim_res[[1]])
  ve_hat_adj <- c(ve_hat_adj,sim_res[[3]])
  ve_hat_adj_var <- c(ve_hat_adj_var, sim_res[[5]])
  ve_hat_unadj <- c(ve_hat_unadj,sim_res[[2]])
  ve_hat_unadj_var <- c(ve_hat_adj_var, sim_res[[4]])


}


########
### Plotting
########

df_for_plot_obs <- data.frame(rho=rho_to_test[1:length(true_ve)],
                              ve_hat_adj = ve_hat_adj,
                              ve_hat_unadj = ve_hat_unadj,
                              true_ve = true_ve)

p1_obs <- ggplot(df_for_plot_obs,aes(x=rho,y=(ve_hat_adj-true_ve)/true_ve,col='Adjusted')) + geom_line() + theme_bw() + ylab("Bias") + xlab("Network Vaccination Correlation") +
  geom_line(aes(x=rho,y=(ve_hat_unadj-true_ve)/true_ve,col='Unadjusted')) + geom_hline(yintercept = 0,linetype="dashed")
p1_obs <- p1_obs +  theme(legend.title=element_blank())
print (p1_obs)
