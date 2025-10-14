################################################################################
# Import source functions and fix some simulation parameters
################################################################################
source("SimulationSetup.R")

N_sim <- 1000
n_cores <- 60

S_split <- 10

# 1D with gam and 5D with xgboost
dim_ml_mat <- rbind(c("1D", "gam"),
                    c("5D", "xgboost"))

instr_setting_vec <- c("linear", "nonlinear")
beta_setting_vec <- c("hom", "het")


################################################################################
# Weak IV case, where we vary the IV strength and fix N and (strong) 
# error correlation
################################################################################

N_vec <- c(500, 1000)
strength_vec <- c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
error_dependence <- "strong"

for(i in 1:nrow(dim_ml_mat)){
  dim_setting <- dim_ml_mat[i, 1]
  ml_method <- dim_ml_mat[i, 2]
  print(paste("Dim: ", dim_setting))
  print(paste("ML: ", ml_method))
  for(beta_setting in beta_setting_vec){
    print(paste("beta: ", beta_setting))
    for(N in N_vec){
      print(paste("N: ", N))
      for(instr_setting in instr_setting_vec){
        print(paste("instr.: ", instr_setting))
        for(strength in strength_vec){
          save_sim_results(N_sim, n_cores, N, dim_setting, beta_setting, instr_setting, strength, error_dependence, ml_method, S_split)
        }
      }
    }
  }
}


################################################################################
# Strong IV case, where we vary N and fix the IV strength and (moderate)
# error correlation
################################################################################

N_vec <- c(250, 500, 1000, 1500, 2000)
strength <- 1
error_dependence <- "moderate"

for(i in 1:nrow(dim_ml_mat)){
  dim_setting <- dim_ml_mat[i, 1]
  ml_method <- dim_ml_mat[i, 2]
  print(paste("Dim: ", dim_setting))
  print(paste("ML: ", ml_method))
  for(beta_setting in beta_setting_vec){
    print(paste("beta: ", beta_setting))
    for(N in N_vec){
      print(paste("N: ", N))
      for(instr_setting in instr_setting_vec){
        print(paste("instr.: ", instr_setting))
        save_sim_results(N_sim, n_cores, N, dim_setting, beta_setting, instr_setting, strength, error_dependence, ml_method, S_split)
      }
    }
  }
}