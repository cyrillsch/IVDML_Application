# to install the IVDML package:
# devtools::install_github("cyrillsch/IVDML")

library(IVDML)
library(parallel)

# N is the sample size
# dim_setting is either "1D" or "5D"
# beta_setting is either "hom" or "het"
# instr_setting is either "linear" or "nonlinear"
# strength is a scalar with which the psi0 gets multiplied
# error_dependence is either "moderate" or "strong"
generate_data <- function(N, dim_setting, beta_setting, instr_setting, strength, error_dependence){
  if(beta_setting == "hom"){
    beta <- function(a){return(1)}
  }
  if(beta_setting == "het"){
    beta <- function(a){return(2 * exp(-0.5 * a^2))}
  }
  if(instr_setting == "linear"){
    psi0 <- function(z){z}
  } 
  if(instr_setting == "nonlinear"){
    psi0 <- function(z){cos(z) + 0.2 * z}
  }
  psi <- function(z){strength * psi0(z)}
  if(dim_setting == "1D"){
    X <- rnorm(N)
    Z <- 0.5 * X + rnorm(N)
    f <- function(Z, X){return(-sin(X) + psi(Z))}
    g <- function(X){return(tanh(X))}
  }
  if(dim_setting == "5D"){
    Sig <- cbind(c(1, 0.5, 0.5, 0.5, 0.5), c(0.5, 1, 0.5, 0.5, 0.5), c(0.5, 0.5, 1, 0.5, 0.5), c(0.5, 0.5, 0.5, 1, 0.5), c(0.5, 0.5, 0.5, 0.5, 1))
    cholSig <- chol(Sig)
    X0 <- matrix(rnorm(5 * N), ncol = 5)
    X <- X0 %*% cholSig
    Z <- 0.5 * X[,1] - 0.5 * X[, 2] + rnorm(N)
    f <- function(Z, X){return(-sin(X[, 1]) + X[, 2] + psi(Z))}
    g <- function(X){return(tanh(X[, 1])- X[, 3])}
  }
  H <- rnorm(N)
  if(error_dependence == "moderate"){
    del <- 0.7 * H + 0.7 * rnorm(N)
    eps <-  sign(H) - 0.5 + 0.5 * rnorm(N)
  }
  if(error_dependence == "strong"){
    del <- 0.7 * H + 0.1 * rnorm(N)
    eps <- 0.7 * H + 0.1 * rnorm(N)
  }
  D <- f(Z, X) + del
  Y <- beta(as.matrix(X)[, 1]) * D + g(X) + eps
  return(list(Y = Y, D = D, X = X, Z = Z, A = as.matrix(X)[, 1]))
}


one_sim <- function(N, dim_setting, beta_setting, instr_setting, strength, error_dependence, ml_method, S_split){
  data_sim <- generate_data(N, dim_setting = dim_setting, beta_setting = beta_setting, instr_setting = instr_setting, strength = strength, error_dependence = error_dependence)
  dml_result <- fit_IVDML(Y = data_sim$Y, D = data_sim$D, Z = data_sim$Z, X = data_sim$X, A = data_sim$A, 
                          A_deterministic_X = TRUE, ml_method = ml_method, iv_method = c("linearIV", "mlIV"), S_split = S_split)
  if(beta_setting == "hom"){
    beta <- function(a){return(1)}
  }
  if(beta_setting == "het"){
    beta <- function(a){return(2 * exp(-0.5 * a^2))}
  }
  # record point estimates, standard deviations, and p-value for robust CI for homogeneous treatment effect
  a <- 0
  coef_hom_linearIV <- coef(dml_result, iv_method = "linearIV")
  coef_hom_mlIV <- coef(dml_result, iv_method = "mlIV")
  se_hom_linearIV <- se(dml_result, iv_method = "linearIV")
  se_hom_mlIV <- se(dml_result, iv_method = "mlIV")
  candidate_value <- beta(a)
  prob_hom_linearIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "linearIV")
  prob_hom_mlIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "mlIV")
  
  # record point estimates, standard deviations, and p-value for robust CI for heterogeneous treatment effect
  bandwidth <- bandwidth_normal(data_sim$A) * length(data_sim$A)^0.2 / length(data_sim$A)^(2/7)
  kernel_name <- "epanechnikov"
  
  # a = 0 has a lot of data nearby
  a <- 0
  coef_het0_linearIV <- coef(dml_result, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  coef_het0_mlIV <- coef(dml_result, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  se_het0_linearIV <- se(dml_result, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  se_het0_mlIV <- se(dml_result, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  candidate_value <- beta(a)
  prob_het0_linearIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  prob_het0_mlIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  
  # a = 1.5 does not have a lot of data nearby
  a <- 1.5
  coef_het1_linearIV <- coef(dml_result, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  coef_het1_mlIV <- coef(dml_result, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  se_het1_linearIV <- se(dml_result, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  se_het1_mlIV <- se(dml_result, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  candidate_value <- beta(a)
  prob_het1_linearIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "linearIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  prob_het1_mlIV <- robust_p_value_aggregated(dml_result, candidate_value = candidate_value, iv_method = "mlIV", a = a, A = data_sim$A, kernel_name = kernel_name, bandwidth = bandwidth)
  
  result <- c(coef_hom_linearIV, coef_hom_mlIV, se_hom_linearIV, se_hom_mlIV, prob_hom_linearIV, prob_hom_mlIV, 
              coef_het0_linearIV, coef_het0_mlIV, se_het0_linearIV, se_het0_mlIV, prob_het0_linearIV, prob_het0_mlIV,
              coef_het1_linearIV, coef_het1_mlIV, se_het1_linearIV, se_het1_mlIV, prob_het1_linearIV, prob_het1_mlIV)
  
  return(result)
}


save_sim_results <- function(N_sim, n_cores, N, dim_setting, beta_setting, instr_setting, strength, error_dependence, ml_method, S_split){
  temp_result <- mclapply(1:N_sim, function(i){one_sim(N, dim_setting, beta_setting, instr_setting, strength, error_dependence, ml_method, S_split)}, mc.cores = n_cores)
  temp_result <- do.call(rbind, temp_result)
  colnames(temp_result) <- c("coef_hom_linearIV", "coef_hom_mlIV", "se_hom_linearIV", "se_hom_mlIV", "prob_hom_linearIV", "prob_hom_mlIV",
                             "coef_het0_linearIV", "coef_het0_mlIV", "se_het0_linearIV", "se_het0_mlIV", "prob_het0_linearIV", "prob_het0_mlIV",
                             "coef_het1_linearIV", "coef_het1_mlIV", "se_het1_linearIV", "se_het1_mlIV", "prob_het1_linearIV", "prob_het1_mlIV")
  filename <- paste0("SimResults/N", N,"_ml", ml_method, "_instr", instr_setting, "_beta",
                    beta_setting, "_dim", dim_setting, "_strength", strength,
                    "_errordep", error_dependence, ".RData")
  save(temp_result, file = filename)
}
