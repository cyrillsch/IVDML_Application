## R script to reproduce Figure 2, Figure 3, Figure 6, Figure 7 and Figure 8
##################################################################

# to install the IVDML package:
# devtools::install_github("cyrillsch/IVDML")

library(IVDML)
library(parallel)

# N is the sample size
# dim_setting is either "1D" or "5D"
# beta_setting is either "hom" or "het"
# instr_setting is either "linear" or "nonlinear"
# strength_setting is either "weak" or "strong"
generate_data <- function(N, dim_setting, beta_setting, instr_setting, strength_setting){
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
  if(strength_setting == "strong"){
    psi <- function(z){psi0(z)}
  }
  if(strength_setting == "weak"){
    psi <- function(z){0.1 * psi0(z)}
  }
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
  del <- 0.7 * H + 0.7 * rnorm(N)
  eps <-  sign(H) - 0.5 + 0.5 * rnorm(N)
  D <- f(Z, X) + del
  Y <- beta(as.matrix(X)[, 1]) * D + g(X) + eps
  return(list(Y = Y, D = D, X = X, Z = Z, A = as.matrix(X)[, 1]))
}

strength_setting_vec <- c("strong", "weak")
dim_setting_vec <- c("1D", "5D")
instr_setting_vec <- c("linear", "nonlinear")
beta_setting_vec <- c("hom", "het")

N_vec <- c(250, 500, 1000, 1500, 2000)

S_split <- 10

N_sim <- 500

n_cores <- 20

ml_method_vec <- c("gam", "xgboost")



one_sim <- function(N, dim_setting, beta_setting, instr_setting, strength_setting, ml_method, S_split){
  data_sim <- generate_data(N, dim_setting = dim_setting, beta_setting = beta_setting, instr_setting = instr_setting, strength_setting = strength_setting)
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


# Dim is either "1D" or "5D"
RNGkind("L'Ecuyer-CMRG")
set.seed(1214)
strength_setting <- "strong"
for(ml_method in ml_method_vec){
  print(paste("ML: ", ml_method))
  for(dim_setting in dim_setting_vec){
    print(paste("Dim: ", dim_setting))
    for(beta_setting in beta_setting_vec){
      print(paste("beta : ", beta_setting))
      for(instr_setting in instr_setting_vec){
        print(paste("Instrument : ", instr_setting))
        list_result <- list()
        for(N in N_vec){
          temp_result <- mclapply(1:N_sim, function(i){one_sim(N, dim_setting, beta_setting, instr_setting, strength_setting, ml_method, S_split)}, mc.cores = n_cores)
          temp_result <- do.call(rbind, temp_result)
          colnames(temp_result) <- c("coef_hom_linearIV", "coef_hom_mlIV", "se_hom_linearIV", "se_hom_mlIV", "prob_hom_linearIV", "prob_hom_mlIV",
                                     "coef_het0_linearIV", "coef_het0_mlIV", "se_het0_linearIV", "se_het0_mlIV", "prob_het0_linearIV", "prob_het0_mlIV",
                                     "coef_het1_linearIV", "coef_het1_mlIV", "se_het1_linearIV", "se_het1_mlIV", "prob_het1_linearIV", "prob_het1_mlIV")
          list_result[[paste("N",N, sep = "")]] <- temp_result
          print(paste("N = ", N, " done!"))
        }
        filename <- paste("SimulationsMainResults/","ml", ml_method, "_instr", instr_setting, "_beta",
                          beta_setting, "_dim", dim_setting, "_strength", strength_setting, ".RData", sep = "")
        save(list_result, file = filename)
      }
    }
  }
}









## Generate Plots

generate_plots <- function(ml_method, instr_setting, beta_setting, dim_setting, strength_setting){
  if(beta_setting == "hom"){
    beta <- function(a){return(1)}
  }
  if(beta_setting == "het"){
    beta <- function(a){return(2 * exp(-0.5 * a^2))}
  }
  filename <- paste("SimulationsMainResults/","ml", ml_method, "_instr", instr_setting, "_beta",
                    beta_setting, "_dim", dim_setting, "_strength", strength_setting, ".RData", sep = "")
  load(filename)
  za <- qnorm(0.975)
  calc_statistics <- function(NX){
    NX <- as.data.frame(NX)
    # record MSE, standard deviations and coverages of CIs for hom. treatment effect
    a <- 0
    MSE_hom_linearIV <- mean((NX$coef_hom_linearIV - beta(a))^2)
    MSE_hom_mlIV <- mean((NX$coef_hom_mlIV- beta(a))^2)
    st_cov_hom_linearIV <- mean(abs(NX$coef_hom_linearIV - beta(a)) <= za * NX$se_hom_linearIV)
    st_cov_hom_mlIV <- mean(abs(NX$coef_hom_mlIV - beta(a)) <= za * NX$se_hom_mlIV)
    rob_cov_hom_linearIV <- mean(NX$prob_hom_linearIV >= 0.05)
    rob_cov_hom_mlIV <- mean(NX$prob_hom_mlIV >= 0.05)
    # record MSE, standard deviations and coverages of CIs for het. treatment effect
    a <- 0
    MSE_het0_linearIV <- mean((NX$coef_het0_linearIV - beta(a))^2)
    MSE_het0_mlIV <- mean((NX$coef_het0_mlIV- beta(a))^2)
    st_cov_het0_linearIV <- mean(abs(NX$coef_het0_linearIV - beta(a)) <= za * NX$se_het0_linearIV)
    st_cov_het0_mlIV <- mean(abs(NX$coef_het0_mlIV - beta(a)) <= za * NX$se_het0_mlIV)
    rob_cov_het0_linearIV <- mean(NX$prob_het0_linearIV >= 0.05)
    rob_cov_het0_mlIV <- mean(NX$prob_het0_mlIV >= 0.05)

    a <- 1.5
    MSE_het1_linearIV <- mean((NX$coef_het1_linearIV - beta(a))^2)
    MSE_het1_mlIV <- mean((NX$coef_het1_mlIV- beta(a))^2)
    st_cov_het1_linearIV <- mean(abs(NX$coef_het1_linearIV - beta(a)) <= za * NX$se_het1_linearIV)
    st_cov_het1_mlIV <- mean(abs(NX$coef_het1_mlIV - beta(a)) <= za * NX$se_het1_mlIV)
    rob_cov_het1_linearIV <- mean(NX$prob_het1_linearIV >= 0.05)
    rob_cov_het1_mlIV <- mean(NX$prob_het1_mlIV >= 0.05)
    return(data.frame(MSE_hom_linearIV = MSE_hom_linearIV, MSE_hom_mlIV = MSE_hom_mlIV, st_cov_hom_linearIV = st_cov_hom_linearIV,
                      st_cov_hom_mlIV = st_cov_hom_mlIV, rob_cov_hom_linearIV = rob_cov_hom_linearIV, rob_cov_hom_mlIV = rob_cov_hom_mlIV,
                      MSE_het0_linearIV = MSE_het0_linearIV, MSE_het0_mlIV = MSE_het0_mlIV, st_cov_het0_linearIV = st_cov_het0_linearIV,
                      st_cov_het0_mlIV = st_cov_het0_mlIV, rob_cov_het0_linearIV = rob_cov_het0_linearIV, rob_cov_het0_mlIV = rob_cov_het0_mlIV,
                      MSE_het1_linearIV = MSE_het1_linearIV, MSE_het1_mlIV = MSE_het1_mlIV, st_cov_het1_linearIV = st_cov_het1_linearIV,
                      st_cov_het1_mlIV = st_cov_het1_mlIV, rob_cov_het1_linearIV = rob_cov_het1_linearIV, rob_cov_het1_mlIV = rob_cov_het1_mlIV))
  }
  results_agg <- do.call(rbind, lapply(list_result, calc_statistics))
  results_agg$N <- N_vec
  with(results_agg,
       {setting_title <- paste("(", beta_setting, ".)/(Z ", ifelse(instr_setting == "linear", "lin.", "nonlin."), ")", sep = "")

        plot(N, MSE_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.0005, max(c(MSE_het1_mlIV, MSE_het1_linearIV))),
             log = "y", xlab = "N", ylab = "MSE", main = paste("MSE \n", setting_title), cex.lab = 1.2)
        lines(N, MSE_het1_linearIV, col = "red", lty = 3)
        lines(N, MSE_het0_linearIV, col = "red", lty = 2)
        lines(N, MSE_het0_mlIV, col = "blue", lty = 2)
        if(beta_setting == "hom"){
          lines(N, MSE_hom_linearIV, col = "red")
          lines(N, MSE_hom_mlIV, col = "blue")
        }

        plot(N, st_cov_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.8, 1),
             xlab = "N", ylab = "Coverage", main = paste("Coverage of st. CI \n", setting_title), cex.lab = 1.2)
        lines(N, st_cov_het1_linearIV, col = "red", lty = 3)
        lines(N, st_cov_het0_linearIV, col = "red", lty = 2)
        lines(N, st_cov_het0_mlIV, col = "blue", lty = 2)
        if(beta_setting == "hom"){
          lines(N, st_cov_hom_mlIV, col = "blue")
          lines(N, st_cov_hom_linearIV, col = "red")
        }
        abline(h = 0.95, col = "grey")
        if (beta_setting == "hom") {
          legend("bottomright",
                 legend = c("linearIV", "mlIV", "hom. T.E.", "het. T.E., a = 0", "het. T.E., a = 1.5"),
                 col = c("red", "blue", "black", "black", "black"), pch = c(20, 20, NA, NA, NA), lty = c(NA, NA, 1, 2, 3))
        } else {
          legend("bottomright",
                 legend = c("linearIV", "mlIV", "het. T.E., a = 0", "het. T.E., a = 1.5"),
                 col = c("red", "blue", "black", "black"), pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 3))
        }

        plot(N, rob_cov_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.8, 1),
             xlab = "N", ylab = "Coverage", main = paste("Coverage of rob. CI \n", setting_title), cex.lab = 1.2)
        lines(N, rob_cov_het1_linearIV, col = "red", lty = 3)
        lines(N, rob_cov_het0_linearIV, col = "red", lty = 2)
        lines(N, rob_cov_het0_mlIV, col = "blue", lty = 2)
        if(beta_setting == "hom"){
          lines(N, rob_cov_hom_mlIV, col = "blue")
          lines(N, rob_cov_hom_linearIV, col = "red")
        }
        abline(h = 0.95, col = "grey")
    })
}




## Plots for Section 4
strength_setting <- "strong"
dim_setting <- "1D"
ml_method <- "gam"


for(beta_setting in beta_setting_vec){
  pdf(paste("Plots/SimGam1D_beta", beta_setting, ".pdf", sep = ""), width = 6.7, height = 0.75 * 6.7)
  par(mfrow = c(2, 3), mar = c(4.1, 3.1, 3.1, 1.6))
  par(mgp = c(2, 1, 0))
  par(bg = "white")
  for(instr_setting in instr_setting_vec){
    generate_plots(ml_method, instr_setting, beta_setting, dim_setting, strength_setting)
  }
  dev.off()
}


## Plots for Appendix E

strength_setting <- "strong"
dim_setting <- "5D"
ml_method <- "gam"

for(beta_setting in beta_setting_vec){
  pdf(paste("Plots/SimGam5D_beta", beta_setting, ".pdf", sep = ""), width = 6.7, height = 0.75 * 6.7)
  par(mfrow = c(2, 3), mar = c(4.1, 3.1, 3.1, 1.6))
  par(mgp = c(2, 1, 0))
  par(bg = "white")
  for(instr_setting in instr_setting_vec){
    generate_plots(ml_method, instr_setting, beta_setting, dim_setting, strength_setting)
  }
}

strength_setting <- "strong"
dim_setting <- "5D"
ml_method <- "xgboost"

for(beta_setting in beta_setting_vec){
  pdf(paste("Plots/SimXgboost5D_beta", beta_setting, ".pdf", sep = ""), width = 6.7, height = 0.75 * 6.7)
  par(mfrow = c(2, 3), mar = c(4.1, 3.1, 3.1, 1.6))
  par(mgp = c(2, 1, 0))
  par(bg = "white")
  for(instr_setting in instr_setting_vec){
    generate_plots(ml_method, instr_setting, beta_setting, dim_setting, strength_setting)
  }
}






