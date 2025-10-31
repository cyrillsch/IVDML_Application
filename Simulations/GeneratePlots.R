################################################################################
# The script "RunSimulations.R" needs to be run first!
################################################################################

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

generate_plots_strength <- function(ml_method, instr_setting, beta_setting, dim_setting, error_dependence, N, strength_vec){
  if(beta_setting == "hom"){
    beta <- function(a){return(1)}
  } else {
    beta <- function(a){return(2 * exp(-0.5 * a^2))}
  }
  za <- qnorm(0.975)
  results_list <- list()
  for(strength in strength_vec){
    filename <- paste0("SimResults/N", N,"_ml", ml_method, "_instr", instr_setting, "_beta",
                      beta_setting, "_dim", dim_setting, "_strength", strength,
                      "_errordep", error_dependence, ".RData")
    load(filename)
    results_list[[as.character(strength)]] <- temp_result
  }
  calc_statistics <- function(NX){
    NX <- as.data.frame(NX)
    # hom.
    a <- 0
    MSE_hom_linearIV <- mean((NX$coef_hom_linearIV - beta(a))^2)
    MSE_hom_mlIV     <- mean((NX$coef_hom_mlIV - beta(a))^2)
    st_cov_hom_linearIV <- mean(abs(NX$coef_hom_linearIV - beta(a)) <= za * NX$se_hom_linearIV)
    st_cov_hom_mlIV     <- mean(abs(NX$coef_hom_mlIV - beta(a)) <= za * NX$se_hom_mlIV)
    rob_cov_hom_linearIV <- mean(NX$prob_hom_linearIV >= 0.05)
    rob_cov_hom_mlIV     <- mean(NX$prob_hom_mlIV >= 0.05)
    
    # het, a = 0
    a <- 0
    MSE_het0_linearIV <- mean((NX$coef_het0_linearIV - beta(a))^2)
    MSE_het0_mlIV     <- mean((NX$coef_het0_mlIV - beta(a))^2)
    st_cov_het0_linearIV <- mean(abs(NX$coef_het0_linearIV - beta(a)) <= za * NX$se_het0_linearIV)
    st_cov_het0_mlIV     <- mean(abs(NX$coef_het0_mlIV - beta(a)) <= za * NX$se_het0_mlIV)
    rob_cov_het0_linearIV <- mean(NX$prob_het0_linearIV >= 0.05)
    rob_cov_het0_mlIV     <- mean(NX$prob_het0_mlIV >= 0.05)
    
    # het, a = 1.5
    a <- 1.5
    MSE_het1_linearIV <- mean((NX$coef_het1_linearIV - beta(a))^2)
    MSE_het1_mlIV     <- mean((NX$coef_het1_mlIV - beta(a))^2)
    st_cov_het1_linearIV <- mean(abs(NX$coef_het1_linearIV - beta(a)) <= za * NX$se_het1_linearIV)
    st_cov_het1_mlIV     <- mean(abs(NX$coef_het1_mlIV - beta(a)) <= za * NX$se_het1_mlIV)
    rob_cov_het1_linearIV <- mean(NX$prob_het1_linearIV >= 0.05)
    rob_cov_het1_mlIV     <- mean(NX$prob_het1_mlIV >= 0.05)
    
    data.frame(MSE_hom_linearIV = MSE_hom_linearIV, MSE_hom_mlIV = MSE_hom_mlIV,
               st_cov_hom_linearIV = st_cov_hom_linearIV, st_cov_hom_mlIV = st_cov_hom_mlIV,
               rob_cov_hom_linearIV = rob_cov_hom_linearIV, rob_cov_hom_mlIV = rob_cov_hom_mlIV,
               MSE_het0_linearIV = MSE_het0_linearIV, MSE_het0_mlIV = MSE_het0_mlIV,
               st_cov_het0_linearIV = st_cov_het0_linearIV, st_cov_het0_mlIV = st_cov_het0_mlIV,
               rob_cov_het0_linearIV = rob_cov_het0_linearIV, rob_cov_het0_mlIV = rob_cov_het0_mlIV,
               MSE_het1_linearIV = MSE_het1_linearIV, MSE_het1_mlIV = MSE_het1_mlIV,
               st_cov_het1_linearIV = st_cov_het1_linearIV, st_cov_het1_mlIV = st_cov_het1_mlIV,
               rob_cov_het1_linearIV = rob_cov_het1_linearIV, rob_cov_het1_mlIV = rob_cov_het1_mlIV)
  }
  
  # aggregate
  results_agg <- do.call(rbind, lapply(results_list, calc_statistics))
  results_agg$strength <- strength_vec
  with(results_agg, {
    setting_title <- paste("(", beta_setting, ".)/(Z ", ifelse(instr_setting == "linear", "lin.", "nonlin."), ")", sep = "")
    
    plot(strength, MSE_het1_mlIV, col = "blue", type = "l", lty = 3,
         ylim = c(0.0005, max(c(MSE_het1_mlIV, MSE_het1_linearIV))),
         log = "y", xlab = "Strength", ylab = "MSE", main = paste("MSE \n", setting_title), cex.lab = 1.2)
    lines(strength, MSE_het1_linearIV, col = "red", lty = 3)
    lines(strength, MSE_het0_linearIV, col = "red", lty = 2)
    lines(strength, MSE_het0_mlIV, col = "blue", lty = 2)
    if(beta_setting == "hom"){
      lines(strength, MSE_hom_linearIV, col = "red")
      lines(strength, MSE_hom_mlIV, col = "blue")
    }
    
    # coverage standard CI
    plot(strength, st_cov_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.5, 1),
         xlab = "Strength", ylab = "Coverage", main = paste("Coverage of st. CI \n", setting_title), cex.lab = 1.2)
    lines(strength, st_cov_het1_linearIV, col = "red", lty = 3)
    lines(strength, st_cov_het0_linearIV, col = "red", lty = 2)
    lines(strength, st_cov_het0_mlIV, col = "blue", lty = 2)
    if(beta_setting == "hom"){
      lines(strength, st_cov_hom_mlIV, col = "blue")
      lines(strength, st_cov_hom_linearIV, col = "red")
    }
    abline(h = 0.95, col = "grey")
    if (beta_setting == "hom") {
      legend("bottomright",
             legend = c("linearIV", "mlIV", "hom. T.E.", "het. T.E., v = 0", "het. T.E., v = 1.5"),
             col = c("red", "blue", "black", "black", "black"), pch = c(20, 20, NA, NA, NA), lty = c(NA, NA, 1, 2, 3))
    } else {
      legend("bottomright",
             legend = c("linearIV", "mlIV", "het. T.E., v = 0", "het. T.E., v = 1.5"),
             col = c("red", "blue", "black", "black"), pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 3))
    }
    
    # coverage robust CI
    plot(strength, rob_cov_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.5, 1),
         xlab = "Strength", ylab = "Coverage", main = paste("Coverage of rob. CI \n", setting_title), cex.lab = 1.2)
    lines(strength, rob_cov_het1_linearIV, col = "red", lty = 3)
    lines(strength, rob_cov_het0_linearIV, col = "red", lty = 2)
    lines(strength, rob_cov_het0_mlIV, col = "blue", lty = 2)
    if(beta_setting == "hom"){
      lines(strength, rob_cov_hom_mlIV, col = "blue")
      lines(strength, rob_cov_hom_linearIV, col = "red")
    }
    abline(h = 0.95, col = "grey")
  })
}


for(i in 1:nrow(dim_ml_mat)){
  dim_setting <- dim_ml_mat[i, 1]
  ml_method <- dim_ml_mat[i, 2]
  for(beta_setting in beta_setting_vec){
    for(N in N_vec){
      plottitle <- paste0("Plots/IVSTRENGTH_N", N, "_dim", dim_setting,
                          "_ml", ml_method, "_beta", beta_setting, ".pdf")
      pdf(plottitle, width = 6.7, height = 0.75 * 6.7)
      par(mfrow = c(2, 3), mar = c(4.1, 3.1, 3.1, 1.6))
      par(mgp = c(2, 1, 0))
      par(bg = "white")
      for(instr_setting in instr_setting_vec){
        generate_plots_strength(ml_method, instr_setting, beta_setting, dim_setting, error_dependence, N, strength_vec)
      }
      dev.off()
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

generate_plots_N <- function(ml_method, instr_setting, beta_setting, dim_setting, error_dependence, N_vec, strength){
  if(beta_setting == "hom"){
    beta <- function(a){return(1)}
  } else {
    beta <- function(a){return(2 * exp(-0.5 * a^2))}
  }
  za <- qnorm(0.975)
  results_list <- list()
  for(N in N_vec){
    filename <- paste0("SimResults/N", N,"_ml", ml_method, "_instr", instr_setting, "_beta",
                       beta_setting, "_dim", dim_setting, "_strength", strength,
                       "_errordep", error_dependence, ".RData")
    load(filename)
    results_list[[as.character(N)]] <- temp_result
  }
  calc_statistics <- function(NX){
    NX <- as.data.frame(NX)
    # hom.
    a <- 0
    MSE_hom_linearIV <- mean((NX$coef_hom_linearIV - beta(a))^2)
    MSE_hom_mlIV     <- mean((NX$coef_hom_mlIV - beta(a))^2)
    st_cov_hom_linearIV <- mean(abs(NX$coef_hom_linearIV - beta(a)) <= za * NX$se_hom_linearIV)
    st_cov_hom_mlIV     <- mean(abs(NX$coef_hom_mlIV - beta(a)) <= za * NX$se_hom_mlIV)
    rob_cov_hom_linearIV <- mean(NX$prob_hom_linearIV >= 0.05)
    rob_cov_hom_mlIV     <- mean(NX$prob_hom_mlIV >= 0.05)
    
    # het, a = 0
    a <- 0
    MSE_het0_linearIV <- mean((NX$coef_het0_linearIV - beta(a))^2)
    MSE_het0_mlIV     <- mean((NX$coef_het0_mlIV - beta(a))^2)
    st_cov_het0_linearIV <- mean(abs(NX$coef_het0_linearIV - beta(a)) <= za * NX$se_het0_linearIV)
    st_cov_het0_mlIV     <- mean(abs(NX$coef_het0_mlIV - beta(a)) <= za * NX$se_het0_mlIV)
    rob_cov_het0_linearIV <- mean(NX$prob_het0_linearIV >= 0.05)
    rob_cov_het0_mlIV     <- mean(NX$prob_het0_mlIV >= 0.05)
    
    # het, a = 1.5
    a <- 1.5
    MSE_het1_linearIV <- mean((NX$coef_het1_linearIV - beta(a))^2)
    MSE_het1_mlIV     <- mean((NX$coef_het1_mlIV - beta(a))^2)
    st_cov_het1_linearIV <- mean(abs(NX$coef_het1_linearIV - beta(a)) <= za * NX$se_het1_linearIV)
    st_cov_het1_mlIV     <- mean(abs(NX$coef_het1_mlIV - beta(a)) <= za * NX$se_het1_mlIV)
    rob_cov_het1_linearIV <- mean(NX$prob_het1_linearIV >= 0.05)
    rob_cov_het1_mlIV     <- mean(NX$prob_het1_mlIV >= 0.05)
    
    data.frame(MSE_hom_linearIV = MSE_hom_linearIV, MSE_hom_mlIV = MSE_hom_mlIV,
               st_cov_hom_linearIV = st_cov_hom_linearIV, st_cov_hom_mlIV = st_cov_hom_mlIV,
               rob_cov_hom_linearIV = rob_cov_hom_linearIV, rob_cov_hom_mlIV = rob_cov_hom_mlIV,
               MSE_het0_linearIV = MSE_het0_linearIV, MSE_het0_mlIV = MSE_het0_mlIV,
               st_cov_het0_linearIV = st_cov_het0_linearIV, st_cov_het0_mlIV = st_cov_het0_mlIV,
               rob_cov_het0_linearIV = rob_cov_het0_linearIV, rob_cov_het0_mlIV = rob_cov_het0_mlIV,
               MSE_het1_linearIV = MSE_het1_linearIV, MSE_het1_mlIV = MSE_het1_mlIV,
               st_cov_het1_linearIV = st_cov_het1_linearIV, st_cov_het1_mlIV = st_cov_het1_mlIV,
               rob_cov_het1_linearIV = rob_cov_het1_linearIV, rob_cov_het1_mlIV = rob_cov_het1_mlIV)
  }
  
  # aggregate
  results_agg <- do.call(rbind, lapply(results_list, calc_statistics))
  results_agg$N <- N_vec
  with(results_agg, {
    setting_title <- paste("(", beta_setting, ".)/(Z ", ifelse(instr_setting == "linear", "lin.", "nonlin."), ")", sep = "")
    
    plot(N, MSE_het1_mlIV, col = "blue", type = "l", lty = 3,
         ylim = c(0.0005, max(c(MSE_het1_mlIV, MSE_het1_linearIV))),
         log = "y", xlab = "N", ylab = "MSE", main = paste("MSE \n", setting_title), cex.lab = 1.2)
    lines(N, MSE_het1_linearIV, col = "red", lty = 3)
    lines(N, MSE_het0_linearIV, col = "red", lty = 2)
    lines(N, MSE_het0_mlIV, col = "blue", lty = 2)
    if(beta_setting == "hom"){
      lines(N, MSE_hom_linearIV, col = "red")
      lines(N, MSE_hom_mlIV, col = "blue")
    }
    
    # coverage standard CI
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
             legend = c("linearIV", "mlIV", "hom. T.E.", "het. T.E., v = 0", "het. T.E., v = 1.5"),
             col = c("red", "blue", "black", "black", "black"), pch = c(20, 20, NA, NA, NA), lty = c(NA, NA, 1, 2, 3))
    } else {
      legend("bottomright",
             legend = c("linearIV", "mlIV", "het. T.E., v = 0", "het. T.E., v = 1.5"),
             col = c("red", "blue", "black", "black"), pch = c(20, 20, NA, NA), lty = c(NA, NA, 2, 3))
    }
    
    # coverage robust CI
    plot(N, rob_cov_het1_mlIV, col = "blue", type = "l", lty = 3, ylim = c(0.8, 1),
         xlab = "Strength", ylab = "Coverage", main = paste("Coverage of rob. CI \n", setting_title), cex.lab = 1.2)
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

for(i in 1:nrow(dim_ml_mat)){
  dim_setting <- dim_ml_mat[i, 1]
  ml_method <- dim_ml_mat[i, 2]
  for(beta_setting in beta_setting_vec){
    plottitle <- paste0("Plots/VARYN", "_dim", dim_setting,
                        "_ml", ml_method, "_beta", beta_setting, "strength", strength, ".pdf")
    pdf(plottitle, width = 6.7, height = 0.75 * 6.7)
    par(mfrow = c(2, 3), mar = c(4.1, 3.1, 3.1, 1.6))
    par(mgp = c(2, 1, 0))
    par(bg = "white")
    for(instr_setting in instr_setting_vec){
      generate_plots_N(ml_method, instr_setting, beta_setting, dim_setting, error_dependence, N_vec, strength)
    }
    dev.off()
  }
}