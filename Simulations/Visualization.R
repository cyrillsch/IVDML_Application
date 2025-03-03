## R script to reproduce Figure 1
##################################################################

# to install the IVDML package:
# devtools::install_github("cyrillsch/IVDML")

library(IVDML)

# generate data
beta <- function(a){return(2 * exp(-0.5 * a^2))}
generate_data <- function(N, instr_setting){
  if(instr_setting == "linear"){
    psi <- function(z){z}
  } 
  if(instr_setting == "nonlinear"){
    psi <- function(z){cos(z) + 0.2 * z}
  }
  X <- rnorm(N)
  Z <- 0.5 * X + rnorm(N)
  f <- function(Z, X){return(-sin(X) + psi(Z))}
  g <- function(X){return(tanh(X))}
  H <- rnorm(N)
  del <- 0.7 * H + 0.7 * rnorm(N)
  eps <-  sign(H) - 0.5 + 0.5 * rnorm(N)
  D <- f(Z, X) + del
  Y <- beta(X) * D + g(X) + eps
  A <- X
  return(list(Y = Y, D = D, X = X, Z = Z, A = X))
}


N <- 1000
S_split <- 10

ml_method <- "gam"

set.seed(1214)

for(instr_setting in c("linear", "nonlinear")){
  data_sim <- generate_data(N, instr_setting)
  dml_result <- fit_IVDML(Y = data_sim$Y, D = data_sim$D, Z = data_sim$Z, X = data_sim$X, A = data_sim$A, 
                          ml_method = ml_method, A_deterministic_X = TRUE, iv_method = c("linearIV", "mlIV"), S_split = S_split)
  filename <- paste("VisualizationResults/", "instr", instr_setting, ".RData", sep = "")
  save(dml_result, file = filename)
}



pdf("Plots/Visualization.pdf", width = 8, height = 6)
par(mfrow = c(2,2))
par(mgp = c(2, 1, 0))
par(bg = "white")
aa <- seq(-3.8, 3.8, 0.1)
iv_vec <- c("linearIV", "mlIV")
col_vec <- c("red", "blue")
transparent_col_vec <- c(rgb(1, 0, 0, 0.1), rgb(0, 0, 1, 0.1))

for(instr_setting in c("linear", "nonlinear")){
  filename <- paste("VisualizationResults/", "instr", instr_setting, ".RData", sep = "")
  load(filename)
  h_normal <- bandwidth_normal(dml_result$A)
  h_undersmooth <- h_normal * N^0.2 / N^(2/7)
  for(h in c(h_normal, h_undersmooth)){
    if(instr_setting == "linear"){
      plot_title <- paste("(het.)/(Z lin.) \n h = ", round(h, 2), sep = "")
    }
    if(instr_setting == "nonlinear"){
      plot_title <- paste("(het.)/(Z nonlin.) \n h = ", round(h, 2), sep = "")
    }
    plot(aa, Vectorize(beta)(aa), type = "l", ylim = c(-1.3, 3), main = plot_title, xlab = "a", ylab = "beta(a)", cex.lab = 1.2)
    rug(dml_result$A)
    for(i in 1:length(iv_vec)){
      iv <- iv_vec[i]
      beta_hat <- Vectorize(function(a){
        return(coef(dml_result, iv, a, dml_result$A, "epanechnikov", h))
      })
      se_hat <- Vectorize(function(a){
        return(se(dml_result, iv, a, dml_result$A, "epanechnikov", h))
      })
      lines(aa, beta_hat(aa), col = col_vec[i])
      za <- qnorm(0.975)
      lines(aa, beta_hat(aa) + za * se_hat(aa), col = col_vec[i], lty = 3)
      lines(aa, beta_hat(aa) - za * se_hat(aa), col = col_vec[i], lty = 3)
      yy <- seq(-1.6, 3.2, length.out = 55) + 0.1 * i
      ay <- expand.grid(aa = aa, yy = yy)
      shaded <- function(a, y){
        return(robust_p_value_aggregated(dml_result, y, iv, a, dml_result$A, "epanechnikov", h) >= 0.05)
      }
      in_int <- mapply(shaded, ay$aa, ay$yy)
      in_int <- ay[in_int, ]
      points(in_int$aa, in_int$yy, col = transparent_col_vec[i], pch = 19, cex = 0.3)
    }
    legend("topright", legend = iv_vec, col = col_vec, lty = 1)
  }
}
dev.off()

