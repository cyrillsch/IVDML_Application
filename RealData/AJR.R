## R script to reproduce the tables and figures from Section 4.3.1
##################################################################


# to install the IVDML package:
# devtools::install_github("cyrillsch/IVDML")

library(IVDML)
library(hdm)
data(AJR)

Y <- AJR$GDP
D <- AJR$Exprop
Z <- AJR$logMort

X <- cbind(AJR$Latitude, AJR$Africa, AJR$Asia, AJR$Namer, AJR$Samer)
A <- AJR$Latitude


ml_vec <- c("gam", "xgboost", "randomForest")
iv_method <- c("linearIV", "mlIV", "mlIV_direct")

S_split <- 200

set.seed(1013)

list_res <- list()
for(ml in ml_vec){
  tic <- Sys.time()
  list_res[[ml]] <- fit_IVDML(Y = Y, D = D, Z = Z, X = X, A = A, ml_method = ml,
                            ml_par = list(), A_deterministic_X = TRUE, K_dml = 5, iv_method = iv_method, S_split = S_split)
  toc <- Sys.time()
  print(paste(ml, " took ", as.numeric(toc - tic, units = "mins"), " minutes.", sep = ""))
}

save(list_res, file = "AJRResult.RData")

load(file = "AJRResult.RData")


# results for the homogeneous treatment effect

for(ml in ml_vec){
  print(paste("Coefficient (standard deviation) for ", ml, sep = ""))
  print("---------------------------------------------")
  for(iv in iv_method){
    print(paste("IV method ", iv, ": ", round(coef(list_res[[ml]], iv), 2), " (", round(se(list_res[[ml]], iv), 2), ")", sep = ""))
  }
  print("_____________________________________________")
}

# standard confidence intervals
for(ml in ml_vec){
  print(paste("Standard confidence intervals for ", ml, sep = ""))
  print("---------------------------------------------")
  for(iv in iv_method){
    print(paste("IV method ", iv, sep = ""))
    print(paste("Standard Confidence Interval: ", "[", paste(round(standard_confint(list_res[[ml]], iv)$CI, 2), collapse = ", "), "]", sep = ""))
  }
  print("_____________________________________________")
}

# robust confidence intervals
for(ml in ml_vec){
  print(paste("Robust confidence intervals for ", ml, sep = ""))
  print("---------------------------------------------")
  for(iv in iv_method){
    print(paste("IV method ", iv, sep = ""))
    print(robust_confint(list_res[[ml]], iv, CI_range = c(-50, 50)))
  }
  print("_____________________________________________")
}


pdf("Plots/AJR.pdf", width = 8, height = 6)
par(mfrow = c(2,2))
par(mgp = c(2, 1, 0))
par(bg = "white")
aa <- seq(0, 0.7, length.out = 80)
iv_vec <- c("linearIV", "mlIV")
col_vec <- c("red", "blue")
transparent_col_vec <- c(rgb(1, 0, 0, 0.1), rgb(0, 0, 1, 0.1))

h_normal <- bandwidth_normal(A)
h_undersmooth <- h_normal * length(A)^0.2 / length(A)^(2/7)

for(ml_method in c("gam", "xgboost")){
  dml_result <- list_res[[ml_method]]
  for(h in c(h_normal, h_undersmooth)){
    plot_title <- paste("ML method: ", ml_method, "\n h = ", round(h, 2), sep = "")
    plot(1, 1, type = "n", xlim = c(0, 0.7), ylim = c(-5, 5), 
         xlab = "latitude", ylab = "beta", main = plot_title, cex.lab = 1.2)
    for(i in 1:length(iv_vec)){
      iv <- iv_vec[i]
      beta_hat <- Vectorize(function(a){
        return(coef(dml_result, iv, a, dml_result$A, "epanechnikov", h))
      })
      se_hat <- Vectorize(function(a){
        return(se(dml_result, iv, a, dml_result$A, "epanechnikov", h))
      })
      lines(aa, beta_hat(aa), col = col_vec[i], type = "l")
      za <- qnorm(0.975)
      lines(aa, beta_hat(aa) + za * se_hat(aa), col = col_vec[i], lty = 3)
      lines(aa, beta_hat(aa) - za * se_hat(aa), col = col_vec[i], lty = 3)
      yy <- seq(-5.2, 5.2, length.out = 60)
      ay <- expand.grid(aa = aa, yy = yy)
      shaded <- function(a, y){
        return(robust_p_value_aggregated(dml_result, y, iv, a, dml_result$A, "epanechnikov", h) >= 0.05)
      }
      in_int <- mapply(shaded, ay$aa, ay$yy)
      in_int <- ay[in_int, ]
      points(in_int$aa, in_int$yy, col = transparent_col_vec[i], pch = 19, cex = 0.3)
    }
    rug(dml_result$A)
    legend("topright", legend = iv_vec, col = col_vec, lty = 1)
  }
}
dev.off()

