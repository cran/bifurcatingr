#' Bias-Corrected Confidence intervals of Least Squares Estimators for Bifurcating
#' Autoregressive Models
#'
#' This function performs bias correction confidence intervals on the least squares
#' estimators of the autoregressive coefficients in a BAR(p) model using single,
#' fast-double, and the Bias-corrected and accelerated bootstrapping as described
#' in Elbayoumi and Mostafa (2023).
#'
#' @param z a numeric vector containing the tree data
#' @param p an integer determining the order of bifurcating autoregressive model
#'   to be fit to the data
#' @param method method of bias correction. Currently, "boot1", "boot2fast", and
#'  "BCa" are supported and they implement single bootstrap, fast-double bootstrap,
#'   and bias-corrected and accelerated bootstrap, respectively. Defaults to "BCa".
#' @param conf_int type of the confidence interval. Currently, "standard_normal_bc",
#'  "percentile", and "percentile_bc" are supported and they implement corrected
#'  standard normal bootstrap CI, uncorrected percentile bootstrap CI, and corrected
#'  percentile bootstrap CI, respectively. If "boot1" method is selected, the
#'  "standard_normal_bc", "percentile", "percentile_bc" confidence intervals can
#'  be obtained. If "boot2fast" method is selected, the "standard_normal_bc" and
#'  "percentile_bc" confidence intervals can be obtained. No effect for conf_int,
#'  the "BCa" method is selected. Defaults to standard_normal_bc".
#' @param conf_level confidence level to be used in computing confidence intervals
#'   for model coefficients. Defaults to \code{0.95}.
#' @param B number of bootstrap samples (replicates).
#' @param burn number of tree generations to discard before starting the
#'  bootstrap sample (replicate). Defaults to 5.
#' @return \item{Bias_corrected_coef}{a matrix containing the
#'  bias-correction least squares estimates of the autoregressive coefficients}
#'  \item{BCa_ci}{a matrix containing the lower and upper limits of corrected BCa
#'  confidence intervals,if \code{method="BCa"}}
#'  \item{standard_normal_bc_ci}{a matrix containing the lower and upper limits of
#'  corrected confidence intervals, if \code{method="boot1"} and \code{conf_int="standard_normal_bc"}
#'  or \code{conf_int="All"}}
#'  \item{percentile_ci}{a matrix containing the lower and upper limits of uncorrected
#'  percentile confidence intervals, if \code{method="boot1"} and \code{conf_int="percentile"}
#'  or \code{conf_int="All"}}
#'  \item{percentile_bc_ci}{a matrix containing the lower and upper limits of corrected
#'  percentile confidence intervals, if \code{method="boot1"} and \code{conf_int="percentile_bc"}
#'  or \code{conf_int="All"}}
#'  \item{standard_normal_bc_ci}{a matrix containing the lower and upper limits of corrected
#'  confidence intervals, if \code{method="boot2fast"} and \code{conf_int="standard_normal_bc"}
#'  or \code{conf_int="All"}}
#'  \item{percentile_bc_ci}{a matrix containing the lower and upper limits of corrected percentile
#'  confidence intervals, if \code{method="boot2fast"} and \code{conf_int="percentile_bc"} or \code{conf_int="All"}}
#' @references Elbayoumi, T. M. & Mostafa, S. A. (2023). Impact of Bias Correction
#'  of the Least Squares Estimation on Bootstrap Confidence Intervals for Bifurcating
#'  Autoregressive Models. \emph{Journal of Data Science}, 1-20, doi.org/10.6339/23-JDS1092.
#' @export
#' @examples
#' # Generating Non-contaminated normal BAR(1) tree and calculating the bias corrected
#' # standard normal CI for the autoregressive coefficients of the BAR(1) model
#' # Note that in this example (B=2) for speeding up the calculations.
#' # B must be set to 499 or more for calculation accuracy.
#' z <- bfa_tree_gen(15, 1, 1, 1, -0.9, -0.9, 0, 10, c(-0.5))
#' bfa_ls_bc_ci(z, p=1, method="boot1", B=2)
bfa_ls_bc_ci <- function(z, p, method="BCa", conf_int="standard_normal_bc", conf_level=0.95, B=5, burn=5){
  if (method=="boot1"){
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    boot_ci<- matrix(NA, nrow=p,ncol=2)
    perc_boot1_ci<- matrix(NA, nrow=p,ncol=2)
    percbc_ci<- matrix(NA, nrow=p,ncol=2)
    a1_ls_boot1_star <- bfa_boot_ls_bc(z, p, method=method, burn = burn, B=B, boot_est=TRUE)$boot_bcest
    coef_ls_bc <- 2*coef_ls$coef[-1] - colMeans(a1_ls_boot1_star)
    if (conf_int=="All"){
      a1_ls_star <- bfa_boot1_ls(z, p, burn = burn, B=B, boot_est = TRUE, boot_data = FALSE)
    for (i in 1:p) {
      boot_ci[i,] <- bfa_boot_ci(coef_ls_bc[i], a1_ls_boot1_star[,i], conf_level)
      perc_boot1_ci[i,] <- bfa_perc_ci(a1_ls_boot1_star[,i], conf_level)
      percbc_ci[i,] <- bfa_perc_bc_ci(z, coef_ls$coef[1,(i+1)], a1_ls_star$boot_est[,(i+1)], conf_level)
    }
    out <- list(Bias_corrected_coef = as.matrix(coef_ls_bc), standard_normal_bc_ci=boot_ci,percentile_ci=perc_boot1_ci, percentile_bc_ci=percbc_ci)
    colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
    rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$standard_normal_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$standard_normal_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$percentile_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$percentile_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$percentile_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$percentile_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf_int=="standard_normal_bc"){
      for (i in 1:p) {
        boot_ci[i,] <- bfa_boot_ci(coef_ls_bc[i], a1_ls_boot1_star[,i], conf_level)
      }
      out <- list(Bias_corrected_coef = as.matrix(coef_ls_bc),standard_normal_ci=boot_ci)
      colnames(out$standard_normal_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$standard_normal_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf_int=="percentile"){
      for (i in 1:p) {
        perc_boot1_ci[i,] <- bfa_perc_ci(a1_ls_boot1_star[,i], conf_level)
      }
      out <- list(Bias_corrected_coef = as.matrix(coef_ls_bc),percentile_ci=perc_boot1_ci)
      colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
      rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$percentile_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf_int=="percentile_bc"){
      a1_ls_star <- bfa_boot1_ls(z, p, burn = burn, B=B, boot_est = TRUE, boot_data = FALSE)
      for (i in 1:p) {
        percbc_ci[i,] <- bfa_perc_bc_ci(z, coef_ls$coef[1,(i+1)], a1_ls_star$boot_est[,(i+1)], conf_level)
      }
      out <- list(Bias_corrected_coef = as.matrix(coef_ls_bc), percentile_bc_ci=percbc_ci)
      colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
      rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$percentile_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
  }
  if (method=="boot2fast"){
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    boot2fast_ci<- matrix(NA, nrow=p,ncol=2)
    perc_boot2fast_ci<- matrix(NA, nrow=p,ncol=2)
    a1_ls_boot2fast_star <- bfa_boot_ls_bc(z, p, method=method, burn = burn, B=B, boot_est=TRUE)$boot_bcest
    bt2 <- bfa_boot2fast_ls(z, p, burn = burn, B=B)
    bias1 <- colMeans(bt2[[1]])[-1]-coef_ls$coef[-1]
    bt2_1 <- matrix(unlist(lapply(bt2[[2]], FUN=unlist)), nrow=B, byrow=TRUE)
    bias2 <- bias1 - colMeans(bt2_1)[-1] + colMeans(bt2[[1]])[-1]
    coef_ls_bc <- coef_ls$coef[-1] - bias1 - bias2
    if (conf_int=="All"){
      a1_ls_star <- bfa_boot2_ls(z, p, burn = burn, B1=B, B2=1)
      for (i in 1:p) {
        boot2fast_ci[i,] <- bfa_boot_ci(coef_ls_bc[i], a1_ls_boot2fast_star[,i], conf_level)
        perc_boot2fast_ci[i,] <- bfa_perc_ci(a1_ls_boot2fast_star[,i], conf_level)
      }
      out <- list(Bias_corrected_coef= as.matrix(coef_ls_bc),standard_normal_bc_ci=boot2fast_ci,percentile_bc_ci=perc_boot2fast_ci)
      colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
      rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$standard_normal_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$standard_normal_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$percentile_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf_int=="standard_normal_bc"){
      for (i in 1:p) {
        boot2fast_ci[i,] <- bfa_boot_ci(coef_ls_bc[i], a1_ls_boot2fast_star[,i], conf_level)
      }
      out <- list(Bias_corrected_coef=as.matrix(coef_ls_bc), standard_normal_bc_ci=boot2fast_ci)
      colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
      rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$standard_normal_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$standard_normal_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
    if (conf_int=="percentile_bc"){
      for (i in 1:p) {
        perc_boot2fast_ci[i,] <- bfa_perc_ci(a1_ls_boot2fast_star[,i], conf_level)
      }
      out <- list(Bias_corrected_coef=as.matrix(coef_ls_bc),percentile_bc_ci=perc_boot2fast_ci)
      colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
      rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
      colnames(out$percentile_bc_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
      rownames(out$percentile_bc_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
    }
  }
  if (method=="BCa"){
    coef_ls <- bfa_ls(z, p, error_cor = FALSE)
    a1_ls_boot1_star <- bfa_boot_ls_bc(z, p, burn = burn, B=B, boot_est=TRUE)$boot_bcest
    coef_ls_bc <- 2*coef_ls$coef[-1] - colMeans(a1_ls_boot1_star)
    bca_ci <- boot_bca_ci(z,p,B,J=4,conf_level)$BCa
    out <- list(Bias_corrected_coef=as.matrix(coef_ls_bc),BCa_ci=bca_ci)
    colnames(out$Bias_corrected_coef) <- c(paste0("Estimates"))
    rownames(out$Bias_corrected_coef) <- c(paste0('X_[t/',2^c(1:p),"]"))
    colnames(out$BCa_ci) <- c(paste0(((1-conf_level)/2)*100,"%"), paste0((1-(1-conf_level)/2)*100,"%"))
    rownames(out$BCa_ci) <- c(paste0('X_[t/',2^c(1:p),"]"))
  }
  out
}
