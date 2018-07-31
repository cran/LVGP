#' @title The Prediction Function of \code{LVGP} Package
#' @description Predicts the output and associated uncertainties of the GP model fitted by \code{\link[LVGP]{LVGP_fit}}.
#'
#' @param X_new Matrix or vector containing the input(s) where the predictions are to be made. Each row is an input vector.
#' @param model The LVGP model fitted by \code{\link[LVGP]{LVGP_fit}}.
#' @param MSE_on A scalar indicating whether the uncertainty (i.e., mean squared error \code{MSE}) is calculated.
#'     Set to a non-zero value to calculate \code{MSE}.
#'
#' @return A prediction list containing the following components:
#' \itemize{
#' \item{\code{Y_hat}} {A vector containing the mean prediction values}
#' \item{\code{MSE}} {A vector containing the prediction uncertainty (i.e., the covariance or covariance matrix for the output(s) at prediction location(s)) }
#' }
#'
#' @export
#'
#' @references
#' \enumerate{
#' \item "A Latent Variable Approach to Gaussian Process Modeling with Qualitative and Quantitative Factors", Yichi Zhang, Siyu Tao, Wei Chen, and Daniel W. Apley (\href{https://arxiv.org/abs/1806.07504}{arXiv})
#' }
#'
#' @seealso
#' \code{\link[LVGP]{LVGP_fit}} to fit LVGP model to the datasets.\cr
#' \code{\link[LVGP]{LVGP_plot}} to plot the features of the fitted model.
#'
#' @examples
#' # see the examples in the documentation of the function LVGP_fit.

LVGP_predict <-  function(X_new, model, MSE_on = 0){

  ## import model and check inputs
  if (class(model) != "LVGP model"){
    stop('    The 2nd input should be a model of class "LVGP model".')
  }
  if (length(MSE_on)!=1){
    stop('    MSE_on should be a scalar flag. Set it to 1 to turn it "on".')
  }
  if (is.list(X_new)) {
    stop('    The package do not yet support "list" type for X_new.')
  }
  p_all <- model$data$p_all
  if (is.null(nrow(X_new))) {
    if (length(X_new) != p_all){
      stop('    The dimension of X_new is not correct!')
    }
  } else {
    if (ncol(X_new) != p_all){
      stop('    The dimension of X_new is not correct!')
    }
  }
  p_qual <- model$data$p_qual
  p_quant <- model$data$p_quant

  X_quant_min <- model$data$X_quant_min
  X_quant_max <- model$data$X_quant_max
  Y_min <- model$data$Y_min
  Y_max <- model$data$Y_max

  X_old_full <- model$data$X_full
  lvs_qual <- model$data$lvs_qual
  n_lvs_qual <- model$data$n_lvs_qual
  ind_qual <- model$data$ind_qual

  phi <- model$quant_param$phi
  dim_z <- model$qual_param$dim_z
  z_vec <- model$qual_param$z_vec

  beta_hat <- as.numeric(model$fit_detail$beta_hat)
  RinvPYminusMbetaP <- model$fit_detail$RinvPYminusMbetaP

  ## process X_new
  if (is.null(nrow(X_new))) {
    m <- 1
    X_new <- matrix(X_new, 1, p_all)
    X_new_2 <- as.data.frame(rbind(X_new, X_new)) # convert to dataframe
    X_new <- X_new_2[1,]
  } else {
    m <- nrow(X_new)
  }

  if (p_qual == 0) { # no qualitative/categorical variables
    X_new_qual <- NULL
    X_new_quant <- X_new
    if (is.matrix(X_new_quant) == FALSE) {
      X_new_quant <- as.matrix(X_new_quant)
    }
    if (!all(is.finite(X_new_quant)) || !is.numeric(X_new_quant)){
      stop('    All the elements of X_new_quant must be finite numbers.')
    }
    X_new_quant <- t((t(X_new_quant)-X_quant_min)/(X_quant_max-X_quant_min))

    R_old_new <- corr_mat(X_old_full, X_new_quant, phi)
    R_new_new <- corr_mat(X_new_quant, X_new_quant, phi)
  } else {
    X_new_qual <- X_new[, ind_qual]
    if (is.data.frame(X_new_qual) == FALSE) {
      X_new_qual <- as.data.frame(X_new_qual)
    }

    if (p_qual == p_all) { # all the variables are qualitative/categorical
      X_new_quant <- NULL
    } else {
      X_new_quant <- X_new[, -ind_qual] # data of quant variables
      if (is.matrix(X_new_quant) == FALSE) {
        X_new_quant <- as.matrix(X_new_quant)
      }
      if (!all(is.finite(X_new_quant)) || !is.numeric(X_new_quant)){
        stop('    All the elements of X_new_quant must be finite numbers.')
      }
      X_new_quant <- t((t(X_new_quant)-X_quant_min)/(X_quant_max-X_quant_min))
    }

    X_new_qual_la <- to_latent(X_new_qual, lvs_qual, n_lvs_qual, p_qual, z_vec, dim_z, m)
    X_new_full <- cbind(X_new_quant, X_new_qual_la)
    phi_full <- c(phi, rep(0, p_qual*dim_z))
    R_old_new <- corr_mat(X_old_full, X_new_full, phi_full)
    R_new_new <- corr_mat(X_new_full, X_new_full, phi_full)
  }
  R_new_new <- (R_new_new + t(R_new_new))/2

  ## calc predictions
  Y_hat <- beta_hat+t(R_old_new)%*%RinvPYminusMbetaP
  Y_hat = t(t(Y_hat)*(Y_max-Y_min) + Y_min)

  prediction <- list()
  prediction$Y_hat <- Y_hat

  if (MSE_on) {
    # import relevant quantities
    sigma2 <- model$fit_detail$sigma2
    LTinv <- model$fit_detail$LTinv
    MTRinvM <- model$fit_detail$MTRinvM
    MTLinv <- model$fit_detail$MTLinv

    temp <- LTinv%*%R_old_new
    W <- 1-MTLinv%*%temp
    MSE <- sigma2*(R_new_new-t(temp)%*%temp+t(W)%*%W/MTRinvM)*((Y_max-Y_min)^2)
    prediction$MSE <- MSE
  }

  return(prediction)
}
