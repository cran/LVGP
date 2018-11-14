#' @title The Fitting Function of \code{LVGP} Package
#'
#' @description Fits a latent-variable Gaussian process (LVGP) model to a dataset as described in \code{reference 1}.
#'     The input variables can be quantitative or qualitative/categorical or mixed.
#'     The output variable is quantitative and scalar.
#'
#' @param X Matrix or data frame containing the inputs of training data points. Each row is a data point.
#' @param Y Vector containing the outputs of training data points
#' @param ind_qual Vector containing the indices of columns of qualitative/categorical variables
#' @param dim_z Dimensionality of latent space, usually 1 or 2 but can be higher
#' @param eps The vector of smallest eigen values that the correlation matrix is allowed to have, which determines the nugget added to the correlation matrix.
#' @param lb_phi_ini,ub_phi_ini The initial lower and upper search bounds of the scale/roughness parameters (\code{phi}) of quantitative variables
#' @param lb_phi_lat,ub_phi_lat The later lower and upper search bounds of the scale/roughness parameters (\code{phi}) of quantitative variables
#' @param lb_z,ub_z The lower and upper search bounds of the latent parameters (\code{z}) of qualitative variables
#' @param n_opt The number of times the log-likelihood function is optimized
#' @param max_iter_ini The maximum number of iterations for each optimization run for largest (first) eps case
#' @param max_iter_lat The maximum number of iterations for each optimization run for after first eps cases
#' @param seed An integer for the random number generator. Use this to make the results reproducible.
#' @param progress The switch determining whether to print function run details
#' @param parallel The switch determining whether to use parallel computing
#' @param noise The switch for whether the data are assumed noisy
#'
#' @import lhs
#' @import randtoolbox
#' @import parallel
#'
#' @return A model of class "LVGP model" list of the following items:
#' \itemize{
#' \item{\code{quant_param}} {A list containing the estimated parameter \code{phi} and its search bounds for quantitative variables}
#' \item{\code{qual_param}} {A list containing the estimated parameter \code{z} and its dimensionality, vectorized form and search bounds for qualitative variables}
#' \item{\code{data}} {A list containing the fitted dataset in verbose format}
#' \item{\code{fit_detail}} {A list of more detailed variables for fitting and prediction process}
#' \item{\code{optim_hist}} {Optimization history}
#' \item{\code{setting}} {Settings for the optimization and fitting process}
#' }
#'
#' @references
#' \enumerate{
#' \item "A Latent Variable Approach to Gaussian Process Modeling with Qualitative and Quantitative Factors", Yichi Zhang, Siyu Tao, Wei Chen, and Daniel W. Apley (\href{https://arxiv.org/abs/1806.07504}{arXiv})
#' }
#'
#' @export
#'
#' @seealso
#' \code{\link[stats]{optim}} for the details on \code{L-BFGS-B} algorithm used in optimization.\cr
#' \code{\link[LVGP]{LVGP_predict}} to use the fitted LVGP model for prediction.\cr
#' \code{\link[LVGP]{LVGP_plot}} to plot the features of the fitted model.
#'
#' @examples
#' # Math example with 2 quantitative and 1 qualitative variables (dataset included in the package):
#' #     Fit a model (with default settings) and evaluate the performance
#' #     by computing the root mean squared error (RMSE) in prediction.
#' #     Also, plot the latent variable parameters.
#' X_tr <- math_example$X_tr
#' Y_tr <- math_example$Y_tr
#' X_te <- math_example$X_te
#' Y_te <- math_example$Y_te
#' n_te <- nrow(X_te)
#' model <- LVGP_fit(X_tr, Y_tr, ind_qual = c(3))
#' output <- LVGP_predict(X_te, model)
#' Y_hat <- output$Y_hat
#' RRMSE <- sqrt(sum((Y_hat-Y_te)^2)/n_te)/(max(Y_te)-min(Y_te))
#' LVGP_plot(model)
#'

LVGP_fit <- function(X, Y, ind_qual = NULL, dim_z = 2, eps = 10^(seq(-1, -8, length.out = 15)),
                     lb_phi_ini = -2, ub_phi_ini = 2, lb_phi_lat = -8, ub_phi_lat = 3, lb_z = -3, ub_z = 3, n_opt = 8,
                     max_iter_ini = 100, max_iter_lat = 20, seed = 123, progress = FALSE,
                     parallel = FALSE, noise = FALSE) {

  set.seed(seed)

  ## function input checking & preprocessing
  if (progress){
    cat('** Checking & preprocessing the inputs...\n')
  }
  if (missing(X)){
    stop('    X must be provided.')
  }
  if (!is.matrix(X) && !is.data.frame(X)){
    stop('    X must be a matrix or a data frame.')
  }
  k <- nrow(X) # number of data points
  p_all <- ncol(X) # total number of input variables
  # check qualitative and quantitative parts of dataset
  if (is.null(ind_qual)) { # no qualitative/categorical variables
    p_qual <- as.integer(0)
    X_qual <- NULL
    X_quant <- X # data of quant variables
    if (!all(is.finite(X_quant)) || !is.numeric(X_quant)){
      stop('    All the elements of X_quant must be finite numbers.')
    }
    if (is.matrix(X_quant) == FALSE) {
      X_quant <- as.matrix(X_quant)
    }
    lvs_qual <- NULL
    n_lvs_qual <- NULL
    n_z <- 0
  } else {
    if (!all(is.finite(ind_qual)) || !is.numeric(ind_qual)){
      stop('    All the elements of ind_qual must be finite numbers.')
    }
    ind_qual <- as.vector(as.integer(ind_qual)) # indices of qual variables in X mat
    p_qual <- length(ind_qual) # number of qual variables
    X_qual <- X[, ind_qual] # data of qual variables
    if (is.data.frame(X_qual) == FALSE) {
      X_qual <- as.data.frame(X_qual)
    }
    if (p_qual == p_all) { # all the variables are qualitative/categorical
      X_quant <- NULL
    } else {
      X_quant <- X[, -ind_qual]
      if (!all(is.finite(X_quant)) || !is.numeric(X_quant)){
        stop('    All the elements of X_quant must be finite numbers.')
      }
      if (is.matrix(X_quant) == FALSE) {
        X_quant <- as.matrix(X_quant)
      }
    }

    lvs_qual <- list() # levels of each qual variable
    n_lvs_qual <- rep(0, p_qual) # number of levels of each qual variable
    for (i in 1:p_qual) {
      lvs_qual[[i]] <- sort(unique(X_qual[, i]))
      n_lvs_qual[i] <- length(lvs_qual[[i]])
    }
    n_z <- dim_z*(sum(n_lvs_qual)-p_qual) # number of latent parameters
  }


  if (missing(Y)){
    stop('    Y must be provided.')
  }
  if (!all(is.finite(Y)) || !is.numeric(Y)){
    stop('    All the elements of Y must be finite numbers.')
  }
  if (is.matrix(Y) == FALSE) {
    Y <- as.matrix(Y)
  }
  if (k != nrow(Y)){
    stop('    The number of rows (i.e., observations) in X and Y should match!')
  }

  dim_z <- as.integer(dim_z)
  if (dim_z != 1 && dim_z != 2) {
    warning('    The dimensionality of latent space is uncommon!')
  }

  eps <- as.numeric(eps)
  if (any(eps < (10^(round(log10(.Machine$double.eps))+3)))){
    stop(paste('    Increase the smallest member of eps. The minimum allowable is ', toString(10^(round(log10(.Machine$double.eps))+3))))
  }
  if (any(diff(eps) > 0)){
    stop('The elements of eps should be in a descending order.')
  }

  lb_phi_ini <- as.numeric(lb_phi_ini)
  ub_phi_ini <- as.numeric(ub_phi_ini)
  lb_phi_lat <- as.numeric(lb_phi_lat)
  ub_phi_lat <- as.numeric(ub_phi_lat)
  lb_z <- as.numeric(lb_z)
  ub_z <- as.numeric(ub_z)
  n_opt <- as.integer(n_opt) # number of starting points for optimization
  max_iter_ini <- as.integer(max_iter_ini)
  max_iter_lat <- as.integer(max_iter_lat)

  p_quant <- p_all - p_qual

  ## normalization of X_quant and Y
  if (p_quant > 0) {
    X_quant_min <- apply(X_quant, 2, min)
    X_quant_max <- apply(X_quant, 2, max)
    X_quant <- t((t(X_quant)-X_quant_min)/(X_quant_max-X_quant_min))
  } else {
    X_quant_min <- NULL
    X_quant_max <- NULL
  }

  Y_min <- apply(Y, 2, min)
  Y_max <- apply(Y, 2, max)
  Y <- t((t(Y)-Y_min)/(Y_max-Y_min))

  ## initiation for optimization
  n_hyper <- p_quant+n_z
  lb_ini <- c(rep(lb_phi_ini, p_quant), rep(lb_z, n_z))
  ub_ini <- c(rep(ub_phi_ini, p_quant), rep(ub_z, n_z))
  lb_lat <- c(rep(lb_phi_lat, p_quant), rep(lb_z, n_z))
  ub_lat <- c(rep(ub_phi_lat, p_quant), rep(ub_z, n_z))
  if (dim_z == 2 && p_qual != 0) {
    temp_ind <- p_quant
    for (i in 1:p_qual) {
      n_lvs <- n_lvs_qual[i]
      lb_ini[temp_ind+dim_z] <- -1e-4
      ub_ini[temp_ind+dim_z] <- 1e-4
      lb_lat[temp_ind+dim_z] <- -1e-4
      ub_lat[temp_ind+dim_z] <- 1e-4
      temp_ind <- temp_ind + dim_z*(n_lvs-1)
    }
  }

  opt_ctrl_ini <- c(trace=progress, maxit = max_iter_ini,  REPORT = 10, lmm = 100)
  opt_ctrl_lat <- c(trace=progress, maxit = max_iter_lat,  REPORT = 10, lmm = 100)

  setting <- list(max_iter_ini = max_iter_ini, max_iter_lat = max_iter_lat, seed = seed, n_opt = n_opt,
                  lb_phi_ini = lb_phi_ini, ub_phi_ini = ub_phi_ini,
                  lb_phi_lat = lb_phi_lat, ub_phi_lat = ub_phi_lat,
                  lb_z = lb_z, ub_z = ub_z, trace=progress,
                  parallel = parallel, eps = eps)

  start_time <- proc.time()[3]

  A <- lhs::maximinLHS(n_opt, n_hyper)
  hyper0 <- t(t(A)*(ub_ini - lb_ini) + lb_ini) # starting points (of hyperparameters) for optimization

  M <- matrix(1, nrow = k, ncol = 1)

  ## optimization runs
  # prepare parallel start
  if (parallel) {
    n_cores <- detectCores()
    optimfun <- function(x0_ele, objfun, addinput) {
      p_quant <- addinput$p_quant
      p_qual <- addinput$p_qual
      lvs_qual <- addinput$lvs_qual
      n_lvs_qual <- addinput$n_lvs_qual
      dim_z <- addinput$dim_z
      X_quant <- addinput$X_quant
      X_qual <- addinput$X_qual
      Y <- addinput$Y
      eps_i <- addinput$eps_i
      k <- addinput$k
      M <- addinput$M
      lb <- addinput$lb
      ub <- addinput$ub
      opt_ctrl <- addinput$opt_ctrl
      x_sol_ele <- stats::optim(x0_ele, objfun, gr = NULL, p_quant, p_qual, lvs_qual, n_lvs_qual,
                                dim_z, X_quant, X_qual, Y, eps_i, k, M,
                                method = 'L-BFGS-B', lower = lb, upper = ub, control = opt_ctrl)
      return(x_sol_ele)
    }
    addinput <- list(p_quant = p_quant, p_qual = p_qual, lvs_qual = lvs_qual, n_lvs_qual = n_lvs_qual,
                     dim_z = dim_z, X_quant = X_quant, X_qual = X_qual, Y = Y, eps_i = NULL,
                     k = k, M = M, lb = lb_ini, ub = ub_ini, opt_ctrl = opt_ctrl_ini)
  }
  # prepare parallel end
  n_try <- length(eps)
  optim_hist <- list()
  #optim_hist$hyper0[[1]] <- hyper0
  for (i_try in 1:n_try) {
    eps_i <- eps[i_try]
    n_opt_i <- nrow(hyper0)
    hyper_sol <- matrix(0, n_opt_i, n_hyper)
    obj_sol <- matrix(0, n_opt_i, 1)
    flag_sol <- matrix(0, n_opt_i, 1)
    # note: could shrink the search space at first step, then expand
    if (parallel) {
      if (i_try == 2) {
        addinput$opt_ctrl <- opt_ctrl_lat
        addinput$lb <- lb_lat
        addinput$ub <- ub_lat
      }
      clust <- makeCluster(min(n_cores,n_opt_i))
      hyper0_list <- list()
      for (j in 1: n_opt_i) {
        hyper0_list[[j]] <- hyper0[j, ]
      }
      addinput$eps_i <- eps_i
      temp_list <- parLapply(clust, hyper0_list, optimfun, neg_log_l, addinput)
      stopCluster(clust)
      for (j in 1: n_opt_i){
        hyper_sol[j, ] <- temp_list[[j]]$par
        obj_sol[j] <- temp_list[[j]]$value
        flag_sol[j] <- temp_list[[j]]$convergence
      }
    } else {
      for (j in 1: n_opt_i){
        if (i_try == 1) {
          temp <- stats::optim(hyper0[j, ], neg_log_l, gr = NULL, p_quant, p_qual, lvs_qual, n_lvs_qual,
                               dim_z, X_quant, X_qual, Y, eps_i, k, M,
                               method = 'L-BFGS-B', lower = lb_ini, upper = ub_ini, control = opt_ctrl_ini)
        } else {
          temp <- stats::optim(hyper0[j, ], neg_log_l, gr = NULL, p_quant, p_qual, lvs_qual, n_lvs_qual,
                               dim_z, X_quant, X_qual, Y, eps_i, k, M,
                               method = 'L-BFGS-B', lower = lb_lat, upper = ub_lat, control = opt_ctrl_lat)
        }
        hyper_sol[j, ] <- temp$par
        obj_sol[j] <- temp$value
        flag_sol[j] <- temp$convergence
      }
    }
    ID = sort(obj_sol, index.return=TRUE)$ix
    obj_sol <- as.matrix(obj_sol[ID, drop = FALSE])
    flag_sol <- as.matrix(flag_sol[ID, drop = FALSE])
    hyper_sol <- as.matrix(hyper_sol[ID, , drop = FALSE])
    hyper0 <- as.matrix(hyper0[ID, , drop = FALSE])

    rownames(obj_sol) <- 1:n_opt_i
    obj_sol_unique  <-  as.matrix(unique(round(obj_sol, 2)))
    flag_sol_unique <- as.matrix(as.data.frame(flag_sol)[c(rownames(obj_sol_unique)), ])
    hyper_sol_unique <- as.matrix(as.data.frame(hyper_sol)[c(rownames(obj_sol_unique)), ])

    #optim_hist$nug_opt[[i]] <- nug_opt
    #optim_hist$raw_min_eig[[i]] <-  raw_min_eig
    optim_hist$hyper0[[i_try]] <- hyper0
    optim_hist$hyper_sol[[i_try]] <-  hyper_sol
    optim_hist$hyper_sol_unique[[i_try]] <-  hyper_sol_unique
    optim_hist$obj_sol[[i_try]] <-  obj_sol
    optim_hist$obj_sol_unique[[i_try]] <- obj_sol_unique
    optim_hist$obj_sol_best[[i_try]] <-  obj_sol[1]
    optim_hist$flag_sol_unique[[i_try]] <- flag_sol_unique

    if (i_try < n_try) {
      hyper0 <- hyper_sol_unique
    }

  }

  fit_time <- proc.time()[3] - start_time


  ## post-processing
  if (!noise) {
    id_best_try <- n_try
  } else {
    id_best_try <- which.min(optim_hist$obj_sol_best)
  }
  hyper_full <- as.matrix(optim_hist$hyper_sol[[id_best_try]][1, ])
  min_n_log_l <- min(optim_hist$obj_sol_best[[id_best_try]])

  if (p_quant == 0) {
    phi <- NULL
  } else {
    phi <- hyper_full[1:p_quant]
  }
  if (p_qual == 0) {
    z_vec <- NULL
    z <- NULL
  } else {
    z_vec <- hyper_full[-c(1:p_quant)]
    z <- list() # z is a list of matrices corresponding to each qualitative variable
    ind_temp <- as.integer(0)
    for (i in 1:p_qual) {
      n_lvs <- n_lvs_qual[i]
      z_i <- z_vec[(dim_z*(ind_temp)+1) : (dim_z*(ind_temp+n_lvs-1))]
      ind_temp <- ind_temp+n_lvs-1
      mat_tmp <- matrix(0, nrow = n_lvs, ncol = dim_z)
      mat_tmp[-1, ] <- matrix(z_i, nrow = n_lvs-1, byrow = TRUE)
      z[[i]] <- mat_tmp
    }
  }


  # calc convenient quantities (stored in model$fit_detail)
  if (p_qual == 0) { # no qualitative/categorical variables
    R <- corr_mat(X_quant, X_quant, phi)
    X_full <- X_quant
  } else {
    X_qual_trans <- to_latent(X_qual, lvs_qual, n_lvs_qual, p_qual, z_vec, dim_z, k)
    X_full <- cbind(X_quant, X_qual_trans)
    omega <- c(phi, rep(0, p_qual*dim_z))
    R <- corr_mat(X_full, X_full, omega)
  }
  R <- (R + t(R))/2

  raw_min_eig = min(eigen(R, symmetric = TRUE, only.values = TRUE)$values)
  if (raw_min_eig < eps[id_best_try]) {
    nug_opt <- eps[id_best_try] - raw_min_eig
    R <- R + diag(x = 1, k, k)*nug_opt
  } else {
    nug_opt <- 0
  }

  L <- t(chol(R)) # R = L%*%t(L)
  Linv <- solve(L)
  LinvM <- solve(L, M)
  MTRinvM <- sum(LinvM^2)
  beta_hat <- t(LinvM)%*%solve(L, Y)/MTRinvM
  beta_hat <- as.numeric(beta_hat)
  temp <- solve(L, Y - M*beta_hat)
  sigma2 <- sum(temp^2)/k
  if ( sigma2 < 1e-300) {
    sigma2 <- 1e-300
  }
  RinvPYminusMbetaP <- solve(t(L), temp)

  ## Save the fitted model
  model <- NULL
  model$quant_param <- list('phi' = phi, 'lb_phi_ini' = lb_phi_ini, 'ub_phi_ini' = ub_phi_ini,
                            'lb_phi_lat' = lb_phi_lat, 'ub_phi_lat' = ub_phi_lat)
  model$qual_param <-  list('dim_z' = dim_z, 'z' = z, 'z_vec' = z_vec,
                            'lb_z' = lb_z, 'ub_z' = ub_z)
  model$data <- list('X' = X, 'X_quant' = X_quant, 'X_qual' = X_qual, 'X_full' = X_full, 'Y' = Y,
                     'X_quant_min' = X_quant_min, 'X_quant_max' = X_quant_max,
                     'Y_min' = Y_min, 'Y_max' = Y_max,
                     'ind_qual' = ind_qual, 'lvs_qual' = lvs_qual, 'n_lvs_qual' = n_lvs_qual,
                     'p_all' = p_all, 'p_quant' = p_quant, 'p_qual' = p_qual)
  model$fit_detail <- list('beta_hat' = beta_hat, 'sigma2'=sigma2, 'MTRinvM' = MTRinvM,
                           'Linv' = Linv, 'LinvM' = LinvM, 'RinvPYminusMbetaP'=RinvPYminusMbetaP,
                           'raw_min_eig' = raw_min_eig,
                           'nug_opt' = nug_opt, 'min_n_log_l' = min_n_log_l, 'fit_time' = fit_time)
  model$optim_hist <- optim_hist
  model$setting <- setting
  class(model) <- 'LVGP model'

  return(model)
}
