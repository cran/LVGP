#' Dataset for the example in function 'LVGP_fit'
#'
#' Data are sampled from the modified math function based on the first example in the paper listed in code{references}.
#' There are still 2 quantitative and 1 qualitative variables, but the qualitative variable has only 3 levels.
#' For each level, there are 8 training data points and 30 testing data points, all generated with Latin hypercube sampling.
#' In total, there are 24 training data points and 90 testing data points.
#'
#' @docType data
#'
#' @usage data(math_example)
#'
#' @format A named list containing training and test data:
#' \describe{
#'   \item{"X_tr"}{24-by-3 matrix for 24 training data inputs, 3rd column being the qualitative variable}
#'   \item{"Y_tr"}{24-by-1 matrix for 24 training data outputs}
#'   \item{"X_te"}{90-by-3 matrix for 90 testing data inputs, 3rd column being the qualitative variable}
#'   \item{"Y_te"}{90-by-1 matrix for 90 testing data outputs}
#' }
#'
#' @keywords example dataset
#'
#' @references
#' \enumerate{
#' \item "A Latent Variable Approach to Gaussian Process Modeling with Qualitative and Quantitative Factors", Yichi Zhang, Siyu Tao, Wei Chen, and Daniel W. Apley (\href{https://arxiv.org/abs/1806.07504}{arXiv})
#' }
#'
#' @source The dataset can be generated with the code at the end of this description file.
#'
#' @examples
#' data(math_example)
#' X_tr <- math_example$X_tr
#' Y_tr <- math_example$Y_tr
#' X_te <- math_example$X_te
#' Y_te <- math_example$Y_te
"math_example"
#
# n_lv <- 3 # number of levels for the qualitative variable
# lb <- c(0, 0); ub <- c(1, 1)
# n_tr_lv <- 8 # number of training samples for each level
# n_te_lv <- 30 # number of testing samples for each level
# n_tr <- n_lv*n_tr_lv # total number of training samples
# n_te <- n_lv*n_te_lv # total number of training samples
# X_tr <- matrix(0, nrow = n_tr, ncol = 3)
# Y_tr <- matrix(0, nrow = n_tr, ncol = 1)
# X_te <- matrix(0, nrow = n_te, ncol = 3)
# Y_te <- matrix(0, nrow = n_te, ncol = 1)
# var_fctr <- c(0, 12.0, 6.0)
# for (i in 1:n_lv) {
#   A <- lhs::maximinLHS(n_tr_lv, 2)
#   lhs_temp <- t(t(A)*(ub - lb) + lb)
#   X_tr[((i-1)*n_tr_lv+1):(i*n_tr_lv),1:2] <- lhs_temp
#   X_tr[((i-1)*n_tr_lv+1):(i*n_tr_lv),3] <- rep(i, n_tr_lv)
#   Y_tr[((i-1)*n_tr_lv+1):(i*n_tr_lv),] <- sin(2*pi*lhs_temp[,2]-pi) +
#     7*sin(2*pi*lhs_temp[,1]-pi)+var_fctr[i]*sin(2*pi*lhs_temp[,2]-pi)
#
#   A <- lhs::maximinLHS(n_te_lv, 2)
#   lhs_temp <- t(t(A)*(ub - lb) + lb)
#   X_te[((i-1)*n_te_lv+1):(i*n_te_lv),1:2] <- lhs_temp
#   X_te[((i-1)*n_te_lv+1):(i*n_te_lv),3] <- rep(i, n_te_lv)
#   Y_te[((i-1)*n_te_lv+1):(i*n_te_lv),] <- sin(2*pi*lhs_temp[,2]-pi) +
#     7*sin(2*pi*lhs_temp[,1]-pi)+var_fctr[i]*sin(2*pi*lhs_temp[,2]-pi)
# }
#
# model <- LVGP_fit(X_tr, Y_tr, ind_qual = c(3))
# output <- LVGP_predict(X_te, model)
# Y_hat <- output$Y_hat
# RRMSE <- sqrt(sum((Y_hat-Y_te)^2)/n_te)/(max(Y_te)-min(Y_te))
# LVGP_plot(model)
#
# math_example = list(X_tr=X_tr,Y_tr=Y_tr,X_te=X_te,Y_te=Y_te)
