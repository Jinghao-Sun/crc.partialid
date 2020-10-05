#' Profile likelihood confidence interval
#'
#' Use profile likelihood to compute the confidence interval for the population size M
#'
#' @importFrom stats qchisq rpois
#' @param mem_list A data.frame that contains the aggregated capture histories. For a
#' capture-recapture study with k samples, mem_list contains one row per observed capture
#' history followed by its frequency. mem_list has 2^k - 1 rows and k+1 columns, where the first k columns
#' can only have 0 or 1 (1 indicates captured in this sample), and the last column is the frequency.
#' See ./data/PWID_dev.RData for an example.
#' @param restriction_type A character string, either "pairwise" (for restrictions on the pairwise dependence),
#' or "highest" (for restrictions on the highest order interaction term lambda_c in a log-linear model).
#' @param alpha A real number between 0 and 1. Specify the confidence level 1-alpha.
#' @param r_cond A vector. Suppose there are restrictions on omega pairs of samples. r_cond is a vector of length omega.
#' It contains the indices of one end of the pairs.Ignored if restriction_type != "pairwise".
#' @param t_cond A vector. t_cond is a vector of length omega. It contains the indices of the other end of the pairs.
#' Ignored if restriction_type != "pairwise".
#' @param eta_cond A vector of length omega. Lower bound of the pairwise restrictions, corresponding to r_cond, t_cond.
#' Ignored if restriction_type != "pairwise".
#' @param xi_cond A vector. Upper bound of the pairwise restrictions, corresponding to r_cond, t_cond. Ignored if
#' restriction_type != "pairwise".
#' @param gamma A positive number. Restrictions on the  highest order interaction term lambda_c in a
#' log-linear model, such that lambda_c >= -gamma, lambda_c >= -gamma <= gamma. Ignored if restriction_type != "highest".
#' @param lb A positive number. lb times the totals counts of the observed units will define the lower bound
#' the searching algorithm. Vary lb and ub to avoid local minimum.
#' @param ub A positive number. ub times the totals counts of the observed units will define the upper bound
#' the searching algorithm. Vary lb and ub to avoid local minimum.
#' @param tol A positive number. Tolerance/precision of the optimization algorithm when finding the endpoints of the confidence
#' interval.
#' @param verbose A logical. If TRUE, print intermediate information of the computation.
#' @param xtol_rel A real number. Control the tolerance when finding the optimum of log-likelihood.
#' @param maxeval An integer. Control the max iterations when finding the optimum of log-likelihood.
#' @return (1 - alpha) confidence interval for the population size
#' @export
#' @examples

#' CI = ci.pl(mem_list, restriction_type = "highest", gamma = 0.05, lb = 0.5, ub = 60, tol = 1e-3)
#'
#' CI = ci.pl(mem_list, restriction_type = "pairwise", eta_cond = c(0.5, 0.5, 0.5),
#' xi_cond = c(5, 5, 5), lb = 0.5, ub = 30, tol = 1e-3)


ci.pl= function(mem_list, restriction_type = "pairwise",
                         alpha = 0.05, gamma = 0.3, r_cond = c(1, 1, 2), t_cond = c(2, 3, 3),
                         eta_cond = c(0.1, 0.1, 0.1), xi_cond = c(10, 10, 10),
                         lb = 0.3, ub = 100, tol = 1e-6,
                         verbose = F, xtol_rel = 1e-3, maxeval = 300, ...){
  # need to vary lb and ub to avoid local minimum, to let ub not be too large.
  # pair, ub = 20; high, ub = 60
  # the run time is very long (5 minutes for one CI, especially for pairwise)
  N = mem_list$N
  k = length(mem_list)-1
  c = 2^k -1

  eval_l_pair <- function(M) {

    eval_f <- function(m) {
      return( list( "objective" = sum (m - N * log(m)),
                    "gradient" = 1 - N/m  ) )
    }

    eval_g_ineq <- function(m){
      constr = rep(0, 2*omega+1)
      grad = matrix(0, ncol = 2*omega+1, nrow = c)
      for (l in 1:omega){
        r = r_cond[l]
        t = t_cond[l]
        eta = eta_cond[l]
        xi = xi_cond[l]
        m_10 = sum(m[(mem_list[r]==1) & (mem_list[t]==0)])
        m_01 = sum(m[(mem_list[r]==0) & (mem_list[t]==1)])
        m_11 = sum(m[(mem_list[r]==1) & (mem_list[t]==1)])
        m_00 = sum(m[(mem_list[r]==0) & (mem_list[t]==0)])
        constr[2*l-1] = m_11 * (M - m_01 - m_10 - m_11) - xi * m_10 * m_01
        constr[2*l] = -m_11 * (M - m_01 - m_10 - m_11) + eta * m_10 * m_01

        dg_dmp = matrix(c(0, 0,
                          -m_11 - xi * m_10, m_11 + eta * m_10,
                          -m_11 - xi * m_01, m_11 + eta * m_01,
                          M - m_01 - m_10 - 2*m_11, -(M - m_01 - m_10 - 2*m_11)), nrow = 2)
        # print(dg_dmp)
        dmp_dm = matrix(0, nrow = 4, ncol = c)
        count = 0
        for (i in c(0, 1)){
          for (j in c(0,1)){
            count = count+1
            for (w in 1:c){
              dmp_dm[count, w] = as.numeric(((mem_list[r]==i) & (mem_list[t]==j))[w])
            }
          }
        }
        dg_dm = dg_dmp %*% dmp_dm
        grad[, 2*l-1] = dg_dm[1,]
        grad[, 2*l] = dg_dm[2,]
      }
      constr[ 2*omega+1] = sum(m) - M
      grad[, 2*omega+1] = matrix(rep(1, 4), nrow = 1) %*% dmp_dm
      # return(constr)
      return(list( "constraints"=constr, "jacobian"=t(grad) ) )
    }

    m_init <- N
    lb <- rep(0.1, c)
    ub <- rep(100*sum(N), c)
    local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                        "xtol_rel" = xtol_rel ) # 1.0e-6

    opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", # only one that can deal with init value out of bound
                  "xtol_rel" = xtol_rel, # 1.0e-6
                  "maxeval" = maxeval, # good amount of iterations, to make sure point-id parameters will be estimated with precision, e.g. 2000.
                  "local_opts" = local_opts, "print_level" = 0, "check_derivatives" = F )

    res <- nloptr::nloptr( x0=m_init,
                   eval_f=eval_f,
                   lb=lb,
                   ub=ub,
                   eval_g_ineq=eval_g_ineq,
                   eval_g_eq=NULL,
                   opts=opts)
    if (all(eval_g_ineq(res$solution)$constraints <= 1e-4)) {
      l = res$objective
    } else{
      l = NA
      print("Infeasible!")
    }
    return(l)
  }

  eval_l_high <- function(M) {

    eval_f <- function(m) {
      return( list( "objective" = sum (m - N * log(m)),
                    "gradient" = 1 - N/m  ) )
    }

    eval_g_ineq <- function(m){
      constr = rep(0, 3)
      grad = matrix(0, ncol = 3, nrow = c)

      constr[1] = (prod(m[I]) - exp(gamma) * prod(m[J])*(M - sum(m)))/M^(2^(k-2))
      constr[2] = (prod(m[J])*(M - sum(m)) - exp(gamma) * prod(m[I]))/M^(2^(k-2))
      constr[3] = (sum(m) - M)/M

      for (i in 1:c){
        if (I[i]){
          I_tilde = I
          I_tilde[i] = FALSE
          grad[i, 1] = (prod(m[I_tilde]) + exp(gamma) * prod(m[J]))/M^(2^(k-2))
          grad[i, 2] = (-exp(gamma) * prod(m[I_tilde]) - prod(m[J]))/M^(2^(k-2))
          grad[i, 3] = 1/M
        } else{
          J_tilde = J
          J_tilde[i] = FALSE
          grad[i, 1] = -exp(gamma) * (  prod(m[J_tilde]) * (M - sum(m)) - prod(m[J]) )/M^(2^(k-2))
          grad[i, 2] = (  prod(m[J_tilde]) * (M - sum(m)) - prod(m[J]) )/M^(2^(k-2))
          grad[i, 3] = 1/M
        }
      }

      return( list( "constraints"=constr, "jacobian"=t(grad) ) )
    }

    m_init <- N
    lb <- rep(0.1, c)
    ub <- rep(100*sum(N), c)
    local_opts <- list( "algorithm" = "NLOPT_LD_MMA",
                        "xtol_rel" = xtol_rel ) # 1e-6

    opts <- list( "algorithm" = "NLOPT_LD_AUGLAG", # only one that can deal with init value out of bound
                  "xtol_rel" = xtol_rel,
                  "maxeval" = maxeval,
                  "local_opts" = local_opts, "print_level" = 0, "check_derivatives" = F )

    res <- nloptr::nloptr( x0=m_init,
                   eval_f=eval_f,
                   lb=lb,
                   ub=ub,
                   eval_g_ineq=eval_g_ineq,
                   eval_g_eq=NULL,
                   opts=opts)
    if (all(eval_g_ineq(res$solution)$constraints <= 1e-4)) {
      l = res$objective
    } else{
      l = NA
      print("Infeasible!")
    }
    return(l)
  }

  if (restriction_type == "pairwise"){
    eval_l = eval_l_pair
    omega = length(r_cond)
  } else if (restriction_type == "highest"){
    eval_l = eval_l_high
    I = rowSums(mem_list[,1:k])%%2 == 1
    J = rowSums(mem_list[,1:k])%%2 == 0
  }

  lb <- lb * sum(N)
  ub <- ub * sum(N)

  # determine a better init
  M_init_vec = seq(lb, ub, length.out = ub/sum(N))
  l_vec = sapply(M_init_vec, eval_l)
  M_init = M_init_vec[which.min(l_vec)]

  m0_init <- M_init
  local_opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                      "xtol_rel" = 1.0e-4 )

  opts <- list( "algorithm" = "NLOPT_LN_COBYLA",
                "xtol_rel" = 1.0e-4,
                "maxeval" = 100, "print_level" = 0, "check_derivatives" = FALSE )

  res <- nloptr::nloptr( x0=m0_init,
                 eval_f=eval_l,
                 lb=lb,
                 ub=ub,
                 eval_g_eq=NULL,
                 opts=opts)
  print(res)
  center = res$solution
  obj_real = res$objective

  close_val = function(x, y, tol = 1e-2){
    tmp = abs(x - y)
    if (tmp < tol) {return(TRUE)} else {return(FALSE)}
  }


  # CI lower bound

  x0 = lb
  x1 = center
  # tol = 1e-1
  print("Finding left endpoint of confidence interval using bisection...")
  while (abs(x0 - x1) > tol){
    x = (x0 + x1)/2
    obj_tmp = eval_l(x)
    if (verbose) {cat("Current Value, Current, Objective: ", x, obj_tmp, "\n")}
    if (obj_tmp <= obj_real + qchisq(1 - alpha, 1)/2){
      x1 = x
    } else{
      x0 = x
    }
  }
  CI_lb = x0

  # CI upper bound


  x0 = center
  x1 = ub
  # tol = 1e-1
  print("Finding right endpoint of confidence interval using bisection...")
  while (abs(x0 - x1) > tol){
    x = (x0 + x1)/2
    obj_tmp = eval_l(x)
    if (verbose) {cat("Current Value, Current, Objective: ", x, obj_tmp, "\n")}
    if (obj_tmp <= obj_real + qchisq(1 - alpha, 1)/2){
      x0 = x
    } else{
      x1 = x
    }
  }
  CI_ub = x1
  CI = c(CI_lb, CI_ub)
  cat("1-",alpha, " confidence interval:", CI, "\n")
  if (close_val(CI_lb, lb)) {warning("Decrease lower bound (CI)!")}
  if (close_val(CI_ub, ub)) {warning("Increase upper bound (CI)!")}
  # if (compute_ID){return(c(ID, CI))} else {return(CI)}
  return(CI)
}

