#' Test Inversion Bootstrap Confidence Interval
#'
#' Use bootstrap resampling techniques to compute the confidence interval for the population size M.
#' @importFrom stats qchisq rpois
#'
#' @param mem_list A data.frame that contains the aggregated capture histories. For a
#' capture-recapture study with k samples, mem_list contains one row per observed capture
#' history followed by its frequency. mem_list has 2^k - 1 rows and k+1 columns, where the first k columns
#' can only have 0 or 1 (1 indicates captured in this sample), and the last column is the frequency.
#' See ./data/PWID_dev.RData for an example.
#' @param restriction_type A character string, either "pairwise" (for restrictions on the pairwise dependence),
#' or "highest" (for restrictions on the highest order interaction term lambda_c in a log-linear model).
#' @param alpha A real number between 0 and 1. Specify the confidence level 1-alpha.
#' @param multi_correction An integer. In the two step algorithm, beta is alpha/multi_correction.
#' @param B An integer. Number of bootstrap resamples.
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
#' @param tol A positive number. Tolerance/precision of the optimization algorithm when finding the endpoints of the
#' confidence interval.
#' @param verbose A logical. If TRUE, print intermediate information of the computation.
#' @param search_step A psitive number. The length of interval of a grid on the positive real line, when testing all M's.
#' @param tasks An integer. Number of parallel jobs in bootstrap resampling. tasks = -1 will use all the available cores.
#'
#'
#' @return (1 - alpha) confidence interval for the population size
#' @export
#' @examples

#' CI = ci.tib(mem_list, restriction_type = "pairwise",eta_cond = c(0.5, 0.5, 0.5),
#' xi_cond = c(4, 4, 4), B = 100,
#' lb = 0.5, ub = 10, tol = 1e-2,
#' search_step = 30, tasks = 2, verbose = TRUE)
#'
#' CI = ci.tib(mem_list, restriction_type = "highest", B = 500, gamma = 0.05,
#' lb = 0.5, ub = 20, tol = 1e-3,
#' search_step = 30, tasks = 2, verbose = TRUE)



#'
#' @return 1 - alpha confidence interval for the population size
#' @export
#' @examples
#' print("123")

# library("tidyverse")
# library(nloptr)
# library(parallel)
# library(MASS)
# library(foreach)
# library(doParallel)
# numCores <- detectCores()
# print(numCores)
# registerDoParallel(numCores)  # use multicore, set to the number of our cores
# # load("./rdat/PWID_dev.RData")

ci.tib = function(mem_list, restriction_type = "pairwise", alpha = 0.05,
                       multi_correction = 10, B = 5000,
                       r_cond = c(1, 1, 2), t_cond = c(2, 3, 3),
                       eta_cond = c(0.1, 0.1, 0.1), xi_cond = c(10, 10, 10),
                       gamma = 0.1, verbose = F, lb = 0.8, ub = 3, tol = 0.1, search_step = 10, tasks = -1, ...){
  # lb and ub are used to find acceptance point (phi = 0); then one should ajust lb and ub to find endpoints for bisection such that phi = 1.
  # increase B to increase precision of the endpoints.
  total_cores = parallel::detectCores()
  if (tasks == -1){numCores <- total_cores} else {numCores = min(tasks, total_cores)}
  doParallel::registerDoParallel(numCores)

  # load helper functions
  test_stat_max = function(vec){
    t = max(vec)
    return(t)
  }
  test_stat = test_stat_max

  g_high_unbiased = function(mem_list, M, r_cond =NULL, t_cond=NULL, eta_cond=NULL, xi_cond=NULL, gamma = 1,...){
    g = c()
    var_g = c()
    k = length(mem_list)-1
    # omega = length(r_cond)
    N = mem_list$N
    I = rowSums(mem_list[,1:k])%%2 == 1
    J = rowSums(mem_list[,1:k])%%2 == 0
    g_1 = (prod(N[I]) + exp(gamma) * prod(N[J])*(sum(N) - 2^(k-1) + 1 - M))/M^(2^(k-2))
    g_2 = -(prod(N[I]) + exp(-gamma) * prod(N[J])*(sum(N) - 2^(k-1) + 1 - M))/M^(2^(k-2))
    g = c(g_1, g_2)

    # use an unbiased estimator for variance

    var_g_1 = 1/M^(2^(k-1)) * (prod(N[I]^2) - prod(N[I]^2 - N[I]) +
                                 exp(2*gamma) * ((sum(N[I])^2) * prod(N[J]^2) -
                                                   prod(N[J]^2-N[J]) * (sum(N[I])^2 - sum(N[I])))+
                                 exp(2*gamma) * ((prod(N[J]^2)*(sum(N[J]) - 2^(k-1) + 1 - M)^2) -
                                                   prod(N[J]^2 - N[J])* (sum((N[J]-2)*(N[J]-3)) -
                                                                           (2*M*sum(N[J]-2)) + (M^2))) +
                                 2*exp(gamma) * ( (2^(k-1)-1) * prod(N[I]) * prod(N[J])) +
                                 2*exp(2*gamma) * sum(N[I])*(prod(N[J]^2)*(sum(N[J]) + 1 - M - 2^(k-1)) +
                                                               prod(N[J]^2 - N[J])*(M - sum(N[J]-2)))
    )
    var_g_2 = 1/M^(2^(k-1)) * (prod(N[I]^2) - prod(N[I]^2 - N[I]) +
                                 exp(2*(-gamma)) * ((sum(N[I])^2) * prod(N[J]^2) -
                                                      prod(N[J]^2-N[J]) * (sum(N[I])^2 - sum(N[I])))+
                                 exp(2*(-gamma)) * ((prod(N[J]^2)*(sum(N[J]) - 2^(k-1) + 1 - M)^2) -
                                                      prod(N[J]^2 - N[J])* (sum((N[J]-2)*(N[J]-3)) -
                                                                              (2*M*sum(N[J]-2)) + (M^2))) +
                                 2*exp((-gamma)) * ( (2^(k-1)-1) * prod(N[I]) * prod(N[J])) +
                                 2*exp(2*(-gamma)) * sum(N[I])*(prod(N[J]^2)*(sum(N[J]) + 1 - M - 2^(k-1)) +
                                                                  prod(N[J]^2 - N[J])*(M - sum(N[J]-2)))
    )
    g_dat = tibble::tibble(g = g, std_g = c(sqrt(var_g_1), sqrt(var_g_2)))
    return(g_dat)
  }


  g_pair = function(mem_list, M, r_cond =NULL, t_cond=NULL, eta_cond=NULL, xi_cond=NULL, gamma = 1,...){
    # check conditions bound
    g = c()
    var_g = c()
    omega = length(r_cond)
    N = mem_list$N
    for (l in 1:omega){
      N_10 = sum(N[(mem_list[r_cond[l]]==1) & (mem_list[t_cond[l]]==0)])
      N_01 = sum(N[(mem_list[r_cond[l]]==0) & (mem_list[t_cond[l]]==1)])
      N_11 = sum(N[(mem_list[r_cond[l]]==1) & (mem_list[t_cond[l]]==1)])
      N_00_obs = sum(N[(mem_list[r_cond[l]]==0) & (mem_list[t_cond[l]]==0)])
      g_tmp = c((-N_11^2 - N_10*N_11 - N_01*N_11-xi_cond[l]*N_10*N_01 + N_11+M*N_11)/M^2 , -(-N_11^2 - N_10*N_11 - N_01*N_11-eta_cond[l]*N_10*N_01 + N_11+M*N_11)/M^2)

      var_g_1 = (N_11*N_10*(N_10 + N_11 - 1) + N_11*N_01*(N_01 + N_11 - 1) + 4*N_11^3 - 6*N_11^2 + (M^2 + 2*M+4)*N_11 + xi_cond[l]^2*N_10*N_01*(N_10 + N_01 - 1) + (2 + 4*xi_cond[l])*N_11*N_10*N_01 + 2*(N_10 + N_01 - M - 1)*(2*N_11^2-N_11) - 2*(M+1)*N_11*(N_10 + N_01))/M^4
      var_g_2 = (N_11*N_10*(N_10 + N_11 - 1) + N_11*N_01*(N_01 + N_11 - 1) + 4*N_11^3 - 6*N_11^2 + (M^2 + 2*M+4)*N_11 + eta_cond[l]^2*N_10*N_01*(N_10 + N_01 - 1) + (2 + 4*eta_cond[l])*N_11*N_10*N_01 + 2*(N_10 + N_01 - M - 1)*(2*N_11^2-N_11) - 2*(M+1)*N_11*(N_10 + N_01))/M^4
      g = c(g, g_tmp)
      var_g = c(var_g, var_g_1, var_g_2)
    }
    g_dat = tibble::tibble(g = g, var_g = var_g, std_g = sqrt(var_g))
    return(g_dat)
  }
  set_est_high_sa = function(mem_list, gamma = 1, ...){
    # sample analog estimator
    k = length(mem_list)-1
    N = mem_list$N
    I = rowSums(mem_list[,1:k])%%2 == 1
    J = rowSums(mem_list[,1:k])%%2 == 0
    M_max = sum(N) + prod(N[I])/prod((N[J]+1))*exp(gamma)
    M_min = sum(N) + prod(N[I])/prod((N[J]+1))*exp(-gamma)
    return(c(M_min, M_max))
  }

  set_est_pair_sa = function(mem_list,r_cond = NULL, t_cond = NULL, eta_cond = NULL, xi_cond = NULL, gamma = 1, ...){
    # sample analog estimator
    N = mem_list$N
    omega = length(r_cond)
    M_min = c()
    M_max = c()
    for (l in 1:omega){
      N_10 = sum(N[(mem_list[r_cond[l]]==1) & (mem_list[t_cond[l]]==0)])
      N_01 = sum(N[(mem_list[r_cond[l]]==0) & (mem_list[t_cond[l]]==1)])
      N_11 = sum(N[(mem_list[r_cond[l]]==1) & (mem_list[t_cond[l]]==1)])
      N_00_obs = sum(N[(mem_list[r_cond[l]]==0) & (mem_list[t_cond[l]]==0)])

      min_tmp = eta_cond[l]*N_10*N_01/(N_11 + 1) + N_10 + N_01 + N_11
      max_tmp = xi_cond[l]*N_10*N_01/(N_11 + 1) + N_10 + N_01 + N_11
      M_min = c(M_min, min_tmp)
      M_max = c(M_max, max_tmp)
    }
    # print(c(M_min, M_max))
    # if (max(max(M_min), sum(N)) > min(M_max)) warning("Estimated upper bound is less than lower bound.")
    return(c(max(max(M_min), sum(N)), min(M_max)))
  }

  dist_sort = function(set_est, M_vec){
    left_dist = abs(M_vec - set_est[1])
    center_dist = abs(M_vec - (set_est[1] + set_est[2])/2)
    right_dist = abs(M_vec - set_est[2])
    dist = pmin(left_dist, right_dist, center_dist)
    return(M_vec[order(dist)])
  }


  if (restriction_type == "pairwise"){
    g_func = g_pair
    set_est = set_est_pair_sa  # not a good set estimate, but can be a good starting point for CI , if N is in Theta.
  } else if (restriction_type == "highest"){
    g_func = g_high_unbiased
    set_est = set_est_high_sa  # "sa" means sample analogue estimator
  }
  k = length(mem_list)-1
  # print(k)
  a = 2^k -1
  beta = alpha/multi_correction
  set_estimate = set_est(mem_list, r_cond =r_cond, t_cond=t_cond, eta_cond=eta_cond, xi_cond=xi_cond, gamma = gamma)
  obs_sum = sum(mem_list$N)

  phi_function = function(M){
    if (verbose) {print(M)}
    g_dat = g_func(mem_list, M = M, r_cond =r_cond, t_cond=t_cond, eta_cond=eta_cond, xi_cond=xi_cond, gamma = gamma)
    mu_K = g_dat$g
    S_K = g_dat$std_g
    # n_boot = t(matrix(rpois(n = a*B, mem_list$N), nrow = a))
    func = foreach::`%dopar%`
    K_n = foreach::`%dopar%`(foreach::foreach(iterators::icount(B), .combine=c), {
      mem_list_tmp = mem_list
      mem_list_tmp$N = rpois(n = a, mem_list$N)
      g_dat_tmp = g_func(mem_list_tmp, M = M, r_cond =r_cond, t_cond=t_cond, eta_cond=eta_cond, xi_cond=xi_cond, gamma = gamma)
      max((mu_K - g_dat_tmp$g)/g_dat_tmp$std_g)
    })
    K_n_crit = sort(K_n)[ceiling(B*(1-beta))]
    Q_n = K_n_crit * S_K + mu_K
    lam_star = pmin(Q_n, 0)
    lam_star

    T_stat = test_stat(mu_K/S_K)
    A_n = foreach::`%dopar%`(foreach::foreach(iterators::icount(B), .combine=c) , {
      mem_list_tmp = mem_list
      mem_list_tmp$N = rpois(n = a, mem_list$N)
      g_dat_tmp = g_func(mem_list_tmp, M = M, r_cond =r_cond, t_cond=t_cond, eta_cond=eta_cond, xi_cond=xi_cond, gamma = gamma)
      test_stat((-mu_K + g_dat_tmp$g + lam_star)/g_dat_tmp$std_g)
    })

    A_n_crit = sort(A_n)[ceiling(B*(1-alpha+beta))]
    phi = 1 - (all(Q_n <= 0) | T_stat <= A_n_crit)
    return(phi)
  }
  # find a point that is in CI
  phi_tmp = 1
  M_vec = seq(lb*obs_sum, ub*obs_sum, by = search_step)
  M_vec = dist_sort(set_est = set_estimate, M_vec = M_vec)
  M_vec = c(set_estimate, M_vec) # can modify to optimize
  # M_middle = mean(set_estimate)
  # should not set M_middle arbitrarily, since it is possible that, a given model is contradictary to the data, evidenced by that all the possible M are rejected under a model. should give warning and stop the program.
  phi_vec = c()
  for (M in M_vec){
    phi_tmp = phi_function(M)
    phi_vec = c(phi_vec, phi_tmp)
    if (phi_tmp == 0){
      M_middle = M
      break
    }
  }
  if (prod(phi_vec)) {
    warning("The model is very likely to be contradictory to the data.")
    return(c(NA, NA))
  }
  # find proper endpoints for bisections
  while (!phi_function(obs_sum * lb)){lb = lb/2}
  if (restriction_type == "pairwise") {while (!phi_function(obs_sum * ub)){ub = ub*2}} # will this work?
  # get left end of CI
  print("Finding left endpoint of confidence interval using bisection...")
  lower_left = obs_sum * lb
  # upper_left = set_estimate[1]
  upper_left = M_middle
  while ((upper_left - lower_left) > tol){
    # if (verbose) print(c(lower_left, upper_left))
    M = mean(c(upper_left, lower_left))
    # if (verbose) print(M)
    phi = phi_function(M)
    if (phi == 1){
      lower_left = M
    } else{
      upper_left = M
    }

  }

  CI_left = lower_left

  # get right CI
  if (restriction_type == "highest") {
    CI_right = Inf
  } else {
    print("Finding right endpoint of confidence interval using bisection...")
    lower_right =  M_middle
    # upper_left = set_estimate[1]
    upper_right = obs_sum * ub
    while ((upper_right - lower_right) > tol){
      # if (verbose) print(c(lower_right, upper_right))
      M = mean(c(upper_right, lower_right))
      # if (verbose) print(M)
      phi = phi_function(M)
      if (phi == 1){
        upper_right = M
      } else{
        lower_right = M
      }
    }
    CI_right = upper_right
  }

  CI = c(CI_left, CI_right)
  # if (verbose) cat("Set estimate: ", set_estimate, "\n")
  # if (verbose) cat("CI", CI, "\n")
  cat("1 -",alpha, " confidence interval:", CI, "\n")
  if (CI_left > CI_right){
    CI = c(-1, -1)
  }
  doParallel::stopImplicitCluster()
  return(CI)
}
