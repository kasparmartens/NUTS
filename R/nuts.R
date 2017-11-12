#' No-U-Turn sampler
#'
#' @param theta Initial value for the parameters
#' @param f log-likelihood function (up to a constant)
#' @param grad_f the gradient of the log-likelihood function
#' @param n_iter Number of MCMC iterations
#' @param M_diag Diagonal elements of the mass matrix in HMC. Defaults to ones.
#' @param M_adapt Parameter M_adapt in algorithm 6 in the NUTS paper
#' @param delta Target acceptance ratio, defaults to 0.5
#' @param max_treedepth Maximum depth of the binary trees constructed by NUTS
#' @param eps Starting guess for epsilon
#' @param verbose logical. Message diagnostic info each iteration? Default \code{TRUE}.
#' @return Matrix with the trace of sampled parameters. Each mcmc iteration in rows and parameters in columns.
#' @export
NUTS <- function(theta, f, grad_f, n_iter, M_diag = NULL, M_adapt = 50, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE) {
  theta_trace   <- matrix(0, n_iter, length(theta))
  epsilon_trace <- rep(NA, n_iter)
  depth_trace   <- rep(NA, n_iter)
  par_list <- list(M_adapt = M_adapt)
  for(iter in 1:n_iter) {
    nuts <- NUTS_one_step(theta, iter, f, grad_f, par_list, delta = delta, max_treedepth = max_treedepth, eps = eps, verbose = verbose)
    theta <- nuts$theta
    par_list <- nuts$pars
    theta_trace[iter, ] <- theta
    epsilon_trace[iter] <- par_list$eps
    depth_trace[iter]   <- par_list$depth
    
    if(verbose) {message(paste("n:", format(iter, width=log(n_iter, 10) + 1), "j:", format(par_list$depth, width=2), "e:", par_list$eps))}
  }
  
  attributes(theta_trace)$epsilon <- epsilon_trace
  attributes(theta_trace)$depth   <- depth_trace
  
  theta_trace
}


NUTS_one_step <- function(theta, iter, f, grad_f, par_list, delta = 0.5, max_treedepth = 10, eps = 1, verbose = TRUE) {
  kappa <- 0.75
  t0 <- 10
  gamma <- 0.05
  M_adapt <- par_list$M_adapt
  M_diag  <- if(is.null(par_list$M_diag)) {rep(1, length(theta))} else {par_list$M_diag}
  
  if(iter == 1) {
    eps     <- find_reasonable_epsilon(theta, f, grad_f, M_diag, eps = eps, verbose = verbose)
    mu      <- log(10*eps)
    H       <- 0
    eps_bar <- 1
  } else {
    eps     <- par_list$eps
    eps_bar <- par_list$eps_bar
    H       <- par_list$H
    mu      <- par_list$mu
  }
  
  eps <- if(iter > M_adapt) {stats::runif(1, 0.9 * eps_bar, 1.1 * eps_bar)} else {eps}
  
  r0  <- stats::rnorm(length(theta), 0, sqrt(M_diag))
  
  ## take log of u, because on meaningful-sized problems u itself will usually underflow, leading to kablooie
  ## u ~ Uniform(0, p) => -log(u) + p ~ Exponential(1)
  log_u <- joint_log_density(theta, r0, f, M_diag) - stats::rexp(1)
  if(!is.finite(log_u)) {
    warning("NUTS: sampled slice u is not a number")
    log_u <- log(stats::runif(1, 0, 1e5))
  }
  
  theta_minus <- theta
  theta_plus  <- theta
  r_minus     <- r0
  r_plus      <- r0
  
  j <- 0
  n <- 1
  s <- 1
  
  while(s == 1) {
    # choose direction {-1, 1}
    direction <- sample(c(-1, 1), 1)
    if(direction == -1) {
      temp <- build_tree(theta_minus, r_minus, log_u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_minus <- temp$theta_minus
      r_minus <- temp$r_minus
    } else {
      temp <- build_tree(theta_plus, r_plus, log_u, direction, j, eps, theta, r0, f, grad_f, M_diag)
      theta_plus <- temp$theta_plus
      r_plus <- temp$r_plus
    }
    if(!is.finite(temp$s)) {temp$s <- 0}
    if(temp$s == 1) {
      if(stats::runif(1) < temp$n / n) {theta <- temp$theta}
    }
    n <- n + temp$n
    s <- check_NUTS(temp$s, theta_plus, theta_minus, r_plus, r_minus)
    j <- j + 1
    if(j > max_treedepth) {
      warning("NUTS: Reached max tree depth")
      break
    }
  }
  if(iter <= M_adapt) {
    H <- (1 - 1/(iter + t0))*H + 1/(iter + t0) * (delta - temp$alpha / temp$n_alpha)
    log_eps <- mu - sqrt(iter)/gamma * H
    eps_bar <- exp(iter**(-kappa) * log_eps + (1 - iter**(-kappa)) * log(eps_bar))
    eps <- exp(log_eps)
  } else {
    eps <- eps_bar
  }
  
  return(list(theta = theta,
              pars = list(eps = eps, eps_bar = eps_bar, H = H, mu = mu, M_adapt = M_adapt, M_diag = M_diag, depth = j)))
}
