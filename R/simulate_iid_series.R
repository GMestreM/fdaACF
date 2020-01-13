simulate_iid_brownian_motion <- function(N, v = seq(from = 0, to = 1, length.out = 100), sig = 1){

  #' Simulate a FTS from a brownian motion process
  #'
  #' @description Generate a functional time series from a Brownian Motion process.
  #' Each functional observation is discretized in the points
  #' indicated in \code{v}. The series obtained is i.i.d.
  #' and does not exhibit any kind of serial correlation
  #'
  #' @param N The number of observations of the simulated data.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param sig Standard deviation of the Brownian Motion process,
  #' by default \code{1}.
  #' @return Return the simulated functional time series as a matrix.
  #' @examples
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 20)
  #' sig <- 2
  #' bmotion <- simulate_iid_brownian_motion(N, v, sig)
  #' matplot(v,t(bmotion), type = "l", xlab = "v", ylab = "Value")
  #' @export simulate_iid_brownian_motion
  dv <- length(v)

  # Initialize
  b_motion <- matrix(0, ncol = dv, nrow = N)
  b_motion[,2:dv] <- t(apply(matrix(stats::rnorm( (dv-1) *N, sd = sig), ncol = dv-1, nrow = N), 1, cumsum))

  return(b_motion)
}

simulate_iid_brownian_bridge <- function(N, v = seq(from = 0, to = 1, length.out = 100), sig = 1){

  #' Simulate a FTS from a brownian bridge process
  #'
  #' @description Generate a functional time series from a Brownian Bridge process.
  #' If \eqn{W(t)} is a Wiener process, the Brownian Bridge is
  #' defined as \eqn{W(t) - tW(1)}.
  #' Each functional observation is discretized in the points
  #' indicated in \code{v}. The series obtained is i.i.d.
  #' and does not exhibit any kind of serial correlation.
  #'
  #' @param N The number of observations of the simulated data.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param sig Standard deviation of the Brownian Motion process,
  #' by default \code{1}.
  #' @return Return the simulated functional time series as a matrix.
  #' @examples
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 20)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' matplot(v,t(bbridge), type = "l", xlab = "v", ylab = "Value")
  #' @export simulate_iid_brownian_bridge
  dv <- length(v)

  # Initialize
  b_bridge <- matrix(0, ncol = dv, nrow = N)
  b_bridge[,2:dv] <- simulate_iid_brownian_motion(N, v = v[1:dv-1], sig)
  b_bridge <- b_bridge -  b_bridge[,dv]*matrix(rep(v, times = N),ncol = dv, nrow = N, byrow = T)

  return(b_bridge)
}



