estimate_iid_distr_MC <- function(Y,v,autocovSurface,matindex,nsimul= 10000,figure = FALSE, ...){

  #' Estimate distribution of the fACF under the iid. hypothesis using MC method
  #'
  #' @description Estimate the distribution of the autocorrelation function
  #' under the hypothesis of strong functional white noise. This
  #' function uses Montecarlo's method to estimate the distribution.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param autocovSurface An \eqn{(m x m)} matrix with the discretized
  #' values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained
  #' by calling the function \code{obtain_autocovariance}.
  #' The value \eqn{m} indicates the number of points observed
  #' in each curve.
  #' @param matindex A vector containing the L2 norm of
  #' the autocovariance function. It can be obtained by calling
  #' function \code{obtain_suface_L2_norm}.
  #' @param nsimul Positive integer indicating the number of
  #' MC simulations that will be used to estimate the distribution
  #' of the statistic. Increasing the number of simulations will
  #' improve the estimation, but it will increase the computational
  #' time.
  #' By default, \code{nsimul = 10000}.
  #' @param figure Logical. If \code{TRUE}, plots the
  #' estimated distribution.
  #' @param ... Further arguments passed to the  \code{plot}
  #' function.
  #' @return Return a list with:
  #' \itemize{
  #'     \item \code{ex}: Knots where the
  #'     distribution has been estimated
  #'     \item \code{ef}: Discretized values of
  #'     the estimated distribution.
  #'     \item \code{Reig}: Raw values of the i.i.d.
  #'     statistic for each MC simulation.
  #' }
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' autocovSurface <- obtain_autocovariance(Y,nlags)
  #' matindex <- obtain_suface_L2_norm (v,autocovSurface)
  #' # Remove lag 0
  #' matindex <- matindex[-1]
  #' MC_dist <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex)
  #' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
  #' grid()
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' N <- 400
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 20
  #' autocovSurface <- obtain_autocovariance(Y,nlags)
  #' matindex <- obtain_suface_L2_norm (v,autocovSurface)
  #' # Remove lag 0
  #' matindex <- matindex[-1]
  #' MC_dist <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex)
  #' plot(MC_dist$ex,MC_dist$ef,type = "l",main = "ecdf obtained by MC simulation")
  #' grid()
  #' }
  #' @export estimate_iid_distr_MC

  mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
  l <- obtain_autocov_eigenvalues(v,Y - mat.means)


  neig <- length(l)
  Reig <- rep(0,nsimul)

  for(jj in 1:neig){
    for(kk in 1:neig){
      Reig <- Reig + stats::rchisq(n = nsimul,df = 1)*l[jj]*l[kk]
    }
  }
  Reig=Reig/nrow(Y)

  ecdf.aux <- stats::ecdf(Reig)
  ex <- stats::knots(ecdf.aux)
  ef <- ecdf.aux(ex)
  if(figure){
    # Check if any additional plotting parameters are present
    arguments <- list(...)
    if("main" %in% names(arguments)){
      graphics::plot(ex,ef,type = "l",...)
    }else{
      graphics::plot(ex,ef,type = "l",main = "ecdf obtained by MC simulation",...)
    }
    graphics::grid()
  }
  iid.distribution <- list()
  iid.distribution$ex <- ex
  iid.distribution$ef <- ef
  iid.distribution$Reig <- Reig
  return(iid.distribution)
}


estimate_iid_distr_Imhof <- function(Y,v,autocovSurface,matindex,figure = FALSE,...){

  #' Estimate distribution of the fACF under the iid. hypothesis using Imhof's method
  #'
  #' @description Estimate the distribution of the autocorrelation function
  #' under the hypothesis of strong functional white noise. This
  #' function uses Imhof's method to estimate the distribution.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param autocovSurface An \eqn{(m x m)} matrix with the discretized
  #' values of the autocovariance operator \eqn{\hat{C}_{0}}, obtained
  #' by calling the function \code{obtain_autocovariance}.
  #' The value \eqn{m} indicates the number of points observed
  #' in each curve.
  #' @param matindex A vector containing the L2 norm of
  #' the autocovariance function. It can be obtained by calling
  #' function \code{obtain_suface_L2_norm}.
  #' @param figure Logical. If \code{TRUE}, plots the
  #' estimated distribution.
  #' @param ... Further arguments passed to the  \code{plot}
  #' function.
  #' @return Return a list with:
  #' \itemize{
  #'     \item \code{ex}: Knots where the
  #'     distribution has been estimated
  #'     \item \code{ef}: Discretized values of
  #'     the estimated distribution.
  #' }
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' autocovSurface <- obtain_autocovariance(Y,nlags)
  #' matindex <- obtain_suface_L2_norm (v,autocovSurface)
  #' # Remove lag 0
  #' matindex <- matindex[-1]
  #' Imhof_dist <- estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex)
  #' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
  #' grid()
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' N <- 400
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' autocovSurface <- obtain_autocovariance(Y,nlags)
  #' matindex <- obtain_suface_L2_norm (v,autocovSurface)
  #' # Remove lag 0
  #' matindex <- matindex[-1]
  #' Imhof_dist <- estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex)
  #' plot(Imhof_dist$ex,Imhof_dist$ef,type = "l",main = "ecdf obtained by Imhof's method")
  #' grid()
  #' }
  #' @export estimate_iid_distr_Imhof

  mat.means <- matrix(rep(colMeans(Y),nrow(Y)),ncol=ncol(Y),byrow = TRUE)
  l <- obtain_autocov_eigenvalues(v,Y - mat.means)

  nl <- length(l)
  Reig <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex,figure = FALSE)$Reig
  w <- seq(from = 0, to = max(Reig), length.out = 250)
  x <- w * nrow(Y)
  nx <- length(x)
  # Compute products of eigenvalues
  #Calculamos serie de productos de autovalores
  L <- rep(0,nl*nl)
  k <- 1;
  for (i in l) {
    for (j in l){
      L[k] <- i*j;
      k <- k+1
    }
  }

  # Obtain ecdf using Imhof function
  cdf <- rep(0,nx)
  for (i in 1:nx){
    cdf[i]=1-CompQuadForm::imhof(x[i],L)$Qq
  }

  ex <- w
  ef <- t(cdf)
  if(figure){
    # Check if any additional plotting parameters are present
    arguments <- list(...)
    if("main" %in% names(arguments)){
      graphics::plot(ex,ef,type = "l",...)
    }else{
      graphics::plot(ex,ef,type = "l",main = "ecdf obtained by Imhof estimation")
    }
    graphics::grid()
  }
  iid.distribution <- list()
  iid.distribution$ex <- ex
  iid.distribution$ef <- ef
  return(iid.distribution)
}




obtain_autocov_eigenvalues <- function(v,Y,epsilon = 0.0001){

  #' Estimate eigenvalues of the autocovariance function
  #'
  #' @description Estimate the eigenvalues of the sample autocovariance
  #' function \eqn{\hat{C}_{0}}. This functions returns the
  #' eigenvalues which are greater than the value \code{epsilon}.
  #'
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param epsilon Value used to determine how many
  #' eigenvalues will be returned. The eigenvalues
  #' \eqn{\lambda_{j}/\lambda_{1} > \code{epsilon}} 
  #' will be returned.
  #' By default \code{epsilon = 0.0001}.
  #' @return A vector containing the \eqn{k} eigenvalues
  #' greater than \code{epsilon}.
  #' @examples
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' lambda <- obtain_autocov_eigenvalues(v = v, Y = Y)
  #' 
  #'@export obtain_autocov_eigenvalues
  # Requires package pracma (needs trapz function)
  nt <- nrow(Y)
  nv <- ncol(Y)

  # Obtain matrix W[nt x nt] where w(i,j) = integral of the product of the curves i and j

  W <- matrix(0,nrow = nt,ncol = nt)

  for(ii in 1:nt){
    mat.aux <- matrix(rep(Y[ii,],each = nt),nrow = nt, ncol = nv)*Y
    for(jj in 1:nt){
      W[ii,jj] = pracma::trapz(v,mat.aux[jj,])
    }
  }

  # Obtain eigenvalues and eigenfunctions of W/n
  eigenlist <- eigen(W/nt)

  if(F){
    # OLD 
    # Select the k eigenvalues higher than epsilon
    k <- which(eigenlist$values>epsilon)
  }else{
    # NEW
    # Obtain all first $m$ eigenvalues such that
    # $$
    # \lambda_{m}/\lambda_{1} > \varepsilon
    # $$
    lambd1 <- eigenlist$values[1]
    frac <- eigenlist$values/lambd1
    k <- which(frac > epsilon)
  }
  

  lambda <- eigenlist$values[k]
  lambda
}


obtain_suface_L2_norm <- function(v,autocovSurface){

  #' Obtain L2 norm of the autocovariance functions
  #'
  #' @description Returns the L2 norm of the lagged autocovariance
  #' functions \eqn{\hat{C}_{h}}. The L2 norm of these
  #' functions is defined as
  #' \deqn{\sqrt(\int \int \hat{C}^{2}_{h}(u,v)du dv)}
  #'
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param autocovSurface An \eqn{(m x m)} matrix with
  #' the discretized values of the autocovariance operator
  #' \eqn{\hat{C}_{0}}, obtained by calling the
  #' function \code{obtain_autocovariance}. The value
  #' \eqn{m} indicates the number of points observed in
  #' each curve.
  #' @return A vector containing the L2 norm of the
  #' lagged autocovariance functions \code{autocovSurface}.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' autocovSurface <- obtain_autocovariance(Y=Y,nlags = nlags)
  #' norms <- obtain_suface_L2_norm(v = v,autocovSurface = autocovSurface)
  #' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 1)
  #' title(sub = paste0("Lag ",1," - L2 Norm: ",norms[2]))
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' N <- 400
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 2
  #' autocovSurface <- obtain_autocovariance(Y=Y,nlags = nlags)
  #' norms <- obtain_suface_L2_norm(v = v,autocovSurface = autocovSurface)
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,3))
  #' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 0)
  #' title(sub = paste0("Lag ",0," - L2 Norm: ",norms[1]))
  #' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 1)
  #' title(sub = paste0("Lag ",1," - L2 Norm: ",norms[2]))
  #' plot_autocovariance(fun.autocovariance = autocovSurface,lag = 2)
  #' title(sub = paste0("Lag ",2," - L2 Norm: ",norms[3]))
  #' par(opar)
  #' }
  #' @export obtain_suface_L2_norm
  repsize <- length(autocovSurface)
  matindex <- rep(NA,repsize)

  nt <- nrow(autocovSurface$Lag0)

  # Obtain L2 norm of squares of autocov surfaces
  for(rr in 1:repsize){
    surf.aux <- autocovSurface[[paste("Lag",rr-1,sep="")]]^2

    norm.vec <- matrix(NA,nrow = nt, ncol = 1)
    for(ii in 1:nrow(surf.aux)){
      norm.vec[ii] <- pracma::trapz(v,surf.aux[ii,])
    }

    norm.aux <- pracma::trapz(v,norm.vec)
    matindex[rr] <- norm.aux
  }
  return(matindex)
}
