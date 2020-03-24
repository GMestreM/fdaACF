obtain_FACF <- function(Y,v,nlags,ci=0.95,estimation = "MC",figure = TRUE,...){

  #' Obtain the autocorrelation function for a given functional time series.
  #'
  #' @description Estimate the lagged autocorrelation function for a given
  #' functional time series and its distribution under the
  #' hypothesis of strong functional white noise. This graphic tool
  #' can be used to identify seasonal patterns in the functional
  #' data as well as auto-regressive or moving average terms.
  #' i.i.d. bounds are included to test the presence of serial
  #' correlation in the data.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param nlags Number of lagged covariance operators
  #' of the functional time series that will be used
  #' to estimate the autocorrelation function.
  #' @param ci A value between 0 and 1 that indicates
  #' the confidence interval for the i.i.d. bounds
  #' of the autocorrelation function. By default
  #' \code{ci = 0.95}.
  #' @param estimation Character specifying the
  #' method to be used when estimating the distribution
  #' under the hypothesis of functional white noise.
  #' Accepted values are:
  #' \itemize{
  #'    \item "MC": Monte-Carlo estimation.
  #'    \item "Imhof": Estimation using Imhof's method.
  #' }
  #' By default, \code{estimation = "MC"}.
  #' @param figure Logical. If \code{TRUE}, plots the
  #' estimated autocorrelation function with the
  #' specified i.i.d. bound.
  #' @param ... Further arguments passed to the \code{plot_FACF}
  #' function.
  #' @return Return a list with:
  #' \itemize{
  #'     \item \code{Blueline}: The upper prediction
  #'     bound for the i.i.d. distribution.
  #'     \item \code{rho}: Autocorrelation values for
  #'     each lag of the functional time series.
  #' }
  #' @references \url{https://www.sciencedirect.com/science/article/pii/S0047259X17303512}
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 5)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' obtain_FACF(Y,v,20)
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' data(elec_prices)
  #' v <- seq(from = 1, to = 24)
  #' nlags <- 50
  #' obtain_FACF(Y = as.matrix(elec_prices), 
  #' v = v,
  #' nlags = nlags,
  #' ci = 0.95,
  #' figure = TRUE)
  #' }
  #' @export obtain_FACF


  # Obtain autocovariance surfaces
  autocovSurface <- obtain_autocovariance(Y,nlags)

  # Obtain L2 norm of autocov surface
  matindex <- obtain_suface_L2_norm (v,autocovSurface)

  # Remove lag 0
  matindex <- matindex[-1]

  # Estimate distribution of iid bound
  if(estimation == "MC"){
    iid.distribution <- estimate_iid_distr_MC(Y,v,autocovSurface,matindex,figure = FALSE)
  }else{
    iid.distribution <- estimate_iid_distr_Imhof(Y,v,autocovSurface,matindex,figure = FALSE)
  }

  # Obtain autocorrelation
  normalization.value <- pracma::trapz(v, diag(autocovSurface$Lag0))
  rho = sqrt(matindex)/normalization.value

  # Obtain iid bound for specified confidence value
  Blueline <- rep(NA,length(ci))
  for(ii in 1:length(ci)){
    Blueline[ii] <- sqrt(iid.distribution$ex[min(which(iid.distribution$ef >= ci[ii])) ])/normalization.value
  }

  if(figure){
    plot_FACF(rho,Blueline,ci,...)
  }


  output <- list()
  output$Blueline <- Blueline
  output$rho      <- rho
  return(output)
}



obtain_autocovariance <- function(Y,nlags){

  #' Estimate the autocovariance function of the series
  #'
  #' @description Obtain the empirical autocovariance function for
  #' lags \eqn{= 0,...,}\code{nlags} of the functional time
  #' series. Given \eqn{Y_{1},...,Y_{T}} a functional time
  #' series, the sample autocovariance functions
  #' \eqn{\hat{C}_{h}(u,v)} are given by:
  #' \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
  #' where
  #' \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
  #' denotes the sample mean function.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param nlags Number of lagged covariance operators
  #' of the functional time series that will be used
  #' to estimate the autocorrelation function.
  #' @return Return a list with the lagged autocovariance
  #' functions estimated from the data. Each function is given
  #' by a \eqn{(m x m)} matrix, where \eqn{m} is the
  #' number of points observed in each curve.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,
  #'                                         nlags = nlags)
  #' image(x = v, y = v, z = lagged_autocov$Lag0)
  #' 
  #' \donttest{
  #' # Example 2
  #'
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 10
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,
  #'                                         nlags = nlags)
  #' image(x = v, y = v, z = lagged_autocov$Lag0)
  #' image(x = v, y = v, z = lagged_autocov$Lag10)
  #'
  #' # Example 3
  #'
  #' require(fields)
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 4
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,
  #'                                         nlags = nlags)
  #' z_lims <- range(lagged_autocov$Lag0)
  #' colors <- heat.colors(12)
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,5))
  #' par(oma=c( 0,0,0,6)) 
  #' for(k in 0:nlags){
  #'    image(x=v,
  #'          y=v,
  #'          z = lagged_autocov[[paste0("Lag",k)]],
  #'          main = paste("Lag",k),
  #'          col = colors,
  #'          xlab = "u",
  #'          ylab = "v")
  #' }
  #' par(oma=c( 0,0,0,2.5)) # reset margin to be much smaller.
  #' image.plot( legend.only=TRUE, legend.width = 2,zlim=z_lims, col = colors)
  #' par(opar)
  #' }
  #' @export obtain_autocovariance

  nt <- nrow(Y)
  nv <- ncol(Y)

  # Obtain mean function of Y
  mean.Y <- colMeans(Y)

  # Compute autocovariance functions
  fun.autocovariance <- list()
  for(kk in 0:nlags){
    fun.autocovariance[[paste("Lag",kk,sep = "")]] <- matrix(0,nrow = nv,ncol = nv)
    for(ind_curve in (1+kk):nt){
      fun.autocovariance[[paste("Lag",kk,sep = "")]] <- fun.autocovariance[[paste("Lag",kk,sep = "")]] + (as.numeric(Y[ind_curve - kk,] - mean.Y)) %*% t(as.numeric(Y[ind_curve,] - mean.Y))
    }
    fun.autocovariance[[paste("Lag",kk,sep = "")]] <- fun.autocovariance[[paste("Lag",kk,sep = "")]]/(nt-1)
  }
  return(fun.autocovariance)
}

obtain_autocorrelation <- function(Y,v = seq(from = 0, to = 1, length.out = ncol(Y)),nlags){

  #' Estimate the autocorrelation function of the series
  #'
  #' @description Obtain the empirical autocorrelation function for
  #' lags \eqn{= 0,...,}\code{nlags} of the functional time
  #' series. Given \eqn{Y_{1},...,Y_{T}} a functional time
  #' series, the sample autocovariance functions
  #' \eqn{\hat{C}_{h}(u,v)} are given by:
  #' \deqn{\hat{C}_{h}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(Y_{i+h}(v) - \overline{Y}_{T}(v))}
  #' where
  #' \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
  #' denotes the sample mean function. By normalizing these
  #' functions using the normalizing factor
  #' \eqn{\int\hat{C}_{0}(u,u)du}, the range of the
  #' autocovariance functions becomes \eqn{(0,1)}; thus
  #' defining the autocorrelation functions of the series
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Discretization points of the curves, by default
  #' \code{seq(from = 0, to = 1, length.out = 100)}.
  #' @param nlags Number of lagged covariance operators
  #' of the functional time series that will be used
  #' to estimate the autocorrelation function.
  #' @return Return a list with the lagged autocorrelation
  #' functions estimated from the data. Each function is given
  #' by a \eqn{(m x m)} matrix, where \eqn{m} is the
  #' number of points observed in each curve.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' lagged_autocor <- obtain_autocorrelation(Y = bbridge,
  #'                                         nlags = nlags)
  #' image(x = v, y = v, z = lagged_autocor$Lag0)
  #' 
  #' \donttest{
  #' # Example 2
  #' require(fields)
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 4
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,
  #'                                         nlags = nlags)
  #' lagged_autocor <- obtain_autocorrelation(Y = bbridge,
  #'                                          v = v,
  #'                                          nlags = nlags)
  #'
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,2))
  #' z_lims <- range(lagged_autocov$Lag1)
  #' colors <- heat.colors(12)
  #' image.plot(x = v, 
  #'            y = v,
  #'            z = lagged_autocov$Lag1,
  #'            legend.width = 2,
  #'            zlim = z_lims,
  #'            col = colors,
  #'            xlab = "u",
  #'            ylab = "v",
  #'            main = "Autocovariance")
  #' z_lims <- range(lagged_autocor$Lag1)
  #' image.plot(x = v, 
  #'            y = v,
  #'            z = lagged_autocor$Lag1,
  #'            legend.width = 2,
  #'            zlim = z_lims,
  #'            col = colors,
  #'            xlab = "u",
  #'            ylab = "v",
  #'            main = "Autocorrelation")
  #' par(opar)
  #' }
  #' @export obtain_autocorrelation

  fun.autocovariance <- obtain_autocovariance(Y,nlags)
  normalization.value <- pracma::trapz(v,diag(fun.autocovariance$Lag0))
  fun.autocorrelation <- list()
  for(kk in 0:nlags){
    fun.autocorrelation[[paste("Lag",kk,sep = "")]] <- fun.autocovariance[[paste("Lag",kk,sep = "")]]/normalization.value;
  }
  return(fun.autocorrelation)
}


obtain_covariance <- function(Y,X,nlags = 0){
  #' Estimate the covariance function of two FTS
  #'
  #' @description Obtain the empirical covariance function 
  #' of the functional time
  #' series \eqn{Y_{t}} and \eqn{X_{t}}. Given 
  #' \eqn{Y_{1},...,Y_{T}} and \eqn{X_{1},...,X_{T}} two
  #' functional time
  #' series, the sample autocovariance functions
  #' \eqn{\hat{\Sigma}_{Y,X}(u,v)} are given by:
  #' \deqn{\hat{\Sigma}_{Y,X}(u,v) =  \frac{1}{T} \sum_{i=1}^{T-h}(Y_{i}(u) - \overline{Y}_{T}(u))(X_{i}(v) - \overline{X}_{T}(v))}
  #' where
  #' \eqn{ \overline{Y}_{T}(u) = \frac{1}{T} \sum_{i = 1}^{T} Y_{i}(t)}
  #' denotes the sample mean function for functional time series \eqn{Y_{t}}.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the first functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param X Matrix containing the discretized values
  #' of the second functional time series. The dimension of the
  #' matrix is \eqn{(n x q)}, where \eqn{n} is the
  #' number of curves and \eqn{q} is the number of points
  #' observed in each curve.
  #' @param nlags Number of lagged covariance operators
  #' of the functional time series that will be estimated.
  #' By default \code{nlags = 0}.
  #' @return Return a list that contains the lagged covariance
  #' functions of time series \eqn{Y} and \eqn{X}. Each 
  #' element of the list contains a matrix with the 
  #' discretized values of the lagged covariance function
  #' estimated from the data. The covariance function is 
  #' given by a \eqn{(m x q)} matrix, where \eqn{m} is the
  #' number of points observed in curve \eqn{Y_{t}} and
  #' \eqn{q} is the number of points observed in curve 
  #' \eqn{X_{t}}.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' X <- simulate_iid_brownian_bridge(N, v, 0.5)
  #' nlags <- 1
  #' cov <- obtain_covariance(Y = Y,
  #'                          X = X,
  #'                          nlags = nlags)
  #' image(x = v, y = v, z = cov$Lag0)
  #' 
  #' \donttest{
  #' # Example 2
  #'
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' Y1 <- Y[1:N-1,]
  #' Y2 <- Y[2:N,]
  #' nlags <- 10
  #' cov <- obtain_covariance(Y = Y1, X = Y2, nlags = nlags)
  #' image(x = v, y = v, z = cov$Lag0, main = "Lag 0")
  #' image(x = v, y = v, z = cov$Lag1, main = "Lag 1")
  #' image(x = v, y = v, z = cov$Lag2, main = "Lag 2")
  #' # Example 3
  #'
  #' require(fields)
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 3
  #' Y1 <- Y[1:N-1,]
  #' Y2 <- Y[2:N,]
  #' cov <- obtain_covariance(Y = Y1,X = Y2, nlags = nlags)
  #' z_lims <- range(cov$Lag0)
  #' colors <- heat.colors(12)
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,5))
  #' par(oma=c( 0,0,0,6)) 
  #' for(k in 0:nlags){
  #'    image(x=v,
  #'          y=v,
  #'          z = cov[[paste0("Lag",k)]],
  #'          main = paste("Lag",k),
  #'          col = colors,
  #'          xlab = "u",
  #'          ylab = "v")
  #' }
  #' par(oma=c( 0,0,0,2.5)) # reset margin to be much smaller.
  #' image.plot( legend.only=TRUE, legend.width = 2,zlim=z_lims, col = colors)
  #' par(opar)
  #' 
  #' # Perspective Plots
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,4))
  #' for(k in 0:nlags){
  #'   persp(x=v,
  #'         y=v,
  #'         z = cov[[paste0("Lag",k)]],
  #'         main = paste("Lag",k),
  #'         xlab = "u",
  #'         ylab = "v",
  #'         theta = 330,
  #'         phi = 30,
  #'         ticktype = "detailed",
  #'         zlab = "")
  #' }
  #' par(opar)
  #' }
  #' @export obtain_covariance
  
  nt <- nrow(Y)
  nu <- ncol(Y)
  nv <- ncol(X)
  
  # Obtain mean function of Y
  mean.Y <- colMeans(Y)
  # Obtain mean function of X
  mean.X <- colMeans(X)
  
  # Compute covariance function
  fun.covariance <- list()
  for(kk in 0:nlags){
    fun.covariance[[paste("Lag",kk,sep = "")]] <- matrix(0,nrow = nu,ncol = nv)
    for(ind_curve in (1+kk):nt){
      fun.covariance[[paste("Lag",kk,sep = "")]] <- fun.covariance[[paste("Lag",kk,sep = "")]] + (as.numeric(Y[ind_curve,] - mean.Y)) %*% t(as.numeric(X[ind_curve - kk,] - mean.X))
    }
    fun.covariance[[paste("Lag",kk,sep = "")]] <- fun.covariance[[paste("Lag",kk,sep = "")]]/(nt-1)
  }
  return(fun.covariance)
}
