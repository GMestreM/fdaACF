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
  #' @references
  #' Mestre, G., Portela, J., Muñoz-San Roque, A., Alonso, E. (2020).
  #' \emph{Forecasting hourly supply curves in the Italian Day-Ahead 
  #' electricity market with a double-seasonal SARMAHX model.}
  #' International Journal of Electrical Power & Energy Systems, 
  #' 121, 1--13. \url{https://www.sciencedirect.com/science/article/pii/S0142061519337135}
  #' 
  #' Kokoszka, P., Rice, G., Shang, H.L. (2017).
  #' \emph{Inference for the autocovariance of a functional
  #'  time series under conditional heteroscedasticity}
  #' Journal of Multivariate Analysis, 
  #' 162, 32--50. \url{https://www.sciencedirect.com/science/article/pii/S0142061519337135}
  #' 
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
  #' nlags <- 30
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

obtain_FPACF <- function(Y,v,nlags,n_harm,ci=0.95,estimation = "MC",figure = TRUE,...){
  #' Obtain the partial autocorrelation function for a given FTS.
  #'
  #' @description Estimate the partial autocorrelation
  #' function for a given functional time series and its
  #' distribution under the hypothesis of strong functional
  #' white noise.
  #'
  #' @param Y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Discretization points of the curves.
  #' @param nlags Number of lagged covariance operators
  #' of the functional time series that will be used
  #' to estimate the partial autocorrelation function.
  #' @param n_harm Number of principal components
  #' that will be used to fit the ARH(p) models.
  #' @param ci A value between 0 and 1 that indicates
  #' the confidence interval for the i.i.d. bounds
  #' of the partial autocorrelation function. By default
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
  #' estimated partial autocorrelation function with the
  #' specified i.i.d. bound.
  #' @param ... Further arguments passed to the \code{plot_FACF}
  #' function.
  #' @return Return a list with:
  #' \itemize{
  #'     \item \code{Blueline}: The upper prediction
  #'     bound for the i.i.d. distribution.
  #'     \item \code{rho}: Partial autocorrelation
  #'     coefficients for
  #'     each lag of the functional time series.
  #' }
  #' 
  #' 
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 5)
  #' sig <- 2
  #' Y <- simulate_iid_brownian_bridge(N, v, sig)
  #' obtain_FPACF(Y,v,20, n_harm = 2)
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' data(elec_prices)
  #' v <- seq(from = 1, to = 24)
  #' nlags <- 30
  #' obtain_FPACF(Y = as.matrix(elec_prices), 
  #' v = v,
  #' nlags = nlags,
  #' n_harm = 5, 
  #' ci = 0.95,
  #' figure = TRUE)
  #' }
  #' @export obtain_FPACF
  
  y <- Y
  dv <- length(v)
  dt <- nrow(y)
  
  # Initialize FPACF vector
  FPACF <- rep(NA, nlags)
  
  FACF <- fdaACF::obtain_FACF(Y = y,
                              v = v,
                              nlags = 1,
                              ci = ci,
                              figure = F)
  
  Blueline <- FACF$Blueline
  FPACF[1] <- FACF$rho[1]
  
  vector_PACF <- FPACF
  
  # Start loop for fitting ARH(p-1)
  for(pp in 2:nlags){
    lag_PACF <- pp
    
    # 1 - Fit ARH(1) to the series
    Yest_ARIMA <- fit_ARHp_FPCA(y = y, v = v, p = lag_PACF-1, n_harm = n_harm)$y_est
    if(F){
      kkk <- 30
      graphics::plot(v,y[kkk,],type = "b")
      graphics::lines(v,Yest_ARIMA[kkk,],col = "red")
    }
    
    # 2 - Fit ARH(1) to the REVERSED series
    y_rev <- y[seq(from = nrow(y), to = 1, by = -1),]
    
    # y_rev[1,] == y[dt,]
    
    Yest_ARIMA_REV <- fit_ARHp_FPCA(y = y_rev, v = v, p = lag_PACF-1, n_harm = n_harm)$y_est
    
    # 3 - Estimate covariance surface for PACF 
    Yest_1 <- Yest_ARIMA;
    Yest_2 <- Yest_ARIMA_REV[seq(from = nrow(y), to = 1, by = -1),]
    
    res_filt_1 = y - Yest_1;
    res_filt_2 = y - Yest_2;
    
    # Cross-covariance surface
    sup_cov = matrix(0,dv,dv);
    ini_serie = max(which(is.na(Yest_1[,1])))+2 #max(find(isnan(Yest_1(:,1))))+1;
    fin_serie = min(which(is.na(Yest_2[,1])))-2 #min(find(isnan(Yest_2(:,1))))-1;
    
    count <- 0;
    for(jj in ini_serie:fin_serie){
      epsilon_1 <- as.matrix(res_filt_1[jj,])
      epsilon_2 <- as.matrix(res_filt_2[jj-lag_PACF,])
      
      sup_cov <- sup_cov + (epsilon_1 %*% t(epsilon_2))
      
      count <- count + 1;
    }
    
    sup_cov = (1/count)*sup_cov
    
    if(F){
      graphics::persp(v,v,sup_cov,theta = 330, phi = 20, 
                      main = paste0("Covariance surface PACF lag ",lag_PACF),
                      ticktype = "detailed")
    }
    
    # L2 norm of covariance surface
    if(F){
      sqrt(fdaACF::obtain_suface_L2_norm(v, list(Lag0 = sup_cov)))
    }
    
    
    # Estimate traces
    var_1 <- matrix(0,dv,dv);
    count <- 0;
    for (jj in ini_serie:dt){
      epsilon_1 = as.matrix(res_filt_1[jj,]);
      # epsilon_2 = as.matrix(res_filt_1[jj,]);
      
      var_1 <- var_1 + (epsilon_1 %*% t(epsilon_1))
      
      count <- count + 1
      
    }
    var_1 = (1/count)*var_1
    
    
    
    if(F){
      graphics::persp(v,v,var_1,theta = 330, phi = 20, ticktype = "detailed", main = "var 1") 
    }
    traza_1 <- pracma::trapz(v,diag(var_1));
    
    var_2 <- matrix(0,dv,dv);
    count <- 0;
    for (jj in 1:fin_serie){
      epsilon_1 <- as.matrix(res_filt_2[jj,]);
      epsilon_2 <- as.matrix(res_filt_2[jj,]);
      
      var_2 <- var_2 + (epsilon_1 %*% t(epsilon_2))
      
      count <- count + 1;
      
      
    }
    var_2 = (1/count)*var_2
    if(F){
      graphics::persp(v,v,var_2,theta = 330, phi = 20, ticktype = "detailed", main = "var 2") 
    }
    
    traza_2 = pracma::trapz(v,diag(var_2));
    
    
    sup_corr = sup_cov/(sqrt(traza_1)*sqrt(traza_2));
    
    # L2 norm of cross-correlation surface
    if(F){
      sqrt(fdaACF::obtain_suface_L2_norm(v, list(Lag0 = sup_corr)))    
    }
    
    vector_PACF[lag_PACF] <- sqrt(fdaACF::obtain_suface_L2_norm(v, list(Lag0 = sup_corr)))
  }
  
  if(figure){
    plot_FACF(vector_PACF,Blueline,ci,...)
  }
  
  # Prepare output
  output <- list(Blueline = Blueline,
                 rho = vector_PACF)
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



