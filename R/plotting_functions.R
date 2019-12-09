plot_autocovariance <- function(fun.autocovariance,lag = 0){

  #' Generate a 3D plot of the autocovariance surface of a given FTS
  #'
  #' Obtain a 3D plot of the autocovariance surfaces of a
  #' given functional time series. This visualization is
  #' useful to detec any kind of dependency between
  #' the discretization points of the series.
  #'
  #' @param fun.autocovariance A list obtained by
  #' calling the function \code{obtain_autocovariance}.
  #' @param lag An integer between 0 and \code{nlags}, indicating
  #' the lagged autocovariance function to be plotted.
  #' By default 0.
  #' @examples
  #' \dontrun{
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 4
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,nlags = nlags)
  #' set.panel()
  #' set.panel(1,5)
  #' for(k in 0:nlags){
  #'    plot_autocovariance(lagged_autocov,k)
  #' }
  #' set.panel()
  #' }
  # Color palette
  col.pal<-colorRampPalette(c("blue", "red"))
  colors<-col.pal(100)

  # Select the lagged autocovariance surface
  z <- fun.autocovariance[[paste("Lag",lag,sep="")]]
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)

  persp(z,theta = 360-45, phi = 30,
        col=colors[z.facet.range], shade = NA,
        ticktype='detailed',
        expand = 0.6,
        xlab = "u",
        ylab = "v",
        zlab = "",
        main = paste("Lag",lag,sep=" "))
}




plot_FACF <- function(rho,Blueline,ci,...){
  #' Plot the autocorrelation function of a given FTS
  #'
  #' Plot a visual representation of the autocorrelation function
  #' of a given functional time series, including the upper i.i.d.
  #' bound.
  #'
  #' @param rho Autocorrelation values for each lag of
  #' the functional time series obtained by calling the
  #' function \code{obtain_FACF}.
  #' @param Blueline The upper prediction bound for the
  #' i.i.d. distribution obtained by calling the
  #' function \code{obtain_FACF}.
  #' @param ci Value between 0 and 1 that was used
  #' when calling the function \code{obtain_FACF}.
  #' This value is only used to display information
  #' in the figure.
  #' @examples
  #' \dontrun{
  #' N <- 200
  #' v <- seq(from = 0, to = 1, length.out = 30)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 15
  #' upper_bound <- 0.95
  #' fACF <- obtain_FACF(Y = bbridge,v = v,nlags = nlags,ci=upper_bound,figure = FALSE)
  #' plot_FACF(rho = fACF$rho,Blueline = fACF$Blueline,ci = upper_bound)
  #' }
  plot(x = seq(1,length(rho),by = 1),
       y= rho,
       type="h",
       lwd = 2,
       xlab = "Lag",
       ylab = "Autocorrelation",
       ylim = c(0,max(rho)*1.5),...)
  abline(h = Blueline,col = "blue",lwd = 2, lty = 2)
  legend(x = "topleft",
         legend = c(paste("i.i.d. bound (",ci*100," % conf.)",sep="")),
         col = "blue",
         lty = 2,
         lwd = 2)
}
