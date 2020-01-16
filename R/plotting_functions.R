plot_autocovariance <- function(fun.autocovariance,lag = 0,...){

  #' Generate a 3D plot of the autocovariance surface of a given FTS
  #'
  #' @description Obtain a 3D plot of the autocovariance surfaces of a
  #' given functional time series. This visualization is
  #' useful to detect any kind of dependency between
  #' the discretization points of the series.
  #'
  #' @param fun.autocovariance A list obtained by
  #' calling the function \code{obtain_autocovariance}.
  #' @param lag An integer between 0 and \code{nlags}, indicating
  #' the lagged autocovariance function to be plotted.
  #' By default 0.
  #' @param ... Further arguments passed to the  \code{persp}
  #' function.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 1
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,nlags = nlags)
  #' plot_autocovariance(lagged_autocov,1)
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' N <- 500
  #' v <- seq(from = 0, to = 1, length.out = 50)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 4
  #' lagged_autocov <- obtain_autocovariance(Y = bbridge,nlags = nlags)
  #' opar <- par(no.readonly = TRUE)
  #' par(mfrow = c(1,5))
  #' for(k in 0:nlags){
  #'    plot_autocovariance(lagged_autocov,k)
  #' }
  #' par(opar)
  #' }
  #' @export plot_autocovariance
  # Color palette
  col.pal<-grDevices::colorRampPalette(c("blue", "red"))
  colors<-col.pal(100)

  # Select the lagged autocovariance surface
  z <- fun.autocovariance[[paste("Lag",lag,sep="")]]
  z.facet.center <- (z[-1, -1] + z[-1, -ncol(z)] + z[-nrow(z), -1] + z[-nrow(z), -ncol(z)])/4
  z.facet.range<-cut(z.facet.center, 100)

  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))     arguments$xlab <- "u"
  if(!"ylab" %in% names(arguments))     arguments$ylab <- "v"
  if(!"zlab" %in% names(arguments))     arguments$zlab <- ""
  if(!"main" %in% names(arguments))     arguments$main  <- paste("Lag",lag,sep=" ")
  if(!"expand" %in% names(arguments))   arguments$expand  <- 0.7
  if(!"ticktype" %in% names(arguments)) arguments$ticktype  <-'detailed'
  if(!"shade" %in% names(arguments))    arguments$shade  <- NA
  if(!"col" %in% names(arguments))      arguments$col <- colors[z.facet.range]
  if(!"theta" %in% names(arguments))    arguments$theta  <- 315
  if(!"phi" %in% names(arguments))      arguments$phi  <- 30
  arguments$z = z

  do.call(graphics::persp,arguments)

}




plot_FACF <- function(rho,Blueline,ci,...){
  #' Plot the autocorrelation function of a given FTS
  #'
  #' @description Plot a visual representation of the autocorrelation function
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
  #' @param ... Further arguments passed to the  \code{plot}
  #' function.
  #' @examples
  #' # Example 1
  #' 
  #' N <- 100
  #' v <- seq(from = 0, to = 1, length.out = 10)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 15
  #' upper_bound <- 0.95
  #' fACF <- obtain_FACF(Y = bbridge,v = v,nlags = nlags,ci=upper_bound,figure = FALSE)
  #' plot_FACF(rho = fACF$rho,Blueline = fACF$Blueline,ci = upper_bound)
  #' 
  #' \donttest{
  #' # Example 2
  #' 
  #' N <- 200
  #' v <- seq(from = 0, to = 1, length.out = 30)
  #' sig <- 2
  #' bbridge <- simulate_iid_brownian_bridge(N, v, sig)
  #' nlags <- 15
  #' upper_bound <- 0.95
  #' fACF <- obtain_FACF(Y = bbridge,v = v,nlags = nlags,ci=upper_bound,figure = FALSE)
  #' plot_FACF(rho = fACF$rho,Blueline = fACF$Blueline,ci = upper_bound)
  #' }
  #' @export plot_FACF
  # Check if any additional plotting parameters are present
  arguments <- list(...)
  if(!"xlab" %in% names(arguments))  arguments$xlab <- "Lag"
  if(!"ylab" %in% names(arguments))  arguments$ylab <- "Autocorrelation"
  if(!"ylim" %in% names(arguments))  arguments$ylim <- c(0,max(rho)*1.5)
  if(!"lwd"  %in% names(arguments))   arguments$lwd  <- 2
  if(!"main" %in% names(arguments))  arguments$main  <- ""
  if(!"xlim" %in% names(arguments))  arguments$xlim  <- c(0,length(rho))
  arguments$x = seq(1,length(rho),by = 1)
  arguments$y = rho
  arguments$type = "h"
  
  do.call(graphics::plot,arguments)
  graphics::abline(h = Blueline,col = "blue",lwd = 2, lty = 2)
  graphics::legend(x = "topleft",
                   legend = c(paste("i.i.d. bound (",ci*100," % conf.)",sep="")),
                   col = "blue",
                   lty = 2,
                   lwd = 2)
}
