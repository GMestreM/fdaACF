mat2fd <- function(mat_obj, range_val = c(0,1),argvals = NULL){
  #' Obtain a fd object from a matrix
  #'
  #' @description This function returns a \code{fd} object
  #' obtained from the discretized functional observations
  #' contained in \code{mat_obj}.
  #'
  #' It is assumed that the functional observations
  #' contained in \code{mat_obj} are real observations,
  #' hence a poligonal base will be used to obtain
  #' the functional object.
  #'
  #' @param mat_obj A matrix that contains the
  #' discretized functional observations.
  #' @param range_val A numeric vector of length 2
  #' that contains the range of the observed functional
  #' data. By default \code{range_val = c(0,1)}.
  #' @param argvals Optinal argument that contains
  #' a strictly increasing vector
  #' of argument values at which line segments join
  #' to form a polygonal line. If \code{argvals = NULL},
  #' it is assumed a equidistant discretization vector.
  #' By default \code{argvals = NULL}.
  #' @return A fd object obtained from the functional
  #' observations in \code{mat_obj}.
  #'
  #' @examples
  #' # Example 1
  #' N <- 100
  #' dv <- 30
  #' v <- seq(from = 0, to = 1, length.out = dv)
  #' set.seed(150) # For replication
  #' mat_func_obs <- matrix(rnorm(N*dv), 
  #'                        nrow = N, 
  #'                        ncol = dv)
  #' fd_func_obs <- mat2fd(mat_obj = mat_func_obs,
  #'                       range_val = range(v),
  #'                       argvals = v)
  #' plot(fd_func_obs)
  #' @export mat2fd
  
  # Create functional basis
  num.discr <- ncol(mat_obj)
  if (is.null(argvals)){
    arg.vals  <- seq(from = range_val[1],
                     to = range_val[2],
                     length.out = num.discr)
  }else{
    arg.vals <- argvals
  }
  
  basis.obj <- fda::create.polygonal.basis(rangeval = range_val,
                                           argvals  = arg.vals)
  
  # Obtain fd object
  fd.obj <- fda::Data2fd(y = t(mat_obj),
                         basisobj = basis.obj,
                         argvals  = arg.vals)
  
  return(fd.obj)
}


reconstruct_fd_from_PCA <- function(pca_struct,scores,centerfns = T){
  #' Obtain the reconstructed curves after PCA
  #'
  #' @description This function reconstructs the
  #' functional curves from a score vector
  #' using the basis obtained after applying
  #' functional PCA. This allows the user to
  #' draw estimations from the joint density of
  #' the FPCA scores and reconstruct the curves
  #' for those new scores.
  #'
  #' @param pca_struct List obtained after calling
  #' function \code{pca.fd}.
  #' @param scores Numerical vector that contains the
  #' scores of the fPCA decomposition for one
  #' functional observation.
  #' @param centerfns Logical value specifying
  #' wheter the FPCA performed used \code{centerfns = T}
  #' or \code{centerfns = F}. By default
  #' \code{centerfns = T}.
  #'
  #' @return Returns a object of type \code{fd}
  #' that contains the reconstructed curve.
  #'
  #' @examples
  #' # Example 1
  #'
  #' # Simulate fd
  #' nobs <- 200
  #' dv <- 10
  #' basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=10)
  #' set.seed(5)
  #' C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
  #' fd_sim <- fda::fd(coef=t(C),basis)
  #'
  #' # Perform FPCA
  #' pca_struct <- fda::pca.fd(fd_sim,nharm = 6)
  #'
  #' # Reconstruct first curve
  #' fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
  #' plot(fd_sim[1])
  #' plot(fd_rec, add = TRUE, col = "red")
  #' legend("topright",
  #'        legend = c("Real Curve", "PCA Reconstructed"),
  #'        col = c("black","red"),
  #'        lty = 1)
  #'
  #' # Example 2 (Perfect reconstruction)
  #'
  #' # Simulate fd
  #' nobs <- 200
  #' dv <- 7
  #' basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=dv)
  #' set.seed(5)
  #' C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
  #' fd_sim <- fda::fd(coef=t(C),basis)
  #'
  #' # Perform FPCA
  #' pca_struct <- fda::pca.fd(fd_sim,nharm = dv)
  #'
  #' # Reconstruct first curve
  #' fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
  #' plot(fd_sim[1])
  #' plot(fd_rec, add = TRUE, col = "red")
  #' legend("topright",
  #'        legend = c("Real Curve", "PCA Reconstructed"),
  #'        col = c("black","red"),
  #'        lty = 1)
  #' @export reconstruct_fd_from_PCA
  
  if (centerfns){
    out_curve <- pca_struct$meanfd
    for(rr in 1:ncol(pca_struct$scores)){
      out_curve <- out_curve + pca_struct$harmonics[rr]*scores[rr]
    }
    return(out_curve)
  }else{
    out_curve <- pca_struct$harmonics[1]*scores[1]
    for(rr in 2:ncol(pca_struct$scores)){
      out_curve <- out_curve + pca_struct$harmonics[rr]*scores[rr]
    }
    return(out_curve)
  }
  
}

integral_operator <- function(operator_kernel, curve, v){
  #' Integral transformation of a curve using an integral operator
  #'
  #' @description Compute the integral transform of the
  #' curve \eqn{Y_i} with respect to a given integral
  #' operator \eqn{\Psi}. The transformation is given by
  #' \deqn{\Psi(Y_{i})(v) = \int \psi(u,v)Y_{i}(u)du}
  #'
  #' @param operator_kernel Matrix with the values
  #' of the kernel surface of the integral operator.
  #' The dimension of the matrix is \eqn{(g x m)}, where
  #' \eqn{g} is the number of discretization points of the
  #' input curve and \eqn{m} is the number of discretization
  #' points of the output curve.
  #' @param curve Vector containing the discretized values
  #' of a functional observation. The dimension of the
  #' matrix is \eqn{(1 x m)}, where \eqn{m} is the
  #' number of points observed in the curve.
  #' @param v Numerical vector specifying the
  #' discretization points of the curves.
  #' @return Returns a matrix the same size as
  #' \code{curve} with the transformed values.
  #' @examples
  #' # Example 1
  #'
  #' v <- seq(from = 0, to = 1, length.out = 20)
  #' set.seed(10)
  #' curve <- sin(v) + rnorm(length(v))
  #' operator_kernel <- 0.6*(v %*% t(v))
  #' hat_curve <- integral_operator(operator_kernel,curve,v)
  #'
  #' @export integral_operator
  
  # Initialize output
  yhat <- rep(0, times = ncol(operator_kernel))
  
  # Perform the integral transformation
  for(ind_v in 1:ncol(operator_kernel)){
    yhat[ind_v] <- pracma::trapz(y = operator_kernel[,ind_v]*curve,
                                 x = v)
  }
  return(yhat)
}


fit_ARHp_FPCA <- function(y, v, p, n_harm, show_varprop = T){
  #' Fit an ARH(p) to a given functional time series
  #' 
  #' @description Fit an \eqn{ARH(p)} model to a given functional
  #' time series. The fitted model is based on the model proposed
  #' in (Aue et al, 2015), first decomposing the original
  #' functional observations into a vector time series of \code{n_harm}
  #' FPCA scores, and then fitting a vector autoregressive
  #' model of order \eqn{p} (\eqn{VAR(p)}) to the time series
  #' of the scores. Once fitted, the Karhunen-Loève expansion
  #' is used to re-transform the fitted values into functional
  #' observations.
  #' 
  #' @param y Matrix containing the discretized values
  #' of the functional time series. The dimension of the
  #' matrix is \eqn{(n x m)}, where \eqn{n} is the
  #' number of curves and \eqn{m} is the number of points
  #' observed in each curve.
  #' @param v Numeric vector that contains the 
  #' discretization points of the curves.
  #' @param p Numeric value specifying the order 
  #' of the functional autoregressive 
  #' model to be fitted.
  #' @param n_harm Numeric value specifying the number
  #' of functional principal components to be used when fitting
  #' the \eqn{ARH(p)} model.
  #' @param show_varprop Logical. If \code{show_varprop = TRUE},
  #' a plot of the proportion of variance explained by the first
  #' \code{n_harm} functional principal components will be shown.
  #' By default \code{show_varprop = TRUE}.
  #' 
  #' @examples 
  #' # Example 1
  #' 
  #' # Simulate an ARH(1) process
  #' N <- 250
  #' dv <- 20
  #' v <- seq(from = 0, to = 1, length.out = 20)
  #' 
  #' phi <- 1.3 * ((v) %*% t(v))
  #' 
  #' persp(v,v,phi, 
  #'       ticktype = "detailed", 
  #'       main = "Integral operator")
  #' 
  #' set.seed(3)
  #' white_noise <-  simulate_iid_brownian_bridge(N, v = v)
  #' 
  #' y <- matrix(nrow = N, ncol = dv)
  #' y[1,] <- white_noise[1,]
  #' for(jj in 2:N){
  #'     y[jj,] <- white_noise[jj,];
  #'     
  #'     y[jj,]=y[jj,]+integral_operator(operator_kernel = phi,
  #'                                     v = v,
  #'                                     curve = y[jj-1,])
  #' }
  #' 
  #' # Fit an ARH(1) model
  #' mod <- fit_ARHp_FPCA(y = y,
  #'                      v = v,
  #'                      p = 1,
  #'                      n_harm = 5)
  #' 
  #' # Plot results
  #' plot(v, y[50,], type = "l", lty = 1, ylab = "")
  #' lines(v, mod$y_est[50,], col = "red")
  #' legend("bottomleft", legend = c("real","est"),
  #'        lty = 1, col = c(1,2))
  #' 
  #' @references Aue, A., Norinho, D. D., Hormann, S. (2015).
  #' \emph{On the Prediction of Stationary Functional 
  #' Time Series}
  #' Journal of the American Statistical Association, 
  #' 110, 378--392. \url{https://doi.org/10.1080/01621459.2014.909317}
  #' @export fit_ARHp_FPCA
  dt <- nrow(y)
  dv <- length(v)
  
  # Step 1: FPCA decomposition of the curves
  y_fd <- mat2fd(mat_obj = y,argvals = v, range_val = range(v))
  pca <- fda::pca.fd(fdobj = y_fd, nharm = n_harm)
  
  if(show_varprop){
    graphics::plot(1:n_harm,cumsum(pca$varprop),
                   type = "b",
                   pch = 20,
                   main = "% Variance explained by FPCA",
                   xlab = "Number of components",
                   ylab = "% Var. Expl")
  }
  
  y_scores <- pca$scores
  
  
  # Step 2: Fit a VAR(p) model to the scores series
  y_scores_aux <- as.data.frame(y_scores)
  names(y_scores_aux) <- paste0("y",1:ncol(y_scores))
  mod <- vars::VAR(y_scores_aux,p = p)
  if(show_varprop){
    summary(mod)
  }
  
  fitted_vals <- stats::fitted(mod)
  
  # Fill with NA
  fitted_vals <- rbind(matrix(NA,nrow = nrow(y_scores)-nrow(fitted_vals),
                              ncol = ncol(y_scores)),
                       fitted_vals)
  
  if(F){
    graphics::plot(y_scores[,1], type = "l")
    graphics::lines(fitted_vals[,1], col = "red")
  }
  
  
  # Step 3: reconstruct the curves
  y_rec <- matrix(nrow = dt, ncol = dv)
  for(ii in 1:nrow(y)){
    # fd_aux <- reconstruct_fd_from_PCA(pca_struct = pca,scores = y_scores[ii,],centerfns = T)
    fd_aux <- reconstruct_fd_from_PCA(pca_struct = pca,scores = fitted_vals[ii,],centerfns = T)
    
    fd_aux_mat <- fda::eval.fd(evalarg = v, fdobj = fd_aux)
    
    if(F){
      graphics::plot(v,y[ii,], type = "b")
      graphics::lines(v,fd_aux_mat,col = "red")
    }
    
    y_rec[ii,] <- fd_aux_mat
    
  }
  
  
  # Prepare output
  output <- list(y_est = y_rec,
                 mod = mod,
                 fpca = pca,
                 fitted_vals = fitted_vals,
                 y = y,
                 v = v)
  
  return(output)
}