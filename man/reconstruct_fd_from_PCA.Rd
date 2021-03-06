% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misc.R
\name{reconstruct_fd_from_PCA}
\alias{reconstruct_fd_from_PCA}
\title{Obtain the reconstructed curves after PCA}
\usage{
reconstruct_fd_from_PCA(pca_struct, scores, centerfns = T)
}
\arguments{
\item{pca_struct}{List obtained after calling
function \code{pca.fd}.}

\item{scores}{Numerical vector that contains the
scores of the fPCA decomposition for one
functional observation.}

\item{centerfns}{Logical value specifying
wheter the FPCA performed used \code{centerfns = T}
or \code{centerfns = F}. By default
\code{centerfns = T}.}
}
\value{
Returns a object of type \code{fd}
that contains the reconstructed curve.
}
\description{
This function reconstructs the
functional curves from a score vector
using the basis obtained after applying
functional PCA. This allows the user to
draw estimations from the joint density of
the FPCA scores and reconstruct the curves
for those new scores.
}
\examples{
# Example 1

# Simulate fd
nobs <- 200
dv <- 10
basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=10)
set.seed(5)
C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
fd_sim <- fda::fd(coef=t(C),basis)

# Perform FPCA
pca_struct <- fda::pca.fd(fd_sim,nharm = 6)

# Reconstruct first curve
fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
plot(fd_sim[1])
plot(fd_rec, add = TRUE, col = "red")
legend("topright",
       legend = c("Real Curve", "PCA Reconstructed"),
       col = c("black","red"),
       lty = 1)

# Example 2 (Perfect reconstruction)

# Simulate fd
nobs <- 200
dv <- 7
basis<-fda::create.bspline.basis(rangeval=c(0,1),nbasis=dv)
set.seed(5)
C <- matrix(rnorm(nobs*dv), ncol =  dv, nrow = nobs)
fd_sim <- fda::fd(coef=t(C),basis)

# Perform FPCA
pca_struct <- fda::pca.fd(fd_sim,nharm = dv)

# Reconstruct first curve
fd_rec <- reconstruct_fd_from_PCA(pca_struct = pca_struct, scores = pca_struct$scores[1,])
plot(fd_sim[1])
plot(fd_rec, add = TRUE, col = "red")
legend("topright",
       legend = c("Real Curve", "PCA Reconstructed"),
       col = c("black","red"),
       lty = 1)
}
