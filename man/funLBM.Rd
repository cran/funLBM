\name{funLBM}
\alias{funLBM}
\title{
The functional latent block model
}
\description{
The funLBM algorithm allows to simultaneously cluster the rows and the columns of a data matrix where each entry of the matrix is a function or a time series.
}
\usage{
funLBM(X, K, L, maxit = 50, burn = 25, basis.name = "fourier", nbasis = 15,
        nbinit = 1, gibbs.it = 3, display = FALSE, init = "kmeans", mc.cores = 1, ...)
}
\arguments{
  \item{X}{
The data array (n x p x T) where each entry corresponds to the measure of one individual i, i=1,...,n, for a functional variable j, j=1,...,p, at point t, t=1,...,T,
}
  \item{K}{
The number of row clusters,
}
  \item{L}{
The number of column clusters,
}
  \item{maxit}{
The maximum number of iterations of the SEM-Gibbs algorithm (default is 100),
}
  \item{burn}{
The number of of iterations of the burn-in period (default is 50),
}
  \item{basis.name}{
The name ('fourier' or 'spline') of the basis functions used for the decomposition of the functions (default is 'fourier'),
}
  \item{nbasis}{
Number of the basis functions used for the decomposition of the functions (default is 15),
}
  \item{nbinit}{
Number of initializations (default is 3),
}
  \item{gibbs.it}{
Number of Gibbs iterations (default is 3),
}
  \item{display}{
Binary value. If TRUE, information about the iterations is displayed,
}
  \item{init}{
The type of initialization: 'random', 'kmeans' or 'funFEM'. Default is 'kmeans',
}
  \item{mc.cores}{
The number of cores for parallel computing (default is 1),
}
  \item{\dots}{
Additional parameters to provide to sub-functions.
}
}
\value{
The resulting object contains, in addition to call information: 
\item{prms}{A list containing all fited parameters for the best model (according to ICL)}
\item{Z}{The dummy matrix of row clustering}
\item{W}{The dummy matrix of column clustering}
\item{row_clust}{The group memberships of rows}
\item{col_clust}{The group memberships of columns}
\item{allPrms}{A list containing the fited parameters for all tested models}
\item{loglik}{The log-likelihood of the best model}
\item{icl}{The value of ICL for the best model}
}
\references{
C. Bouveyron, L. Bozzi, J. Jacques and F.-X. Jollois, The Functional Latent Block Model for the Co-Clustering of Electricity Consumption Curves, Journal of the Royal Statistical Society, Series C, 2018 (https://doi.org/10.1111/rssc.12260).
}
\examples{
## Co-clustering on simulated data with K=4 and L=3
set.seed(12345)
X = simulateData(n = 30, p = 30, t = 15)
out = funLBM(X$data,K=4,L=3)

# Visualization of results
plot(out,type='blocks')
plot(out,type='proportions')
plot(out,type='means')

# Evaluating clustering results
ari(out$col_clust,X$col_clust)
ari(out$row_clust,X$row_clust)

\donttest{
## Co-clustering on simulated data with parallel model selection
## To enable parallel computation (on Linux and Mac systems), 
## set mc.cores to the number of available CPUs
## and increase nbinit to allow the best possible fit
set.seed(12345)
X = simulateData(n = 30, p = 30, t = 15)
out = funLBM(X$data,K=2:4,L=2:4,init = "funFEM",mc.cores = 1,nbinit = 1)

# Evaluating clustering results
ari(out$col_clust,X$col_clust)
ari(out$row_clust,X$row_clust)

## Co-clustering of Velib data
data(Velib)
out = funLBM(Velib$data,K=4,L=2)

# Visualization of results
plot(out,type='blocks')
plot(out,type='proportions')
plot(out,type='means')
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Clustering}
\keyword{Functional data}
