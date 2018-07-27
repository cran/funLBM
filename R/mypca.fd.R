mypca.fd <-
function(fdobj, nharm = 2, harmfdPar=fdPar(fdobj), centerfns = TRUE){
  #  Carry out a functional PCA with regularization
  #  Arguments:
  #  FDOBJ      ... Functional data object
  #  NHARM     ... Number of principal components or harmonics to be kept
  #  HARMFDPAR ... Functional parameter object for the harmonics
  #  CENTERFNS ... If TRUE, the mean function is first subtracted from each function
  #
  #  Returns:  An object PCAFD of class "pca.fd" with these named entries:
  #  harmonics  ... A functional data object for the harmonics or eigenfunctions
  #  values     ... The complete set of eigenvalues
  #  scores     ... A matrix of scores on the principal components or harmonics
  #  varprop    ... A vector giving the proportion of variance explained
  #                 by each eigenfunction
  #  meanfd     ... A functional data object giving the mean function
  Ti = rep(1,ncol(fdobj$coefs))
  
  #  Check FDOBJ
  if (!(inherits(fdobj, "fd"))) stop(
    "Argument FD  not a functional data object.")
  
  #  compute mean function and center if required
  # browser()
  meanfd <- mean.fd(fdobj)
  # if (centerfns) fdobj <- center.fd(fdobj)
  
  if (centerfns){
    coefmean <- apply(t(as.matrix(Ti) %*% matrix(1,1,nrow(fdobj$coefs))) * fdobj$coefs, 1, sum) / sum(Ti)
    fdobj$coefs <- sweep(fdobj$coefs, 1, coefmean)
    meanfd$coefs = as.matrix(data.frame(mean=coefmean))
  }
  
  #  get coefficient matrix and its dimensions
  coef  <- fdobj$coefs
  coefd <- dim(coef)
  ndim  <- length(coefd)
  nrep  <- coefd[2]
  coefnames <- dimnames(coef)
  if (nrep < 2) stop("PCA not possible without replications.")
  
  basisobj <- fdobj$basis
  nbasis   <- basisobj$nbasis
  type     <- basisobj$type
  
  #  set up HARMBASIS
  #  currently this is required to be BASISOBJ
  harmbasis <- basisobj
  
  #  set up LFDOBJ and LAMBDA
  Lfdobj <- harmfdPar$Lfd
  lambda <- harmfdPar$lambda
  
  #  compute CTEMP whose cross product is needed
  ctemp <- coef
  
  #  set up cross product and penalty matrices
  #   Cmat <- crossprod(t(ctemp))/nrep
  Cmat = (Ti * ctemp) %*% t(ctemp) / nrep
  
  Jmat <- eval.penalty(basisobj, 0)
  if(lambda > 0) {
    Kmat <- eval.penalty(basisobj, Lfdobj)
    Wmat <- Jmat + lambda * Kmat
  } else {    Wmat <- Jmat  }
  Wmat <- (Wmat + t(Wmat))/2
  
  #  compute the Choleski factor of Wmat
  Lmat    <- chol(Wmat)
  Lmat.inv <- solve(Lmat)
  
  #  set up matrix for eigenanalysis
  if(lambda > 0) { Cmat <- t(Lmat.inv) %*% Jmat %*% Cmat %*% Jmat %*% Lmat.inv  }
  else { Cmat <- Lmat %*% Cmat %*% t(Lmat)  }
  
  #  eigenalysis
  Cmat    <- (Cmat + t(Cmat))/2
  result  <- eigen(Cmat)
  eigvalc <- result$values
  eigvecc <- as.matrix(result$vectors[, 1:nharm])
  sumvecc <- apply(eigvecc, 2, sum)
  eigvecc[,sumvecc < 0] <-  - eigvecc[, sumvecc < 0]
  
  varprop <- eigvalc[1:nharm]/sum(eigvalc)
  
  
  harmcoef <- Lmat.inv %*% eigvecc
  U = t(Lmat) %*% eigvecc
  harmscr  <- t(ctemp) %*% U
  
  harmnames <- rep("", nharm)
  for(i in 1:nharm)
    harmnames[i] <- paste("PC", i, sep = "")
  harmnames <- list(coefnames[[1]], harmnames,"values")
  harmfd   <- fd(harmcoef, basisobj, harmnames)
  
  pcafd  <- list(harmonics=harmfd,values=eigvalc,scores=harmscr,U=U,varprop=varprop,meanfd=meanfd,W=Wmat)
  class(pcafd) <- "pca.fd"
  return(pcafd)
}
