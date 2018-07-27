funLBM.main <-
function(X,K,L,maxit=100,burn=50,basis.name='fourier',nbasis=15,gibbs.it=5,display=FALSE,init='kmeans',...){
  # This function implements a SEM-Gibbs algorithm for the LFBM
  
  N = dim(X)[1]; p = dim(X)[2]; T = dim(X)[3];
  
  # Computing functional coefs
  if (basis.name == 'spline') basis <- create.bspline.basis(c(0,T),nbasis)
  else if (basis.name == 'fourier') basis <- create.fourier.basis(c(0,T),nbasis)
  else stop('Unavailable basis functions!')
  nbasis = basis$nbasis
  Xcol = apply(X,3,cbind)
  obj = smooth.basis(1:T,t(Xcol),basis)$fd
  Coefs = t(obj$coefs)
  C = array(NA,c(N,p,nbasis))
  for(j in 1:p) C[,j,] = as.matrix(Coefs[((j-1)*N+1):(j*N),])
  
  # Parameter initialization
  alpha = rep(1/K,K)
  beta = rep(1/L,L)
  if (init=='funFEM'){
    if (display) cat('Initialization: FunFEM...\n')
    if (basis.name == 'spline') basisfunFEM <- create.bspline.basis(c(0,p*T),nbasis)
    else if (basis.name == 'fourier') basisfunFEM <- create.fourier.basis(c(0,p*T),nbasis)
    Z = t(apply(funFEM(smooth.basis(1:(p*T),apply(X,1,cbind),basisfunFEM)$fd,K,maxit=10)$P,1,function(x) (x>=max(x))+0))
  }
  else if (init=='kmeans'){
    if (display) cat('Initialization: kmeans...\n'); 
    cls = kmeans(t(apply(X,1,cbind)),K)$cluster
    Z = dummy(cls,K)
    }
  else {cat('Initialization: random...\n');Z = t(rmultinom(N,1,alpha))}
  Z = empty.class.check(Z)
  if (display){cat('Z:'); print(colSums(Z))}
  if (init=='funFEM'){
    if (basis.name == 'spline') basisfunFEM <- create.bspline.basis(c(0,N*T),nbasis)
    else if (basis.name == 'fourier') basisfunFEM <- create.fourier.basis(c(0,N*T),nbasis)
    W = t(apply(funFEM(smooth.basis(1:(N*T),apply(X,2,cbind),basisfunFEM)$fd,L,maxit=10)$P,1,function(x) (x>=max(x))+0))
  }
  else if (init=='kmeans'){cls = kmeans(t(apply(X,2,cbind)),L)$cluster; W = dummy(cls,L)}
  else {W = t(rmultinom(p,1,beta))}
  W = empty.class.check(W)
  if (display){cat('W:'); print(colSums(W))}
  Q = rep(list(list()),K)
  mu = array(NA,c(K,L,nbasis))
  for (l in 1:L) mu[,l,] = matrix(rep(colMeans(Coefs),K)+rnorm(K*nbasis,0,0.01),nrow=K,byrow = TRUE)
  a = b = d = matrix(NA,K,L)
  
  # Parameter storage
  lik = rep(NA,maxit)
  Alphas = matrix(NA,maxit,K); Betas = matrix(NA,maxit,L)
  Alphas[1,] = alpha; Betas[1,] = beta
  Zs = matrix(NA,N,maxit); Ws = matrix(NA,p,maxit)
  Zs[,1] = max.col(Z); Ws[,1] = max.col(W)
  
  for (it in 1:maxit){
    if (display) cat('.')
    # M step for a specific block
    for (k in 1:K){
      for (l in 1:L){
        if (sum(W[,l])==1 | sum(Z[,k])==1){x = C[Z[,k]==1,W[,l]==1,]}
        else {x = apply(C[Z[,k]==1,W[,l]==1,],3,cbind)}
        if (is.vector(x)) x = t(x)
        alpha[k] = sum(Z[,k]==1) / N
        alpha[alpha<0.05] = 0.05
        beta[l] = sum(W[,l]==1) / p
        beta[beta<0.05] = 0.05
        mu[k,l,] = colMeans(x)
        obj$coefs = t(x)
        dc = mypca.fd(obj,nharm = nbasis)
        dc$values[dc$values < 0] = min(dc$values[dc$values > 0]) 
        d[k,l] = cattell(dc$values)
        a[k,l] = mean(dc$values[1:d[k,l]])
        b[k,l] = mean(dc$values[(d[k,l]+1):nbasis])
        Q[[k]][[l]] = dc$U[,1:d[k,l]]
      }
    }
    
    # SE step
    P = array(NA,c(N,K,L))
    for (r in 1:gibbs.it){
      # Z part
      Pz = estep.Z(C,alpha,beta,mu,a,b,d,Q,W)
      Z = t(apply(Pz,1,function(x){rmultinom(1,1,x)}))
      Z = empty.class.check(Z,Pz)
      if (display){cat('Z:'); print(colSums(Z))}
      
      # W part
      Pw = estep.W(C,alpha,beta,mu,a,b,d,Q,Z)
      W = t(apply(Pw,1,function(x){rmultinom(1,1,x)}))
      W = empty.class.check(W,Pw)
      if (display){cat('W:'); print(colSums(W))}
    }
    
    # Computing complete likelihood and parameter storage
    lik[it] = compute.CompleteLikelihood(C,alpha,beta,mu,a,b,d,Q,Z,W)
    Alphas[it,] = alpha; Betas[it,] = beta
    Zs[,it] = max.col(Z); Ws[,it] = max.col(W)
    
    # Test for early ending
    if (burn > 5 & it > burn) if (sum(abs(diff(lik[it:(it-5)]))) < 1e-6){
      burn = it - 5
      break
    }
  }
  if (display) cat('\n')
  
  # Averaging and computing MAP parameters
  alpha = colMeans(Alphas[burn:it,])
  beta = colMeans(Betas[burn:it,])
  Z = dummy(apply(Zs[,burn:it],1,Mode),K)
  Z = empty.class.check(Z,Pz)
  W = dummy(apply(Ws[,burn:it],1,Mode),L)
  W = empty.class.check(W,Pw)
  for (k in 1:K){
    for (l in 1:L){
      if (sum(W[,l])==1 | sum(Z[,k])==1){x = C[Z[,k]==1,W[,l]==1,]}
      else {x = apply(C[Z[,k]==1,W[,l]==1,],3,cbind)}
      mu[k,l,] = colMeans(x)
      obj$coefs = t(x)
      dc = mypca.fd(obj,nharm = nbasis)
      dc$values[dc$values < 0] = min(dc$values[dc$values > 0]) 
      d[k,l] = cattell(dc$values)
      a[k,l] = mean(dc$values[1:d[k,l]])
      b[k,l] = mean(dc$values[(d[k,l]+1):nbasis])
      Q[[k]][[l]] = dc$U[,1:d[k,l]]
    }
  }
  
  # ICL-BIC criterion
  nu = L*(K*nbasis +  K-1 + K + K) + sum(d*(nbasis-(d+1)/2))
  crit = compute.CompleteLikelihood(C,alpha,beta,mu,a,b,d,Q,Z,W) - (K-1)/2*log(N) - (L-1)/2*log(p) - nu/2*log(N*p)
  
  # Return results
  prms = list(alpha=alpha,beta=beta,mu=mu,a=a,b=b,d=d,Q=Q)
  allPrms = list(Alphas=Alphas[1:it,],Betas=Betas[1:it,],Zs=Zs[,1:it],Ws=Ws[,1:it])
  out = list(basisName=basis.name,nbasis=nbasis,T=T,K=K,L=L,prms=prms,Z=Z,W=W,
             row_clust=max.col(Z),col_clust=max.col(W),allPrms=allPrms,loglik=lik,icl=crit)
  class(out) = "funLBM"
  out
}
