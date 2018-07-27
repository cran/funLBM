plot.funLBM <-
function(x,type='blocks',...){
  colors = c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00",
             "#CAB2D6","#6A3D9A","#FFFF99","#B15928")
  palette(colors)
  X = eval.parent(x$call[[2]])
  type = tolower(as.character(type))
  allTypes = c("means","evolution","blocks","criterion","likelihood","proportions")
  type = match.arg(type, allTypes)
  old.par <- par(no.readonly = TRUE)
  if (type=='means'){
    T = x$T; K = x$K; L = x$L
    if (x$basisName == 'spline') basis <- create.bspline.basis(c(0,T),x$nbasis)
    else if (x$basisName == 'fourier') basis <- create.fourier.basis(c(0,T),x$nbasis)
    obj = list(basis = basis,coefs=c(),fdnames=list(time=1:T,reps=c(),values=c()))
    class(obj) = 'fd'
    par(mfrow=c(L,1))
    for (l in 1:L){
      obj$coefs = t(x$prms$mu[,l,])
      plot(obj,col=1:6,lwd=2,main=paste('Column group',l),...)
    }
  }
  else if (type=='evolution'){
    par(mfrow=c(1,2))
    matplot(x$allPrms$Alphas,type='l',xlab='Iterations',
            main=expression(alpha),lwd=2,...)
    matplot(x$allPrms$Betas,type='l',xlab='Iterations',
            main=expression(beta),lwd=2,...)
  }
  else if (type=='blocks'){
    T = x$T; K = x$K; L = x$L
    if (x$basisName == 'spline') basis <- create.bspline.basis(c(0,T),x$nbasis)
    else if (x$basisName == 'fourier') basis <- create.fourier.basis(c(0,T),x$nbasis)
    obj = list(basis = basis,coefs=c(),fdnames=list(time=1:T,reps=c(),values=c()))
    class(obj) = 'fd'
    #par(mfrow=c(K,L))
    op <- par(mfrow = c(K,L),mar=c(0.5,0.5,0.5,0.5),oma=c(1,1,1,1))
    coord = which(x$prms$mu == max(x$prms$mu),arr.ind=TRUE)
    mu.max = x$prms$mu[coord[1],coord[2],]
    coord = which(x$prms$mu == min(x$prms$mu),arr.ind=TRUE)
    mu.min = x$prms$mu[coord[1],coord[2],]
    for (k in 1:K)
      for (l in 1:L){
        obj$coefs = t(rbind(mu.min,mu.min))
        plot(obj,col=0,xaxt='n',yaxt='n',xlab='',ylab='',...)
        obj$coefs = replicate(2,x$prms$mu[k,l,])
        lines(obj,col=k*(K-1)+l,lwd=2,...)
      }
  }
  else if (type=='criterion'){
    stop('Not yet implemented!')
  }
  else if (type=='likelihood'){
    plot(x$lik,type='b',col=2,xlab='Iterations',ylab='Complete log-likelihood')
  }
  else if (type=='proportions'){
    par(mfrow=c(1,2))
    barplot(x$prms$alpha,col=1:x$K,main=expression(alpha),names.arg=1:x$K)
    barplot(x$prms$beta,col=1:x$L,main=expression(beta),names.arg=1:x$L)
  }
  par(old.par)
}
