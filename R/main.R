dlasso  <-  function(x, y, s  =  10^-5,
                     intercept = FALSE ,
                     lambda=-1         ,
                     split  = 50       ,
                     maxIter = 300     ,
                     eps = 10^-9       ,
                     adj = 1.05        ,
                     lowlambda=10^-3   ,
                     digit =4          ,
                     trace = FALSE     )
{
  if (intercept)
    x = cbind(1,x)
  n      =  length(y)
  ncx    =  ncol(x)
  if(all(lambda<=0)){
    lmax   =  max(abs(cov(y,x)))*adj
    lambda =  lseq( lowlambda , lmax  , length.out  =  split)
  }else{
    split=length(lambda)
  }
  ls     =  length(s)
  b      =  MASS::ginv( t(x) %*% x) %*% t(x) %*% y
  outputMat  =  matrix(0, ncol  =  ncx+5 , nrow =  split*ls)
  counter = 1
  if (trace)
    cat('\r\n\r\n Ind|Max \t lambda \t s \t BIC \t\t GCV \t Converge \t Error')

  for(ss in s){
    active     =  as.vector(abs(b)>10^-digit)
    tmp = b*0
    for(rep in seq(split) ){
      l          =  lambda[rep]
      iter = err = 1
      while(any(active)  && (err    >  eps)  ){
        newx       =  x[,active,drop=FALSE]
        inV        =  t(newx) %*% newx + n * l * Sigma( b  =  b[active] , s  =  ss )
        xInv       =  chol2inv(chol(inV))
        b[active]  =  xInv %*% t(newx) %*% y
        #################
        active     =  as.vector(abs(b)>10^-digit)
        b[!active] =  0
        err        =  ncx*max(abs(tmp-b))
        tmp        =  b
        iter       =  iter + 1
        if(iter >=  maxIter) break
      }

      e         =  y - x %*% b
      df        =  sum(diag(newx %*% xInv %*% t(newx)))

      ################
      if(ncx > n){
        core    =  sum(dnorm(e, mean(e),sd(y),log = TRUE))
      }else{
        core    =  sum(dnorm(e, mean(e),sd(e),log = TRUE))
      }
      denGCV    = (1 - df/n)^2
      aic       = -2*core +      2*df
      bic       = -2*core + log(n)*df
      gcv       = (-2/n*core)/denGCV
      #gcv       = 1/n * sum(e^2)/denGCV

      if (trace){
        cat('\r\n',
            format(counter,digits=4),
            '|',
            format(split*ls,digits = 3,nsmall = 0),
            '\t ',
            format(l,digits = 4,nsmall = 2),
            '\t',
            format(ss,digits = 4,nsmall = 2),
            '\t',
            format(bic,digits = 4,nsmall = 2),
            '\t',
            format(gcv,digits = 4,nsmall = 2),
            '\t',
            format(iter,digits = 4,nsmall = 0),
            ' \t ',
            format(err,digits = 4,nsmall = 2)
        )
      }else{
        cat('\r ',counter,'|',split*ls,rep(' ',20))
      }
      outputMat[counter,]  =  c(l,  ss, aic, bic, gcv, b)

      colnames(outputMat)  =  c('lambda','s','AIC','BIC','GCV',colnames(x,do.NULL = FALSE,prefix = 'X.'))
      counter              = counter + 1
    }
  }
  output = outputMat
  class(output) = 'dlasso'
  return(output)
}

plot.dlasso <- function (x,label=FALSE,cex.lab=1,all=TRUE,...){
  #par(mfrow=c(3,1))
  object = x
  label=abs(label)
  if (all){
    plot(object[,1,drop=FALSE],object[,3,drop=FALSE],type = 'l',ylab='AIC',xlab='Lambda')
    abline(v=object[,1,drop=FALSE],col='gray87',lty=2)
    abline(v= object[mean(which(object[,3,drop=FALSE] == min(object[,3,drop=FALSE]))),1] ,col=2,lty=2,lwd=2)

    plot(object[,1,drop=FALSE],object[,4,drop=FALSE],type = 'l',ylab='BIC',xlab='Lambda')
    abline(v=object[,1,drop=FALSE],col='gray87',lty=2)
    abline(v= object[mean(which(object[,4,drop=FALSE] == min(object[,4,drop=FALSE]))),1] ,col=3,lty=3,lwd=2)

    plot(object[,1,drop=FALSE],object[,5,drop=FALSE],type = 'l',ylab='GCV',xlab='Lambda')
    abline(v=object[,1,drop=FALSE],col='gray87',lty=2)
    abline(v= object[mean(which(object[,5,drop=FALSE] == min(object[,5,drop=FALSE]))),1] ,col=4,lty=4,lwd=2)
    #par(mfrow=c(1,1))
  }
  coef1=(object[,-c(1:5)])
  s1 <- apply(abs(coef1), 1, sum)
  if(max(s1) != 0)
    ms1 = max(s1)
  else
    ms1 = 1
  s1=s1/ms1
  matplot(s1,coef1,type = 'l',
          ylab='Coefficients',xlab='|beta|/max|beta|',
          xlim=c(0,1+label) ,
          lwd=1,col=c(1:(ncol(object)-5)),
          ...)
  abline(h=0,lty=1,col='gray')
  if(label)
    text(x=1+label,y=object[1,-c(1:5)],labels = colnames(object)[-c(1:5)],cex=cex.lab)
  abline(v= s1[which(object[,3,drop=FALSE] == min(object[,3,drop=FALSE]))],col=2,lty=2,lwd=2)
  abline(v= s1[which(object[,4,drop=FALSE] == min(object[,4,drop=FALSE]))],col=3,lty=3,lwd=2)
  abline(v= s1[which(object[,5,drop=FALSE] == min(object[,5,drop=FALSE]))],col=4,lty=4,lwd=2)
  legend('topleft',legend = c('min.AIC','min.BIC','min.GCV'),lty=2:4,col=2:4,lwd=2)
  invisible()
}

coef.dlasso <-function(object , ...){
  B_AIC  = object[mean(which(object[,3,drop=FALSE] == min(object[,3,drop=FALSE])))+0  , -c(1:5)]
  B_BIC  = object[mean(which(object[,4,drop=FALSE] == min(object[,4,drop=FALSE])))+0  , -c(1:5)]
  B_GCV  = object[mean(which(object[,5,drop=FALSE] == min(object[,5,drop=FALSE])))+0  , -c(1:5)]
  print(apply(object[,3:5],2,function(x){min(x)[1]}))
  result = cbind(coef.AIC=B_AIC,coef.BIC=B_BIC,coef.GCV=B_GCV)

  return(result)
}