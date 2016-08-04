fmax = function(x){
  if(x>= 0){
    r = max(x)
  }  else{
    r=min(x)
  }
  return(r)
}
lseq <- function(from=1, to=5, length.out=6 , adj = 1) {
  r=exp(seq(log(from)/adj, log(to), length.out = length.out))
  return(r)
}
kfold <- function(Nobs,K=5){
  rs <- runif(Nobs)
  id <- seq(Nobs)[order(rs)]
  k <- as.integer(Nobs*seq(1,K-1)/K)
  k <- matrix(c(0,rep(k,each=2),Nobs),ncol=2,byrow=TRUE)
  k[,1] <- k[,1]+1
  l <- lapply(seq.int(K),function(x,k,d)
    list(train=d[!(seq(d) %in% seq(k[x,1],k[x,2]))],
         test=d[seq(k[x,1],k[x,2])]),k=k,d=id)
  return(l)
}
erfa3<-function(x=1,s=1){
  x=x/s
  y=tanh(39*x/sqrt(4*pi)-111/2*atan(35*x/111/sqrt(pi)))
  return(y)
}


Sigma  =  function(b,s){
  r =  ( 2*pnorm(b/s,0,1/sqrt(2)) - 1 + 2* b/s * dnorm( b/s,0,1/sqrt(2)) )/b
  if( length(r) > 1 ){
    r = diag(as.vector(r))
  }
  return(r)
}

