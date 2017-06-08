\name{dlasso}
\alias{dlasso}
\title{
An implementation of dlasso using iterative ridge algorithm
}
\description{
This function allows implementing differentiable lasso (dlasso) for arbitrary values of \eqn{\lambda} and s.
}
\usage{
dlasso (x,
       y,
       s  =  1           ,
       intercept = FALSE ,
       c  = 1            ,
       adp = TRUE        ,
       lambda = NULL     ,
       split  = 50       ,
       maxIter = 500     ,
       adj = 1.1         ,
       lowlambda = 10^-3 ,
       digit = 5         ,
       cauchy = FALSE    ,
       force = 'auto'    ,
       trace = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  Matrix of predictors
}
  \item{y}{
  Response vector
}
  \item{s}{
  A single or a vector of precision value, s, given adp=FALSE. Default is 1. See "adp" parameter.
}
  \item{intercept}{
  Logical flag. If TRUE, an intercept is included in the model, otherwise no intercept is included. Default is FALSE.
}
  \item{c}{
  Choose between dlasso (c=1) and dSCAD (c=-1). Default is dlasso. See futher "force" parameter.
}
  \item{adp}{
  Logical flag. TRUE to use adaptive adjustment for s. If TRUE then the function ignores the initial s.
  }
  \item{lambda}{
  Optional values for the tuning parameter. A single value or a sequence of values. Useful for manually searching   over user defined set of tuning values. Set to any negative value to activate the automatic setting for  \eqn{\lambda.max} and  \eqn{\lambda.min}.
  }
    \item{split}{
    The number of splits between \eqn{\lambda.min} and \eqn{\lambda.max}.
}
  \item{maxIter}{
  The maximum iterations for the algorithm. Default is 500.
}
  \item{adj}{
Positive value. This value adjusts the upper value for the penalty term, \eqn{adj*\lambda.max}.
}
  \item{lowlambda}{
    The lower limit for the tuning parameter. Default is 10^-3.
}
  \item{digit}{
  The maximum number of digits before setting an estimation to zero. The default is 5 digits.
}
  \item{cauchy}{
  Logical flag. Set to TRUE to use Cauchy CDF instead of Gaussian one in the penalty function. The default is Gaussian.
}
  \item{force}{
  Logical flag. Set to TRUE to let only a decrease in absolute estimation of the parameters. Default is 'auto' that is only applied if sqrt(n)>2*log(p) for n the number of observations and p the number of parameters.
  }
  \item{trace}{
    Logical flag. If TRUE, output contains some information about the steps. Default is FALSE.
}
}

\value{
  A "dlasso" object of the form of a matrix ( \eqn{\lambda} | s | AICc | GIC | BIC | GCV | estimations).
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}



\seealso{
\code{\link{coef.dlasso}},\code{\link{plot.dlasso}}
}
\examples{
    # dLASSO
    r = 5
    zr= 95
    n = 50
    b = c(1:r,rep(0,zr))
    x = matrix(rnorm((r+zr)*n),ncol=r+zr)
    y = x \%*\% b +rnorm(n)
    dLasso = dlasso(x=x,y=y,adp=TRUE)
    plot(dLasso,label=.1,cex=.80,all = 1)
    coef(dLasso)
    # dSCAD
    dscad = dlasso(x=x,y=y,c=-1,adj=50)
    plot(dscad,label=.1,cex=.80,all = 1)
    coef(dscad)

}
\keyword{ "dlasso" }% use one of  RShowDoc("KEYWORDS")