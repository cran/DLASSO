\name{dlasso}
\alias{dlasso}
\title{
An implementation of dlasso using iterative ridge algorithm
}
\description{
This function allows implementing differentiable lasso (dlasso) for arbitrary values of \eqn{\lambda} and s.
}
\usage{
dlasso(x, y, s = 10^-5, intercept = FALSE,
       lambda = -1, split = 50,
       maxIter = 300, eps = 10^-9,
       adj = 1.05, lowlambda = 10^-3,
       digit = 4, trace = FALSE)
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
  dlasso precision. Default is 10^-5.
}
  \item{intercept}{
  Logical flag. If TRUE, an intercept is included in the model, otherwise no intercept is included. Default is FALSE.
}
  \item{eps}{
  The maximum error between two consecutive iterations before the loop breaks. Default is 10^-9.
}
  \item{maxIter}{
  The maximum iterations for the algorithm. Default is 300.
}
  \item{adj}{
Positive value. This value adjusts the upper value for the penalty term, \eqn{adj*\lambda.max}.
}
  \item{lambda}{
  Optional values for the tuning parameter. A single value or a sequence of values. Useful for manually searching over user defined set of tuning values. Set to any negative value to activate the automatic setting for  \eqn{\lambda.max} and  \eqn{\lambda.min}.
  }
  \item{lowlambda}{
    The lower limit for the tuning parameter. Default is 10^-3.
}
  \item{split}{
    The number of splits between \eqn{\lambda.min} and \eqn{\lambda.max}.
}
  \item{digit}{
  The maximum number of digits before setting an estimation to zero. The default is 4 digits.
}
  \item{trace}{
    Logical flag. If TRUE, output contains some information about the steps. Default is FALSE.
}
}

\value{
  A "dlasso" object of the form of a matrix ( \eqn{\lambda} | s | AIC | BIC | GCV | estimations).
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}



\seealso{
\code{\link{coef.dlasso}},\code{\link{plot.dlasso}}
}
\examples{
  b = c(1:5,rep(0,5))
  x = matrix(rnorm(10*50),ncol=10)
  y = x \%*\% b +rnorm(50)
  object = dlasso(x=x,y=y)
  plot(object,label=.1,cex=.80)
  coef(object)
}
\keyword{ "dlasso" }% use one of  RShowDoc("KEYWORDS")