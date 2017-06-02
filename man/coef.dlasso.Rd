\name{coef.dlasso}
\alias{coef.dlasso}
\title{
  Extract coefficients from a fitted dlasso model
}
\description{
While dlasso() produces the entire path of solutions, coef.dlasso extracts a particular point along the path corresponded to the minimum AIC, BIC, GIC or GCV.
}
\usage{
\method{coef}{dlasso}(object , ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
  dlasso object.
}
  \item{...}{
  Not working.
}
}


\value{
  A matrix of extracted coefficients.
}

\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}



\seealso{
\code{\link{dlasso}},\code{\link{plot.dlasso}}
}
\examples{
  example(dlasso)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ "dlasso" }% use one of  RShowDoc("KEYWORDS")