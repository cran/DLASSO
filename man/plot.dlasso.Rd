\name{plot.dlasso}
\alias{plot.dlasso}
\title{
  Plot method for dlasso objects
}
\description{
  Produce a plot from a dlasso fit.
}
\usage{
\method{plot}{dlasso}(x, label=FALSE , cex.lab=1 , all=TRUE ,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
  dlasso object
}
  \item{label}{
  Show covariate labels on the right hand side of the plot.
}
  \item{cex.lab}{
  See "cex.lab"" in "par" function.
}
  \item{all}{
  Logical flag. If TRUE shows all plots including AIC, BIC, CGV and complete coefficient path. Otherwise the output is a single plot for coefficient path. Default is TRUE.
}
  \item{\dots}{
  Additional arguments for generic plot.
}
}
\author{
  Hamed Haselimashhadi <hamedhaseli@gmail.com>
}

\seealso{
\code{\link{dlasso}},\code{\link{coef.dlasso}}
}
\examples{
  example(dlasso)
}

