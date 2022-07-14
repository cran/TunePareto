\name{print.TuneParetoResult}
\alias{print.TuneParetoResult}
\alias{print.TuneParetoClassifier}
\alias{print.TuneParetoModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print method for objects used in TunePareto
}
\description{
Customized printing methods for several objects used in TunePareto: For \code{TuneParetoResult} objects, the Pareto-optimal parameter configurations are printed. For \code{TuneParetoClassifier} and \code{TuneParetoModel} objects, information on the classifier and its parameters is printed.
}
\usage{
\method{print}{TuneParetoResult}(x, ...)
\method{print}{TuneParetoClassifier}(x, ...)
\method{print}{TuneParetoModel}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{TuneParetoResult}, \code{TuneParetoClassifier} or \code{TuneParetoModel} to be printed.
}

  \item{\dots}{Further parameters (currently unused).
}
}
\value{
Invisibly returns the printed object.
}

\seealso{
\code{\link{tunePareto}}, \code{\link{tuneParetoClassifier}}, \code{\link{trainTuneParetoClassifier}}
}
\keyword{print}
