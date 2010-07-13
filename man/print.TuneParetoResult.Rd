\name{print.TuneParetoResult}
\alias{print.TuneParetoResult}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Print function for the result of a parameter tuning
}
\description{
Prints the non-dominated parameter configurations determined in a call of \code{tunePareto}.
}
\usage{
print.TuneParetoResult(x, ...)
}
\arguments{
  \item{x}{An object of class \code{TuneParetoResult} to be printed.
}

  \item{\dots}{Further parameters (currently unused).
}
}
\value{
Invisibly returns the printed object.
}

\seealso{
\code{\link{tunePareto}}
}
\examples{
# call tunePareto.knn with different k and l on the 'iris' data set
x <- tunePareto.knn(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 k = c(5,7,9),
                 l = c(1,2,3),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa")))
print(x)                 
}
\keyword{print}
