\name{recalculateParetoSet}
\alias{recalculateParetoSet}
\title{
Recalculate Pareto-optimal solutions
}
\description{
Recalculates the Pareto-optimal solutions in a \code{TuneParetoResult} according to the specified objectives only, and returns another \code{TuneParetoResult} object reduced to these objectives. This avoids time-consuming recalculations of objective values if only a subset of objectives should be considered for a previously evaluated set of parameter combinations.
}
\usage{
recalculateParetoSet(tuneParetoResult, 
                     objectives)
}

\arguments{
  \item{tuneParetoResult}{
The \code{TuneParetoResult} object containing the parameter configurations to be examined
}
  \item{objectives}{
A vector of objective function indices. The Pareto set is recalculated according to these objectives, i.e. omitting other objectives. If this argument is not supplied, all objectives are used, which usually returns a copy of the input.
}
}
\value{
Returns a reduced \code{TuneParetoResult} object. For more details on the object structure, refer to \code{\link{tunePareto}}.
}

\seealso{
\code{\link{tunePareto}}, \code{\link{mergeTuneParetoResults}}
}
\examples{
# optimize the 'cost' parameter of an SVM on
# the 'iris' data set
res <- tunePareto(classifier = tunePareto.svm(),
                  data = iris[, -ncol(iris)], 
                  labels = iris[, ncol(iris)],
                  cost=seq(0.01,0.1,0.01),
                  objectiveFunctions=list(cvWeightedError(10, 10),
                                          cvSensitivity(10, 10, caseClass="setosa"),
                                          cvSpecificity(10, 10, caseClass="setosa")))
print(res)

# select only specificity and sensitivity
print(recalculateParetoSet(res, 2:3))                                
}
\keyword{subset of objectives}
