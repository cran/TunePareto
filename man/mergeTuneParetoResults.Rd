\name{mergeTuneParetoResults}
\alias{mergeTuneParetoResults}
\title{
Calculate optimal solutions from several calls of tunePareto
}
\description{
Merges the results of multiple \code{TuneParetoResult} objects as returned by \code{\link{tunePareto}}, and recalculates the optimal solutions for the merged solution set. All supplied \code{TuneParetoResult} objects must use the same objective functions.
}
\usage{
mergeTuneParetoResults(...)
}
\arguments{
  \item{\dots}{
A set of \code{TuneParetoResult} objects to be merged.
}
}

\value{
A \code{TuneParetoResult} object containing the parameter configurations of all objects in the \code{\dots} argument and selecting the Pareto-optimal solutions among all these configurations. For more details on the object structure, refer to \code{\link{tunePareto}}.
}

\seealso{
\code{\link{tunePareto}}, \code{\link{recalculateParetoSet}}
}
\examples{
# optimize an SVM with small costs on 
# the 'iris' data set
r1 <- tunePareto(classifier = tunePareto.svm(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 cost=seq(0.01,0.1,0.01),
                 objectiveFunctions=list(cvWeightedError(10, 10),
                                         cvSensitivity(10, 10, caseClass="setosa")))
print(r1)
                                         
# another call to tunePareto with higher costs
r2 <- tunePareto(classifier = tunePareto.svm(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 cost=seq(0.5,10,0.5),
                 objectiveFunctions=list(cvWeightedError(10, 10),
                                         cvSensitivity(10, 10, caseClass="setosa")))
print(r2)

# merge the results
print(mergeTuneParetoResults(r1,r2))
}

\keyword{merge results}
