\name{rankByDesirability}
\alias{rankByDesirability}

\title{
Rank results according to their desirabilities
}
\description{
Calculates the desirability index for each Pareto-optimal combination (or for all combinations), and ranks the combinations according to this value. The desirability index was introduced by Harrington in 1965 for multicriteria optimization. Desirability functions specify the desired values of each objective and are aggregated in a single desirability index.
}
\usage{
rankByDesirability(tuneParetoResult, 
                   desirabilityIndex, 
                   optimalOnly = TRUE)
}

\arguments{
  \item{tuneParetoResult}{
A \code{TuneParetoResult} object containing the parameter configurations to be examined
}
  \item{desirabilityIndex}{
A function accepting a vector of objective values and returning a desirability index in [0,1].
}
  \item{optimalOnly}{
If set to true, only the Pareto-optimal solutions are ranked. Otherwise, all tested solutions are ranked. Defaults to TRUE.
}
}
\value{
A matrix of objective values with an additional column for the desirability index. The rows of the matrix are sorted according to the index.
}


\examples{
\donttest{


# optimize the 'cost' parameter of an SVM on
# the 'iris' data set
res <- tunePareto(classifier = tunePareto.svm(),
                  data = iris[, -ncol(iris)], 
                  labels = iris[, ncol(iris)],
                  cost=c(0.01,0.05,0.1,0.5,1,5,10,50,100),
                  objectiveFunctions=list(cvWeightedError(10, 10),
                                          cvSensitivity(10, 10, caseClass="setosa"),
                                          cvSpecificity(10, 10, caseClass="setosa")))
             
# create desirability functions 
# aggregate functions in desirability index (e.g. harrington desirability function)
# here, for the sake of simplicity a random number generator
di <- function(x) {runif(1)}

# rank all tuning results according to their desirabilities
print(rankByDesirability(res,di,optimalOnly=FALSE))
}
}

\keyword{desirability function}
\keyword{desirability index}
