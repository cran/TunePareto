\name{allCombinations}
\alias{allCombinations}
\title{
Build a list of all possible combinations of parameter values
}
\description{
Builds a list of all possible combinations of parameter values from supplied ranges of parameter values. That is, each of the specified values is combined with all specified values for other parameters. The resulting lists can be used in the \code{classifierParameterCombinations} and \code{predictorParameterCombinations} parameters of \code{\link{tunePareto}}.
}
\usage{
allCombinations(parameterRanges)
}

\arguments{
  \item{parameterRanges}{
A list of lists of parameter ranges. That is, each element of the list specifies the values of a single parameter to be tested and is named according to this parameter. It is also possible to set parameters to fixed values by specifying only one value. 
}
}
\value{
Returns a list of lists, where each of the inner lists represents one parameter combination and consists of named elements for the parameters.
}

\seealso{
\code{\link{tunePareto}}
}
\examples{

library(class)
# Combine only valid combinations of 'k' and 'l'
# for the k-NN classifier:
comb <- c(allCombinations(list(k=1,l=0)),
          allCombinations(list(k=3,l=0:2)),
          allCombinations(list(k=5,l=0:4)),
          allCombinations(list(k=7,l=0:6)))
print(comb)

print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.knn(),
                 parameterCombinations = comb,
                 objectiveFunctions = list(cvError(10, 10),
                                           reclassError())))
}

\keyword{parameter combinations}
