\name{TunePareto-package}
\alias{TunePareto-package}
\alias{TunePareto}
\docType{package}
\title{
Multi-objective parameter tuning for classifiers
}
\description{
Generic methods for parameter tuning of classification algorithms using multiple scoring functions
}
\details{
\tabular{ll}{
Package: \tab TunePareto\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-07-13\cr
License: \tab Artistic-2.0 \cr
LazyLoad: \tab yes\cr
}
The methods of this package allow to assess the performance of classifiers with respect to certain parameter values and multiple scoring functions, such as the cross-validation error or the sensitivity. It provides a generic function \code{\link{tunePareto}} which can be configured to run most common classification methods implemented in R. The method tests combinations of the supplied classifier parameters according to the supplied scoring functions and calculates the Pareto front of optimal parameter configurations. The Pareto fronts cal be visualized using \code{\link{plotDominationGraph}}.


A number of predefined scoring functions are provided (see \code{\link{predefinedObjectiveFunctions}}), but the user is free to implement own scores (see \code{\link{createObjective}}).
For certain classifiers, the package includes specialized versions of \code{tunePareto} that are easier to use (see \code{\link{tunePareto.knn}}, \code{\link{tunePareto.tree}}, \code{\link{tunePareto.svm}}). 
}
\author{
Christoph Muessel, Hans Kestler

Maintainer: Hans Kestler <hans.kestler@uni-ulm.de>
}

\keyword{package
         multi-objective parameter tuning}

\examples{

# optimize the 'cost' and 'kernel' parameters according
# to CV error and CV Specificity on the 'iris' data set
r <- tunePareto.svm(data = iris[, -ncol(iris)], 
                   labels = iris[, ncol(iris)],
                   cost=c(0.001,0.01,0.1,1,10), 
                   kernel=c("linear", "polynomial", "radial", "sigmoid"),
                   objectiveFunctions=list(cvError(10, 10), 
                                           cvSpecificity(10, 10, caseClass="setosa")))

# the same example using parallel execution with snowfall
\donttest{library(snowfall)}
\donttest{sfInit(parallel=TRUE, cpus=2, type="SOCK")}
\donttest{r <- tunePareto.svm(data = iris[, -ncol(iris)], 
                    labels = iris[, ncol(iris)],
                    cost=c(0.001,0.01,0.1,1,10), 
                    kernel=c("linear", "polynomial", "radial", "sigmoid"),
                    objectiveFunctions=list(cvError(10, 10),
                                            cvSpecificity(10, 10, caseClass="setosa")),
                    useSnowfall=TRUE)}
\donttest{sfStop()}

# print Pareto-optimal solutions
print(r)

# plot the Pareto fronts
\donttest{plotDominationGraph(r, legend.x="topright")}
}
