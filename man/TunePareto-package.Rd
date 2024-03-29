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

The methods of this package allow to assess the performance of classifiers with respect to certain parameter values and multiple scoring functions, such as the cross-validation error or the sensitivity. It provides the \code{\link{tunePareto}} function which can be configured to run most common classification methods implemented in R. Several sampling strategies for parameters are supplied, including Latin Hypercube sampling, quasi-random sequences, and evolutionary algorithms.

Classifiers are wrapped in generic \code{TuneParetoClassifier} objects which can be created using \code{\link{tuneParetoClassifier}}. For state-of-the-art classifiers, the package includes the corresponding wrapper objects (see \code{\link{tunePareto.knn}}, \code{\link{tunePareto.tree}}, \code{\link{tunePareto.randomForest}}, \code{\link{tunePareto.svm}}, \code{\link{tunePareto.NaiveBayes}}). 

The method tests combinations of the supplied classifier parameters according to the supplied scoring functions and calculates the Pareto front of optimal parameter configurations. The Pareto fronts can be visualized using \code{\link{plotDominationGraph}}, \code{\link{plotParetoFronts2D}} and \code{\link{plotObjectivePairs}}.

A number of predefined scoring functions are provided (see \code{\link{predefinedObjectiveFunctions}}), but the user is free to implement own scores (see \code{\link{createObjective}}).
}
\author{
Christoph Müssel, Ludwig Lausser, Hans Kestler

Maintainer: Hans Kestler <hans.kestler@uni-ulm.de>
}

\references{
Christoph Müssel, Ludwig Lausser, Markus Maucher, Hans A. Kestler (2012). Multi-Objective Parameter Selection for Classifiers. Journal of Statistical Software, 46(5), 1-27. DOI https://doi.org/10.18637/jss.v046.i05.
}

\keyword{multi-objective parameter tuning}

\examples{
\donttest{
# optimize the 'cost' and 'kernel' parameters of an SVM according
# to CV error and CV Specificity on the 'iris' data set
# using several predefined values for the cost
r <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier=tunePareto.svm(),
                cost=c(0.001,0.01,0.1,1,10), 
                kernel=c("linear", "polynomial", 
                         "radial", "sigmoid"),
                objectiveFunctions=list(cvError(10, 10), 
                                        cvSpecificity(10, 10, caseClass="setosa")))

# print Pareto-optimal solutions
print(r)

# use a continuous interval for the 'cost' parameter 
# and optimize it using evolutionary algorithms and
# parallel execution with snowfall
library(snowfall)
sfInit(parallel=TRUE, cpus=2, type="SOCK")
r <- tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.svm(), 
                 cost = as.interval(0.001,10), 
                 kernel = c("linear", "polynomial",
                            "radial", "sigmoid"),
                 sampleType="evolution",
                 numCombinations=20,
                 numIterations=20,                      
                 objectiveFunctions = list(cvError(10, 10),
                                           cvSensitivity(10, 10, caseClass="setosa"),
                                           cvSpecificity(10, 10, caseClass="setosa")),
                useSnowfall=TRUE)
sfStop()

# print Pareto-optimal solutions
print(r)

# plot the Pareto fronts
plotDominationGraph(r, legend.x="topright")}
}
