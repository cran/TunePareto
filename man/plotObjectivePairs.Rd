\name{plotObjectivePairs}
\alias{plotObjectivePairs}
\title{
Plot a matrix of Pareto front panels
}
\description{
Plots a matrix of Pareto front panels for each pair of objectives. The plot for \code{n} objectives consists of \code{n x n} panels, where the panel in row \code{i} and column \code{j} depicts the Pareto fronts of the \code{i}-th and the \code{j}-th objective. Each of the panels is drawn in the same way as \code{\link{plotParetoFronts2D}}. 
}
\usage{
plotObjectivePairs(tuneParetoResult, 
                   drawLabels = TRUE, 
                   drawBoundaries = TRUE, 
                   labelPos = 4, 
                   ...)
}
\arguments{
  \item{tuneParetoResult}{
An object of class \code{TuneParetoResult} as returned by \code{\link{tunePareto}}.
}
  \item{drawLabels}{
If set to true, the descriptions of the configurations are printed next to the points in the plot.
}

  \item{drawBoundaries}{
If set to true, the upper or lower objective limits supplied in the \code{objectiveBoundaries} parameter of \code{\link{tunePareto}} are drawn as horizontal and vertical lines.
}

  \item{labelPos}{
The position of the configuration labels in the plot. Values of 1, 2, 3 and 4 denote positions below, to the left of, above and to the right of the points on the Pareto fronts.
}  

  \item{...}{
Further graphical parameters to be passed to the \code{\link{plot}} function.
}
}
\value{
This function does not have a return value.
}

\seealso{
\code{\link{tunePareto}}, \code{\link{plotParetoFronts2D}}, \code{\link{plotDominationGraph}}
}
\examples{
# optimize the 'cost' parameter of an SVM according
# to CV error, CV error variance, and CV Specificity 
# on two classes of the 'iris' data set
r <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier = tunePareto.svm(),
                cost=c(0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50), 
                objectiveFunctions=list(cvError(10, 10),
                                        cvErrorVariance(10, 10),
                                        cvSpecificity(10, 10, caseClass="virginica")))

# plot the matrix of Pareto fronts
plotObjectivePairs(r)         
}

\keyword{multi-objective optimization
         MOO
         Pareto front}
