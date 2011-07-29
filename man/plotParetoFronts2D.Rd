\name{plotParetoFronts2D}
\alias{plotParetoFronts2D}
\title{
A classical 2-dimensional plot of Pareto fronts
}
\description{
Draws a classical Pareto front plot of 2 objectives of a parameter tuning. The first objective is on the x axis of the plot, and the second objective is on the y axis. Points on a Pareto front are connected by lines. Each Pareto front is drawn in a different color.
}
\usage{
plotParetoFronts2D(tuneParetoResult, 
                   objectives, 
                   drawLabels = TRUE, 
                   drawBoundaries = TRUE, 
                   labelPos = 4,
                   fitLabels=TRUE,
                   cex.conf=0.5, 
                   lty.fronts=1, 
                   pch.fronts=8, 
                   col.fronts,
                   ...)
}
\arguments{
  \item{tuneParetoResult}{
An object of class \code{TuneParetoResult} as returned by \code{\link{tunePareto}}.
}
  \item{objectives}{
The names or indices of the two objectives to plot. Pareto-optimality is determined only on the basis of these two objectives. Optional if the parameter tuning has exactly two objectives.
}
  \item{drawLabels}{
If set to true, the descriptions of the configurations are printed next to the points in the plot.
}

  \item{drawBoundaries}{
If set to true, the upper or lower objective limits supplied in the \code{objectiveBoundaries} parameter of \code{\link{tunePareto}} are drawn as horizontal and vertical lines.
}

  \item{labelPos}{
The position of the configuration labels in the plot (if \code{drawLabels} is true). Values of 1, 2, 3 and 4 denote positions below, to the left of, above and to the right of the points on the Pareto fronts.
}

  \item{fitLabels}{
If this parameter is true (and \code{drawLabels} is true), overlapping or partially hidden configuration labels are removed from the plot to improve the readability of the remaining labels.
}

  \item{cex.conf}{
The size of the configuration labels in the plots (if \code{drawLabels} is true).  
}

  \item{lty.fronts}{
A vector of line types to use for the Pareto fronts. By default, straight lines are drawn for all fronts.  
}

  \item{pch.fronts}{
A vector of symbols to use for points on the Pareto fronts. All points on the same front will have the same symbol. By default, an asterisk is used.
}

  \item{col.fronts}{
A vector of colors to use for the Pareto fronts. By default, a predefined set of colors is used.  
}

  \item{...}{
Further graphical parameters to be passed to the \code{\link{plot}} function.
}
}

\value{
This function does not have a return value.
}

\seealso{
\code{\link{tunePareto}}, \code{\link{plotDominationGraph}}
}
\examples{
# optimize the 'cost' parameter according
# to CV error and CV Specificity on the 'iris' data set
r <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier = tunePareto.svm(),
                cost=c(0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50), 
                objectiveFunctions=list(cvError(10, 10),
                                        cvSpecificity(10, 10, caseClass="setosa")))
                   
# plot the Pareto graph
plotParetoFronts2D(r)
}
\keyword{multi-objective optimization
         MOO
         Pareto front}
