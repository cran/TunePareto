\name{plotDominationGraph}
\alias{plotDominationGraph}
\title{
Visualize the Pareto fronts of parameter configuration scores
}
\description{
Draws the Pareto fronts and domination relations of tested parameter configurations in a graph. Here, the leftmost column of nodes represents the non-dominated configurations (i.e. the first Pareto front). The second column contains the second Pareto front, i.e. the configurations that are only dominated by the first Pareto front, and so on. An edge between two configurations indicate that the first configuration is dominated by the second.
}
\usage{
plotDominationGraph(tuneParetoResult, 
                    transitiveReduction = TRUE, 
                    drawDominatedObjectives = TRUE, 
                    drawLabels = TRUE, 
                    drawLegend = TRUE,
                    x.legend = "topleft",
                    cex.legend = 0.7,
                    col.indicator, 
                    pch.indicator = 15, 
                    cex.indicator = 0.8,
                    ...)
}

\arguments{
  \item{tuneParetoResult}{
An object of class \code{TuneParetoResult} as returned by \code{\link{tunePareto}}.
}
  \item{transitiveReduction}{
If this is true, transitive edges in the graph are removed to enhance readability. That is, if configuration \code{c1} dominates configuration \code{c2} and \code{c2} dominates \code{c3}, no edge from \code{c3} to \code{c1} is drawn.
}
  \item{drawDominatedObjectives}{
If set to true, color indicators are drawn next to the nodes. Here, each color corresponds to one objective. The color is drawn next to a node if this node has the best score in this objectives among all solutions of the same Pareto front (i.e., column of the graph).
}
  \item{drawLabels}{
Specifies whether the parameter configurations should be printed next to the corresponding edges.
}
  \item{drawLegend}{
If \code{drawDominatedObjectives=TRUE}, this specifies whether a legend with the objective colors should be drawn.
}

  \item{x.legend}{
The position of the legend. For details, refer to the \code{x} parameter of \code{\link{legend}}.
}

  \item{cex.legend}{
Specifies the size of the text in the legend if \code{drawLegend} is true.  
}

  \item{col.indicator}{
Specifies an optional list of colors, one for each objective function. These colors will be used for the indicators if \code{drawDominatedObjectives} is true. By default, a predefined set of colors is used. 
}

  \item{pch.indicator}{
Specifies a single plotting character or a list of plotting characters for the objective functions in the indicators which is used for the indicators if \code{drawDominatedObjectives} is true.  
}

  \item{cex.indicator}{
Specifies the size of the symbols in the indicators which is be used for the indicators if \code{drawDominatedObjectives} is true. This can also be a vector of sizes for the symbols of the objectives. 
}
  
  \item{\dots}{
Further graphical parameters for \code{\link[igraph:plot.graph]{plot.igraph}}.
}
}

\value{
Invisibly returns the \code{igraph} object representing the graph.
}

\seealso{
\code{\link{tunePareto}}
}
\examples{
# call tunePareto using a k-NN classifier 
# with different 'k' and 'l' on the 'iris' data set
x <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier = tunePareto.knn(),
                k = c(5,7,9),
                l = c(1,2,3),
                objectiveFunctions=list(cvError(10, 10),
                                        cvSpecificity(10, 10, caseClass="setosa")))

# plot the graph                 
\donttest{plotDominationGraph(x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pareto front domination graph objective function}
