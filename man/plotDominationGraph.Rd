\name{plotDominationGraph}
\alias{plotDominationGraph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the Pareto fronts of parameter configuration scores
}
\description{
Draws the Pareto fronts and domination relations of tested parameter configurations in a graph. Here, the leftmost column of nodes represents the non-dominated configurations (i.e. the first Pareto front). The second column contains the second Pareto front, i.e. the configurations that are only dominated by the first Pareto front, and so on. An edge between two configurations indicate that the first configuration dominates the second.
}
\usage{
plotDominationGraph(tuneParetoResult, 
                    transitiveReduction = TRUE, 
                    drawDominatedObjectives = TRUE, 
                    drawLabels = TRUE, 
                    drawLegend = TRUE,
                    legend.x = "topleft", 
                    ...)
}

\arguments{
  \item{tuneParetoResult}{
An object of class \code{TuneParetoResult} as returned by \code{\link{tunePareto}}.
}
  \item{transitiveReduction}{
If this is true, transitive edges in the graph are removed to enhance readability. That is, if configuration \code{c1} dominates configuration \code{c2} and \code{c2} dominates \code{c3}, no edge from \code{c1} to \code{c3} is drawn.
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

  \item{legend.x}{
The position of the legend. For details, refer to the \code{x} parameter of \code{\link{legend}}.
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
# call tunePareto.knn with different k and l on the 'iris' data set
x <- tunePareto.knn(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 k = c(5,7,9),
                 l = c(1,2,3),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa")))

# plot the graph                 
\donttest{plotDominationGraph(x)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Pareto front
         domination graph
         objective function}