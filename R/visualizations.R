#################################################################
# Visualization routines
#################################################################

# Retrieves a list of Pareto fronts from 
# a domination matrix <mat>.
#
# Returns a list of vectors.
calculateParetoFronts <- function(mat)
{
  paretoFronts <-list()
  indexList <- 1:nrow(mat)
  
  # calculate the Pareto fronts
  while (length(indexList) != 0)
  {
    nonDominated <- which(!apply(mat,2,any))
    
    paretoFronts[[length(paretoFronts) + 1]] <- indexList[nonDominated]
   
    # remove the Pareto front from the list of points to consider 
    indexList <- indexList[-nonDominated]
    mat <- mat[-nonDominated, -nonDominated, drop=FALSE]
  }
  return(paretoFronts)
}


# Plots a graph in which each column of nodes represents one Pareto front
# and edges represent a domination relation.
# <tuneParetoResult> is a TuneParetoResult object to be plotted.
# If <transitiveReduction> is true, transitive edges are removed from the graph.
# If <drawDominatedObjectives> is true, color bars that indicate the objectives in
# which a configuration is the best in its Pareto front are drawn.
# <drawLabels> specifies whether the configuration descriptions should be drawn.
# <drawLegend> specifies whether a legend for the color bars should be plotted.
# <legend.x> specifies the position of this legend.
#
# Invisibly returns the igraph object.
plotDominationGraph <- function(tuneParetoResult, transitiveReduction=TRUE, 
                                drawDominatedObjectives=TRUE,
                                drawLabels=TRUE,
                                drawLegend=TRUE,
                                legend.x="topleft", ...)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")
  
  require(igraph)
  
  colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
      "purple","orange","seagreen","tomato","darkgray","chocolate",
      "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
      "burlywood2")
  
  mat <- tuneParetoResult$dominationMatrix
  edges <- tuneParetoResult$dominationMatrix
  
  if (transitiveReduction)
  # remove transitive edges
  {
    transitiveGraph <- edges
    # first, build the complete transitive graph
    for (i in 1:nrow(edges))
    {
      for (j in 1:nrow(edges))
      {
        for (k in 1:nrow(edges))
        {
          if (edges[i,j] && edges[j,k])
            transitiveGraph[i,k] <- TRUE
        }
      }
    }

    # now check if there are transitive dependencies
    # and remove the corresponding edges
    for (i in 1:nrow(transitiveGraph))
    {
      for (j in 1:nrow(transitiveGraph))
      {
        for (k in 1:nrow(transitiveGraph))
        {
          if (transitiveGraph[i,j] && transitiveGraph[j,k])
            edges[i,k] <- FALSE
        }
      }
    }

  }
  
  paretoFronts <- calculateParetoFronts(mat)
  
  # create igraph object
  g <- graph.adjacency(t(edges), mode="directed")
  
  dim_y <- 10
  
  positions <- matrix(ncol=2, nrow=nrow(edges))
  
  # calculate graph layout
  for (i in 1:length(paretoFronts))
  {

    indices <- paretoFronts[[i]]
    
    if (i %% 2 == 0)
    {
      distance <- dim_y/length(indices)
      offset <- 0
    }
    else
    {
      distance <- (dim_y - 1)/length(indices)
      offset <- 0.5
    }  
    positions[indices, ] <- t(sapply(1:length(indices),function(j)
                              {
                                c(i + 1, offset + (j-1)*distance+distance/2)
                              }))
  }
  
  args <- list(...)
  
  # check for certain graphical parameters in ... 
  # that have different default values in this plot
  if (is.null(args$vertex.size))
    args$vertex.size <- 2
    
  if (is.null(args$vertex.color))
    args$vertex.color <- "grey"
    
  if (is.null(args$vertex.label.cex))
    args$vertex.label.cex <- 0.7
    
  if (is.null(args$vertex.label.dist))
    args$vertex.label.dist <- 0.2
    
  if (is.null(args$edge.arrow.size))
    args$edge.arrow.size <- 0.3
    
  if (is.null(args$edge.curved))
    args$edge.curved <- !transitiveReduction
    
  if (drawLabels)
    args$vertex.label <- rownames(edges)
  else
    args$vertex.label <- NA
  
  plot(g, vertex.label=args$vertex.label, vertex.color=args$vertex.color,
      vertex.label.cex = args$vertex.label.cex, vertex.label.dist=args$vertex.label.dist, 
      vertex.size=args$vertex.size, edge.arrow.size=args$edge.arrow.size, 
      layout=positions[1:nrow(edges),,drop=FALSE], ...)
  
  # normalize layout
  positions <- layout.norm(positions,-1, 1, -1, 1)

  cols <- c()  
  if (drawDominatedObjectives)
  {
    # calculate the dominated objectives for the color indicators
    dominatedObjectives <- lapply(paretoFronts,function(front)
    {
     r <- lapply(1:ncol(tuneParetoResult$testedObjectiveValues),function(j)
            {
              vec <- tuneParetoResult$testedObjectiveValues[front,j]
              if (tuneParetoResult$minimizeObjectives[j])
                front[which(vec == min(vec))]
              else
                front[which(vec == max(vec))]
            })
      matrix(sapply(front,function(x)
            {
              sapply(r,function(y)
              {
                x %in% y
              })                
            }), ncol=length(front))
    })
        
    # calculate the positions and colors of the color indicators      
    for (i in 1:length(paretoFronts))
    {
      for (c in 1:ncol(dominatedObjectives[[i]]))
      {
        sumDom <- 1          
        for (r in 1:nrow(dominatedObjectives[[i]]))
        {       
          if (dominatedObjectives[[i]][r,c])
          {
            positions <- rbind(positions,
                               c(positions[paretoFronts[[i]][c],1] - sumDom*0.02 - 0.02, 
                                 positions[paretoFronts[[i]][c],2]))
            cols <- c(cols, colorSet[r])
            sumDom <- sumDom + 1
          }         
        }
      }
    }
    
    # draw the indicators     
    pts <- positions[-(1:nrow(edges)),,drop=FALSE]
    points(pts[,1], pts[,2], pch=15, col=cols,
           cex=0.8)
           
    if (drawLegend)
      # draw the legend       
      legend(x=legend.x, pch=15, ncol=1,
             col=colorSet[1:ncol(tuneParetoResult$testedObjectiveValues)],
             legend = paste("Dominates front in",colnames(tuneParetoResult$testedObjectiveValues)),
             cex=0.75, box.lty=0)
  } 
     
  # return the graph
  return(invisible(g))
}

# Internal function called from plotParetoFronts2D and plotObjectivePairs
# to plot the Pareto fronts of two objectives.
# <objectiveVals> is a matrix of objective values.
# <minimizeObjectives> is a Boolean vector specifying which objectives are minimized.
# If <drawLabels> is TRUE, the parameter configurations are printed.
# If <plotNew> is true, a new plot is started, otherwise lines are added to existing plots.
# <xlim>, <ylim>, <xlab> and <ylab> correspond to the parameters in the generic plot routine.
# <labelPos> specifies the position of the labels with respect to the points
internal.plotParetoFronts2D <- function(objectiveVals, minimizeObjectives, boundaries,
                                        drawLabels=TRUE, drawBoundaries=TRUE, plotNew=TRUE, 
                                        xlim, ylim, xlab, ylab, labelPos=4, ...)
{
  colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
      "purple","orange","seagreen","tomato","darkgray","chocolate",
      "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
      "burlywood2")

  # calculate domination matrix of the 2 chosen objectives
  domination <- calculateDominationMatrix(objectiveVals, minimizeObjectives)
  
  # calculate the Pareto fronts
  paretoFronts <- calculateParetoFronts(domination)
  
  args <- list(...)
  
  if (missing(xlim))
    xlim <- c(min(objectiveVals[,1]), max(objectiveVals[,1]))
  
  if (missing(ylim))
    ylim <- c(min(objectiveVals[,2]), max(objectiveVals[,2]))
    
  if (missing(xlab))
    xlab <- colnames(objectiveVals)[1]
  
  if (missing(ylab))
    ylab <- colnames(objectiveVals)[2]  
    
  if (plotNew)
  {
        plot(1, type="n", 
           xlim=xlim,
           ylim=ylim,
           xlab=xlab,
           ylab=ylab, ...)
  }
  
  if (drawBoundaries && !is.null(boundaries))
  {
    if (!is.na(boundaries[1]) && !is.null(boundaries[1]))
      abline(v=boundaries[1], col="darkgrey", lty=2)
    if (!is.na(boundaries[2]) && !is.null(boundaries[2]))
      abline(h=boundaries[2], col="darkgrey", lty=2)
  }
  
  for (i in 1:length(paretoFronts))
  {
    # extract points and order them by their x value
    pts <- objectiveVals[paretoFronts[[i]],,drop=FALSE]
    pts <- pts[order(pts[,1]),,drop=FALSE]
    
    # add new Pareto front
    lines(pts[,1], pts[,2], col=colorSet[i %% length(colorSet)], pch=8, type="o")
 
    if (drawLabels)
      # add the configuration descriptions
      text(pts[,1], pts[,2], rownames(pts), cex=0.5, pos=labelPos)
  }
 
}

# Plots a 2D plot of two objectives in an optimization.
# <tuneParetoResult> is an object of class TuneParetoResult.
# <objectives> is a vector of indices or names of the objectives to plot.
# If <drawLabels> is true, the descriptions of the configurations are written next to the plot
# <labelPos> specifies the position of the labels with respect to the points
plotParetoFronts2D <- function(tuneParetoResult, objectives, drawLabels=TRUE, drawBoundaries=TRUE, labelPos=4,...)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
  {
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")
  }
  
  if (missing(objectives))
  {
    if (ncol(tuneParetoResult$testedObjectiveValues) == 2)
      objectives = c(1,2)
    else
      stop("Please supply exactly 2 objectives!")
  }
  
  
  if (length(objectives) != 2)
  {
    stop("Please supply exactly 2 objectives!")
  }
  
  colorSet <- c("blue","green","red","darkgoldenrod","gold","brown","cyan",
      "purple","orange","seagreen","tomato","darkgray","chocolate",
      "maroon","darkgreen","gray12","blue4","cadetblue","darkgoldenrod4",
      "burlywood2")
  
  objectiveVals <- tuneParetoResult$testedObjectiveValues[,objectives,drop=FALSE]
  minimizeObjectives <- tuneParetoResult$minimizeObjectives[objectives]
  boundaries <- tuneParetoResult$objectiveBoundaries[objectives]
  
  internal.plotParetoFronts2D(objectiveVals, minimizeObjectives, boundaries, 
                              drawLabels=drawLabels, drawBoundaries=drawBoundaries,
                              labelPos=labelPos,...)
}

# Plots the Pareto fronts of pairs of objectives in <tuneParetoResult>
# in a matrix-like plot. If <drawLabels> is true, the parameter configurations
# are printed in the plot.
plotObjectivePairs <- function(tuneParetoResult, drawLabels=TRUE, drawBoundaries=TRUE, labelPos=4, ...)
{
  xpos <- 0
  ypos <- 1
  
  pairs(data.frame(tuneParetoResult$testedObjectiveValues), 
        panel=function(x,y)
              {
                # calculate the current x and y position in the plot
                xpos <<- xpos + 1
                if (xpos == ypos)
                  xpos <<- xpos + 1
                  
                if (xpos > length(tuneParetoResult$minimizeObjectives))
                {
                  ypos <<- ypos + 1
                  xpos <<- 1
                }
                
                # build output matrix
                data <- data.frame(x,y)
                colnames(data) <- names(tuneParetoResult$minimizeObjectives)[c(xpos,ypos)]
                rownames(data) <- rownames(tuneParetoResult$testedObjectiveValues)
                
                # plot Pareto fronts
                internal.plotParetoFronts2D(data, 
                                            tuneParetoResult$minimizeObjectives[c(xpos,ypos)],
                                            tuneParetoResult$objectiveBoundaries[c(xpos,ypos)], 
                                            drawLabels=drawLabels, drawBoundaries=drawBoundaries,
                                            labelPos=labelPos,
                                            plotNew=FALSE, ...)
              })
}

