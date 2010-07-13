\name{tunePareto.knn}
\alias{tunePareto.knn}
\alias{tunePareto.svm}
\alias{tunePareto.tree}
\alias{tunePareto.randomForest}
\title{
Wrappers of tunePareto for certain classifiers
}
\description{
Simplified wrappers of \code{\link{tunePareto}} for the k-Nearest Neighbour classifier, support vector machines, and trees.
}
\usage{
tunePareto.knn(data, labels, 
               k, l, use.all, 
               numCombinations, 
               objectiveFunctions, 
               objectiveBoundaries,
               keepSeed = TRUE, 
               useSnowfall = FALSE,
               verbose=TRUE)
               
tunePareto.svm(data, labels, 
               kernel, degree, gamma, 
               coef0, cost, nu, 
               class.weights, cachesize, 
               tolerance, epsilon, subset, na.action, 
               numCombinations, 
               objectiveFunctions, 
               objectiveBoundaries,               
               keepSeed = TRUE, 
               useSnowfall = FALSE,
               verbose=TRUE)
               
tunePareto.tree(data, labels, 
                weights, subset, 
                na.action, method,
                split, mincut, minsize, mindev, 
                numCombinations, 
                objectiveFunctions,
                objectiveBoundaries,                
                keepSeed = TRUE, 
                useSnowfall = FALSE,
                verbose=TRUE)
                
tunePareto.randomForest(data, labels,
                        subset,
                        na.action,
                        ntree,
                        mtry,
                        replace, 
                        classwt, 
                        cutoff, 
                        strata,
                        sampsize,
                        nodesize,
                        maxnodes,
                        numCombinations,                 
                        objectiveFunctions,
                        objectiveBoundaries,
                        keepSeed = TRUE,
                        useSnowfall = FALSE,
                        verbose=TRUE)                                           
}
\arguments{
  \item{data}{
The data set to be used for the parameter tuning. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{labels}{
A vector of class labels for the samples in \code{data}.
}
  \item{k, l, use.all}{
The corresponding parameters of the \code{\link[class]{knn}} function. For all these parameters, lists of values can be supplied. The corresponding sets of parameter values are assessed in the parameter tuning.
}

  \item{kernel, degree, gamma, coef0, cost, nu}{
The corresponding parameters of the \code{\link[e1071]{svm}} function. For all these parameters, lists of values can be supplied. The corresponding sets of parameter values are assessed in the parameter tuning.
}

  \item{class.weights, cachesize, tolerance, epsilon}{
Further parameters of \code{\link[e1071]{svm}} as above.
}

  \item{weights, method, split, mincut, minsize, mindev}{The corresponding parameters of the \code{\link[tree]{tree}} function. For all these parameters, lists of values can be supplied. The corresponding sets of parameter values are assessed in the parameter tuning.}
  
\item{ntree, mtry, replace, classwt, cutoff, strata, sampsize, nodesize, maxnodes}{The corresponding parameters of the \code{\link[randomForest]{randomForest}} function. For all these parameters, lists of values can be supplied. The corresponding sets of parameter values are assessed in the parameter tuning.}  
 
  \item{na.action, subset}{Common parameters of \code{\link[e1071]{svm}}, \code{\link[tree]{tree}}, and \code{\link[randomForest]{randomForest}}. These parameters take different values in the functions -- please refer to the corresponding help files. Again, lists of parameter values can be passed to these arguments.}

  \item{numCombinations}{
If this parameter is set, at most \code{numCombinations} randomly chosen parameter configurations are tested. Otherwise, all possible combinations of the supplied parameter ranges are tested.
}
  \item{objectiveFunctions}{
A list of objective functions used to tune the parameters. There are a number of predefined objective functions (see \code{\link{predefinedObjectiveFunctions}}. Custom objective functions can be created using \code{\link{createObjective}}.
}

  \item{objectiveBoundaries}{
If this parameter is set, it specifies boundaries of the objective functions for valid solutions. That is, each element of the supplied vector specifies the upper or lower limit of an objective (depending on whether the objective is maximized or minimized). Parameter combinations that do not meet all these restrictions are not included in the result set, even if they are Pareto-optimal.
}

  \item{keepSeed}{
If this is true, the random seed is reset to the same value for each of the tested parameter configurations. This is an easy way to guarantee comparability in randomized objective functions, e.g. in cross-validation runs of the classifiers.
}
  \item{useSnowfall}{
If this parameter is true, the routine loads the \pkg{snowfall} package and processes the parameter configurations in parallel. Please note that the \pkg{snowfall} cluster has to be initialized properly before running the tuning function and stopped after the run.
}

  \item{verbose}{
If this parameter is true, status messages are printed. In particular, the algorithm prints the currently tested combination.
}

}
\value{
Returns a list of class \code{TuneParetoResult} as described in \code{\link{tunePareto}}
}


\seealso{
\code{\link{tunePareto}}, \code{\link{predefinedObjectiveFunctions}}, \code{\link{createObjective}}
}
\examples{
# call tunePareto.knn with different k and l on 
# the 'iris' data set
print(tunePareto.knn(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 k = c(5,7,9),
                 l = c(1,2,3),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa"))))
                 
# call tunePareto.svm with different costs on 
# the 'iris' data set
print(tunePareto.svm(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 cost=c(0.001,0.01,0.1,1,10),
                 objectiveFunctions=list(cvWeightedError(10, 10),
                                         cvSensitivity(10, 10, caseClass="setosa"))))

# call tunePareto.tree with different splitting criteria on 
# the 'iris' data set
print(tunePareto.tree(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 split=c("deviance","gini"),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvErrorVariance(10, 10))))

# call tunePareto.randomForest with different numbers of trees on
# the 'iris' data set
print(tunePareto.randomForest(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 ntree=seq(50,300,50),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa"))))
                             
}

\keyword{parameter tuning
         multi-objective optimization
         MOO
         knn
         svm
         tree
         randomForest}
