\name{tunePareto}
\alias{tunePareto}
\title{
Generic function for multi-objective parameter tuning of classifiers
}
\description{
This generic function tunes parameters of arbitrary classifiers in a multi-objective setting and returns the Pareto-optimal parameter combinations. 
}
\usage{
tunePareto(..., data, labels, 
           classifier, parameterCombinations,
           sampleType = c("full","uniform",
                          "latin","halton",
                          "niederreiter","sobol",
                          "evolution"), 
           numCombinations, 
           mu=10, lambda=20, numIterations=100,
           objectiveFunctions, objectiveBoundaries,
           keepSeed = TRUE, useSnowfall = FALSE, verbose=TRUE)
}

\arguments{
  \item{data}{
The data set to be used for the parameter tuning. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{labels}{
A vector of class labels for the samples in \code{data}.
}
  \item{classifier}{
A \code{TuneParetoClassifier} wrapper object containing the classifier to tune. A number of state-of-the-art classifiers are included in \pkg{TunePareto} \if{latex}{\cr}(see \code{\link{predefinedClassifiers}}). Custom classifiers can be employed using \code{\link{tuneParetoClassifier}}.}

  \item{parameterCombinations}{If not all combinations of parameter ranges for the classifier are meaningful, you can set this parameter instead of specifying parameter values in the \dots argument. It holds an explicit list of possible combinations, where each element of the list is a named sublist with one entry for each parameter.}

  \item{sampleType}{
Determines the way parameter configurations are sampled. 

If \code{type="full"}, all possible combinations are tried. This is only possible if all supplied parameter ranges are discrete or if the combinations are supplied explicitly in \code{parameterCombinations}. 

If \code{type="uniform"}, \code{numCombinations} combinations are drawn uniformly at random.

If \code{type="latin"}, Latin Hypercube sampling is applied. This is particularly encouraged when tuning using continuous parameters.

If \code{type="halton","niederreiter","sobol"}, \code{numCombinations} parameter combinations are drawn on the basis of the corresponding quasi-random sequences (initialized at a random step to ensure that different values are drawn in repeated runs). This is particularly encouraged when tuning using continuous parameters. \code{type="niederreiter"} and \code{type="sobol"} require the \pkg{gsl} package to be installed.

If \code{type="evolution"}, an evolutionary algorithm is applied. In details, this employs \code{mu+lambda} Evolution Strategies with uncorrelated mutations and non-dominated sorting for survivor selection. This is encouraged when the space of possible parameter configurations is very large. For smaller parameter spaces, the above sampling methods may be faster.
}

  \item{numCombinations}{
If this parameter is set, at most \code{numCombinations} randomly chosen parameter configurations are tested. Otherwise, all possible combinations of the supplied parameter ranges are tested.
}


  \item{mu}{
The number of individuals used in the Evolution Strategies if \code{type="evolution"}.
}

  \item{lambda}{
The number of offspring per generation in the Evolution Strategies if \if{latex}{\cr}\code{type="evolution"}.
}

  \item{numIterations}{
The number of iterations/generations the evolutionary algorithm is run if \if{latex}{\cr}\code{type="evolution"}.
}
  
  \item{objectiveFunctions}{
A list of objective functions used to tune the parameters. There are a number of predefined objective functions (see \code{\link{predefinedObjectiveFunctions}}). Custom objective functions can be created using \code{\link{createObjective}}.
}

  \item{objectiveBoundaries}{
If this parameter is set, it specifies boundaries of the objective functions for valid solutions. That is, each element of the supplied vector specifies the upper or lower limit of an objective (depending on whether the objective is maximized or minimized). Parameter combinations that do not meet all these restrictions are not included in the result set, even if they are Pareto-optimal. If only some of the objectives should have bounds, supply \code{NA} for the remaining objectives.
}

  \item{keepSeed}{
If this is true, the random seed is reset to the same value for each of the tested parameter configurations. This is an easy way to guarantee comparability in randomized objective functions. E.g., cross-validation runs of the classifiers will all start with the same seed, which results in the same partitions. 

\bold{Attention: } If you set this parameter to \code{FALSE}, you must ensure that all configuration are treated equally in the objective functions: There may be randomness in processes such as classifier training, but there should be no random difference in the rating itself. In particular, the choice of subsets for subsampling experiments should always be the same for all configurations. For example, you can provide precalculated fold lists to the cross-validation objectives in the \code{foldList} parameter. If parameter configurations are rated under varying conditions, this may yield over-optimistic or over-pessimistic ratings for some configurations due to outliers.
}

  \item{useSnowfall}{
If this parameter is true, the routine loads the \pkg{snowfall} package and processes the parameter configurations in parallel. Please note that the \pkg{snowfall} cluster has to be initialized properly before running the tuning function and stopped after the run.
}

  \item{verbose}{
If this parameter is true, status messages are printed. In particular, the algorithm prints the currently tested combination.
}

  \item{\dots}{
The parameters of the classifier and predictor functions that should be tuned. The names of the parameters must correspond to the parameters specified in \code{classifierParameterNames} and \code{predictorParameterNames} of \code{\link{tuneParetoClassifier}}.
Each supplied argument describes the possible values of a single parameter. These can be specified in two ways: Discrete parameter ranges are specified as lists of possible values. Continous parameter ranges are specified as intervals using \code{\link{as.interval}}. The algorithm then generates combinations of possible parameter values. Alternatively, the combinations can be defined explicitly using the \code{parameterCombinations} parameter.
}

}
\details{
This is a generic function that allows for parameter tuning of a wide variety of classifiers. You can either specify the values or intervals of tuned parameters in the \code{\dots} argument, or supply selected combinations of parameter values using \code{parameterCombinations}. In the first case, combinations of parameter values specified in the \code{\dots} argument are generated. If \code{sampleType="uniform"}, \code{sampleType="latin"}, \code{sampleType="halton"}, \code{sampleType="niederreiter"} or \code{sampleType="sobol"}, a random subset of the possible combinations is drawn. \if{latex}{\cr}If \code{sampleType="evolution"}, random parameter combinations are generated and optimized using Evolution Strategies.

In the latter case, only the parameter combinations specified explicitly in \code{parameterCombinations} are tested. This is useful if certain parameter combinations are invalid. You can create parameter combinations by concatenating results of calls to \code{\link{allCombinations}}. Only \code{sampleType="full"} is allowed in this mode.

For each of the combinations, the specified objective functions are calculated. This usually involves training and testing a classifier. From the resulting objective values, the non-dominated parameter configurations are calculated and returned.

The \code{\dots} argument is the first argument of \code{tunePareto} for technical reasons (to prevent partial matching of the supplied parameters with argument names of \code{tunePareto}. This requires all arguments to be named.
}
\value{
Returns a list of class \code{TuneParetoResult} with the following components:
\item{bestCombinations}{A list of Pareto-optimal parameter configurations. Each element of the list consists of a sub-list with named elements corresponding to the parameter values.}
\item{bestObjectiveValues}{A matrix containing the objective function values of the Pareto-optimal configurations in \code{bestCombinations}. Each row corresponds to a parameter configuration, and each column corresponds to an objective function.}
\item{testedCombinations}{A list of all tested parameter configurations with the same structure as \if{latex}{\cr}\code{bestCombinations}.}
\item{testedObjectiveValues}{A matrix containing the objective function values of all tested configurations with the same structure as \code{bestObjectiveValues}.}
\item{dominationMatrix}{A Boolean matrix specifying which parameter configurations dominate each other. If a configuration \code{i} dominates a configuration \code{j}, the entry in the \code{i}th row and the \code{j}th column is \code{TRUE}.}
\item{minimizeObjectives}{A Boolean vector specifying which of the objectives are minimization objectives. This is derived from the objective functions supplied to \code{tunePareto}.}
\item{additionalData}{A list containing additional data that may have been returned by the objective functions. The list has one element for each tested parameter configuration, each comprising one sub-element for each objective function that returned additional data. The structure of these sub-elements depends on the corresponding objective function. For example, the predefined objective functions (see \code{\link{predefinedObjectiveFunctions}}) save the trained models here if \code{saveModel} is true.}
}
\seealso{
\code{\link{predefinedClassifiers}}, \code{\link{predefinedObjectiveFunctions}}, \code{\link{createObjective}}, \code{\link{allCombinations}}
}
\examples{
\donttest{

# tune 'k' of a k-NN classifier 
# on two classes of the 'iris' data set --
# see ?knn
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.knn(), 
                 k = c(1,3,5,7,9),
                 objectiveFunctions = list(cvError(10, 10),
                                           reclassError())))
                 
# example using predefined parameter configurations,
# as certain combinations of k and l are invalid:
comb <- c(allCombinations(list(k=1,l=0)),
          allCombinations(list(k=3,l=0:2)),
          allCombinations(list(k=5,l=0:4)),
          allCombinations(list(k=7,l=0:6)))

print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.knn(), 
                 parameterCombinations = comb,
                 objectiveFunctions = list(cvError(10, 10),
                                           reclassError())))
                                           

# tune 'cost' and 'kernel' of an SVM on
# the 'iris' data set using Latin Hypercube sampling --
# see ?svm and ?predict.svm
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.svm(), 
                 cost = as.interval(0.001,10), 
                 kernel = c("linear", "polynomial",
                          "radial", "sigmoid"),
                 sampleType="latin",
                 numCombinations=20,                          
                 objectiveFunctions = list(cvError(10, 10),
                                           cvSensitivity(10, 10, caseClass="setosa"))))

# tune the same parameters using Evolution Strategies
\donttest{print(tunePareto(data = iris[, -ncol(iris)], 
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
                                           cvSpecificity(10, 10, caseClass="setosa"))))}
}
}
\keyword{parameter tuning multi-objective optimization MOO parallel multi-core}
