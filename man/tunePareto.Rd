\name{tunePareto}
\alias{tunePareto}
\title{
Generic function for multi-objective parameter tuning of classifiers
}
\description{
This generic function can be configured to tune parameters of arbitrary classifiers in a multi-objective setting. Predefined versions for knn, tree, and support vector machines exist (see \code{\link{tunePareto.knn}}, \code{\link{tunePareto.tree}}, \code{\link{tunePareto.svm}}).
}
\usage{
tunePareto(data, labels, 
           classifier, classifierParameterRanges, classifierParameterCombinations,
           predictor, predictorParameterRanges, predictorParameterCombinations,
           numCombinations, objectiveFunctions, objectiveBoundaries,
           keepSeed = TRUE, useSnowfall = FALSE, verbose=TRUE,
           useFormula = FALSE, formulaName = "formula", 
           trainDataName = "x", trainLabelName = "y", 
           testDataName = "newdata", modelName = "object")
}

\arguments{
  \item{data}{
The data set to be used for the parameter tuning. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{labels}{
A vector of class labels for the samples in \code{data}.
}
  \item{classifier}{
The classification function to use. If \code{predictor} is \code{NULL}, this function is an all-in-one classification method that receives both training data and test data and returns the predicted labels for the test data. If \code{predictor} is not \code{NULL}, this is the training function of the classifier that builds a model from the training data. This model is then passed to \code{predictor} along with the test data to obtain the predicted labels for the test data.
}

  \item{classifierParameterRanges}{
A list of lists of parameter ranges to try for \code{classifier}. That is, each element of the list specifies the values of a single parameter to be tested and is named according to this parameter. It is also possible to set parameters to fixed values by specifying only one value. 
}

  \item{classifierParameterCombinations}{If not all combinations of parameter ranges are meaningful, you can set this parameter instead of \code{classifierParameterRanges}. It holds an explicit list of possible combinations, where each element of the list is a named sublist with one entry for each parameter.}

  \item{predictor}{
If the classification method consists of separate training and prediction functions, this points to the prediction function that receives a model and the test data as inputs and returns the predicted class labels.
}

  \item{predictorParameterRanges}{
If \code{predictor} is not \code{NULL}, this argument can be used to supply parameters to the predictor function in the same way as \code{classifierParameterRanges}.
}

\item{predictorParameterCombinations}{
If not all combinations of parameter ranges are meaningful, you can set this parameter instead of \code{predictorParameterRanges}. It holds an explicit list of possible combinations, where each element of the list is a named sublist with one entry for each parameter.}

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

  \item{useFormula}{
Set this to true if the classifier expects a formula to describe the relation between features and class labels. The formula itself is built automatically.
}
  \item{formulaName}{
If \code{useFormula} is true, this is the name of the parameter of the classifier's training function that holds the formula.
}
  \item{trainDataName}{
The name of the paramater of the classifier's training function that holds the training data. 
}
  \item{trainLabelName}{
If \code{useFormula=FALSE}, the name of the paramater of the classifier's training function that holds the training labels. Otherwise, the training labels are added to the training data and supplied in parameter \code{trainDataName}.
}
  \item{testDataName}{
If \code{predictor=NULL}, this is the name of the parameter of \code{classifier} that receives the test data. Otherwise, it is the parameter of \code{predictor} that holds the test data.
}
  \item{modelName}{
If \code{predictor} is not NULL, this is the name of the parameter of \code{predictor} that receives the training model (i.e., the return value of \code{classifier}).
}
}
\details{
This is a generic function that allows for parameter tuning of a wide variety of classifiers. You can either specify the parameters \code{classifierParameterRanges} and \code{predictorParameterRanges}, or the parameters \code{classifierParameterCombinations} and \code{predictorParameterCombinations}. In the first case, all combinations of parameter values specified in \code{classifierParameterRanges} and \code{predictorParameterRanges} are generated. If \code{numCombinations} is specified, a random subset of these combinations is drawn.

In the latter case, only the parameter combinations specified explicitly in \code{classifierParameterCombinations} and \code{predictorParameterCombinations} are tested. This is useful if certain parameter combinations are invalid. You can create parameter combinations by concatenating results of calls to \code{\link{allCombinations}}. Again, if \code{numCombinations} is specified, a random subset of these combinations is drawn.

For each of the combinations, the specified objective functions are calculated. This usually involves training and testing a classifier. From the resulting objective values, the non-dominated parameter configurations are calculated and returned.
}
\value{
Returns a list of class \code{TuneParetoResult} with the following components:
\item{bestCombinations}{A list of Pareto-optimal parameter configurations. Each element of the list consists of a sub-list with named elements corresponding to the parameter values.}
\item{bestObjectiveValues}{A matrix containing the objective function values of the Pareto-optimal configurations in \code{bestCombinations}. Each row corresponds to a parameter configuration, and each column corresponds to an objective function.}
\item{testedCombinations}{A list of all tested parameter configurations with the same structure as \code{bestCombinations}.}
\item{testedObjectiveValues}{A matrix containing the objective function values of all tested configurations with the same structure as \code{bestObjectiveValues}.}
\item{dominationMatrix}{A Boolean matrix specifying which parameter configurations dominate each other. If a configuration \code{i} dominates a configuration \code{j}, the entry in the \code{i}th row and the \code{j}th column is \code{TRUE}.}
\item{minimizeObjectives}{A Boolean vector specifying which of the objectives are minimization objectives. This is derived from the objective functions supplied to \code{tunePareto}.}
}
\seealso{
\code{\link{tunePareto.knn}}, \code{\link{tunePareto.tree}}, \code{\link{tunePareto.svm}}, \code{\link{predefinedObjectiveFunctions}}, \code{\link{createObjective}}, \code{\link{allCombinations}}
}
\examples{


# tune k of knn on two classes of
# the 'iris' data set --
# equivalent to tunePareto.knn
# - see ?knn
library(class)
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = knn, 
                 predictor = NULL, 
                 classifierParameterRanges = list(k=c(1,3,5,7,9)),
                 objectiveFunctions = list(cvError(10, 10),
                                           reclassError()),
                 trainDataName = "train",
                 testDataName = "test",
                 trainLabelName = "cl"))
                 
# example using predefined parameter configurations,
# as certain combinations of k and l are invalid:
comb <- c(allCombinations(list(k=1,l=0)),
          allCombinations(list(k=3,l=0:2)),
          allCombinations(list(k=5,l=0:4)),
          allCombinations(list(k=7,l=0:6)))

print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = knn, 
                 predictor = NULL, 
                 classifierParameterCombinations = comb,
                 objectiveFunctions = list(cvError(10, 10),
                                           reclassError()),
                 trainDataName = "train",
                 testDataName = "test",
                 trainLabelName = "cl"))

# tune cost and kernel of svm on
# the 'iris' data set --
# equivalent to tunePareto.svm
# see ?svm and ?predict.svm
library(e1071)
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = svm, 
                 predictor = predict, 
                 classifierParameterRanges = list(cost=c(0.001,0.01,0.1,1,10), 
                                                  kernel=c("linear", "polynomial",
                                                 "radial", "sigmoid")),
                 objectiveFunctions = list(cvError(10, 10),
                                           cvSensitivity(10, 10, caseClass="setosa")),
                 trainDataName = "x",
                 trainLabelName = "y",
                 testDataName = "newdata",
                 modelName = "object"))

}
\keyword{parameter tuning
         multi-objective optimization
         MOO
         parallel
         multi-core}
