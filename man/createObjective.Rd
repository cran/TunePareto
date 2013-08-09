\name{createObjective}
\alias{createObjective}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a new objective function
}
\description{
Creates a new \code{TuneParetoObjective} object. An objective consists of two parts: The precalculation function, which applies the classifier to the data, and the objective itself, which is calculated from the predicted class labels.
}
\usage{
createObjective(precalculationFunction, 
                precalculationParams = NULL, 
                objectiveFunction, 
                objectiveFunctionParams = NULL,
                direction = c("minimize", "maximize"), 
                name)
}
\arguments{
  \item{precalculationFunction}{
  The name of the precalculation function that applies the classifiers to the data. Two predefined precalculation functions are \code{reclassification} and \code{crossValidation}.
}
  \item{precalculationParams}{
A named list of parameters for the precalculation function. 
}
  \item{objectiveFunction}{
The name of the objective function that calculates the objective from the precalculated class labels.
}

  \item{objectiveFunctionParams}{
  A named list of further parameters for the objective function.
}
 
  \item{direction}{
Specifies whether the objective is minimized or maximized.
}
  \item{name}{
A readable name of the objective.
}
}
\details{
The objective calculation is divided into a precalculation step and the objective calculation itself. The main reason for this is the possibility to aggregate precalculation across objectives. For example, if both the specificity and the sensitivity of a cross-validation (with the same parameters) are required, the cross-validation is run only once to save computational time. Afterwards, the results are passed to both objective functions.

A precalculation function has the following parameters:
\describe{
\item{data}{
The data set to be used for the precalculation. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{labels}{
A vector of class labels for the samples in \code{data}.
}
  \item{classifier}{
A \code{TuneParetoClassifier} wrapper object containing the classifier to tune. A number of state-of-the-art classifiers are included in \pkg{TunePareto} (see \code{\link{predefinedClassifiers}}). Custom classifiers can be employed using \code{\link{tuneParetoClassifier}}.
}

  \item{classifierParams}{
A named list of parameter assignments for the classifier.
}

  \item{predictorParams}{
If the classifier has separate training and prediction functions, a named list of parameter assignments for the predictor.
}

}
Additionally, the function can have further parameters which are supplied in \code{precalculationParams}. To train a classifier and obtain predictions, the precalculation function can call the generic \if{latex}{\cr}\code{\link{trainTuneParetoClassifier}} and \code{\link{predict.TuneParetoModel}} functions. 

The precalculation function usually returns the predicted labels, the true labels and the model, but the only requirement of the return value is that it can be processed by the corresponding objective function. Predefined precalculation functions are \code{\link{reclassification}} and \code{\link{crossValidation}}.

The objective function has a single obligatory parameter named \code{result} which supplies the result of the precalculation. Furthermore, optional parameters can be specified. Their values are taken from \code{objectiveFunctionParams}. The function either returns a single number specifying the objective value, or a list with a \code{score} component containing the objective value and a \code{additionalData} component that contains additional information to be stored in the \code{additionalData} component of the \code{TuneParetoResult} object (see \code{\link{tunePareto}}).

}
\value{
Retuns an object of class \code{TuneParetoObjective} with the following components:
\item{precalculationFunction}{The supplied precalculation function}

\item{precalculationParams}{The additional parameters to be passed to the precalculation function}

\item{objectiveFunction}{The objective function}

\item{minimize}{\code{TRUE} if the objective is minimized, \code{FALSE} if it is maximized.}

\item{name}{The readable name of the objective.}
}

\seealso{
\code{\link{predefinedObjectiveFunctions}}, \code{\link{trainTuneParetoClassifier}}, \code{\link{predict.TuneParetoModel}}
}
\examples{

# create new objective minimizing the number of support vectors
# for a support vector machine

reclassSupportVectors <- function (saveModel = FALSE) 
{
    createObjective(precalculationFunction = reclassification, 
        precalculationParams = NULL, objectiveFunction = 
        function(result, saveModel) 
        {
	        if(result$model$classifier$name != "svm")
		        stop("This objective function can only be applied 
		              to classifiers of type tunePareto.svm()")

      		res <- result$model$model$tot.nSV

		      if (saveModel) 
		      # return a list containing the objective value as well as the model
		      {
		         return(list(additionalData = result$model, fitness = res))
		      }
		      else 
		      # only return the objective value
		        return(res)
        }, 
        objectiveFunctionParams = list(saveModel = saveModel), 
        direction = "minimize", 
        name = "Reclass.SupportVectors")
}

# tune error vs. number of support vectors on the 'iris' data set
r <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier = tunePareto.svm(),
                cost=c(0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50),
                objectiveFunctions=list(reclassError(), reclassSupportVectors()))

print(r)
}
\keyword{objective function
         multi-objective optimization
         MOO}
