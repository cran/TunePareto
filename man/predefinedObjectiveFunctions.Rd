\name{predefinedObjectiveFunctions}
\alias{predefinedObjectiveFunctions}
\alias{reclassError}
\alias{reclassWeightedError}
\alias{reclassSensitivity}
\alias{reclassSpecificity}
\alias{reclassConfusion}
\alias{cvError}
\alias{cvErrorVariance}
\alias{cvWeightedError}
\alias{cvSensitivity}
\alias{cvSpecificity}
\alias{cvConfusion}
\title{
Predefined objective functions for parameter tuning
}
\description{
Predefined objective functions that calculate the classification error, the sensitivity or the specificity of reclassification or cross-validation experiments.
}
\usage{
reclassError()
reclassWeightedError()
reclassSensitivity(caseClass)
reclassSpecificity(caseClass)
reclassConfusion(trueClass, predictedClass)
cvError(nfold = 10, ntimes = 10, 
        leaveOneOut = FALSE, stratified=FALSE, 
        foldList=NULL)
cvErrorVariance(nfold = 10, ntimes = 10, 
                leaveOneOut = FALSE, stratified=FALSE, 
                foldList=NULL)
cvWeightedError(nfold = 10, ntimes = 10, 
                leaveOneOut = FALSE, stratified=FALSE, 
                foldList=NULL)
cvSensitivity(nfold = 10, ntimes = 10, 
              leaveOneOut = FALSE, stratified=FALSE, 
              foldList=NULL, caseClass)
cvSpecificity(nfold = 10, ntimes = 10, 
              leaveOneOut = FALSE, stratified=FALSE, 
              foldList=NULL, caseClass)
cvConfusion(nfold = 10, ntimes = 10, 
            leaveOneOut = FALSE, stratified=FALSE,
            foldList=NULL, trueClass, predictedClass)
}
\arguments{
  \item{nfold}{
The number of groups of the cross-validation. Ignored if \code{leaveOneOut=TRUE}.
}
  \item{ntimes}{
The number of repeated runs of the cross-validation.
}
  \item{leaveOneOut}{
If this is true, a leave-one-out cross-validation is performed, i.e. each sample is left out once in the training phase and used as a test sample
}

  \item{stratified}{
If set to true, a stratified cross-validation is carried out. That is, the percentage of samples from different classes in the cross-validation folds corresponds to the class sizes in the complete data set. If set to false, the folds may be unbalanced.
}

  \item{foldList}{
If this parameter is set, the other cross-validation parameters (\code{ntimes}, \code{nfold}, \code{leaveOneOut}, \code{stratified}) are ignored. Instead, the precalculated cross-validation partition supplied in \code{foldList} is used. This allows for using the same cross-validation experiment in multiple \code{\link{tunePareto}} calls. Partitions can be generated using \code{\link{generateCVRuns}}.
}

  \item{caseClass}{The class containing the positive samples for the calculation of specificity and sensitivity. All samples with different class labels are regarded as controls (negative samples).}

  \item{trueClass}{
When calculating the confusion of two classes, the class to which a sample truly belongs.
}

  \item{predictedClass}{
When calculating the confusion of two classes, the class to which a sample is erreneously assigned.
}
}
\details{
The functions do not calculate the objectives, but return a structure of class \code{TuneParetoObjectives} that provides all information on the objective function for later use in \code{\link{tunePareto}}.

The behaviour of the functions in \code{\link{tunePareto}} is as follows:

The reclassification functions train the classifiers with the full data set. Afterwards, the classifiers are applied to the same data set. \code{reclassError} measures the fraction of misclassified samples, \code{reclassWeightedError} calculates the sum of fractions of misclassified samples in each class weighted by the class size. \code{reclassSensitivity} measures the sensitivity, and \code{reclassSpecificity} measures the specificity of the reclassification experiment. \code{reclassConfusion} calculates the fraction of samples in \code{trueClass} that have been confused with \code{predictedClass}. \code{reclassError}, \code{reclassWeightedError} and \code{reclassConfusion} are minimization objectives, whereas \code{reclassSensitivity} and \code{reclassSpecificity} are maximization objectives.

The cross-validation functions partition the samples in the data set into a number of groups (depending on \code{nfold} and \code{leaveOneOut}). Each of these groups is left out once in the training phase and used for prediction. The whole procedure is usually repeated several times (as specified in \code{ntimes}), and the results are averaged. Similar to the reclassification functions, \code{cvError} calculates the average fraction of misclassified samples over the runs, and \code{cvWeightedError} calculates the mean sum of fractions of misclassified samples in each class weighted by the class size. \code{cvErrorVariance} calculates the variance of the cross-validation error. \code{cvSensitivity} calculates the average sensitivity and \code{cvSpecificity} calculates the average specificity. \code{cvConfusion} calculates the average fraction of samples in \code{trueClass} that have been confused with \code{predictedClass}. \code{cvError}, \code{cvWeightedError}, \code{cvErrorVariance} and \code{cvConfusion} are minimization objectives, and \code{cvSensitivity} and \code{cvSpecificity} are maximization objectives.
}
\value{
An object of class \code{TuneParetoObjective} representing the objective function. For more details, see \code{\link{createObjective}}.
}

\seealso{
\code{\link{createObjective}}, \code{\link{tunePareto}}, \code{\link{generateCVRuns}}
}
\examples{

# build a list of objective functions
objectiveFunctions <- list(cvError(10, 10),
                           reclassSpecificity(caseClass="setosa"), 
                           reclassSensitivity(caseClass="setosa"))

# pass them to tunePareto
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.knn(),
                 k = c(3,5,7,9),
                 objectiveFunctions = objectiveFunctions))
}

\keyword{objective function
         multi-objective optimization
         MOO}
