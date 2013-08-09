\name{predefinedObjectiveFunctions}
\alias{predefinedObjectiveFunctions}
\alias{reclassAccuracy}
\alias{reclassError}
\alias{reclassWeightedError}
\alias{reclassSensitivity}
\alias{reclassRecall}
\alias{reclassTruePositive}
\alias{reclassSpecificity}
\alias{reclassTrueNegative}
\alias{reclassFallout}
\alias{reclassFalsePositive}
\alias{reclassMiss}
\alias{reclassFalseNegative}
\alias{reclassPrecision}
\alias{reclassPPV}
\alias{reclassNPV}
\alias{reclassConfusion}
\alias{cvAccuracy}
\alias{cvError}
\alias{cvErrorVariance}
\alias{cvWeightedError}
\alias{cvSensitivity}
\alias{cvRecall}
\alias{cvTruePositive}
\alias{cvSpecificity}
\alias{cvTrueNegative}
\alias{cvFallout}
\alias{cvFalsePositive}
\alias{cvMiss}
\alias{cvFalseNegative}
\alias{cvPrecision}
\alias{cvPPV}
\alias{cvNPV}
\alias{cvConfusion}
\title{
Predefined objective functions for parameter tuning
}
\description{
Predefined objective functions that calculate several performance criteria of reclassification or cross-validation experiments.
}
\usage{
reclassAccuracy(saveModel = FALSE)
reclassError(saveModel = FALSE)
reclassWeightedError(saveModel = FALSE)
reclassSensitivity(caseClass, saveModel = FALSE)
reclassRecall(caseClass, saveModel = FALSE)
reclassTruePositive(caseClass, saveModel = FALSE)
reclassSpecificity(caseClass, saveModel = FALSE)
reclassTrueNegative(caseClass, saveModel = FALSE)
reclassFallout(caseClass, saveModel = FALSE)
reclassFalsePositive(caseClass, saveModel = FALSE)
reclassMiss(caseClass, saveModel = FALSE)
reclassFalseNegative(caseClass, saveModel = FALSE)
reclassPrecision(caseClass, saveModel = FALSE)
reclassPPV(caseClass, saveModel = FALSE)
reclassNPV(caseClass, saveModel = FALSE)
reclassConfusion(trueClass, predictedClass, saveModel = FALSE)

cvAccuracy(nfold = 10, ntimes = 10,
           leaveOneOut = FALSE, stratified = FALSE,
           foldList = NULL,
           saveModel = FALSE) 
cvError(nfold = 10, ntimes = 10, 
        leaveOneOut = FALSE, stratified=FALSE, 
        foldList=NULL,
        saveModel = FALSE)
cvErrorVariance(nfold = 10, ntimes = 10, 
                leaveOneOut = FALSE, stratified=FALSE, 
                foldList=NULL,
                saveModel = FALSE)
cvWeightedError(nfold = 10, ntimes = 10, 
                leaveOneOut = FALSE, stratified=FALSE, 
                foldList=NULL,
                saveModel = FALSE)
cvSensitivity(nfold = 10, ntimes = 10, 
              leaveOneOut = FALSE, stratified=FALSE, 
              foldList=NULL, caseClass,
              saveModel = FALSE)
cvRecall(nfold = 10, ntimes = 10,
         leaveOneOut = FALSE, stratified=FALSE,
         foldList=NULL, caseClass,
         saveModel = FALSE)
cvTruePositive(nfold = 10, ntimes = 10, 
               leaveOneOut = FALSE, stratified=FALSE, 
               foldList=NULL, caseClass,
               saveModel = FALSE)
cvSpecificity(nfold = 10, ntimes = 10, 
              leaveOneOut = FALSE, stratified=FALSE, 
              foldList=NULL, caseClass,
              saveModel = FALSE)
cvTrueNegative(nfold = 10, ntimes = 10, 
               leaveOneOut = FALSE, stratified=FALSE, 
               foldList=NULL, caseClass,
               saveModel = FALSE)
cvFallout(nfold = 10, ntimes = 10, 
          leaveOneOut = FALSE, stratified=FALSE, 
          foldList=NULL, caseClass,
          saveModel = FALSE)
cvFalsePositive(nfold = 10, ntimes = 10, 
                leaveOneOut = FALSE, stratified=FALSE, 
                foldList=NULL, caseClass,
                saveModel = FALSE)          
cvMiss(nfold = 10, ntimes = 10, 
       leaveOneOut = FALSE, stratified=FALSE, 
       foldList=NULL, caseClass,
       saveModel = FALSE)
cvFalseNegative(nfold = 10, ntimes = 10,
                leaveOneOut = FALSE, stratified=FALSE,
                foldList=NULL, caseClass,
                saveModel = FALSE)
cvPrecision(nfold = 10, ntimes = 10,
            leaveOneOut = FALSE, stratified=FALSE,
            foldList=NULL, caseClass,
            saveModel = FALSE)
cvPPV(nfold = 10, ntimes = 10,
      leaveOneOut = FALSE, stratified=FALSE,
      foldList=NULL, caseClass,
      saveModel = FALSE)
cvNPV(nfold = 10, ntimes = 10,
      leaveOneOut = FALSE, stratified=FALSE,
      foldList=NULL, caseClass,
      saveModel = FALSE)
cvConfusion(nfold = 10, ntimes = 10, 
            leaveOneOut = FALSE, stratified=FALSE,
            foldList=NULL, trueClass, predictedClass,
            saveModel = FALSE)
}
\arguments{
  \item{nfold}{
The number of groups of the cross-validation. Ignored if \code{leaveOneOut=TRUE}.
}
  \item{ntimes}{
The number of repeated runs of the cross-validation. Ignored if \code{leaveOneOut=TRUE}.
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

  \item{saveModel}{
If set to true, the trained model(s) are stored to the \code{additionalData} component of the resulting \code{TuneParetoResult} objects (see \code{\link{tunePareto}} for details). In case of a reclassification, a single model is stored. In case of a cross-validation, a list of length \code{nruns}, each containing a sub-list of \code{nfold} models, is stored. If the size of a model is large, setting \code{saveModel = TRUE} can result in a high memory consumption. As the model information is the same for all reclassification objectives or for cross-validation objectives with the same parameters, it is usually sufficient to set \code{saveModel=TRUE} for only one of the objective functions.
}  
}
\details{
The functions do not calculate the objectives directly, but return a structure of class \code{TuneParetoObjectives} that provides all information on the objective function for later use in \code{\link{tunePareto}}.

The behaviour of the functions in \code{\link{tunePareto}} is as follows:

The reclassification functions train the classifiers with the full data set. Afterwards, the classifiers are applied to the same data set. \code{reclassAccuracy} measures the fraction of correctly classified samples, while \code{reclassError} calculates the fraction of misclassified samples. \code{reclassWeightedError} calculates the sum of fractions of misclassified samples in each class weighted by the class size. \code{reclassSensitivity} measures the sensitivity, and \code{reclassSpecificity} measures the specificity of the reclassification experiment. \code{reclassTruePositive} and \code{reclassRecall} are aliases for \code{reclassSensitivity}, and \code{reclassTrueNegative} is an alias for \code{reclassSpecificity}. \code{reclassFallout} and its equivalent alias \code{reclassFalsePositive} give the ratio of false positives to all negative samples, and \code{reclassMiss} and its alias \code{reclassFalseNegative} measure the ratio of false negatives to all positive samples. \code{reclassPrecision} calculates the precision of the reclassification experiment, i.e. the ratio of true positives to all samples classified as positive. This is equivalent to the positive predictive value (\code{reclassPPV}). \code{reclassNPV} measures the negative predictive value, i.e. the ratio of true negatives to all samples classified as negative. \code{reclassConfusion} calculates the fraction of samples in \code{trueClass} that have been confused with \code{predictedClass}.

\code{reclassError}, \code{reclassWeightedError}, \code{reclassFallout}, \code{reclassFalsePositive}, \code{reclassMiss}, \code{reclassFalsePositive}  and \code{reclassConfusion} are minimization objectives, whereas \code{reclassAccuracy}, \code{reclassSensitivity}, \code{reclassTruePositive}, \code{reclassRecall}, \code{reclassSpecificity}, \code{reclassTrueNegative} \code{reclassPrecision}, \code{reclassPPV} and \code{reclassNPV} are maximization objectives.

The cross-validation functions partition the samples in the data set into a number of groups (depending on \code{nfold} and \code{leaveOneOut}). Each of these groups is left out once in the training phase and used for prediction. The whole procedure is usually repeated several times (as specified in \code{ntimes}), and the results are averaged. Similar to the reclassification functions, \code{cvAccuracy} calculates the fraction of correctly classified samples over the runs, \code{cvError} calculates the average fraction of misclassified samples over the runs, and \code{cvWeightedError} calculates the mean sum of fractions of misclassified samples in each class weighted by the class size. \code{cvErrorVariance} calculates the variance of the cross-validation error. \code{cvSensitivity}, \code{cvRecall} and \code{cvTruePositive} calculate the average sensitivity, and \code{cvSpecificity} and \code{cvTrueNegative} calculate the average specificity. \code{cvFallout} and \code{cvFalsePositive} calculate the average false positive rate over the runs. \code{cvMiss} and \code{cvFalseNegative} calculate the average false negative rate over the runs. \code{cvPrecision} and \code{cvPPV} calculate the average precision/positive predictive rate. \code{cvNPV} gives the average negative predictive rate over all runs. \code{cvConfusion} calculates the average fraction of samples in \code{trueClass} that have been confused with \code{predictedClass}.

\code{cvError}, \code{cvWeightedError}, \code{cvErrorVariance}, \code{cvFallout}, \code{cvFalsePositive}, \code{cvMiss}, \code{cvFalseNegative} and \code{cvConfusion} are minimization objectives, and \code{cvAccuracy}, \code{cvSensitivity}, \code{cvRecall}, \code{cvTruePositive}, \code{cvSpecificity}, \code{cvTrueNegative}, \code{cvPrecision}, \code{cvPPV} and \code{cvNPV} are maximization objectives.
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
