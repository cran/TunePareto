\name{precalculation}
\alias{precalculation}
\alias{reclassification}
\alias{crossValidation}
\title{
Predefined precalculation functions for objectives
}
\description{
These predefined precalculation functions can be employed to create own objectives using \code{\link{createObjective}}. They perform a reclassification or a cross-validation and return the true labels and the predictions.
}
\usage{
reclassification(data, labels, 
                 classifier, classifierParams, predictorParams)

crossValidation(data, labels, 
                classifier, classifierParams, predictorParams, 
                ntimes = 10, nfold = 10, 
                leaveOneOut = FALSE, stratified = FALSE,
                foldList = NULL)
}
\arguments{
\item{data}{
The data set to be used for the precalculation. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{labels}{
A vector of class labels for the samples in \code{data}.
}
  \item{classifier}{
A \code{TuneParetoClassifier} wrapper object containing the classifier to tune. A number of state-of-the-art classifiers are included in \pkg{TunePareto} \if{latex}{\cr} (see \code{\link{predefinedClassifiers}}). Custom classifiers can be employed using \code{\link{tuneParetoClassifier}}.}

  \item{classifierParams}{
A named list of parameter assignments for the training routine of the classifier.
}

  \item{predictorParams}{
If the classifier consists of separate training and prediction functions, a named list of parameter assignments for the predictor function.
}

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

}
\details{
\code{reclassification} trains the classifier with the full data set. Afterwards, the classifier is applied to the same data set. 

\code{crossValidate} partitions the samples in the data set into a number of groups (depending on \code{nfold} and \code{leaveOneOut}). Each of these groups is left out once in the training phase and used for prediction. The whole procedure is repeated several times (as specified in \code{ntimes}). 

}
\value{
\code{reclassification} returns a list with the following components:
\describe{
\item{trueLabels}{The original labels of the dataset as supplied in \code{labels}}
\item{predictedLabels}{A vector of predicted labels of the data set}
\item{model}{The \code{TuneParetoModel} object resulting from the classifier training}
}
\code{crossValidation} returns a nested list structure. At the top level, there is one list element for each run of the cross-validation. Each of these elements consists of a list of sub-structures for each fold. The sub-structures have the following components:
\describe{
\item{trueLabels}{The original labels of the test samples in the fold}
\item{predictedLabels}{A vector of predicted labels of the test samples in the fold}
\item{model}{The \code{TuneParetoModel} object resulting from the classifier training in the fold}
}
That is, for a cross-validation with \code{n} runs and \code{m} folds, there are \code{n} top-level lists, each having \code{m} sub-lists comprising the true labels and the predicted labels.
}


\seealso{
\code{\link{createObjective}}, \code{\link{generateCVRuns}}. 
}
\examples{
\donttest{
# create new objective minimizing the
# false positives of a reclassification

cvFalsePositives <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, foldList=NULL, caseClass)
{
  return(createObjective(
            precalculationFunction = "crossValidation",
            precalculationParams = list(nfold=nfold, 
                                        ntimes=ntimes, 
                                        leaveOneOut=leaveOneOut,
                                        foldList=foldList),
            objectiveFunction = 
            function(result, caseClass)
            {
             
              # take mean value over the cv runs
              return(mean(sapply(result,
                    function(run)
                    # iterate over runs of cross-validation
                    {
                      # extract all predicted labels in the folds
                      predictedLabels <- 
                            unlist(lapply(run,
                                         function(fold)fold$predictedLabels))
    
                      # extract all true labels in the folds
                      trueLabels <- 
                            unlist(lapply(run,
                                          function(fold)fold$trueLabels))
                      
                      # calculate number of false positives in the run
                      return(sum(predictedLabels == caseClass & 
                                 trueLabels != caseClass))
                    })))
            },
            objectiveFunctionParams = list(caseClass=caseClass),
            direction = "minimize",        
            name = "CV.FalsePositives"))                  
}

# use the objective in an SVM cost parameter tuning on the 'iris' data set
r <- tunePareto(data = iris[, -ncol(iris)], 
                labels = iris[, ncol(iris)],
                classifier = tunePareto.svm(),
                cost = c(0.001,0.005,0.01,0.05,0.1,0.5,1,5,10,50),
                objectiveFunctions=list(cvFalsePositives(10, 10, caseClass="setosa")))
print(r)
}
}

\keyword{objective function
         multi-objective optimization
         MOO}
