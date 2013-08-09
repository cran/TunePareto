\name{generateCVRuns}
\alias{generateCVRuns}
\title{
Generate cross-validation partitions
}
\description{
This function generates a set of partitions for a cross-validation. It can be employed if the same cross-validation settings should be used in the objective functions of several experiments. The resulting fold list can be passed to the cross-validation objective functions (see \if{latex}{\cr}\code{\link{predefinedObjectiveFunctions}}) and the internal cross-validation precalculation function \code{\link{crossValidation}}. 
}
\usage{
generateCVRuns(labels, 
               ntimes = 10, 
               nfold = 10, 
               leaveOneOut = FALSE, 
               stratified = FALSE)
}
\arguments{
 
  \item{labels}{
A vector of class labels of the data set to be used for the cross-validation.
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

}

\value{
A list with \code{ntimes} elements, each representing a cross-validation run. Each of the runs is a list of \code{nfold} vectors specifying the indices of the samples to be left out in the folds.
}

\seealso{
\code{\link{predefinedObjectiveFunctions}}, \code{\link{crossValidation}}
}
\examples{

# precalculate the cross-validation partitions
foldList <- generateCVRuns(labels = iris[, ncol(iris)],
                           ntimes = 10,
                           nfold = 10,
                           stratified=TRUE)

 # build a list of objective functions
objectiveFunctions <- list(cvError(foldList=foldList),
                           cvSensitivity(foldList=foldList,caseClass="setosa"))

# pass them to tunePareto
print(tunePareto(data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 classifier = tunePareto.knn(),
                 k = c(3,5,7,9),
                 objectiveFunctions = objectiveFunctions))
}
\keyword{cross-validation}
