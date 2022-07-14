\name{trainTuneParetoClassifier}
\alias{trainTuneParetoClassifier}
\title{
Train a TunePareto classifier
}
\description{
Trains a classifier wrapped in a \code{TuneParetoClassifier} object. The trained classifier model can then be passed to \code{\link{predict.TuneParetoModel}}.
}
\usage{
trainTuneParetoClassifier(classifier, trainData, trainLabels, ...)
}
\arguments{
  \item{classifier}{
A \code{TuneParetoClassifier} object as returned by \code{\link{tuneParetoClassifier}} or one of the predefined classification functions (see \code{\link{predefinedClassifiers}}).
}
  \item{trainData}{
The data set to be used for the classifier training. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{trainLabels}{
A vector of class labels for the samples in \code{trainData}.
}
  \item{\dots}{
Further parameters to be passed to the classifier. These must be parameters specified in the \code{classifierParameterNames} parameter of \code{\link{tuneParetoClassifier}} and usually correspond to the tuned parameters.
}
}
\value{
Returns an object of class \code{TuneParetoModel} with the following entries
\item{classifier}{The classifier object supplied in the \code{classifier} parameter}
\item{classifierParams}{The additional parameters supplied to the classifier in the \code{\dots} parameter}
\item{trainData}{If \code{classifier} is an all-in-one classifier without a separate prediction method, this stores the input training data.}
\item{trainLabels}{If \code{classifier} is an all-in-one classifier without a separate prediction method, this stores the input training labels.}
\item{model}{If \code{classifier} consists of separate training and prediction methods, this contains the trained classifier model.}
}

\seealso{
\code{\link{tuneParetoClassifier}}, \code{\link{predefinedClassifiers}}, \code{\link{predict.TuneParetoModel}}
}
\examples{

# train an SVM classifier
cl <- tunePareto.svm()
tr <- trainTuneParetoClassifier(cl, 
                                iris[,-ncol(iris)], 
                                iris[,ncol(iris)], 
                                cost=0.001)

# re-apply the classifier to predict the training data
print(iris[,ncol(iris)])
print(predict(tr, iris[,-ncol(iris)]))
}

\keyword{training}
\keyword{classifier training}
