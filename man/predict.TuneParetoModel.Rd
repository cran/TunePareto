\name{predict.TuneParetoModel}
\alias{predict.TuneParetoModel}
\title{
Prediction method for TuneParetoClassifier objects
}
\description{
S3 method that predicts the labels of unknown samples using a trained \code{TunePareteModel} model of a \code{TuneParetoClassifier} object.
}
\usage{
predict.TuneParetoModel(object, newdata, ...)
}

\arguments{
  \item{object}{
  A \code{TuneParetoTraining} object as returned by \code{\link{trainTuneParetoClassifier}}.
}
  \item{newdata}{
The samples whose class labels are predicted. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{\dots}{
Further parameters for the predictor. These must be defined in the \code{predictorParamNames} argument of \code{\link{tuneParetoClassifier}}.
}
}
\value{
Returns a vector of class labels for the samples in \code{newdata}-
}

\seealso{
\code{\link{tuneParetoClassifier}}, \code{\link{predefinedClassifiers}}, \code{\link{trainTuneParetoClassifier}}
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

\keyword{prediction}
