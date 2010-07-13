\name{callClassifier}
\alias{callClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generic function to call a classifier
}
\description{
This internal helper function can be used in the precalculation function of an objective (see \code{\link{createObjective}}). It calls a classifier and returns the predicted labels of the supplied test samples. You do not have to use it directly if you do not create own objectives.
}
\usage{
callClassifier(trainData, 
               trainLabels, 
               testData, 
               classifier,
               classifierParams, 
               predictor, 
               predictorParams, 
               useFormula = FALSE, 
               formulaName = "formula", 
               trainDataName = "x", 
               trainLabelName = "y", 
               testDataName = "newData", 
               modelName = "object")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{trainData}{
The samples used to train the classifier. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}
  \item{trainLabels}{
A vector of class labels for the training samples.
}
  \item{testData}{
The samples whose labels are predicted by the classifier. This is usually a matrix or data frame with the samples in the rows and the features in the columns.
}   
  \item{classifier}{
The classification function to use. If \code{predictor} is \code{NULL}, this function is an all-in-one classification method that receives both training data and test data and returns the predicted labels for the test data. If \code{predictor} is not \code{NULL}, this is the training function of the classifier that builds a model from the training data. This model is then passed to \code{predictor} along with the test data to obtain the predicted labels for the test data.
}

  \item{classifierParams}{
A named list of parameter assignments for the classifier.
}

  \item{predictor}{
If the classification method consists of separate training and prediction functions, this points to the prediction function that receives a model and the test data as inputs and returns the predicted class labels.
}
  
  \item{predictorParams}{
If \code{predictor != NULL}, a named list of parameter assignments for the predictor.
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

\value{
Returns a vector of predicted class labels for the samples in \code{testData}
}

\seealso{
\code{\link{createObjective}}
}

