\name{tuneParetoClassifier}
\alias{tuneParetoClassifier}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Create a classifier object
}
\description{
Creates a wrapper object mapping all information necessary to call a classifier which can be passed to \code{\link{tunePareto}}.
}
\usage{
tuneParetoClassifier(name,
                     classifier, 
                     classifierParamNames = NULL, 
                     predefinedClassifierParams = NULL, 
                     predictor, 
                     predictorParamNames = NULL, 
                     predefinedPredictorParams = NULL, 
                     useFormula = FALSE, 
                     formulaName = "formula", 
                     trainDataName = "x", 
                     trainLabelName = "y", 
                     testDataName = "newData", 
                     modelName = "object", 
                     requiredPackages)
}
\arguments{

  \item{name}{
A human-readable name of the classifier
}

  \item{classifier}{
The classification function to use. If \code{predictor} is \code{NULL}, this function is an all-in-one classification method that receives both training data and test data and returns the predicted labels for the test data. If \code{predictor} is not \code{NULL}, this is the training function of the classifier that builds a model from the training data. This model is then passed to \code{predictor} along with the test data to obtain the predicted labels for the test data.
}
  \item{classifierParamNames}{
A vector of names of possible arguments for \code{classifier}.
}
  \item{predefinedClassifierParams}{
A named list of default values for the classifier parameters.
}
  \item{predictor}{
If the classification method consists of separate training and prediction functions, this points to the prediction function that receives a model and the test data as inputs and returns the predicted class labels.
}
  \item{predictorParamNames}{
If \code{predictor != NULL}, a vector of names of possible arguments for \code{predictor}.
}
  \item{predefinedPredictorParams}{
If \code{predictor != NULL}, a named list of default values for the parameters of \code{predictor}.
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
  \item{requiredPackages}{
A vector containing the names of packages that are required to run the classifier. These packages are loaded automatically when running the classifier using \code{\link{tunePareto}}. They are also loaded in the \pkg{snowfall} cluster if necessary.
}
}
\details{
TunePareto classifier objects are wrappers containing all information necessary to run the classifier, including the training and prediction function, the required packages, and the names of certain arguments. \pkg{TunePareto} provides a set of predefined objects for state-of-the-art classifiers (see \code{\link{predefinedClassifiers}}).

The main \code{\link{tunePareto}} routine evaluates \code{TuneParetoClassifier} objects to call the training and prediction methods. Furthermore, direct calls to the classifiers are possible using \if{latex}{\cr} \code{\link{trainTuneParetoClassifier}} and \code{\link{predict.TuneParetoModel}}.
}
\value{
An object of class \code{TuneParetoClassifier} with components corresponding to the above parameters.
}

\seealso{
\code{\link{trainTuneParetoClassifier}}, \code{\link{predict.TuneParetoModel}}, \code{\link{tunePareto}}, \code{\link{predefinedClassifiers}}
}
\examples{
  # equivalent to tunePareto.svm()
  cl <- tuneParetoClassifier(name = "svm",
                             classifier = svm, 
                             predictor = predict, 
                             classifierParamNames = c("kernel", "degree", "gamma",
                                                      "coef0", "cost", "nu",
                                                      "class.weights", "cachesize", 
                                                      "tolerance", "epsilon",
                                                      "subset", "na.action"),
                              useFormula = FALSE,
                              trainDataName = "x",
                              trainLabelName = "y",
                              testDataName = "newdata",
                              modelName = "object",
                              requiredPackages="e1071")
  
  # call TunePareto with the classifier
  print(tunePareto(classifier = cl,
                   data = iris[, -ncol(iris)], 
                   labels = iris[, ncol(iris)],
                   cost = c(0.001,0.01,0.1,1,10), 
                   objectiveFunctions=list(cvError(10, 10),
                                           cvSpecificity(10, 10, caseClass="setosa"))))                           
}

\keyword{TuneParetoClassifier}
