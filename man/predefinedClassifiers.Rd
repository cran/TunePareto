\name{predefinedClassifiers}
\alias{predefinedClassifiers}
\alias{tunePareto.knn}
\alias{tunePareto.svm}
\alias{tunePareto.tree}
\alias{tunePareto.randomForest}
\alias{tunePareto.NaiveBayes}
\title{
TunePareto wrappers for certain classifiers
}
\description{
Creates TunePareto classifier objects for the k-Nearest Neighbour classifier, support vector machines, and trees.
}
\usage{
tunePareto.knn()
               
tunePareto.svm()
               
tunePareto.tree()
                
tunePareto.randomForest()

tunePareto.NaiveBayes()
}
\details{

\code{tunePareto.knn} encapsulates a k-Nearest Neighbour classifier as defined in \code{link[class]{knn}} in package \pkg{class}. The classifier allows for supplying and tuning the following parameters of \code{link[class]{knn}}: 

\code{k, l, use.all}

\code{tunePareto.svm} encapsulates the support vector machine \code{\link[e1071]{svm}} classifier in package \pkg{e1071}. The classifier allows for supplying and tuning the following parameters: 

\code{kernel, degree, gamma, 
      coef0, cost, nu, 
      class.weights, cachesize, 
      tolerance, epsilon, 
      scale, shrinking, fitted,
      subset, na.action}

\code{tunePareto.tree} encapsulates the CART classifier \code{\link[tree]{tree}} in package \pkg{tree}. The classifier allows for supplying and tuning the following parameters: 

\code{weights, subset, 
      na.action, method,
      split, mincut, minsize, mindev}
      
as well as the \code{type} parameter of \code{\link[tree]{predict.tree}}.

\code{tunePareto.randomForest} encapsulates the \code{\link[randomForest]{randomForest}} classifier in package \pkg{randomForest}. The classifier allows for supplying and tuning the following parameters:

\code{subset, na.action,
      ntree,  mtry,
      replace, classwt, 
      cutoff, strata,
      sampsize, nodesize,
      maxnodes}
      
\code{tunePareto.NaiveBayes} encapsulates the \code{\link[klaR]{NaiveBayes}} classifier in package \pkg{klaR}. The classifier allows for supplying and tuning the following parameters:

\code{prior, usekernel, fL, subset,
      na.action, bw, adjust, kernel, weights,
      window, width, give.Rkern, n,
      from, to, cut, na.rm}

}

\value{
Returns objects of class \code{TuneParetoClassifier} as described in \code{\link{tuneParetoClassifier}}. These can be passed to functions like \code{\link{tunePareto}} or \code{\link{trainTuneParetoClassifier}}.
}

\seealso{
\code{\link{tuneParetoClassifier}}, \code{\link{tunePareto}}, \code{\link{trainTuneParetoClassifier}}
}
\examples{
\donttest{
# tune a k-NN classifier with different 'k' and 'l' 
# on the 'iris' data set
print(tunePareto(classifier = tunePareto.knn(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 k = c(5,7,9),
                 l = c(1,2,3),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa"))))
                 
# tune an SVM with different costs on 
# the 'iris' data set
# using Halton sequences for sampling
print(tunePareto(classifier = tunePareto.svm(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 cost = as.interval(0.001,10),
                 sampleType = "halton",
                 numCombinations=20,                 
                 objectiveFunctions=list(cvWeightedError(10, 10),
                                         cvSensitivity(10, 10, caseClass="setosa"))))

# tune a CART classifier with different 
# splitting criteria on the 'iris' data set
print(tunePareto(classifier = tunePareto.tree(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 split = c("deviance","gini"),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvErrorVariance(10, 10))))

# tune a Random Forest with different numbers of trees 
# on the 'iris' data set
print(tunePareto(classifier = tunePareto.randomForest(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 ntree = seq(50,300,50),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa"))))

# tune a Naive Bayes classifier with different kernels
# on the 'iris' data set
print(tunePareto(classifier = tunePareto.NaiveBayes(),
                 data = iris[, -ncol(iris)], 
                 labels = iris[, ncol(iris)],
                 kernel = c("gaussian", "epanechnikov", "rectangular",
                            "triangular", "biweight",
                            "cosine", "optcosine"),
                 objectiveFunctions=list(cvError(10, 10),
                                         cvSpecificity(10, 10, caseClass="setosa"))))
                             
}
}

\keyword{TuneParetoClassifier knn svm tree randomForest NaiveBayes}
