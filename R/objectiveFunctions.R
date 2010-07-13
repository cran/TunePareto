#################################################################
# Objective functions and utilities for objective functions
#################################################################

# Creates a new objective from the specified objective
# function <objectiveFunction> and its precalculation
# method <precalculationFunction> with additional parameters
# <precalculationParams>. 
# <direction> specifies whether the objective is 
# minimized or maximized. 
# <name> is a short description of the objective.
#
# Returns an object of class TuneParetoObjective
createObjective <- function(precalculationFunction,
                            precalculationParams = NULL,
                            objectiveFunction,
                            objectiveFunctionParams = NULL,
                            direction=c("minimize","maximize"),
                            name)
{
  res <- list(precalculationFunction = match.fun(precalculationFunction),
              precalculationParams = precalculationParams,
              objectiveFunction = match.fun(objectiveFunction),
              objectiveFunctionParams = objectiveFunctionParams,
              minimize = (match.arg(direction, c("minimize","maximize")) == "minimize"),
              name = name)
  class(res) <- "TuneParetoObjective"
  return(res)
}

# Groups all objectives in <objectiveFunctionList> that
# use the same preprocessing with the same parameters.
# This speeds up calculations.
groupByPrecalculation <- function(objectiveFunctionList)
{
  # determine unique set of precalculation functions
  funcs  <- sapply(objectiveFunctionList,function(obj)obj$precalculationFunction)
  uniqueFuncs <- unique(funcs)
  resultList <- list()
  indices <- c()
  
  for (func in uniqueFuncs)
  {
    # find objectives with the same preprocessing as <func>
    newIndices <- which(sapply(funcs, function(f)identical(f,func)))
    indices <- c(indices,newIndices)
    
    # find unique set of parameters for the preprocessing
    subset <- objectiveFunctionList[newIndices]
    params <- lapply(subset,function(obj)obj$precalculationParams)
    uniqueParams <- unique(params)
    newObjectives <- lapply(uniqueParams,function(param)
                            {
                              # join objectives with the same parameters
                              toJoin <- subset[sapply(params, function(p)identical(p,param))]
                              
                              res <- toJoin[[1]]

                              if (length(toJoin) > 1)
                              {
                                res$objectiveFunctionParams[[1]] <- list(res$objectiveFunctionParams)
                                for (i in 2:length(toJoin))
                                {
                                  res$objectiveFunction <- c(res$objectiveFunction, toJoin[[i]]$objectiveFunction)
                                  res$objectiveFunctionParams[[i]] <- 
                                        list(toJoin[[i]]$objectiveFunctionParams)
                                  res$minimize <- c(res$minimize, toJoin[[i]]$minimize)
                                  res$name <- c(res$name, toJoin[[i]]$name)
                                }
                              }
                              res
                            })
                            
    # add new joint objective to result list
    resultList[(length(resultList) + 1):(length(resultList) + length(newObjectives))] <- newObjectives
  }
  
  # calculate permutation of original indices for subsequent reordering
  permutation <- indices
  for (i in 1:length(permutation))
    permutation[indices[i]] <- i
    
  return(list(grouping=resultList, permutation=permutation))
  
}

# Predefined precalculation function that trains the classifier on the whole
# data set and predicts the same samples with the trained classifier.
# For parameters, see callClassifier.
#
# Returns a list containing a vector of predicted labels and a vector of true labels
reclassification <- function(data, labels, classifier, classifierParams,
                             predictor, predictorParams,
                             useFormula = FALSE, formulaName = "formula",
                             trainDataName = "x", trainLabelName = "y", 
                             testDataName = "newdata", modelName = "object")
{
  predictedLabs <- callClassifier(data, labels, data, classifier, classifierParams,
                                    predictor, predictorParams, useFormula, formulaName, 
                                    trainDataName, trainLabelName, 
                                    testDataName, modelName)
  
  #labels <- as.integer(as.character(labels))

  res <- list(predictedLabels=predictedLabs, trueLabels=labels)
  class(res) <- "ClassificationResult"
                             
  return(res)
}

# Predefined precalculation function that performs a cross-validation on <data>.
# <ntimes> is the number of repetitions of the cross-validation.
# <nfold> is the number of groups in each cross-validation run.
# If <leaveOneOut> is true, a leave-one-out cross-validation is performed.
# For further parameters, see callClassifier.
# If <stratified> is true, a stratified cross-validation is carried out.
#
# Returns a list containing a sub-list for each run. Each of these sub-lists contains
# a vector of true labels and predicted labels for each fold.
crossValidation <- function(data, labels, classifier, classifierParams,
                            predictor, predictorParams,
                            useFormula = FALSE, formulaName = "formula",
                            trainDataName = "x", trainLabelName = "y", 
                            testDataName = "newdata", modelName = "object",
                            ntimes = 10, nfold = 10, leaveOneOut=FALSE, stratified = FALSE)
{
  numSamples <- nrow(data)
  
  res <- lapply(1:ntimes,function(run)
  # for each run
  {
    # calculate folds
	  if (leaveOneOut)
	  {
		  indices <- as.list(1:numSamples)
	  }
	  else
	  {
		if(stratified)
		{
			classes <- unique(labels)
			sing.perm <- lapply(classes, function(cl){
				index <- which(labels == cl)
				sample(index, length(index))
			})
			permut <- unlist(sing.perm)
			indices <- lapply(1:nfold,function(i){c()})
			for(i in 1:numSamples)
			{
				k = i%%nfold
				if(k==0)
				 k = nfold
				
				indices[[k]] <- c(indices[[k]], permut[i])
			}
		}
		else
		{
		  # permute the indices of the samples
		  permut <- sample(1:numSamples, numSamples,replace=FALSE)
		  indices <- lapply(1:nfold, function(i)
		  {
			  # split the samples in nfold groups
			  permut[seq(i, numSamples, nfold)]
		  })
		}
	  }
	  
	  return(lapply(indices, function(fold)
	  # for each fold
	  {
	    # split up data
	    trainData <- data[-fold,,drop=FALSE]
	    trainLabels <- labels[-fold]
	    testData <- data[fold,,drop=FALSE]
	    testLabels <- labels[fold]
	    
	    # predict test data
	    res1 <- callClassifier(trainData, trainLabels, testData, 
	                           classifier, classifierParams, predictor, predictorParams,
                             useFormula, formulaName, trainDataName, trainLabelName, testDataName, modelName)
    
      # return the true labels and the predicted labels
      res1 <- list(predictedLabels=res1, trueLabels=testLabels)
      class(res1) <- "ClassificationResult"
	    return(res1)
	  }))
	})
	return(res)
}

# Predefined objective calculating the sensitivity
# of a reclassification experiment
reclassSensitivity <- function(caseClass)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass)
                                      {
                                        return(sum(result$predictedLabels[result$trueLabels == caseClass] 
                                                    == caseClass,na.rm=TRUE) / 
                                               sum(result$trueLabels == caseClass))
                                      },
                  objectiveFunctionParams = list(caseClass=caseClass),
                  direction="maximize",
                  name="Reclass.Sensitivity")
}

# Predefined objective calculating the specificity
# of a reclassification experiment
reclassSpecificity <- function(caseClass)
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result, caseClass)
                                      {
                                        return(sum(result$predictedLabels[result$trueLabels != caseClass] 
                                                    != caseClass,na.rm=TRUE) / 
                                               sum(result$trueLabels != caseClass))
                                      },
                  objectiveFunctionParams = list(caseClass=caseClass),                                      
                  direction="maximize",
                  name="Reclass.Specificity")
}

# Predefined objective calculating the error percentage
# of a reclassification experiment
reclassError <- function()
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result)
                                      {
                                        sum(is.na(result$predictedLabels) | 
                                            result$predictedLabels != result$trueLabels) / 
                                        length(result$trueLabels)
                                      },
                  direction="minimize",
                  name="Reclass.Error")
}

# Predefined objective calculating the weighted error percentage
# of a reclassification experiment 
reclassWeightedError <- function()
{
  createObjective(precalculationFunction = reclassification,
                  precalculationParams = NULL,
                  objectiveFunction = function(result)
                                      {
                                        classes <- sort(unique(result$trueLabels))
                                        sum(sapply(classes,function(cl)
                                            {
                                              sum(result$trueLabels == cl &
                                                  result$predictedLabels != cl)/
                                              sum(result$trueLabels == cl)
                                            })) / length(classes)
                                      },
                  direction="minimize",
                  name="Reclass.WeightedError")
}

# Predefined objective calculating the error percentage
# of a cross-validation experiment
cvError <- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold=nfold, ntimes=ntimes, leaveOneOut=leaveOneOut,stratified = stratified),
                    objectiveFunction = function(result)
                                        {
                                          numSamples <- sum(sapply(result[[1]], function(fold)length(fold$trueLabels)))
                                          return(mean(sapply(result,function(run)
                                                      {
                                                        sum(unlist(lapply(run,function(fold)
                                                        {
                                                          is.na(fold$predictedLabels) | 
                                                          fold$predictedLabels != fold$trueLabels
                                                        })))
                                                      }))/numSamples)
                                        },
                    direction="minimize",
                    name="CV.Error")
}

cvErrorVariance<- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold=nfold, ntimes=ntimes, leaveOneOut=leaveOneOut,stratified = stratified),
                    objectiveFunction = function(result)
                                        {
                                          numSamples <- sum(sapply(result[[1]], function(fold)length(fold$trueLabels)))
                                          return(var(sapply(result,function(run)
                                                      {
                                                        sum(unlist(lapply(run,function(fold)
                                                        {
                                                          is.na(fold$predictedLabels) | 
                                                          fold$predictedLabels != fold$trueLabels
                                                        })))
                                                      }))/numSamples)
                                        },
                    direction="minimize",
                    name="CV.Error.Var")
}

# Predefined objective calculating the error percentage
# of a cross-validation experiment
cvWeightedError<- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold=nfold, ntimes=ntimes, leaveOneOut=leaveOneOut,stratified = stratified),
                    objectiveFunction = function(result)
                                        {
                                          return(mean(sapply(result,
                                                 function(run)
                                          {
                                            predictedLabels <- unlist(lapply(run,
                                                                      function(fold)fold$predictedLabels))
                                            trueLabels <- unlist(lapply(run,
                                                                 function(fold)fold$trueLabels))
                                            classes <- sort(unique(trueLabels))
                                            
                                            sum(sapply(classes,function(cl)
                                              {
                                                sum(trueLabels == cl &
                                                    predictedLabels != cl) /
                                                sum(trueLabels == cl)                                                    
                                                
                                              })) / length(classes)
                                          })))
                                        },
                    direction="minimize",
                    name="CV.WeightedError")
}

# Predefined objective calculating the sensitivity
# of a cross-validation experiment
cvSensitivity<- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE, caseClass)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold=nfold, ntimes=ntimes, leaveOneOut=leaveOneOut,stratified = stratified),
                    objectiveFunction = function(result, caseClass)
                                        {
                                          return(mean(sapply(result,
                                                 function(run)
                                                      {
                                                        predictedLabels <- unlist(lapply(run,
                                                                                  function(fold)fold$predictedLabels))
                                                        trueLabels <- unlist(lapply(run,
                                                                             function(fold)fold$trueLabels))
                                                        
                                                        return(sum(predictedLabels[trueLabels == caseClass] 
                                                                   == caseClass,na.rm=TRUE) / 
                                                               sum(trueLabels == caseClass))
                                                      })))
                                        },
                    objectiveFunctionParams = list(caseClass=caseClass),                                        
                    direction="maximize",
                    name="CV.Sensitivity")
}

# Predefined objective calculating the specificity
# of a cross-validation experiment
cvSpecificity<- function(nfold=10, ntimes=10, leaveOneOut=FALSE, stratified = FALSE,caseClass)
{
    createObjective(precalculationFunction = crossValidation,
                    precalculationParams = list(nfold=nfold, ntimes=ntimes, leaveOneOut=leaveOneOut,stratified = stratified),
                    objectiveFunction = function(result, caseClass)
                                        {
                                          return(mean(sapply(result,
                                                 function(run)
                                                      {
                                                        predictedLabels <- unlist(lapply(run,
                                                                                  function(fold)fold$predictedLabels))
                                                        trueLabels <- unlist(lapply(run,
                                                                             function(fold)fold$trueLabels))
                                                        return(sum(predictedLabels[trueLabels != caseClass] 
                                                                   != caseClass,na.rm=TRUE) / 
                                                               sum(trueLabels != caseClass))
                                                      })))
                                        },
                    objectiveFunctionParams = list(caseClass=caseClass),                                        
                    direction="maximize",
                    name="CV.Specificity")
}

