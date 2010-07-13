#################################################################
# Main TunePareto routines
#################################################################

# Trains a classifier with the data supplied in <trainData> (one sample per row)
# and <trainLabels>, and returns a vector of predicted labels for <testData>.
# <classifier> is the classifier training function. 
# If <predictor> is empty, <classifier> is a combined training 
# and prediction function. 
# Otherwise, <predictor> specifies the prediction part of the classifier.
# <parameters> are a set of parameters to be passed to <classifier>.
# If <useFormula> is true, the classifier is supplied with a formula specifying the relation-
# ship between the data and the classes in a parameter called <formulaName>, and
# the class labels are automatically included in the training data frame.
# <trainDataName>, <trainLabelName>, <testDataName>, and <modelName> specify the names
# of the parameters of <classifier> and <predictor> for the training data, the training labels,
# the test data, and the classification model respectively.
callClassifier <- function(trainData, trainLabels, testData, classifier, classifierParams, 
                           predictor, predictorParams,
                           useFormula = FALSE, formulaName = "formula",
                           trainDataName = "x", trainLabelName = "y", testDataName = "newData",
                           modelName = "object")
{
  classifierParams <- lapply(classifierParams, unlist, recursive=FALSE)
  predictorParams <- lapply(predictorParams, unlist, recursive=FALSE)
  
  if (missing(predictor) || is.null(predictor))
  # empty predictor => combined train/predict method
  {
    # build parameter list
    
    if (useFormula)
    {
      trainData <- data.frame(trainLabels,trainData)
      colnames(trainData)[1] <- "Class"
      
      testData <- as.data.frame(testData)
      colnames(testData) <- colnames(trainData)[-1]
      
      paramList <- list(as.formula("Class ~ ."), trainData, testData)
      names(paramList) <- c(formulaName, trainDataName, testDataName)
    }
    else
    {
      paramList <- list(trainData, trainLabels, testData)
      names(paramList) <- c(trainDataName, trainLabelName, testDataName)
    }
    paramList <- c(paramList, as.list(classifierParams))
    
    # call classifier and return predicted labels
    return(do.call(classifier, paramList))
  }
  else
  # separate training and prediction methods
  {
    # build parameter list for training
   
    if (useFormula)
    {
      trainData <- data.frame(trainLabels,trainData)
      colnames(trainData)[1] <- "Class"
      
      testData <- as.data.frame(testData)
      colnames(testData) <- colnames(trainData)[-1]
      
      paramList <- list(as.formula("Class ~ ."), trainData)
      names(paramList) <- c(formulaName, trainDataName)
      
    }
    else
    {  
      paramList <- list(trainData, trainLabels)
      names(paramList) <- c(trainDataName, trainLabelName)
    }
    paramList <- c(paramList, classifierParams)
  
    # train the classifier
    model <- do.call(classifier, paramList)
    
    # build parameter list for testing
    paramList <- list(model, testData)
    names(paramList) <- c(modelName, testDataName)
    paramList <- c(paramList, predictorParams)
    
    # predict the unknown samples
    return(do.call(predictor, paramList))
  }
}

# Calculates a matrix specifying which configuration is dominated by which other combination
# based on a matrix of objective values <objectiveValues>.
# <minimizeObjectives> is a Boolean vector specifying which of the objectives
# are minimized.
#
# Returns a Boolean domination matrix.
calculateDominationMatrix <- function(objectiveValues, minimizeObjectives)
{
  domination <- matrix(apply(objectiveValues, 1, function(val1)
                # iterate over configurations
                {
                  apply(objectiveValues, 1, function(val2)
                  # iterate over configurations
                  {
                    oneBetter <- FALSE
                    dom <- sapply(1:length(minimizeObjectives),function(k)
                                  # iterate over objectives
                                 {
                                   if ((minimizeObjectives[k] && 
                                       val1[k] > val2[k]) ||
                                       (!minimizeObjectives[k] && 
                                       val1[k] < val2[k]))
                                   {
                                       oneBetter <<- TRUE
                                       return(TRUE)
                                   }
                                   return((minimizeObjectives[k] && 
                                       val1[k] >= val2[k]) ||
                                       (!minimizeObjectives[k] && 
                                       val1[k] <= val2[k]))
                                 })
                        return (all(dom) && oneBetter)                            
                  })
                }),ncol=nrow(objectiveValues))
  rownames(domination) <- rownames(objectiveValues)              
  colnames(domination) <- rownames(objectiveValues)
  return(domination)                
}

                          
# Tunes the parameters of a classifier for multiple objective functions.
# <data> is a dataset to use for tuning, where each row of the dataset
# corresponds to one sample. 
# <labels> are the corresponding class labels.
# <classifier> is the classifier training function. 
# If <predictor> is empty, <classifier> is a combined training 
# and prediction function. 
# Otherwise, <predictor> specifies the prediction part of the classifier.
# <classifierParameterRanges> and <predictorParameterRanges> are lists of parameter values to test. 
# Each element of the list must name a parameter of <classifier> or <predictor> respectively and consist 
# of a sub-list of possible values for this parameter.
# Alternatively, the complete lists of parameters can be specified in <classifierParameterCombinations>
# and <predictorParameterCombinations>.
# If <numCombinations> is specified, only a subsample of <numCombination> parameter
# combination is tested. 
# <objectiveFunctions> is a list of objective functions, i.e. objects of class TuneParetoObjective.
# If <keepSeed> is true, all parameter configurations are tested with the same random seed.
# If <useSnowfall> is true, parameter configurations are evaluated in parallel using snowfall.
# If <useFormula> is true, the classifier is supplied with a formula specifying the relation-
# ship between the data and the classes in a parameter called <formulaName>, and
# the class labels are automatically included in the training data frame.
# <trainDataName>, <trainLabelName>, <testDataName>, and <modelName> specify the names
# of the parameters of <classifier> and <predictor> for the training data, the training labels,
# the test data, and the classification model respectively.
# 
# Returns a set of Pareto-optimal parameter assignments and the corresponding
# objective function values. If <numCombinations> is specified, a complete list of tested
# combinations is returned as well.
tunePareto <- function(data, labels, 
                      classifier, classifierParameterRanges, classifierParameterCombinations,
                      predictor, predictorParameterRanges, predictorParameterCombinations,
                      numCombinations, objectiveFunctions, objectiveBoundaries,
                      keepSeed = TRUE, useSnowfall=FALSE, verbose=TRUE,
                      useFormula = FALSE, formulaName = "formula",                      
                      trainDataName = "x", trainLabelName = "y", 
                      testDataName = "newdata",modelName = "object")
{

  # standardize labels
  #classes <- sort(unique(labels))
  #newLabels <- rep(NA,length(labels))
  #for (i in 1:length(classes))
  #  newLabels[which(labels == classes[i])] <- i
  #labels <- as.factor(newLabels)
  
  labels <- as.factor(as.character(labels))
 
  if (!missing(objectiveBoundaries))
  {
    if (length(objectiveBoundaries) != length(objectiveFunctions))
        stop("Please supply exactly one boundary for each objective function!")
  }
 
  if (missing(classifierParameterRanges))
  {
    if (!missing(predictorParameterRanges) || missing (classifierParameterCombinations))
      stop(paste("Please provide either classifierParameterRanges/predictorParameterRanges or",
                 "classifierParameterCombinations/predictorParameterCombinations!"))
    
    if (missing(predictorParameterCombinations))
      combinations <- lapply(classifierParameterCombinations,function(x)list(classifierParams=x,
                                                                             predictorParams=NULL))
    else
      combinations <- unlist(lapply(classifierParameterCombinations, function(comb1)
                              {
                                lapply(predictorParameterCombinations,function(comb2)
                                {
                                  list(classifierParams=comb1, predictorParams=comb2)
                                })
                              }), recursive=FALSE)

    if (!missing(numCombinations) && !is.null(numCombinations))
      combinations <- combinations[sample(1:length(combinations), size=numCombinations, replace=FALSE)]
                              
    meaningfulParams <- NULL       
  }
  else
  {
    if (missing(predictorParameterRanges))
      predictorParameterRanges <- NULL
      
      
    # in the description of parameter configurations, use only those parameter that do not have a fixed value
    meaningfulParams <- unlist(c(sapply(classifierParameterRanges,function(param)length(param)>1),
                               sapply(predictorParameterRanges,function(param)length(param)>1)))  
    
    # determine the combinations to be tested
    if (missing(numCombinations) || is.null(numCombinations))
      combinations <- allCombinations(append(classifierParameterRanges, 
                                             predictorParameterRanges))
    else
      combinations <- sampleCombinations(append(classifierParameterRanges, 
                                                predictorParameterRanges), numCombinations)
  }
  
  # group objective functions that have the same precalculation routine
  # to save computational time
  groups <- groupByPrecalculation(objectiveFunctions)
  groupedObjectives <- groups$grouping

  if (verbose)
    cat("Testing parameter combinations...\n")
  
  # store random seed
  runif(n=1)  
  seed <- .Random.seed
  
  calculateObjectiveVals <- function(parameters)
    {
       if (verbose)
       {
         if (is.null(meaningfulParams))
         {
           allParams <- unlist(unname(parameters), recursive=FALSE)
           cat("Evaluating parameter set:",
                paste(paste(names(allParams),"=",allParams), collapse=", "),"\n")
         }
         else
         {
           cat("Evaluating parameter set:",
               paste(paste(names(parameters)[meaningfulParams],"=",parameters[meaningfulParams]), collapse=", "),"\n")
         }
       }
       
        # calculate objective function values for configuration
        res <- unlist(lapply(groupedObjectives,function(objective)
              {                 
                
                if (keepSeed)
                {
                  runif(n=1)
                  .Random.seed <<- seed
                }  
                # build parameter list for precalculation
                
                if (is.null(meaningfulParams))
                {
                  # If there are pre-specified parameter lists,
                  # the parameters for classifiers and predictors are
                  # supplied in two sub-lists
                  classifierParams <- parameters$classifierParams
                  predictorParams <- parameters$predictorParams
                }
                else
                {
                  # if all combinations were generated automatically, parameters for classifiers and
                  # predictors are stored in the same list
                  classifierParams <- parameters[1:length(classifierParameterRanges)]
                  predictorParams <- parameters[-(1:length(classifierParameterRanges))]
                }
                
                params <- list(data, labels, classifier, classifierParams, 
                               predictor, predictorParams,
                               useFormula, formulaName, trainDataName, trainLabelName, 
                               testDataName, modelName)
                names(params) <- c("data", "labels", "classifier", "classifierParams",
                                   "predictor", "predictorParams", 
                                   "useFormula", "formulaName",
                                   "trainDataName", "trainLabelName",
                                   "testDataName", "modelName")
                
                params <- c(params, as.list(objective$precalculationParams))
                
                # call precalculation
                precalc <- do.call(objective$precalculationFunction,params)             

                # call associated objective functions
                if (is.function(objective$objectiveFunction))
                {                                 
                  do.call(objective$objectiveFunction,
                          append(list(result=precalc), objective$objectiveFunctionParams))
		                
		            }
                else
                {
                    
                  mapply(function(func, params)
                         {
                            do.call(func, append(list(result=precalc),unlist(params,recursive=FALSE)))
                         }, objective$objectiveFunction, objective$objectiveFunctionParams)
                 }
                
              }))
         return(res)
      }
  
  if (useSnowfall)
  {
    require(snowfall)
    
     # export objects and functions needed in the cluster
    sfLibrary("snowfall", character.only=TRUE)
    sfLibrary("TunePareto", character.only=TRUE)
    sfExport("data","labels","classifier","predictor",
             "useFormula", "formulaName","trainDataName", "trainLabelName",
             "testDataName", "modelName","keepSeed","seed",
             "groupedObjectives","useSnowfall","verbose","meaningfulParams")
             
    # sfExport("callClassifier", namespace="TunePareto")
    
    # parallel evaluation of combinations in cluster
    objectiveValues <- t(sfSapply(combinations, calculateObjectiveVals))
  }
  else
  {
    # sequential evaluation of combinations
    objectiveValues <- t(sapply(combinations, calculateObjectiveVals))
  }
  
  if (verbose)
    cat("Calculating Pareto-optimal combinations...\n")
 
  objectiveValues <- matrix(objectiveValues, nrow=length(combinations))[, groups$permutation, drop=FALSE]
                         
  colnames(objectiveValues) <- sapply(objectiveFunctions,function(obj)obj$name)                        
  
  rownames(objectiveValues) <- sapply(combinations, 
                                      function(comb)
                                      {
                                          if (!is.null(meaningfulParams))
                                            comb <- comb[meaningfulParams]
                                          else
                                            comb <- unlist(unname(comb), recursive=FALSE)
                                          paste(paste(names(comb),"=",comb), collapse=", ")
                                      })

  # calculate a matrix specifying which configuration is dominated by which other combination
  minimizeObjectives <- sapply(objectiveFunctions,function(x)x$minimize)
  names(minimizeObjectives) <- colnames(objectiveValues)
  domination <- calculateDominationMatrix(objectiveValues, minimizeObjectives)
  
  # determine dominated/non-dominated solutions              
  dominated <- apply(domination,2,any)
  
  if (!missing(objectiveBoundaries))
  {
    dominated <- dominated | 
                 apply(objectiveValues,1,function(val)
                 {
                  any(mapply(function(v, bound, minimize)
                      {
                        if (minimize)
                          v > bound
                        else
                          v < bound
                      }, val, objectiveBoundaries, minimizeObjectives))
                 })
  }

  # build result list                                      
  res <- list(bestCombinations=combinations[!dominated], 
              bestObjectiveValues=objectiveValues[!dominated,,drop=FALSE],
              testedCombinations=combinations,
              testedObjectiveValues=objectiveValues,
              minimizeObjectives = minimizeObjectives,
              objectiveBoundaries =  if (missing(objectiveBoundaries)){NULL}
                                    else {objectiveBoundaries},
              dominationMatrix=domination)
  class(res) <- "TuneParetoResult"              
  return(res)  
}

# Print function for objects of class TuneParetoResult.
# Prints the non-dominated solutions of <x>
print.TuneParetoResult <- function(x, ...)
{ 
  if (is.null(x$objectiveBoundaries))  
    cat("Pareto-optimal parameter sets:\n")
  else
    cat("Pareto-optimal parameter sets matching the objective restrictions:\n")
  print(x$bestObjectiveValues)
  return(invisible(x))
}

# Wrapper for tunePareto that tunes a tree.
# For parameters, see tunePareto and tree.
tunePareto.tree <- function(data, labels,
                            weights, subset,
                            na.action,
                            method,
                            split,
                            mincut, minsize, mindev,
                            numCombinations,                 
                            objectiveFunctions,
                            objectiveBoundaries,
                            keepSeed = TRUE,
                            useSnowfall = FALSE,
                            verbose=TRUE)
{
  if (useSnowfall)
  {
    require(snowfall)  
    sfLibrary("tree", character.only=TRUE)
  }
  else
    require(tree)
    
  paramRanges = list()
  if (!missing(weights))
    paramRanges$weights <- weights
    
  if (!missing(na.action))
    paramRanges$na.action <- na.action
    
  if (!missing(method))
    paramRanges$method <- method
    
  if (!missing(split))
    paramRanges$split <- split

  if (!missing(mincut))
    paramRanges$mincut <- mincut
    
  if (!missing(minsize))
    paramRanges$minsize <- minsize
    
  if (!missing(mindev))
    paramRanges$mindev <- mindev
    
  if (missing(numCombinations))
    numCombinations <- NULL
    
  tunePareto(data = data, 
             labels = labels,
             classifier = tree, 
             predictor = predict, 
             classifierParameterRanges = paramRanges,
             predictorParameterRanges = list(type="class"),
             numCombinations = numCombinations,
             objectiveFunctions = objectiveFunctions,
             objectiveBoundaries = objectiveBoundaries,
             keepSeed = keepSeed,
             useSnowfall = useSnowfall,
             verbose = verbose,             
             useFormula=TRUE,
             formulaName="formula",
             trainDataName="data",
             testDataName="newdata",
             modelName="object")        
}

# Wrapper for tunePareto that tunes a random forest.
# For parameters, see tunePareto and randomForest.
tunePareto.randomForest <- function(data, labels,
                                    subset,
                                    na.action,
                                    ntree,
                                    mtry,
                                    replace, 
                                    classwt, 
                                    cutoff, 
                                    strata,
                                    sampsize,
                                    nodesize,
                                    maxnodes,
                                    numCombinations,                 
                                    objectiveFunctions,
                                    objectiveBoundaries,                                    
                                    keepSeed = TRUE,
                                    useSnowfall = FALSE,
                                    verbose=TRUE)
{
  if (useSnowfall)
  {
    require(snowfall)  
    sfLibrary("randomForest", character.only=TRUE)
  }
  else
    require(randomForest)
    
  paramRanges = list()
  
  if (!missing(subset))
    paramRanges$subset <- subset
    
  if (!missing(na.action))
    paramRanges$subset <- na.action
  
  if (!missing(mtry))
    paramRanges$mtry <- mtry
    
  if (!missing(ntree))
    paramRanges$ntree <- ntree
    
  if (!missing(replace))
    paramRanges$replace <- replace
    
  if (!missing(classwt))
    paramRanges$classwt <- classwt
    
  if (!missing(cutoff))
    paramRanges$cutoff <- cutoff

  if (!missing(strata))
    paramRanges$strata <- strata
    
  if (!missing(sampsize))
    paramRanges$sampsize <- sampsize
    
  if (!missing(nodesize))
    paramRanges$nodesize <- nodesize
    
  if (!missing(maxnodes))
    paramRanges$maxnodes <- maxnodes
    
  if (missing(numCombinations))
    numCombinations <- NULL
    
  tunePareto(data = data, 
             labels = labels,
             classifier = randomForest, 
             predictor = predict, 
             classifierParameterRanges = paramRanges,
             numCombinations = numCombinations,
             objectiveFunctions = objectiveFunctions,
             objectiveBoundaries = objectiveBoundaries,             
             keepSeed = keepSeed,
             useSnowfall = useSnowfall,
             verbose = verbose,             
             useFormula=TRUE,
             trainDataName="data",
             formulaName="formula",
             testDataName="newdata",
             modelName="object")        
}


# Wrapper for tunePareto that tunes a support vector machine.
# For parameters, see tunePareto and svm
tunePareto.svm <- function(data, labels,
                           kernel, degree, gamma,
                           coef0, cost, nu,
                           class.weights, cachesize, 
                           tolerance, epsilon,
                           subset, na.action,
                           numCombinations,                            
                           objectiveFunctions,
                           objectiveBoundaries,                           
                           keepSeed = TRUE,
                           useSnowfall = FALSE,
                           verbose=TRUE)
{
  if (useSnowfall)
  {
    require(snowfall)
    sfLibrary("e1071", character.only=TRUE)
  }
  else
    require(e1071)
    
  paramRanges = list()
  if (!missing(kernel))
    paramRanges$kernel <- kernel
    
  if (!missing(degree))
    paramRanges$degree <- degree
    
  if (!missing(gamma))
    paramRanges$gamma <- gamma
    
  if (!missing(coef0))
    paramRanges$coef0 <- coef0

  if (!missing(cost))
    paramRanges$cost <- cost
    
  if (!missing(nu))
    paramRanges$nu <- nu
    
  if (!missing(class.weights))
    paramRanges$class.weights <- class.weights
    
  if (!missing(cachesize))
    paramRanges$cachesize <- cachesize
    
  if (!missing(tolerance))
    paramRanges$tolerance <- tolerance
    
  if (!missing(epsilon))
    paramRanges$epsilon <- epsilon

  if (!missing(subset))
    paramRanges$subset <- subset
    
  if (!missing(na.action))
    paramRanges$na.action <- na.action
    
  if (missing(numCombinations))
    numCombinations <- NULL    
        
  tunePareto(data = data,   
             labels = labels,
             classifier = svm, 
             predictor = predict, 
             classifierParameterRanges = paramRanges,
             numCombinations = numCombinations,
             objectiveFunctions = objectiveFunctions,
             objectiveBoundaries = objectiveBoundaries,            
             keepSeed = keepSeed,
             useSnowfall = useSnowfall,
             verbose = verbose,
             useFormula = FALSE,
             trainDataName = "x",
             trainLabelName = "y",
             testDataName = "newdata",
             modelName = "object")
}

# Wrapper for tunePareto that tunes a k-Nearest Neighbour classifier.
# For parameters, see tunePareto and knn.
tunePareto.knn <- function(data, labels,
                           k, l, use.all,
                           numCombinations,                            
                           objectiveFunctions,
                           objectiveBoundaries,                           
                           keepSeed = TRUE,
                           useSnowfall = FALSE,
                           verbose=TRUE)
{
  if (useSnowfall)
  {
    require(snowfall)
    sfLibrary("class", character.only=TRUE)
  }
  else
    require(class)
    
  paramRanges = list()
  if (!missing(k))
    paramRanges$k <- k
    
  if (!missing(l))
    paramRanges$l <- l
    
  if (!missing(use.all))
    paramRanges$use.all <- use.all
    
  if (missing(numCombinations))
    numCombinations <- NULL    
        
  tunePareto(data = data, 
             labels = labels,
             classifier = knn, 
             predictor = NULL, 
             classifierParameterRanges = paramRanges,
             numCombinations = numCombinations,
             objectiveFunctions = objectiveFunctions,
             objectiveBoundaries = objectiveBoundaries,
             keepSeed = keepSeed,
             useSnowfall = useSnowfall,
             verbose = verbose,             
             useFormula = FALSE,
             trainDataName = "train",
             testDataName = "test",
             trainLabelName = "cl")
}                           
