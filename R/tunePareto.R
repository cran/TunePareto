#################################################################
# Main TunePareto routines
#################################################################

# Calculates a matrix specifying which configuration is dominated by which other combination
# based on a matrix of objective values <objectiveValues>.
# <minimizeObjectives> is a Boolean vector specifying which of the objectives
# are minimized.
# Returns a Boolean domination matrix.
calculateDominationMatrix <- function(objectiveValues, minimizeObjectives)
{
  domination <- .Call("calculateDominationMatrix",as.numeric(as.matrix(objectiveValues)), as.integer(minimizeObjectives))
  rownames(domination) <- rownames(objectiveValues)              
  colnames(domination) <- rownames(objectiveValues)
  return(domination)  
}

                          
# Tunes the parameters of a classifier for multiple objective functions.
# <data> is a dataset to use for tuning, where each row of the dataset
# corresponds to one sample. 
# <labels> are the corresponding class labels.
# <classifier> is the classifier wrapper object. 
# Parameters to be tuned can be supplied in ...
# Alternatively, the complete lists of parameters can be specified in <parameterCombinations>
# If <sampleType> is not "full" or "evolution", only a subsample of <numCombination> parameter
# combinations is drawn using the specified sampling technique. 
# For <sampleType="evolution">, the number of individuals <mu> and offspring <lambda> and the number
# of generations <numIterations> can be specified.
# <objectiveFunctions> is a list of objective functions, i.e. objects of class TuneParetoObjective.
# <objectiveBoundaries> is a vector of lower/upper bounds for the objectives.
# If <keepSeed> is true, all parameter configurations are tested with the same random seed.
# If <useSnowfall> is true, parameter configurations are evaluated in parallel using snowfall.
# If <verbose> is true, status information is printed
# 
# Returns a set of Pareto-optimal parameter assignments and the corresponding
# objective function values. 
tunePareto <- function(..., data, labels, 
                       classifier, parameterCombinations,
                       sampleType=c("full","uniform","latin","halton","niederreiter","sobol","evolution"), 
                       numCombinations, 
                       mu=10, lambda=20, numIterations=100,
                       objectiveFunctions, objectiveBoundaries,
                       keepSeed = TRUE, useSnowfall=FALSE, verbose=TRUE)
{
  sampleType <- match.arg(sampleType, c("full","uniform","latin","halton","niederreiter","sobol","evolution"))
  
  if ((sampleType %in% c("uniform","latin","halton","niederreiter","sobol")) &&
      (missing(numCombinations) || is.null(numCombinations)))
    stop("For sampleType=\"uniform\",\"latin\",\"halton\",\"niederreiter\" or \"sobol\", you must specify numCombinations!")

  if(!inherits(classifier, "TuneParetoClassifier"))
    stop("\"classifier\" must be a TuneParetoClassifier object!")
  
  if (length(labels) != nrow(data))
    stop("Dimensions of data matrix and class label vector are incompatible!")
  
  if (useSnowfall)
  {
     require(snowfall)
    
    # export libraries needed in the cluster
    sfLibrary("snowfall", character.only=TRUE)
    sfLibrary("TunePareto", character.only=TRUE)
    
    if (length(classifier$requiredPackages) > 0)
      lapply(classifier$requiredPackages,function(package)sfLibrary(package, character.only=TRUE))
  }
  
  labels <- as.factor(as.character(labels))
 
  if (!missing(objectiveBoundaries))
  {
    if (length(objectiveBoundaries) != length(objectiveFunctions))
        stop("Please supply exactly one boundary for each objective function!")
  }


  args <- list(...)

  if (length(args) == 0)
  {
    if (missing (parameterCombinations))
      stop(paste("Please provide either parameterCombinations or",
                 "parameter values in the ... argument!"))
    
    if (sampleType != "full")
      stop(paste("parameterCombinations can only be used with sampleType=\"full\"!"))
    
    
    # determine which parameters belong to the classifier and predictor respectively
    combinations <- parameterCombinations
    classifierParamPos <- intersect(names(parameterCombinations[[1]]),classifier$classifierParamNames)
    predictorParamPos <- intersect(names(parameterCombinations[[1]]),classifier$predictorParamNames)
                             
    meaningfulParams <- NULL
  }
  else
  {
    nonmatch <- setdiff(names(args), c(classifier$classifierParamNames,classifier$predictorParamNames))
    
    if (length(nonmatch) > 0)
      stop("The following unknown parameters have been specified: ",nonmatch)

    if (sampleType == "full" && any(sapply(args,is.interval)))
      stop("For sampleType=\"full\", no intervals can be specified!")

    # determine which parameters belong to the classifier and predictor respectively
    classifierParamPos <- intersect(names(args),classifier$classifierParamNames)
    predictorParamPos <- intersect(names(args),classifier$predictorParamNames)
      
    # in the description of parameter configurations, use only those parameter that do not have a fixed value
    meaningfulParams <- sapply(args,function(param)length(param)>1)
    
    # determine the combinations to be tested
    # using the corresponding sampling technique
    combinations <- 
        switch(sampleType,
               full = allCombinations(args),
               latin = latinHypercube(args, numCombinations),
               evolution = NULL,
               sampleCombinations(args, numCombinations, method=sampleType))
  }
  
  # group objective functions that have the same precalculation routine
  # to save computational time
  groups <- groupByPrecalculation(objectiveFunctions)
  groupedObjectives <- groups$grouping
  
  minimizeObjectives <- sapply(objectiveFunctions,function(x)x$minimize)

  if (verbose)
    cat("Testing parameter combinations...\n")
  
  # store random seed
  runif(n=1)  
  seed <- .Random.seed
  
  # internal function that calculates the objective values
  # for a set of parameter values
  calculateObjectiveVals <- function(parameters)
    {
       if (verbose)
       {
         if (is.null(meaningfulParams))
         {
           cat("Evaluating parameter set:",
                paste(paste(names(parameters),"=",parameters), collapse=", "),"\n")
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
                
                classifierParams <- parameters[classifierParamPos]
                predictorParams <- parameters[predictorParamPos]
                
                params <- list(data, labels, classifier, classifierParams, predictorParams)
                names(params) <- c("data", "labels", "classifier", "classifierParams", "predictorParams")
                
                params <- c(params, as.list(objective$precalculationParams))
                
                tryCatch(
                {
                   # call precalculation
                  precalc <- do.call(objective$precalculationFunction,params)             

                  # call associated objective functions
                  if (is.function(objective$objectiveFunction))
                  {                                 
                    return(do.call(objective$objectiveFunction,
                            append(list(result=precalc), objective$objectiveFunctionParams)))
		                  
		              }
                  else
                  {
                      
                    return(mapply(function(func, params)
                           {
                              do.call(func, append(list(result=precalc),unlist(params,recursive=FALSE)))
                           }, objective$objectiveFunction, objective$objectiveFunctionParams))
                   }
                 },
                 error = function(e)
                 {
                   # an error occurred for the parameter combination
                   allParams <- unlist(parameters, recursive=FALSE)
                   warning(paste("Combination",
                           paste(paste(names(allParams),"=",allParams), collapse=", "),"returned an error:",e))
                   if (is.function(objective$objectiveFunction))
                     NA
                   else
                     rep(NA,length(objective$objectiveFunction))
                 })
              }))
         return(res)
      }
  
  if (sampleType == "evolution")
  # switch to Evolution Strategies
  {
    # initialize individuals and mutation rates
    mutationRates <- sapply(args,function(range)
                            {
                              if (is.interval(range))
                              {
                                (range$upper - range$lower)/100
                              }
                              else
                                NA
                            })
                            
    individuals <- lapply(latinHypercube(args,N=mu), function(ind)
                            list(individual=ind,mutation=mutationRates))
    
    # calculate initial fitness
    if (useSnowfall)
    { 
      # export objects and functions needed in the cluster   
      sfExport("data","labels","classifier","keepSeed","seed",
               "groupedObjectives","useSnowfall","verbose","meaningfulParams")
               
      # parallel evaluation of combinations in cluster
      oldObjectiveValues <- t(sfSapply(individuals, function(cand)calculateObjectiveVals(cand$individual)))
    }
    else
    {
      # sequential evaluation of combinations
      oldObjectiveValues <- t(sapply(individuals, function(cand)calculateObjectiveVals(cand$individual)))
    }
    
    for (i in 1:numIterations)
    # iterate for <numIterations> generations
    {
      if (verbose)
        cat("Iteration ",i,"\n",sep="")
      
      # calculate crowding distances      
      crowd <- apply(oldObjectiveValues,2,function(dim)
         {
          idx <- order(dim)
          dim <- dim[idx]
          
          res <- sapply(1:length(dim),function(i)
                 {
                  if (i == 1)
                    dim[2]-dim[1]
                  else
                  if (i == length(dim))
                    dim[length(dim)]-dim[length(dim) - 1]
                  else
                    mean(c(dim[i]-dim[i-1],dim[i+1]-dim[i]))
                 })
          revidx <- order(idx)
          return(res[revidx])
        })
      crowd <- apply(crowd,1,mean)
      
      # mating probability is based on the rank
      # of the crowding distance
      crowdProbs <- rank(crowd)
      crowdProbs <- crowdProbs/sum(crowdProbs)      
      
      # create offspring  
      candidates <- lapply(1:lambda, function(i)
                    {
                      parents <- sample(1:length(individuals),size=2,replace=FALSE,prob=crowdProbs)
                      mutate(recombine(individuals[[parents[1]]],individuals[[parents[2]]]),
                             args)
                    })
      
      # calculate fitness of offspring
      if (useSnowfall)
      { 
        # export objects and functions needed in the cluster   
        sfExport("data","labels","classifier","keepSeed","seed",
                 "groupedObjectives","useSnowfall","verbose","meaningfulParams")
                 
        # parallel evaluation of combinations in cluster
        objectiveValues <- rbind(oldObjectiveValues,
                                 t(sfSapply(candidates, function(cand)calculateObjectiveVals(cand$individual))))
      }
      else
      {
        # sequential evaluation of combinations
        objectiveValues <- rbind(oldObjectiveValues,
                                 t(sapply(candidates, function(cand)calculateObjectiveVals(cand$individual))))
      }
      
      candidates <- c(individuals, candidates)
      
      # non-dominated sorting on candidates
      dominationMatrix <- calculateDominationMatrix(objectiveValues[, groups$permutation], minimizeObjectives)
      fronts <- calculateParetoFronts(dominationMatrix)
      
      remaining <- mu
      indices <- c()
      for (front in fronts)
      {
        if (remaining < length(front))
        {
          indices <- c(indices, sample(front,size=remaining,replace=FALSE))
          break
        }
        else
        {
          indices <- c(indices, front)
          remaining <- remaining - length(front)
        }
      }
      
      individuals <- candidates[indices]
      oldObjectiveValues <- objectiveValues[indices,,drop=FALSE]
    }
    
    # remove duplicate parameter configurations
    combinations <- lapply(candidates,function(ind)ind$individual)
    dup <- duplicated(combinations)
    combinations <- combinations[!dup]
    
  }
  else
  {
    if (useSnowfall)
    { 
      # export objects and functions needed in the cluster   
      sfExport("data","labels","classifier","keepSeed","seed",
               "groupedObjectives","useSnowfall","verbose","meaningfulParams")
               
      # parallel evaluation of combinations in cluster
      objectiveValues <- t(sfSapply(combinations, calculateObjectiveVals))
    }
    else
    {
      # sequential evaluation of combinations
      objectiveValues <- t(sapply(combinations, calculateObjectiveVals))
    }
  }
  
  # calculate optimal combinations
  if (verbose)
    cat("Calculating Pareto-optimal combinations...\n")
 
  objectiveValues <- matrix(objectiveValues, nrow=length(combinations))[, groups$permutation, drop=FALSE]
  
  invalidEntries <- apply(objectiveValues,1,function(row)any(is.na(row)))
  objectiveValues <- objectiveValues[!invalidEntries,,drop=FALSE]
  combinations <- combinations[!invalidEntries]
                         
  colnames(objectiveValues) <- sapply(objectiveFunctions,function(obj)obj$name)                        
    
  rownames(objectiveValues) <- sapply(combinations, 
                                      function(comb)
                                      {
                                          if (!is.null(meaningfulParams))
                                            comb <- comb[meaningfulParams]
                                          
                                          comb <- sapply(comb,function(x)
                                                        {
                                                          if  (is.numeric(x) && round(x) != x)
                                                            sprintf("%.5g",x)
                                                          else
                                                            as.character(x)
                                                        })  
                                          paste(paste(names(comb),"=",comb), collapse=", ")
                                      })
                                      
  names(minimizeObjectives) <- colnames(objectiveValues)                                      

  # calculate a matrix specifying which configuration is dominated by which other combination
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
                        if (is.na(bound) || is.null(bound))
                          FALSE
                        else
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

# Recombination method for two parents <parent1> and <parent2>
# in the Evolution Strategies
recombine <- function(parent1, parent2)
{
  # recombine parameter values
  child <- mapply(function(gene1, gene2, mut)
                  {
                     
                     if (!is.na(mut))
                     # this is a continuous gene
                     {
                       mean(gene1,gene2)
                     }
                     else
                     # this is a discrete gene
                     {
                       sample(c(gene1,gene2),size=1)
                     }
                  }, parent1$individual, parent2$individual, parent1$mutation,
                  SIMPLIFY=FALSE)
                  
  # recombine mutation rates
  mutation <- mapply(function(mut1, mut2)
                     {
                       if (!is.na(mut1))
                         mean(c(mut1,mut2))
                       else
                         NA
                     },
                     parent1$mutation,parent2$mutation)
  return(list(individual=child,mutation=mutation))
}

# Mutation method for an individual in
# the Evolution Strategies.
# <ranges> specifies the possible values
# for the parameters.
mutate <- function(individual, ranges)
{
  # mutate the mutation rates
  individual$mutation <- sapply(individual$mutation, function(mut)
                               {
                                 if (is.na(mut))
                                   NA
                                 else
                                  mut * exp(rnorm(mean=0,sd=1/sqrt(2*length(individual$mutation)),n=1) +
                                            rnorm(mean=0,sd=1/sqrt(2*sqrt(length(individual$mutation))),n=1))
                               })
  
  # mutate the parameter values
  individual$individual <- mapply(function(gene, mut, range)
                                  {
                                     if (!is.na(mut))
                                     # this is a continuous gene
                                     {
                                       max(range$lower,
                                           min(range$upper,gene + rnorm(mean=0,sd=mut,n=1)))
                                     }
                                     else
                                     if (runif(n=1) < 1/length(ranges))
                                     # perform a discrete mutation only with a small probability
                                     {
                                       if (!is.numeric(range))
                                       # nominally scaled attribute
                                       {
                                         sample(range,size=1)
                                       }
                                       else
                                       # integer attribute => choose from neighbourhood
                                       {
                                         idx <- which(range==gene)
                                         if (idx == 1)
                                           range[2]
                                         else
                                         if (idx == length(range))
                                           range[idx-1]
                                         else
                                          sample(c(range[idx-1],range[idx+1]),size=1)
                                       }
                                     }
                                     else
                                       gene
                                  },
                                  individual$individual, 
                                  individual$mutation,
                                  ranges,
                                  SIMPLIFY=FALSE)
  return(individual)
}

# Recalculate the Pareto-optimal solutions of <tuneParetoResult>
# using only the objectives specified in the index vector <objectives>.
recalculateParetoSet <- function(tuneParetoResult, objectives)
{
  if(!inherits(tuneParetoResult, "TuneParetoResult"))
    stop("\"tuneParetoResult\" must be a TuneParetoResult object!")
  
  if (missing(objectives))
    objectives <- 1:length(tuneParetoResult$minimizeObjectives)
  
  combinations <- tuneParetoResult$testedCombinations
  objectiveValues <- tuneParetoResult$testedObjectiveValues[,objectives,drop=FALSE]
  minimizeObjectives <- tuneParetoResult$minimizeObjectives[objectives]
  objectiveBoundaries <- tuneParetoResult$objectiveBoundaries[objectives]
  
  dominationMatrix <- calculateDominationMatrix(objectiveValues, minimizeObjectives)
  
  # determine dominated/non-dominated solutions              
  dominated <- apply(dominationMatrix,2,any)
  
  if (!is.null(objectiveBoundaries))
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
              objectiveBoundaries =  objectiveBoundaries,
              dominationMatrix=dominationMatrix)
  class(res) <- "TuneParetoResult"  
  return(res)
}

# Merges a list of TuneParetoResult objects
# and calculates the common Pareto-optimal solutions.
mergeTuneParetoResults <- function(...)
{
  l <- list(...)
  
  if (any(sapply(l,function(r)!inherits(r, "TuneParetoResult"))))
    stop("All supplied arguments must be a TuneParetoResult object!")
    
  objectives <- lapply(l,function(r)colnames(r$testedCombinations))
  if (length(unique(objectives)) != 1)
    stop("All supplied TuneParetoResult objects must use the same objective functions!")
    
  res <- l[[1]]
  for (el in l[2:length(l)])
  {
    res$testedCombinations <- c(res$testedCombinations, el$testedCombinations)
    res$testedObjectiveValues <- rbind(res$testedObjectiveValues, el$testedObjectiveValues)
  }
  
  return(recalculateParetoSet(res))  
}

# Print function for objects of class TuneParetoResult.
# Prints the non-dominated solutions of <x>.
print.TuneParetoResult <- function(x, ...)
{ 
  if(!inherits(x, "TuneParetoResult"))
    stop("\"x\" must be a TuneParetoResult object!")
    
  if (is.null(x$objectiveBoundaries))  
    cat("Pareto-optimal parameter sets:\n")
  else
    cat("Pareto-optimal parameter sets matching the objective restrictions:\n")
  print(x$bestObjectiveValues)
  return(invisible(x))
}
