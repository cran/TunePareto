#################################################################
# Helper functions
#################################################################


# Iteratively calculate all combinations of variable assignments,
# where <numValues> is a vector specifying the number of choices
# for each of the variables.
# <currentCombination> is the last used combination from which
# the next new combination is calculated.
#
# Returns the next combination or NULL if all combinations have
# been tested.
nextCombination <- function(numValues, currentCombination)
{	
	while(TRUE)
	{
		for(i in (length(currentCombination)):1)
		{
			currentCombination[i] <- currentCombination[i] + 1	
			if(currentCombination[i] <= numValues[i])
			{
				return(currentCombination)	
			}
			else
			{
				if (i==1)
					return(NULL)	
				else
					currentCombination[i] <- 1
			}		
		}		
	}
}

# Retrieves a data frame of possible parameter combinations
# based on the above method.
# <parameterRanges> is a list of lists, where each of the
# sub-lists specifies the possible values a parameter can take.
allCombinations <- function(parameterRanges)
{
  combination <- rep(1, length(parameterRanges))
  
  # build a data frame from the first combination
  res <- list(mapply(function(parameters, index)
                             {
                              parameters[[index]]
                             },
                             parameterRanges, combination,
                             SIMPLIFY=FALSE))

  numValues <- sapply(parameterRanges,length)

  j <- 1
  while (TRUE)
  {
    j <- j + 1
    combination <- nextCombination(numValues, combination)
    if (is.null(combination))
    # no more combinations
      break
      
    # append new combination to data frame
    res[[j]] <- mapply(function(parameters, index)
                             {
                              parameters[[index]]
                             },
                             parameterRanges, combination,
                             SIMPLIFY=FALSE)
   }
   return(res)
}

# Draws a sample of <N> random parameter combinations from the
# parameter values specified in <parameterRanges>.
#
# Returns the combination list
sampleCombinations <- function(parameterRanges, N)
{
  comb <- allCombinations(parameterRanges)
  return(comb[sample(1:length(comb), size=N, replace=FALSE)])

  #samples <- lapply(parameterRanges, function(param)
  #                  {
  #                    sample(param, size=N, replace=TRUE)
  #                  })
                    
  #res <- lapply(1:N,function(i)
  #              {
  #                r <- lapply(1:length(parameterRanges), function(j)
  #                {
  #                  samples[[j]][[i]]
  #                })
  #                names(r) <- names(parameterRanges)
  #                r
  #              })
  #return(res)                                 
}
