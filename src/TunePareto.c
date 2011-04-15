#include <R.h>
#include <Rinternals.h>
#include <stdbool.h>

SEXP calculateDominationMatrix(SEXP objectiveValues, SEXP minimizeObjectives)
{
  unsigned int numObjectives = length(minimizeObjectives);
  unsigned int numCombinations = length(objectiveValues)/numObjectives;
  
  double * _objectiveValues = REAL(objectiveValues);
  int * _minimizeObjectives = INTEGER(minimizeObjectives);
  
  SEXP res;
  PROTECT(res = allocMatrix(LGLSXP, numCombinations, numCombinations));
  int * _res = LOGICAL(res);

  unsigned int i, j, k;
  
  for (i = 0; i < numCombinations; ++i)
  {
       
    for (j = 0; j < numCombinations; ++j)
    {
      bool anyBetter = false;
      bool domination = true;
      for (k = 0; k < numObjectives; ++k)
      {
        if ((_minimizeObjectives[k] && 
            _objectiveValues[k*numCombinations+i] > _objectiveValues[k*numCombinations+j]) ||
            (!_minimizeObjectives[k] && 
            _objectiveValues[k*numCombinations+i] < _objectiveValues[k*numCombinations+j]))
        {
            anyBetter = true;
        }
        else
          domination = (_minimizeObjectives[k] && 
                        _objectiveValues[k*numCombinations+i] >= _objectiveValues[k*numCombinations+j]) ||
                       (!_minimizeObjectives[k] && 
                        _objectiveValues[k*numCombinations+i] <= _objectiveValues[k*numCombinations+j]);
        if (!domination)
          break;                        
      }
      _res[i*numCombinations+j] = domination && anyBetter;
    }
  }
  
  UNPROTECT(1);
  return res;
}
