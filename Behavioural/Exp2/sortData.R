# Data sorting function
sortData <- function(inData,removeNA=TRUE){
  # Remove NA given args
  if (removeNA == TRUE) { 
    inData <- inData[ which(inData$Response!="NaN"), ]
  }
  # Get data length
  numTrials = dim(inData)[1] 
  # Initialize new vectors
  contextNo = vector(mode="integer",length=numTrials)
  gainVal = vector(mode="integer",length=numTrials)
  lossVal = vector(mode="integer",length=numTrials)
   
  # Iterate through length of data
  for (t in 1:numTrials) {
    # Create new column for context number
    currContextNo = integer()
    if (t == 1) { 
      currContextNo = 1
    } else {
      if (inData$subj_idx[t] !=  inData$subj_idx[t-1]) { 
        currContextNo = 1 
      } else { 
        if  (inData$context[t] !=  inData$context[t-1]) { 
          currContextNo = contextNo[t-1] + 1
        } else { 
          currContextNo = contextNo[t-1]
        }
      }
    }
    contextNo[t] = currContextNo
    # Create gainVal column
    currGainVal = numeric()
    currLossVal = numeric()
    if (inData$context[t] == "OO") { 
      currGainVal = 1.00
      currLossVal = 1.00
    } else if (inData$context[t] == "OF") { 
      currGainVal = 1.00
      currLossVal = 0.50
    } else if (inData$context[t] == "FO") { 
      currGainVal = 0.50
      currLossVal = 1.00
    } else if (inData$context[t] == "FF") {     
      currGainVal = 0.50
      currLossVal = 0.50      
    }
    gainVal[t] = currGainVal
    lossVal[t] = currLossVal
  }
  # Append new and revised vectors to dataframe
  outData = inData
  outData["contextNo"] <- contextNo
  outData["gainVal"] <- gainVal
  outData["lossVal"] <- lossVal
  # Rename outData columns
  colnames(outData)[colnames(outData)=="rt_sec"] <- "rt_sec"
  return(outData)
}


