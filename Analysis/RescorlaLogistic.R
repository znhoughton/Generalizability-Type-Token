RescorlaLogistic = function (cuesOutcomes, traceCue = "h", traceOutcome = "hand", 
    nruns = 1, random = TRUE, randomOrder = NA, alpha = 1, alpha2 = .1,
    lambda = 1, lambda0=0, beta1 = 1, beta2 = 0.1, sigSlope = 1,weightvec=0, ap=FALSE) #Change lambda0 to -1 if using tanh
{
    cues = unique(unlist(strsplit(as.character(cuesOutcomes$Cues), 
        "_")))
    outcomes = unique(unlist(strsplit(as.character(cuesOutcomes$Outcomes), 
        "_")))
    knownCues=NULL
    if (ap==TRUE) {
      ccw = matrix(0, length(cues), length(cues))
      rownames(ccw) <- cues
      colnames(ccw) <- cues
    }
	if (length(weightvec)==1)
	{   
		weightvec = rep(0, length(cues)) 
	}
    names(weightvec) = cues
    res = vector(mode = "numeric")
    dfr = data.frame(cuevector = rep(cuesOutcomes$Cues, cuesOutcomes$Frequency), 
        outcomevector = rep(cuesOutcomes$Outcomes, cuesOutcomes$Frequency), 
        stringsAsFactors = FALSE)
    theOrder = NA
    cnt = 0
    for (run in 1:nruns) {
        if (random) {
            if (is.na(randomOrder[1])) {
                theOrder = sample(1:nrow(dfr))
            }
            else {
                theOrder = randomOrder
            }
            dfr = dfr[theOrder, ]
        }
 
      
        #For each trial
        for (i in 1:dim(dfr)[1]) {
            currentCues = unlist(strsplit(dfr$cuevector[i], "_"))
            currentOutcomes = unlist(strsplit(dfr$outcomevector[i], 
                "_"))
            
            if (ap==T) {
              absentCues=NULL
              for (j in 1:length(currentCues))
              {
                absentCues<-c(absentCues,names(which(ccw[currentCues[j],]==1)))
              }
              absentCues<-absentCues[!(absentCues %in% currentCues)]
              #print(paste("cues:",paste(currentCues,collapse="_"),"outcomes:",paste(currentOutcomes),"absent:",paste(absentCues,collapse="_")))
              
               #keeping track of which cues co-occurred
              if (length(currentCues)>1)
              {
               ccw[currentCues,currentCues]<-1
              }
            }
              
            #Vtotal = sum(weightvec[currentCues])
            #Vtotal = 2*invlogit(sum(weightvec[currentCues]))-1 #this is the tanh function to scale between -1 and 1
            #tanh is needed for second-order retrospective revaluation, to pass on inhibition
            Vtotal = invlogit(sum(weightvec[currentCues]))

            if (traceOutcome %in% currentOutcomes) {
              weightvec[currentCues] = weightvec[currentCues] + 
                alpha * beta1 * (lambda - Vtotal)
              if (ap==T) {
                weightvec[absentCues] = weightvec[absentCues] + 
                  alpha2 * beta1 * (lambda0 - Vtotal)
                #weightvec[absentCues] = weightvec[absentCues] - 
                #alpha2 * beta1 * (lambda - Vtotal)
                
              }
            }
            if (!(traceOutcome %in% currentOutcomes))
            {
              weightvec[currentCues] = weightvec[currentCues] + 
                alpha * beta2 * (lambda0 - Vtotal)
            }
            
            cnt = cnt + 1
            res[cnt] = weightvec[traceCue]
        }
    }
    equilibriumWeights = estimateWeights(cuesOutcomes)
    result <- (list(weightvector = res, 
        traceCue = traceCue, traceOutcome = traceOutcome, randomOrder = theOrder, eqw=0))
    class(result) <- "RescorlaWagner"
    return(result)
}


RescorlaHumble = function (cuesOutcomes, traceCue = "h", traceOutcome = "hand", 
                             nruns = 1, random = TRUE, randomOrder = NA, alpha = 0.1, 
                             lambda = 1, beta1 = 0.1, beta2 = 0.1, sigSlope = 1,weightvec=0,Humble=T) 
{
  cues = unique(unlist(strsplit(as.character(cuesOutcomes$Cues), 
                                "_")))
  outcomes = unique(unlist(strsplit(as.character(cuesOutcomes$Outcomes), 
                                    "_")))
  if (length(weightvec)==1)
  {   
    weightvec = rep(0, length(cues)) 
  }
  names(weightvec) = cues
  res = vector(mode = "numeric")
  dfr = data.frame(cuevector = rep(cuesOutcomes$Cues, cuesOutcomes$Frequency), 
                   outcomevector = rep(cuesOutcomes$Outcomes, cuesOutcomes$Frequency), 
                   stringsAsFactors = FALSE)
  theOrder = NA
  cnt = 0
  for (run in 1:nruns) {
    if (random) {
      if (is.na(randomOrder[1])) {
        theOrder = sample(1:nrow(dfr))
      }
      else {
        theOrder = randomOrder
      }
      dfr = dfr[theOrder, ]
    }
    for (i in 1:nrow(dfr)) {
      currentCues = unlist(strsplit(dfr$cuevector[i], "_"))
      currentOutcomes = unlist(strsplit(dfr$outcomevector[i], 
                                        "_"))
      
      #Modified bit: old Vtotal = sum(weightvec[currentCues])
      Vtotal = sum(weightvec[currentCues])

      if (traceOutcome %in% currentOutcomes) {
        Lambda = lambda
      }
      else {
        Lambda = 0
      }

      if ( (Lambda > 0 & Vtotal < Lambda) | (Lambda==0 & Vtotal > 0) | !Humble  )
      {
      weightvec[currentCues] = weightvec[currentCues] + 
        alpha * beta1 * (Lambda - Vtotal)
      }
      cnt = cnt + 1
      res[cnt] = weightvec[traceCue]
    }
  }
  equilibriumWeights = estimateWeights(cuesOutcomes)
  result <- (list(weightvector = res, 
                  traceCue = traceCue, traceOutcome = traceOutcome, randomOrder = theOrder, eqw=0))
  class(result) <- "RescorlaWagner"
  return(result)
}


RescorlaWagner = function (cuesOutcomes, traceCue = "h", traceOutcome = "hand", 
          nruns = 1, random = TRUE, randomOrder = NA, alpha = 0.1, 
          lambda = 1, beta1 = 0.1, beta2 = 0.1, weightvec=NA) 
{
  cues = unique(unlist(strsplit(as.character(cuesOutcomes$Cues), 
                                "_")))
  outcomes = unique(unlist(strsplit(as.character(cuesOutcomes$Outcomes), 
                                    "_")))
  if(is.na(weightvec)) { weightvec = rep(0, length(cues)) }
  names(weightvec) = cues
  res = vector(mode = "numeric")
  dfr = data.frame(cuevector = rep(cuesOutcomes$Cues, cuesOutcomes$Frequency), 
                   outcomevector = rep(cuesOutcomes$Outcomes, cuesOutcomes$Frequency), 
                   stringsAsFactors = FALSE)
  theOrder = NA
  cnt = 0
  for (run in 1:nruns) {
    Vtotal=rep(0,length(cuesOutcomes$Cues))
    if (random) {
      if (is.na(randomOrder[1])) {
        theOrder = sample(1:nrow(dfr))
      }
      else {
        theOrder = randomOrder
      }
      dfr = dfr[theOrder, ]
    }
    for (i in 1:nrow(dfr)) {
      currentCues = unlist(strsplit(dfr$cuevector[i], "_"))
      currentOutcomes = unlist(strsplit(dfr$outcomevector[i], 
                                        "_"))
      Vtotal[i] = sum(weightvec[currentCues])
      if (traceOutcome %in% currentOutcomes) {
        Lambda = lambda
      }
      else {
        Lambda = 0
      }
      weightvec[currentCues] = weightvec[currentCues] + 
        alpha * beta1 * (Lambda - Vtotal[i])
      cnt = cnt + 1
      res[cnt] = weightvec[traceCue]
    }
  }
  equilibriumWeights = estimateWeights(cuesOutcomes)
  eqw = equilibriumWeights[traceCue, traceOutcome]
  result <- (list(weightvector = res, equilibriumWeight = eqw, 
                  traceCue = traceCue, traceOutcome = traceOutcome, randomOrder = theOrder, Vtotal = Vtotal))
  class(result) <- "RescorlaWagner"
  return(result)
}