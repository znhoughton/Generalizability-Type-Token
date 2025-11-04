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
              if (ap==T) {
              weightvec[absentCues] = weightvec[absentCues] - 
                alpha2 * beta1 * (lambda - Vtotal)
            }}
            
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


RescorlaLinear = function (cuesOutcomes, traceCue = "h", traceOutcome = "hand", 
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
      
      Vtotal = sum(weightvec[currentCues])
      #Vtotal = 2*invlogit(sum(weightvec[currentCues]))-1 #this is the tanh function to scale between -1 and 1
      #tanh is needed for second-order retrospective revaluation, to pass on inhibition
      #Vtotal = invlogit(sum(weightvec[currentCues]))
      
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
        if (ap==T) {
          weightvec[absentCues] = weightvec[absentCues] - 
            alpha2 * beta1 * (lambda - Vtotal)
        }}
      
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

Rescorla = function (cuesOutcomes,  
                     Lambda = 1, 
                     Alpha = 0.1, 
                     Beta = 0.1,logistic=FALSE,sigSlope=1) 
{
  cues = unique(unlist(strsplit(cuesOutcomes$Cues, "_")))
  outcomes = unique(unlist(strsplit(cuesOutcomes$Outcomes, "_")))
  frequency = cuesOutcomes$Frequency
  w = matrix(0, length(cues), length(outcomes))
  rownames(w)=cues
  colnames(w)=outcomes
  theCues = strsplit(cuesOutcomes$Cues, "_")
  theOutcomes = strsplit(cuesOutcomes$Outcomes, "_")
  rate=Alpha*Beta
  
  # new bit
  knownOutcomes = NULL
  x=-0.5
  
  for (i in 1:nrow(cuesOutcomes)) {
    cs = theCues[[i]]
    ou = theOutcomes[[i]] 
    
    # new two lines
    toAdd = ou[!is.element(ou, knownOutcomes)]
    knownOutcomes = c(knownOutcomes, toAdd)
    
    if (logistic==FALSE) {Vtotal = colSums(w[cs,ou,drop=FALSE])}
    else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,ou,drop=FALSE])))}      
    
    w[cs,ou] = t(t(w[cs,ou]) + (Lambda-Vtotal)*rate*frequency[i])
    
    # replaced next line by new line
    # otherou = outcomes[!is.element(outcomes,ou)]
    otherou = knownOutcomes[!is.element(knownOutcomes, ou)]
    # this next line may not work, in which case check if (length(otherou)>0)
    if (!is.null(otherou)) {
      if (logistic==FALSE) {Vtotal = colSums(w[cs,otherou,drop=FALSE])}
      else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,otherou,drop=FALSE])))}
      #1/(1+exp((-1)*sigSlope*sum(weightvec[currentCues])))
      w[cs,otherou] = t(t(w[cs,otherou]) + (0-Vtotal)*rate*frequency[i])
    }
  }
  return(w)
}

BushMosteller = function (cuesOutcomes,  
                     Lambda = 1, 
                     Alpha = 0.1, 
                     Beta = 0.1,logistic=FALSE,sigSlope=1) 
{
  cues = unique(unlist(strsplit(cuesOutcomes$Cues, "_")))
  outcomes = unique(unlist(strsplit(cuesOutcomes$Outcomes, "_")))
  frequency = cuesOutcomes$Frequency
  w = matrix(0, length(cues), length(outcomes))
  rownames(w)=cues
  colnames(w)=outcomes
  theCues = strsplit(cuesOutcomes$Cues, "_")
  theOutcomes = strsplit(cuesOutcomes$Outcomes, "_")
  rate=Alpha*Beta
  
  # new bit
  knownOutcomes = NULL
  x=-0.5
  
  for (i in 1:nrow(cuesOutcomes)) {
    cs = theCues[[i]]
    ou = theOutcomes[[i]] 
    
    # new two lines
    toAdd = ou[!is.element(ou, knownOutcomes)]
    knownOutcomes = c(knownOutcomes, toAdd)
    
    if (logistic==FALSE) {Vtotal = colSums(w[cs,ou,drop=FALSE])}
    else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,ou,drop=FALSE])))}      
    
    w[cs,ou] = t(t(w[cs,ou])) + (Lambda-t(t(w[cs,ou,drop = F])))*rate
    
    # replaced next line by new line
    # otherou = outcomes[!is.element(outcomes,ou)]
    otherou = knownOutcomes[!is.element(knownOutcomes, ou)]
    # this next line may not work, in which case check if (length(otherou)>0)
    if (!is.null(otherou)) {
      if (logistic==FALSE) {Vtotal = colSums(w[cs,otherou,drop=FALSE])}
      else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,otherou,drop=FALSE])))}
      #1/(1+exp((-1)*sigSlope*sum(weightvec[currentCues])))
      w[cs,otherou] = t(t(w[cs,otherou])) + (0-t(t(w[cs,otherou,drop = F])))*rate
    }
  }
  return(w)
}


Rescorla_modified_alpha_for_semantic_cues = function (cuesOutcomes,  
                     Lambda = 1, 
                     Alpha = 0.05,
                     SemAlpha = 0.25,
                     Beta = 0.1,logistic=FALSE,sigSlope=1) 
{
  cues = unique(unlist(strsplit(cuesOutcomes$Cues, "_")))
  outcomes = unique(unlist(strsplit(cuesOutcomes$Outcomes, "_")))
  frequency = cuesOutcomes$Frequency
  w = matrix(0, length(cues), length(outcomes))
  rownames(w)=cues
  colnames(w)=outcomes
  theCues = strsplit(cuesOutcomes$Cues, "_")
  theOutcomes = strsplit(cuesOutcomes$Outcomes, "_")
  rate=Alpha*Beta
  semrate = SemAlpha*Beta
  
  # new bit
  knownOutcomes = NULL
  x=-0.5
  
  for (i in 1:nrow(cuesOutcomes)) {
    cs = theCues[[i]]
    ou = theOutcomes[[i]] 
    
    # new two lines
    toAdd = ou[!is.element(ou, knownOutcomes)]
    knownOutcomes = c(knownOutcomes, toAdd)
    
    if (logistic==FALSE) {Vtotal = colSums(w[cs,ou,drop=FALSE])}
    else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,ou,drop=FALSE])))}      
    
    cue_names = rownames(w[cs,ou, drop=F])
    rates = ifelse(cue_names %in% semantic_cues, semrate, rate)
    
    w[cs,ou] = t(t(w[cs,ou]) + (Lambda-Vtotal)*rates*frequency[i])
    
    # replaced next line by new line
    # otherou = outcomes[!is.element(outcomes,ou)]
    otherou = knownOutcomes[!is.element(knownOutcomes, ou)]
    # this next line may not work, in which case check if (length(otherou)>0)
    if (!is.null(otherou)) {
      if (logistic==FALSE) {Vtotal = colSums(w[cs,otherou,drop=FALSE])}
      else {Vtotal = 1/(1+exp((-1)*sigSlope*colSums(w[cs,otherou,drop=FALSE])))}
      #1/(1+exp((-1)*sigSlope*sum(weightvec[currentCues])))
      
      otherou_cue_names = rownames(w[cs,otherou, drop=F])
      rates_otherou = ifelse(otherou_cue_names %in% semantic_cues, semrate, rate)
      adj_matrix = matrix(rep((0 - Vtotal) * frequency[i], each=length(cs)),
                          nrow=length(cs), byrow=FALSE)
      
      # Multiply each row by its rate using sweep
      adj_matrix = sweep(adj_matrix, 1, rates_otherou, "*")
      w[cs,otherou] = t(t(w[cs,otherou])) + adj_matrix
      
    }
  }
  return(w)
}
