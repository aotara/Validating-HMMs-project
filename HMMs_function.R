################################################################################
#         FUNCTIONS FOR 7 HMMs FITTED TO THE BOAT TRACKING DATA 
#                              &
#         VALIDATION OF BEHAVIOURAL STATES DECODED FROM HMMs
################################################################################
# Load libraries
library(momentuHMM) # HMM analysis
library(dplyr)
# validation
library(caret)
library(mgcv)
library(MLmetrics)
library(cvms)
################################################################################
# BASELINE MODEL - MODEL 0, COMPLETE POOLING EFFECT
HmmVal0 <- function(data,initPara){
  # FIT MODEL
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle  
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                     max = c(initPara[2], initPara[4]))), #sd 
                          # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                           max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                          c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    # fit model
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(mbase) # extract decoded states from base model 
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(mbase)  # extract state probability
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

# MODEL 1 - no pooling effect on state process (transition probability matrix)
HmmVal1<- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit model 1
  Par0_m1 <- getPar0(model = mbase, formula = ~ID) # initial parameters
  m1 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m1$Par,
               stateNames = stateNames, formula = ~ID,estAngleMean = list(angle=T),
               nlmPar=list(hessian=F) )
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m1) # extract decoded states from model 1
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss -  creating dataframe needed
  state_prob <- stateProbs(m1)  # extract state probability from model 1
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

# MODEL 2 - no pooling effect on observed process (step length)
HmmVal2<- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit Model 2
  if (length(initPara) == 8) {
    DM <- list(step = list(mean = ~ID, sd = ~ID))
  } else{
    DM <- list(step = list(mean = ~ID, sd = ~ID , zeromass = ~ID))
  }
  Par0_m2 <- getPar0(model=mbase, DM=DM)   
  m2 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
               DM = DM ,estAngleMean = list(angle=T),stateNames = stateNames,
               nlmPar=list(hessian=F))
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m2) # extract decoded states from model 2
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(m2)  # extract state probability from model 2
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

# MODEL 3 - no pooling effect on observed process (step length) and state process
HmmVal3<- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20    # try different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit Model 3
  if (length(initPara) == 8) {
    DM <- list(step = list(mean = ~ID, sd = ~ID))
  } else{
    DM <- list(step = list(mean = ~ID, sd = ~ID , zeromass = ~ID))
  }
  Par0_m3 <- getPar0(model = mbase,DM=DM,formula=~ID)   
  m3 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
               formula=~ID,DM=DM , estAngleMean = list(angle=T),
               stateNames = stateNames, nlmPar=list(hessian=F))
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m3) # extract decoded states from model 3
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(m3)  # extract state probability from model 3
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

# MODEL 4 - covariate (distance to colony) effect on state process
HmmVal4<- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit model 4
  Par0_m4<- getPar0(model= mbase,formula= ~DIST2COL)#initial parameters
  m4 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m4$Par,
               stateNames = stateNames, formula = ~DIST2COL,nlmPar=list(hessian=F),
               estAngleMean = list(angle=T) )
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m4) # extract decoded states from model 4
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(m4)  # extract state probability from model 4
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD<- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll<- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

# MODEL 5 - covariate (distance to colony) effect on observed process
HmmVal5<- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit model 5
  if (length(initPara) == 8) {
    DM <- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL))
  } else{
    DM <- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL, zeromass = ~DIST2COL))
  }
  Par0_m5 <- getPar0(model = mbase, DM = DM) # initial parameters
  m5 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m5$Par,
              DM = DM ,estAngleMean = list(angle=T),stateNames = stateNames,
              nlmPar=list(hessian=F))
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m5) # extract decoded states from model 5
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(m5)  # extract state probability from model 5
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4))))  
}

# MODEL 6 - covariate (distance to colony) effect on state and observed process
HmmVal6 <- function(data,initPara){
  # Specify parameters
  set.seed(1908)  #reproducibility
  niter <- 20     #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle   
    if (length(initPara) == 8) {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4]))), #sd 
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8]))))
    } else {
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
    }
    stateNames <- c('Foraging', 'Not-foraging')
    allm[[i]] <- momentuHMM::fitHMM(data = data, nbStates = 2, dist = dist,
                                    stateNames = stateNames, Par0 = Par0,
                                    estAngleMean = list(angle = T),
                                    nlmPar = list(hessian = F))
  }
  # Extract likelihoods of fitted model
  allnllk <- unlist(lapply(allm, function(j) j$mod$minimum))
  # Index of best fitting model (smallest negative log-likelihood)
  whichbest <- which.min(allnllk)
  # Best fitting model
  mbase <- allm[[whichbest]] 
  # fit model 6
  if (length(initPara) == 8) {
    DM <- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL))
  } else{
    DM <- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL, zeromass = ~DIST2COL))
  }
  Par0_m6 <- getPar0(model=mbase, DM=DM, formula= ~DIST2COL) # initial parameters
  m6 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m6$Par,
               DM=DM, formula=~DIST2COL,estAngleMean = list(angle=T),
               stateNames = stateNames,nlmPar=list(hessian=F))
  #VALIDATION
  
  decodedS <- momentuHMM::viterbi(m6) # extract decoded states from model 6
  # specify observed and decoded states as factor
  decodedS <- as.factor(decodedS)
  levels(decodedS) <- c('Foraging','Not-foraging')
  # Confusion matrix
  cm <- data.frame(data$observedS, decodedS)
  colnames(cm) <- c("Actual", "Predicted")
  u <- confusionMatrix(cm$Predicted, cm$Actual)    
  conMat<-knitr::kable(u$byClass,'simple',digits=4, caption = "Validation Metrics")
  # Log-loss - creating dataframe needed
  state_prob <- stateProbs(m6)  # extract state probability from model 6
  newD <- data.frame(cbind(state_prob,cm) )
  # rename columns
  colnames(newD)[3] <- 'observed_behaviour'
  colnames(newD)[4] <- 'decoded_behaviour'
  # Specifying 1 for observed foraging, 0 for not-foraging 
  newD <- newD %>%
    mutate(obsForage = if_else(.$observed_behaviour == 'Foraging', 1, 0))
  # calculate log-loss
  ll <- LogLoss(newD$obsForage, newD$Foraging) 
  return(list(conMat,paste("Log-loss is:", round(ll,4)) ))  
}

