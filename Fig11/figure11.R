################################################################################
#              Figure 11 Page 23 
# HMM-fitted state-dependent distribution coloured by decoded behavioural state
# for SouthShian colony
################################################################################

# load libraries
library(momentuHMM)
library(dplyr)

#############################       STEP 1           ###########################
# function for HMM analysis

HmmVal0 <- function(data,initPara){
  # FIT MODEL
  set.seed(1908)  #reproducibility
  niter <- 3    #different starting values
  allm  <- list()  #Save list of fitted models
  for(i in 1:niter) {
    dist <- list(step = 'gamma', angle = 'vm') #specify distribution
    # range of starting values for step length mean and turning angle  
      Par0 <- list(step = c(runif(2, min = c(initPara[1], initPara[3]), 
                                  max = c(initPara[2], initPara[4])),  #mean
                            runif(2, min = c(initPara[1], initPara[3]),
                                  max = c(initPara[2], initPara[4])), #sd 
                            c(initPara[9],initPara[10])),#zeromass
                   # mean, concentration                        
                   angle = c(0, 0, runif(2, min = c(initPara[5], initPara[7]),
                                         max = c(initPara[6], initPara[8])))) 
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
  
  return(plot(mbase,plotTracks=F))  
}

##############################          STEP 2        ##########################

# read data and prepare for HMM analysis
data<-read.csv('SouthShian_Common_CR.csv')
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)
# obtain step length and turning angles
prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.004,0.008,0.01,1,3,3,5,0.014,0.0001) #starting parameters

####################           STEP 3       ####################################
HmmVal0(prep_data,init) #call function  to return figure 11 as output
