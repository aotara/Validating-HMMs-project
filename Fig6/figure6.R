################################################################################
#   Figure 6 Page 18
# Confusion Matrix Metrics for inferred state from HMM and Tern position 
################################################################################

#################    STEPS             #########################################
# (1) MERGING LOCATION AND BEHAVIOURAL DATA OF VISUALLY TRACKED TERNS IN 
# COQUET ISLAD, 2009 DURING CHICK-REARING

# INFORMATION ON BEARING AND DISTANCE OF BIRD FROM BOAT IS CONTAINED IN THE 
# BEHAVIOURAL DATA (dbeh) 

# (2) CALCULATE APPROXIMATE LOCATION DATA OF TERNS

# (3) FIT HMMs TO BOAT AND APPROXIMATE TERN LOCATION DATA

# (4) EXTRACT METRICS TO CREATE DATAFRAME FOR PLOTTING

# (5) PLOT 
################################################################################

# load libraries
library(dplyr)
library(tidyr) 
library("readxl")
library(momentuHMM)
library(NISTunits)

###############       (1)              #########################################

# read location data 
dloc<-read_excel('coquet-2009-tern-tracks.xlsx')
dloc<- dloc %>% 
  arrange(dloc,TRACK_ID, TSECS)%>%
  rename(ID = TRACK_ID)

# read behavioural data with bearing, distance inclusive
dbeh<-read_excel('coquet-2009-behaviour-locs.xlsx')
dbeh<- dbeh %>% 
  arrange(Track_ID,TSecs_GMT) %>%
  filter(!is.na(Bearing) , !is.na(Distance_From_Boat)) %>%
  rename(ID = Track_ID) %>%
  mutate(LonTern = 0) %>%
  mutate(LatTern = 0)

# merge location and behavioural data
l<-merge(dloc, dbeh, by.x=c("ID", 'TSECS','SPECIES'), 
         by.y=c("ID", "TSecs_GMT",'SPECIES'),all.x=T)

new_data<-l %>% 
  group_by(ID) %>% 
  fill(Behaviour_type,Bearing, Distance_From_Boat, .direction = "downup")


########################## (2)  ################################################
# calculate approximate tern location data
for(k in 1:length(new_data$ID)) {
  ll<- geosphere::destPoint(c(new_data$LONGITUDE[k], new_data$LATITUDE[k]),
                            b = new_data$Bearing[k],
                            d = new_data$Distance_From_Boat[k])
  new_data$LonTern[k] <- ll[1]
  new_data$LatTern[k]  <- ll[2]
}

#select columns
data<- subset(new_data, select = c(ID,SPECIES,DATE,TSECS,LATITUDE,LONGITUDE,
                                   DIST2COL, Behaviour_type,LatTern,LonTern))

###########################    (3)   ###########################################
#       subsetting merged data with bearing and distance recorded for 
#          both observed foraging and not foraging behaviours
################################################################################

# ARCTIC TERN species 
dArctic<- data%>% filter(ID==117)
write.csv(dArctic, file ='Arctic-boat-approximate-tern-location.csv',row.names=F) 
dArctic<-read.csv('Arctic-boat-approximate-tern-location.csv')
# subset with boat latitude and longitude columns
d1a<- subset(dArctic, select = -c(LonTern,LatTern) )
# subset with tern latitude and longitude columns 
d1b<- subset(dArctic, select = -c(LONGITUDE,LATITUDE) )
# rename column of tern position data
names(d1b)[names(d1b) == 'LatTern'] <- 'Latitude'
names(d1b)[names(d1b) == 'LonTern'] <- 'Longitude'
# calculate step length and turning angle
p1a<-prepData(d1a,type='LL',coordNames=c("LONGITUDE","LATITUDE"),covNames=c('DIST2COL'))
p1b<-prepData(d1b,type='LL',coordNames=c("Longitude","Latitude"),covNames=c('DIST2COL'))

# checking row corresponding to outliers
nrow<-which(p1b$step>0.02)
nrow
p1b<-p1b[-c(186,237, 365, 480, 500, 564 ,605 ,849 ,1106),]
p1a<-p1a[-c(186,237, 365, 480, 500, 564 ,605 ,849 ,1106),]

# plot histogram of step length and turning angle
par(mfrow=c(1,2))
hist(p1a$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (boat)")
hist(p1b$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (tern)")
hist(p1a$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (boat)")
hist(p1b$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (tern)")

## HMM analysis: Arctic|Track ID=117|Coquet colony 2009 boat data
HmmValArctic <- function(data,initPara){
  # FIT MODEL
  set.seed(1908)  #reproducibility
  niter <- 10     #different starting values
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
  
  #VALIDATION
  M0decodedS <- momentuHMM::viterbi(mbase) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M0decodedS <- as.factor(M0decodedS)
  levels(M0decodedS) <- c('Foraging','Not-foraging')
  
  
  # MODEL 4 covariate (distance to colony) effect on state process
  Par0_m4<- getPar0(model= mbase,formula= ~DIST2COL)#initial parameters: mbase
  # fit model
  m4 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m4$Par,
               stateNames = stateNames, formula = ~DIST2COL, 
               estAngleMean = list(angle=T), nlmPar = list(hessian=F))
  #VALIDATION
  M4decodedS <- momentuHMM::viterbi(m4) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M4decodedS <- as.factor(M4decodedS)
  levels(M4decodedS) <- c('Foraging','Not-foraging')
  
  # 2a. covariate (distance to colony) on the observed process
  DM<- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL,zeromass=~DIST2COL)) 
  Par0_m5<- getPar0(model = mbase, DM = DM) # initial parameters: mbase
  #fit model
  m5<- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m5$Par,
              DM = DM ,estAngleMean = list(angle=T),stateNames = stateNames,
              nlmPar=list(hessian=F))
  #VALIDATION
  M5decodedS <- momentuHMM::viterbi(m5) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M5decodedS <- as.factor(M5decodedS)
  levels(M5decodedS) <- c('Foraging','Not-foraging')
  
  # Model 6 covariate (distance to colony) on state and observed process
  
  DM <- list(step = list(mean = ~DIST2COL,sd = ~DIST2COL,zeromass=~ DIST2COL )) 
  Par0_m6 <- getPar0(model=mbase, DM=DM, formula= ~DIST2COL) 
  #fit model
  m6 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m6$Par,
               DM=DM, formula=~DIST2COL ,estAngleMean = list(angle=T),
               stateNames = stateNames,nlmPar=list(hessian=F))
  #VALIDATION
  M6decodedS <- momentuHMM::viterbi(m6) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M6decodedS <- as.factor(M6decodedS)
  levels(M6decodedS) <- c('Foraging','Not-foraging')
  
  return(list(M0decodedS, M4decodedS, M5decodedS, M6decodedS))  
}

# fit hmm on arctic tern boat and approximate location data - coquet 2009
init<-c(0.001,0.003,0.006,0.009,1,2,4,6,0.017,0.0003) #starting parameters
boatA<-HmmValArctic(p1a,init) #fit HMM to boat location data 
ternA<-HmmValArctic(p1b,init) #fit HMM to approxiate tern location data 

# Confusion matrix for decoded states from boat & approximate tern location data
# resulting matrix is used in generating dataframe in 'generate_dataframe.R' file
confusionMatrix(boatA[[1]],ternA[[1]]) #1st output -model0
confusionMatrix(boatA[[2]],ternA[[2]]) #2nd output -model4
confusionMatrix(boatA[[3]],ternA[[3]]) #3rd output -model5
confusionMatrix(boatA[[4]],ternA[[4]]) #4th output -model6


# COMMON species data with bearing and distance recorded for both observed foraging 
# and not foraging behaviours
dCommon<- subset(data, ID == 41 | ID==107)
write.csv(dCommon, file ='Common-boat-approximate-tern-location.csv',row.names=F) 
dCommon<-read.csv('Common-boat-approximate-tern-location.csv')
# subset with boat latitude and longitude columns
d1a<- subset(dCommon, select = -c(LonTern,LatTern) )
# subset with tern latitude and longitude columns 
d1b<- subset(dCommon, select = -c(LONGITUDE,LATITUDE) )
# rename column of tern position data
names(d1b)[names(d1b) == 'LatTern'] <- 'Latitude'
names(d1b)[names(d1b) == 'LonTern'] <- 'Longitude'
# calculate step length and turning angle
p1a<-prepData(d1a,type='LL',coordNames=c("LONGITUDE","LATITUDE"),covNames=c('DIST2COL'))
p1b<-prepData(d1b,type='LL',coordNames=c("Longitude","Latitude"),covNames=c('DIST2COL'))

# checking row corresponding to outliers
nrow<-which(p1b$step>0.02)
nrow
p1b<-p1b[-c(25,100,125,181,267,291,330,375,626,955,1095,1160,1855,1870, 2135, 
            2160,2445,2456,2496,2636,2675,2681,2976,3001,3126,3216,3448),]
p1a<-p1a[-c(25,100,125,181,267,291,330,375,626,955,1095,1160,1855,1870, 2135, 
            2160,2445,2456,2496,2636,2675,2681,2976,3001,3126,3216,3448),]

# plot histogram of step length and turning angle
par(mfrow=c(1,2))
hist(p1a$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (boat)")
hist(p1b$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (tern)")
hist(p1a$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (boat)")
hist(p1b$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (tern)")

## HMM analysis fucntion for Common and Sanwich tern Coquet colony 2009
HmmValCS <- function(data,initPara){
  # FIT MODEL
  set.seed(1908)  #reproducibility
  niter <- 10     #different starting values
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
  
  #VALIDATION
  M0decodedS <- momentuHMM::viterbi(mbase) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M0decodedS <- as.factor(M0decodedS)
  levels(M0decodedS) <- c('Foraging','Not-foraging')
  
  # MODEL 1 no polling effect on state process
  Par0_m1 <- getPar0(model = mbase, formula = ~ID) # initial parameters
  m1 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m1$Par,
               stateNames = stateNames, formula = ~ID,estAngleMean = list(angle=T),
               nlmPar=list(hessian=F) )
  #VALIDATION
  M1decodedS <- momentuHMM::viterbi(m1) # extract decoded states from model 1
  # specify observed and decoded states as factor
  M1decodedS <- as.factor(M1decodedS)
  levels(M1decodedS) <- c('Foraging','Not-foraging')
  
  # MODEL 2 no pooling effect on the observed process
  DM <- list(step = list(mean = ~ID, sd = ~ID , zeromass = ~ID))
  Par0_m2 <- getPar0(model=mbase, DM=DM)   
  m2 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m2$Par,
               DM = DM ,estAngleMean = list(angle=T),stateNames = stateNames,
               nlmPar=list(hessian=F))
  #VALIDATION
  M2decodedS <- momentuHMM::viterbi(m2) # extract decoded states from model 2
  # specify observed and decoded states as factor
  M2decodedS <- as.factor(M2decodedS)
  levels(M2decodedS) <- c('Foraging','Not-foraging')
  
  # MODEL3 no pooling effect on state and observed process
  DM <- list(step = list(mean = ~ID, sd = ~ID , zeromass = ~ID))
  Par0_m3 <- getPar0(model = mbase,DM=DM,formula=~ID)   
  m3 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m3$Par,
               formula=~ID,DM=DM , estAngleMean = list(angle=T),
               stateNames = stateNames, nlmPar=list(hessian=F))
  #VALIDATION
  M3decodedS <- momentuHMM::viterbi(m3) # extract decoded states from model 3
  # specify observed and decoded states as factor
  M3decodedS <- as.factor(M3decodedS)
  levels(M3decodedS) <- c('Foraging','Not-foraging')
  
  # MODEL 4 covariate (distance to colony) effect on state process
  Par0_m4<- getPar0(model= mbase,formula= ~DIST2COL)#initial parameters: mbase
  # fit model
  m4 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m4$Par,
               stateNames = stateNames, formula = ~DIST2COL, 
               estAngleMean = list(angle=T), nlmPar = list(hessian=F))
  #VALIDATION
  M4decodedS <- momentuHMM::viterbi(m4) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M4decodedS <- as.factor(M4decodedS)
  levels(M4decodedS) <- c('Foraging','Not-foraging')
  
  # 2a. covariate (distance to colony) on the observed process
  DM<- list(step = list(mean= ~DIST2COL,sd= ~DIST2COL,zeromass=~DIST2COL)) 
  Par0_m5<- getPar0(model = mbase, DM = DM) # initial parameters: mbase
  #fit model
  m5<- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m5$Par,
              DM = DM ,estAngleMean = list(angle=T),stateNames = stateNames,
              nlmPar=list(hessian=F))
  #VALIDATION
  M5decodedS <- momentuHMM::viterbi(m5) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M5decodedS <- as.factor(M5decodedS)
  levels(M5decodedS) <- c('Foraging','Not-foraging')
  
  # Model 6 covariate (distance to colony) on state and observed process
  
  DM <- list(step = list(mean = ~DIST2COL,sd = ~DIST2COL,zeromass=~ DIST2COL )) 
  Par0_m6 <- getPar0(model=mbase, DM=DM, formula= ~DIST2COL) 
  #fit model
  m6 <- fitHMM(data = data, nbStates = 2, dist = dist, Par0 = Par0_m6$Par,
               DM=DM, formula=~DIST2COL ,estAngleMean = list(angle=T),
               stateNames = stateNames,nlmPar=list(hessian=F))
  #VALIDATION
  M6decodedS <- momentuHMM::viterbi(m6) # extract decoded states from base model 
  # specify observed and decoded states as factor
  M6decodedS <- as.factor(M6decodedS)
  levels(M6decodedS) <- c('Foraging','Not-foraging')
  
  
  return(list(M0decodedS,M1decodedS,M2decodedS,M3decodedS,M4decodedS,M5decodedS,M6decodedS))  
}
# fit hmm on common tern boat and approximate location data
init<-c(0.001,0.003,0.007,0.011,1,2,3,5,0.0002,0.00008) #starting parameters
boatC<-HmmValCS(p1a,init) #fit HMM to boat location data 
ternC<-HmmValCS(p1b,init) #fit HMM to approximate common tern location data 

# Confusion matrix for decoded states from boat & approximate tern location data
# resulting matrix is used in generating dataframe in 'generate_dataframe.R' file
confusionMatrix(boatC[[1]],ternC[[1]]) #1st output -model0
confusionMatrix(boatC[[2]],ternC[[2]]) #2nd output -model1
confusionMatrix(boatC[[3]],ternC[[3]]) #3rd output -model2
confusionMatrix(boatC[[4]],ternC[[4]]) #4th output -model3
confusionMatrix(boatC[[5]],ternC[[5]]) #5th output -model4
confusionMatrix(boatC[[6]],ternC[[6]]) #6th output -model5
confusionMatrix(boatC[[7]],ternC[[7]]) #7th output -model6

# SANDWICH TERN species data with bearing and distance recorded for both observed 
#foraging and not foraging behaviours 

dSandwich<- subset(data, ID == 20 | ID == 21 | ID==22 | ID == 35|ID == 37)
write.csv(dSandwich, file ='Sandwich-boat-approximate-tern-location.csv',row.names=F) 
dSandwich<-read.csv('Sandwich-boat-approximate-tern-location.csv')
# subset with boat latitude and longitude columns
d1a<- subset(dSandwich, select = -c(LonTern,LatTern) )
# subset with tern latitude and longitude columns 
d1b<- subset(dSandwich, select = -c(LONGITUDE,LATITUDE) )
# rename column of tern position data
names(d1b)[names(d1b) == 'LatTern'] <- 'Latitude'
names(d1b)[names(d1b) == 'LonTern'] <- 'Longitude'
# calculate step length and turning angle
p1a<-prepData(d1a,type='LL',coordNames=c("LONGITUDE","LATITUDE"),covNames=c('DIST2COL'))
p1b<-prepData(d1b,type='LL',coordNames=c("Longitude","Latitude"),covNames=c('DIST2COL'))

# plot histogram of step length and turning angle
par(mfrow=c(1,2))
hist(p1a$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (boat)")
hist(p1b$angle, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Turning angle (tern)")
hist(p1a$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (boat)")
hist(p1b$step, col = blues9[4], border = blues9[6], las = 1, main='',
     ylab = "Frequency",xlab="Step length (tern)")

# checking row corresponding to outliers
nrow<-which(p1b$step>0.02)
nrow
p1b<-p1b[-c( 76,116,153,376,501,761,781,811,851,901,1036,1261,1381,1431,1795,1895,
             1900,2116,2146, 2173, 2446, 2715, 2736, 2877, 2887, 2905, 3041, 3120,
             3169, 3286, 3341, 3396, 3431, 3486, 3541, 3571,3613, 3656, 3811, 3856,
             3996, 4061, 4196, 4261, 4306, 4346, 4398, 4501, 4516, 4556, 4586, 4686,
             4716, 4731,5093,5105, 5170, 5195, 5231, 5269, 5352, 5363, 5503, 5625,
             5655, 5705, 5732, 5845, 5897, 5970, 6062, 6125,6161, 6231, 6383, 6415),]
p1a<-p1a[-c(76,116,153,376,501,761,781,811,851,901,1036,1261,1381,1431,1795,1895,
            1900,2116,2146, 2173, 2446, 2715, 2736, 2877, 2887, 2905, 3041, 3120,
            3169, 3286, 3341, 3396, 3431, 3486, 3541, 3571,3613, 3656, 3811, 3856,
            3996, 4061, 4196, 4261, 4306, 4346, 4398, 4501, 4516, 4556, 4586, 4686,
            4716, 4731,5093,5105, 5170, 5195, 5231, 5269, 5352, 5363, 5503, 5625,
            5655, 5705, 5732, 5845, 5897, 5970, 6062, 6125,6161, 6231, 6383, 6415),]
whichzero <- which(p1a$step==0)
length(whichzero)/nrow(p1a)

# fit hmm on sandwich tern boat and approximate location data-coquet 2009
init<-c(0.001,0.003,0.006,0.012,1,2,3,5,0.016,0.0009) #starting parameters
boatS<-HmmValCS(p1a,init) #fit HMM to boat location data
ternS<-HmmValCS(p1b,init) #fit HMM to approximate sandwich tern location data

# Confusion matrix for decoded states from boat & approximate tern location data
# resulting matrix is used in generating dataframe in 'generate_dataframe.R' file
confusionMatrix(boatS[[1]],ternS[[1]]) #1st output -model0
confusionMatrix(boatS[[2]],ternS[[2]]) #2nd output -model1
confusionMatrix(boatS[[3]],ternS[[3]]) #3rd output -model2
confusionMatrix(boatS[[4]],ternS[[4]]) #4th output -model3
confusionMatrix(boatS[[5]],ternS[[5]]) #5th output -model4
confusionMatrix(boatS[[6]],ternS[[6]]) #6th output -model5
confusionMatrix(boatS[[7]],ternS[[7]]) #7th output -model6

##########################   (4)   #############################################
# GENERATING DATAFARME FOR PLOTTING
################################################################################

## True positive-FF, False negative-FN, False positive-NF, True negative-NN
#metrics
metric<-c('FF','FN','NF','NN')
metrics<-c( (rep(metric,each=4)), (rep(metric,each=7)), (rep(metric,each=7)) )

# species
species<-c((rep('Arctic',16)),(rep('Common',28)),(rep('Sandwich',28)))

#model
mA<-c('Model 0','Model 4','Model 5','Model 6') #Arctic tern
models <- mA %>% 
  c(rep(mA,3))

# common and sandwich tern
mCS<-c('Model 0','Model 1','Model 2','Model 3','Model 4','Model 5','Model 6')
models<-models %>% 
  c(rep(mCS,8))

#row1 FF row2 FN row3 NF row4 NN
# confusion matrix metrics are obtained from STEP 3
value <- c(357,358,352,352,
           0,0  , 0 ,5 ,
           8  ,7  , 4 ,8  , 
           733,733,742,733,
           497 ,491 ,635 ,587 ,500 ,502 ,502 ,
           44  ,24  ,133  ,184 ,43  ,24  ,24  ,  
           26  ,22  ,13  ,12  ,22  ,26  ,26  ,
           2907,2937,2693,2691,2909,2922,2922, 
           2319,2309,2504,2462,2313,2360,2309,
           1  ,0 ,  1 ,  49 ,0  ,9  ,0  ,
           39  ,38  ,44  ,42 ,37  ,30  ,48  ,
           4011,4023,3821,3817,4020,3971,4013
)
# generate dataframe
dj <- data.frame(species,models,metrics,value)
write.csv(dj,file="confusion-matrix-metric_assessing-visual-tracking-data.csv",
          row.names=FALSE)

############################## (5)   ###########################################
# plot figure 6
# read data
dataCMM<-read.csv('confusion-matrix-metric_assessing-visual-tracking-data.csv')
# plot 
dataCMM %>%
  ggplot(aes(x = models , y = value, fill = metrics)) +
  geom_bar(position = "stack", stat = "identity") +
  # axis labels
  xlab("Fitted HMMs")+
  ylab("Total count of metrics")+
  # plot theme
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 16, face = "bold"),
        strip.background.x = element_rect(fill = "white", colour = "black") )+
  theme(axis.text = element_text(size = 16, colour = "black", face = "bold"),
        axis.title = element_text(size = 16, face = "bold")) +
  theme(legend.text = element_text(size = 16, face = "bold"),
        legend.title = element_text(size = 16, face = "bold"))+
  scale_fill_manual(values = c("FF" = "grey50", "FN" = "#E69F00", 
                               "NF" = "#00BFFF", "NN" = '#104E8B'), name = "Metrics")+
  facet_wrap(~ species, scales = "free_y", nrow = 1)
