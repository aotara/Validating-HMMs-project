
################################################################################
# CODE FOR RESULTS IN TABLE 3 (PAGE 19) INCUBATION
# NOTE: RUN ALL FUNCTION IN 'HMMs funtcion.R' before running codes contained here
################################################################################

# load libraries
library(momentuHMM)
library(dplyr)

#COLONY - Cockle ; SPECIES - Arctic IN; MODEL 5 ################################

data<-read.csv('Cockle_Arctic_IN.csv')
View(data)
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST_COL_K,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = DIST_COL_K) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.005,0.012,0.015,1,2,3,5,0.003,0.0009) #starting parameters
HmmVal5(prep_data,init) #call function for model 5
################################################################################

#COLONY - Isle of May ; SPECIES - Arctic IN; MODEL 6 ###########################

data<-read.csv('IsleMay_Arctic_IN.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL_km,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = DIST2COL_km) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
prep_data<-prep_data[-c(6498),]
init<-c(0.002,0.005,0.009,0.012,1,2,3,5,0.003,0.0008) #starting parameters
HmmVal6(prep_data,init) #call function for model 6
################################################################################

#COLONY - Blue Circle ; SPECIES - Roseate IN; MODEL 4 ##############################

data<-read.csv('BlueCircle_Roseate_IN.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST_COL_K,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = DIST_COL_K) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.004,0.006,0.01,1,2,3,5,0.082,0.0001) #starting parameters
HmmVal4(prep_data,init) #call function for model 4
################################################################################

#COLONY - Leith ; SPECIES - Common IN; MODEL 0-base model ######################

data<-read.csv('Leith_Common_IN.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,TSECS_BST,LATITUDE,LONGITUDE,Dist_to_colony_km,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = Dist_to_colony_km) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
prep_data<-prep_data[-c(4074,10643,10644),]
init<-c(0.002,0.004,0.006,0.01,1,2,3,5,0.007,0.0005) #starting parameters
HmmVal0(prep_data,init) #call function for model 0
################################################################################

#COLONY - Cockle ; SPECIES - Sandwich IN; MODEL 6 ##############################

data<-read.csv('Cockle_Sandwich_IN.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST_COL_K,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = DIST_COL_K) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.005,0.008,0.011,1,2,3,5,0.002,0.0003) #starting parameters
HmmVal6(prep_data,init) #call function for model 6




