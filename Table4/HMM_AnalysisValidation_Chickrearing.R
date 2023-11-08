

################################################################################
# CODE FOR RESULTS IN TABLE 4 (PAGE 19) CHICK-REARING
# NOTE: RUN ALL FUNCTION IN 'HMMs funtcion.R' before running codes contained here
################################################################################
# Load libraries
library(momentuHMM) # HMM analysis
library(dplyr)
# validation
library(caret)
library(mgcv)
library(MLmetrics)
library(cvms)

#COLONY - Coquet ; SPECIES - Roseate; MODEL 4    ###############################

data<- read.csv("Coquet_Roseate_CR.csv") 
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)
  
prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
prep_data<-prep_data[-c(739),] #removing outlier
init<-c(0.002,0.006,0.005,0.01,1,3,3,5) #starting parameters
HmmVal4(prep_data,init) #call function for model 4
################################################################################

#COLONY - Coquet ; SPECIES - Common; MODEL0-base model #########################

data<-read.csv('Coquet_Common_CR.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.004,0.007,0.012,1,2,3,5,0.007,0.0009) #starting parameters
HmmVal0(prep_data,init) #call function for model 0
################################################################################

#COLONY - Coquet ; SPECIES - Sandwich; MODEL 3 #################################

data<-read.csv('Coquet_Sandwich_CR.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
prep_data<-prep_data[-c(6011),] #removing outlier
init<-c(0.002,0.005,0.006,0.011,1,2,3,5,0.01,0.0006) #starting parameters
HmmVal3(prep_data,init) #call function for model 3
################################################################################

#COLONY - Coquet ; SPECIES - Arctic; MODEL 1 ###################################

data<-read.csv("Coquet_Arctic_CR.csv")
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.001,0.003,0.007,0.012,1,2,3,5,0.01,0.0003) #starting parameters
HmmVal1(prep_data,init) #call function for model 1
################################################################################

#COLONY - Isle of May ; SPECIES - Arctic; MODEL 1 ###################################

data<-read.csv("IsleMay_Arctic_CR.csv")
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
prep_data<-prep_data[-c(86),] #removing outlier
init<-c(0.002,0.005,0.006,0.009,1,2,3,5,0.0006,0.00001) #starting parameters
HmmVal1(prep_data,init) #call function for model 1
################################################################################

#COLONY - Cemlyn ; SPECIES - Arctic; MODEL 4 ###################################

data<-read.csv('Cemlyn_Arctic_CR.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS_GMT,LATITUDE,LONGITUDE, DIST2COL,
                    Continuous, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.006,0.005,0.009,1,2,3,5,0.008,0.0003) #starting parameters
HmmVal4(prep_data,init) #call function for model 4
################################################################################

#COLONY - Cockle ; SPECIES - Sandwich CR; MODEL 3 ##############################

data<-read.csv('Cockle_Sandwich_CR.csv')
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
init<-c(0.001,0.003,0.007,0.012,1,2,3,5,0.047,0.0001) #starting parameters
HmmVal3(prep_data,init) #call function for model 3
################################################################################

#COLONY - Blue Circle ; SPECIES - Sandwich CR; MODEL 2 #########################
data<-read.csv('BlueCircle_Sandwich_CR.csv')
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
init<-c(0.002,0.005,0.007,0.015,1,2,3,5,0.031,0.0002) #starting parameters
HmmVal2(prep_data,init) #call function for model 2
################################################################################

#COLONY - South Shian ; SPECIES - Common CR; MODEL 4 ##############################

data<-read.csv('SouthShian_Common_CR.csv')
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.004,0.008,0.01,1,3,3,5,0.014,0.0001) #starting parameters
HmmVal4(prep_data,init) #call function for model 4
################################################################################

#COLONY - Cemlyn ; SPECIES - Common CR; MODEL 6 ##############################
data<-read.csv('Cemlyn_Common_CR.csv')
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS_GMT,LATITUDE,LONGITUDE, DIST2COL,
                    Continuous, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.001,0.003,0.008,0.011,1,2,3,5,0.006,0.0007) #starting parameters
HmmVal6(prep_data,init) #call function for model 6
################################################################################

#COLONY - Glas Eileanan ; SPECIES - Common CR; MODEL 0-base model ################

data<-read.csv('GlasEileanan_Common_CR.csv')
dnew<- data %>%
  subset(select= +c(ID,SPECIES,TSECS,LATITUDE,LONGITUDE, DIST2COL,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
prep_data<-prep_data[-c(10759),]#removing outlier
init<-c(0.001,0.003,0.008,0.01,1,2,3,5,0.009,0.0006) #starting parameters
HmmVal0(prep_data,init) #call function for model 0
################################################################################

#COLONY - Leith ; SPECIES - Common CR; MODEL 2 ##############################

data<-read.csv('Leith_Common_CR.csv')
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,TSECS,LATITUDE,LONGITUDE,DIST2COL_K,
                    CONTBEH, observed_state)) %>%
  rename(observedS = observed_state) %>%
  rename(DIST2COL = DIST2COL_K) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.001,0.003,0.007,0.011,1,2,3,5,0.014,0.0001) #starting parameters
HmmVal2(prep_data,init) #call function for model 2
################################################################################

#COLONY - Forvie ; SPECIES - Sandwich CR; MODEL 3 ##############################

data<- read.csv("Forvie_Sandwich_CR.csv")
# data preparation for HMM analysis and validation
dnew<- data %>%
  subset(select= +c(ID,BSTTSECS,SPECIES,LATITUDE,LONGITUDE,DIST2COL,
                    CONTBEH, Observed_2State)) %>%
  rename(observedS = Observed_2State) %>%
  mutate(observedS = replace(observedS, observedS == 'Travelling', 
                             'Not-foraging'))
dnew$observedS <- as.factor(dnew$observedS)

prep_data<-prepData(dnew,type='LL',coordNames=c("LONGITUDE","LATITUDE"),
                    covNames=c('DIST2COL'))
init<-c(0.002,0.004,0.008,0.012,1,2,3,5,0.015,0.0002) #starting parameters
HmmVal3(prep_data,init) #call function for model 3

