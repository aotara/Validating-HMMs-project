################################################################################
# The codes are for

# Figures 7- Page 20 and Figure 8 Pages 21
# Proportion of correctly inferred foraging states within each observed foraging
# event across study sites during incubation
#                   &
# Figure 9 - Page 21: proportion of missed observed foraging events
################################################################################

#load library
library(dplyr)
library(ggplot2)

###################        FUNCTION             ################################
# write function to calculate percentage of identified and missed foraging events
################################################################################

# getting the count of observed foraging behaviour 
FE <- function(data) {
  nLe <- data %>% 
  group_by(factor(ID)) %>% 
  # set counter 
  mutate(counter = unlist(sapply(rle(observedState)$lengths, seq_len)))%>%
  mutate(observedState =factor(observedState)) %>%
  mutate(decodedState =factor(decodedState)) %>%
  filter(observedState == 'Foraging')
  # obtain observed foraging events based on counter specified
  Cc<-do.call(rbind,by(nLe,cumsum(c(0,diff(nLe$counter)!=1)), function(g)
              data.frame(imin = min(g$counter),obsFcount = max(g$counter), 
              decFcount = 0, percentage = 0) ))
  #percentage of decoded foraging state identified within observed foraging events
  for (p in 1:length(Cc$obsFcount)) {
    ndata <-  nLe[1:Cc$obsFcount[p],] # segment data
    cdf <- ndata %>%  filter(decodedState=='Foraging') %>% nrow()
    Cc$decFcount[p] <- cdf
    Cc$percentage[p] <- (Cc$decFcount[p] / Cc$obsFcount[p])*100
  }
  uq <- unique(Cc$obsFcount[Cc$decFcount == 0])
 
  return (list( length(which(Cc$percentage > 0 & Cc$percentage <25)) ,
                length(which(Cc$percentage >= 25 & Cc$percentage <50)) ,
                length(which(Cc$percentage >= 50 & Cc$percentage <75)) ,
                length(which(Cc$percentage >= 75)) , 
                uq   
         )    )
}

#######################     FIGURE 7 Page 20             #######################
# proportion of correctly identified foraging states from observed foraging events
# DURING INCUBATION

### STEP 1 : READ DATA AND CALL FUNCTION 

# read data with observed and decoded state result obatained from HMMs deemed optimal
cockleArcticIN<- read.csv('Cockle-arctic-IN-M5-result.csv') 
cockleSandwichIN<- read.csv('Cockle-sandwich-IN-M6-result.csv') 
bluecircleRoseateIN<- read.csv('BlueCircle-roseate-IN-M4-result.csv') 
islemayArcticIN<- read.csv('IsleMay-arctic-IN-M6-result.csv') 
leithCommonIN<- read.csv('Leith-common-IN-M0-result.csv') 
# apply function to data
A<-FE(cockleArcticIN) ; B<-FE(cockleSandwichIN); C<-FE(bluecircleRoseateIN)
D<-FE(islemayArcticIN) ; E<-FE(leithCommonIN)

### STEP 2 : CREATE DATAFRAME BASED ON OUTPUT FROM CALLED FUNCTION
dh <- data.frame(proportion = c(rep("< 25%", 5), rep("25-49%", 5),
                                rep("50-74%", 5), rep("75-100%", 5)),
                 colonies= c(rep(c("Cockle","Cockle","Blue Circle","Isle of May","Leith"),4)), 
                 species= c(rep(c("Arctic" , "Sandwich","Roseate", "Arctic" ,"Common"),4)), 
                 # values- each column 1- 5 are results from lines 52 - 56
                 #           
                 value = c(A[[1]], B[[1]], C[[1]], D[[1]], E[[1]],  # <25%
                           A[[2]], B[[2]], C[[2]], D[[2]], E[[2]],  # 25-49%
                           A[[3]], B[[3]], C[[3]], D[[3]], E[[3]],  # 50-74%
                           A[[4]], B[[4]], C[[4]], D[[4]], E[[4]])) # 75-100%

### STEP 3 : PLOT FIGURE 7 page 20

dh %>%
  ggplot(aes(x = proportion , y = value, fill = species)) +
  geom_bar(position = "stack", stat = "identity") +
  # axis labels
  xlab("Percentage of correct inferred foraging state within each observed foraging event")+
  ylab("No. of foraging events")+
  # plot theme
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 15, face = "bold"),
        strip.background.x = element_rect(fill = "white", colour = "black"))+
  theme(axis.text = element_text(size = 15, colour = 'black', face = "bold"),
        axis.title = element_text(size = 15, face = "bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.position = "top", legend.direction = "horizontal")+
  scale_fill_manual(values = c("Arctic" = "#104E8B", "Common" = "#6495ED",
                               "Roseate" = "#838B8B", "Sandwich" = '#56B4E9'),
                    name = "Tern species")+
  facet_wrap(~ colonies, scales = "free_y", nrow = 1)
################################################################################

##############           Figure 8 Page 21                #######################
# Proportion of correctly inferred foraging states within each observed foraging
# event across study sites DURING CHICK-REARING

### STEP 1 : READ DATA AND CALL FUNCTION 

# read data with observed and decoded state result obatained from HMMs deemed optimal
cemlynCommonCR<-read.csv('Cemlyn-common-CR-M6-result.csv')
cemlynArcticCR<-read.csv('Cemlyn-arctic-CR-M4-result.csv')
cockleSandwichCR<-read.csv('Cockle-sandwich-CR-M6-result.csv')
coquetArcticCR<-read.csv('Coquet-arctic-CR-M1-result.csv')
coquetCommonCR<-read.csv('Coquet-common-CR-M0-result.csv')
coquetRoseateCR<-read.csv('Coquet-roseate-CR-M4-result.csv')
coquetSandwichCR<-read.csv('Coquet-sandwich-CR-M3-result.csv')
glaseileananCommonCR<-read.csv('GlasEileanan-common-CR-M0-result.csv')
leithCommonCR<-read.csv('Leith-common-CR-M2-result.csv')
southshianCommonCR<-read.csv('SouthShian-common-CR-M4-result.csv')
islemayArcticCR<-read.csv('IsleMay-arctic-CR-M1-result.csv')
bluecircleSandwichCR<-read.csv('BlueCircle-sandwich-CR-M2-result.csv')
forvieSandwichCR<-read.csv('Forvie-sandwich-CR-M3-result.csv')

#call function for obtaining foraging events and proportion of decoded
# foraging states identified within observed foraging event
a<-FE(cemlynCommonCR) ; b<-FE(cemlynArcticCR) ; c<-FE(cockleSandwichCR)
d<-FE(coquetArcticCR) ; e<-FE(coquetCommonCR) ; f<-FE(coquetRoseateCR)
g<-FE(coquetSandwichCR) ; h<-FE(glaseileananCommonCR) ; i<-FE(leithCommonCR)
j<-FE(southshianCommonCR) ; k<-FE(islemayArcticCR) 
m<-FE(bluecircleSandwichCR) ; n<-FE(forvieSandwichCR)

### STEP 2 : CREATE DATAFRAME BASED ON OUTPUT FROM CALLED FUNCTION

d <- data.frame(proportion =rep(c("< 25%","25-49%","50-74%","75-100%"),each=13),
                
                colonies= rep(c('Cemlyn','Cemlyn','Cockle',"Coquet","Coquet",
                                "Coquet","Coquet","Glas Eilean" , "Leith",
                                "South Shian","Isle of May","Blue Circle","Forvie"),4),
                
                species= rep(c("Common","Arctic","Sandwich","Arctic","Common",
                               "Roseate", "Sandwich","Common","Common","Common",
                               "Arctic" , "Sandwich","Sandwich"),4),
                
 value =c(a[[1]],b[[1]],c[[1]],d[[1]],e[[1]],f[[1]],g[[1]],h[[1]],i[[1]],j[[1]],k[[1]],m[[1]],n[[1]],# <25%
          a[[2]],b[[2]],c[[2]],d[[2]],e[[2]],f[[2]],g[[2]],h[[2]],i[[2]],j[[2]],k[[2]],m[[2]],n[[2]],# 25-49%
          a[[3]],b[[3]],c[[3]],d[[3]],e[[3]],f[[3]],g[[3]],h[[3]],i[[3]],j[[3]],k[[3]],m[[3]],n[[3]],# 50-74%
          a[[4]],b[[4]],c[[4]],d[[4]],e[[4]],f[[4]],g[[4]],h[[4]],i[[4]],j[[4]],k[[4]],m[[4]],n[[4]] # 75-100% 
          )) 

### STEP 3 : PLOT FIGURE 8 on page 21

d %>% 
  ggplot(aes(x = proportion , y = value, fill = species)) +
  geom_bar(position = "stack", stat = "identity") +
  # axis labels
  xlab("Percentage of correct inferred foraging state within each observed foraging event")+
  ylab("No. of foraging events")+
  # plot theme
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(strip.text = element_text(size = 11.5, face = "bold"),
        strip.background.x = element_rect(fill = "white", colour = "black"))+
  theme(axis.text = element_text(size = 15, colour = 'black', face = "bold"),
        axis.title = element_text(size = 15, face = "bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15, face = "bold"),
        legend.position = "top", legend.direction = "horizontal")+                     
  scale_fill_manual(values = c("Arctic" = "#6495ED", "Common" = "#104E8B",
                               "Roseate" = "#56B4E9", "Sandwich" = '#838B8B'),
                    name = "Tern species")+
  facet_wrap(~ colonies, scales = "free_y", nrow = 1)
################################################################################

##################       Figure 9 Page 21      #################################
# Missed foraging event in chick rearing and incubation
################################################################################

### STEP 1 : CREATE DATAFRAME BASED ON OUTPUT FROM CALLED FUNCTION

dataMFE<-data.frame(site =c(rep("Cemlyn (A), CR", 2), rep("Glas Eilean (C), CR",5),
                            rep("Leith (C), CR", 5), rep("Leith (C), IN", 11), 
                            "Blue Circle (S), CR",rep("Coquet (S), CR",3),
                           rep("Forvie (S), CR", 6), rep("Cemlyn (C), CR", 32)),
                    # total time for missed foraging event
                    time = c(b[[5]],h[[5]],i[[5]], E[[5]], m[[5]],g[[5]], 
                             n[[5]],a[[5]])
                    )

### STEP 2 : PLOT FIGURE 9 on page 21

dataMFE %>%
  ggplot(aes(x = time, y = site))+
  geom_point(shape = 'triangle', colour = 'orangered3', size = 2)+
  # axis labels
  ylab("Study site (tern species), season")+
  xlab("Time of missed foraging events (secs)")+
  # plot theme
  theme(panel.background = element_rect(fill = "white", colour = "black"))+
  theme(strip.text = element_text(size = 15, face = "bold"),
        strip.background.x = element_rect(fill = "white", colour = "black"))+
  theme(axis.text = element_text(size = 15, colour = 'black', face = "bold"),
        axis.title = element_text(size = 15, face = "bold")) +
  theme(legend.text = element_text(size = 15, face = "bold"),
        legend.title = element_text(size = 15,face = "bold"))
