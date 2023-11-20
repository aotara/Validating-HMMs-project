################################################################################
#     CODE FOR FIGURE 4 ON PAGE 15
################################################################################
# load libraries
library(NISTunits)
library(momentuHMM)
library(dplyr)

## SET WORKING DIRECTORY
setwd(" ")

#####################   STEP 1     #############################################

# read behavioural data with bearing, distance inclusive
dbeh<-read_excel('Fig4/coquet-2009-behaviour-locs.xlsx')
db<- dbeh %>% 
  arrange(Track_ID,TSecs_GMT) %>%
  filter(!is.na(Bearing) , !is.na(Distance_From_Boat)) %>%
  rename(ID = Track_ID) %>%
  mutate(Lon_T = 0) %>%
  mutate(Lat_T = 0)

#####################     STEP 2          ######################################

# calculate approximate tern location (longitude,latitude)
for(i in 1:length(db$ID)) {
  db$Lat_T[i] <- NISTradianTOdeg( asin( 
     sin(NISTdegTOradian(db$Latitude[i]))* cos(db$Distance_From_Boat[i]/6371000)
     + cos(NISTdegTOradian(db$Latitude[i]))* sin(db$Distance_From_Boat[i]/6371000)*cos(NISTdegTOradian(db$Bearing[i]))   )   )
  y <- sin(NISTdegTOradian(db$Bearing[i])) * sin(db$Distance_From_Boat[i]/6371000)* cos(NISTdegTOradian(db$Latitude[i]))
  x <- cos(db$Distance_From_Boat[i]/6371000) - sin(NISTdegTOradian(db$Latitude[i]))*sin(NISTdegTOradian(db$Lat_T[i]))
  db$Lon_T[i] <- db$Longitude[i] + atan2(y,x)
}

######################  STEP 3:  FIGURE 4 Page 15 ##############################
# ID 84
# filter data based on ID
id84<- db%>% filter(ID==84)
write.csv(id84, file ='id84-boat-approximate-tern-location.csv',row.names=F) 
id84<-read.csv('id84-boat-approximate-tern-location.csv')
# subset with boat latitude and longitude columns
d1a<- subset(id84, select = -c(Lon_T,Lat_T) )
# subset with approximate tern latitude and longitude columns 
d1b<- subset(id84, select = -c(Longitude,Latitude) )
# rename column of tern position data
names(d1b)[names(d1b) == 'Lat_T'] <- 'Latitude'
names(d1b)[names(d1b) == 'Lon_T'] <- 'Longitude'
# calculate step length and turning angle for boat and tern tracks
p1a<-prepData(d1a,type='LL',coordNames=c("Longitude","Latitude")) #boat
p1b<-prepData(d1b,type='LL',coordNames=c("Longitude","Latitude")) #tern

# set plot layout 
mat <- matrix(c(1, 2,3,  # First, second, third
                1, 4,5), # first and fourth, fifth plot
              nrow = 2, ncol = 3, byrow = T)
layout(mat = mat,widths = c(2.5,2, 2)) # 1st, 2nd and 3rd column relative widths

# plot figure 4
# plot boat and approximate tern track
plot(id84$Longitude,id84$Latitude, type="o", pch="o",lty=1,col='darkblue',
     main='ID 84', cex=1.5, lwd=1.8,xlab='Longitude', ylab='Latitude',
     cex.lab=1.8, cex.axis=1.7)                                    #boat
points(id84$Lon_T,id84$Lat_T, col="steelblue1",pch="+",cex=1.7)    #tern
lines(id84$Lon_T,id84$Lat_T, col='steelblue1',lty=2, lwd=2)
legend("bottomright",legend=c("boat track","approximate tern track"), 
       col=c("darkblue","steelblue1"),
       pch=c("o","+"),lty=c(1,2), ncol=1,lwd=1.5,cex=1.5)
# plot histogram of step length for boat and approximate tern track
hist(p1a$step, col = "grey60", border = "black", las = 1, main='',
    xlab="Step length (boat)",prob=T,cex.lab=1.2, cex.axis=1.2) #boat
hist(p1b$step, col = "grey60", border = "black", las = 1, main='',
    xlab="Step length (tern)",prob=T,cex.lab=1.2, cex.axis=1.2) #tern
# plot histogram of turning angles for boat and approximate tern track
hist(p1a$angle, col = "grey60", border = "black", las = 1, main='',
    xlab="Turning angle (boat)",prob=T,cex.lab=1.2, cex.axis=1.2) #boat
hist(p1b$angle, col = "grey60", border = "black", las = 1, main='',
     xlab="Turning angle (tern)",prob=T,cex.lab=1.2, cex.axis=1.2) #tern

