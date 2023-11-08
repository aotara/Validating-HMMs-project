
################################################################################
# Load Libraries 
library(ggplot2)
library(tidyverse)
library(marmap)
library(ggsn)
################################################################################
# Figure 10 (a) and (b) Page 22
# Visual tracks from South Shian Colony coloured by observed & decoded states
################################################################################
data<-read.csv('SouthShian_mee_result.csv')  #read data
data[data == 'Travelling'] <- 'Not-foraging' #replace Travelling with Not-foraging
colnames(data)[6] <- 'Observed states'       #rename column name
colnames(data)[7] <- 'Decoded states'        #rename column name

# data for legend 
dd <- data.frame(lat = c(56.58, 56.58), lng = c(-5.65, -5.63))
dc <- data.frame(lat = c(56.56, 56.56), lng = c(-5.65, -5.63))

#get bathymetry
ssh <- getNOAA.bathy(lon1 = -5.7, lon2 = -5.351, lat1 = 56.35, lat2 = 56.6, 
                     resolution = 0.3)
# set breaks and colours
breaks <- seq(0, 500, by = 250)
cols <- c('white','white','grey60','grey60','grey60','grey60','grey60')

#plot
plot<-function(data,n){
  plt<-autoplot(ssh, geom=c('r')) +
    scale_fill_stepsn(colours = cols, breaks = breaks, guide = NULL) +
    geom_path(data = data, aes(x = LONGITUDE,y = LATITUDE, group = ID,
                               color = !! rlang::ensym(n)), size = 0.9)+
    scale_colour_manual(values = c(`Foraging` = 'darkorange2', 
                                   `Not-foraging` = 'steelblue1'))+
    # site location
    geom_text(aes(x = -5.57, y = 56.60, label =  n), size = 9,
              colour = 'black')+
    # legend
    geom_path(data = dd, aes(x = lng, y = lat), size = 0.9,
              colour = 'darkorange2')+
    geom_text(aes(x = -5.6, y = 56.58, label = "Foraging"), size = 7,
              colour = 'black')+
    geom_path(data = dc, aes(x = lng, y = lat), size = 0.9,
              colour = 'steelblue1')+
    geom_text(aes(x = -5.59, y = 56.56,label = "Not-foraging"), size = 7, 
              colour = 'black')+
    geom_point(aes(x = -5.425, y = 56.52), shape = '\u2605', colour = 'black', 
               size = 8) +
    geom_text(aes(x = -5.43, y = 56.505, label = "South Shian"), size = 5, 
              colour = 'black')+
    # axis labels and plot theme
    xlab("Longitude (\u00B0 W)") + ylab("Latitude (\u00B0 N)") +
    theme(axis.text = element_text(size = 24, colour = "black", face = 'bold'),
          axis.title = element_text(size = 25, face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          legend.position = "none")+
    # scalebar and north arrow
    scalebar(x.min = -5.7, x.max = -5.351, y.min = 56.35, y.max = 56.6,  dist = 2,
             dist_unit = "km", st.bottom = T, st.dist = 0.02,  st.size = 4, 
             transform = T, model = "WGS84", anchor = c(x = -5.41, y = 56.4), 
             height = 0.03, box.fill = c("black", "white","black",'white'))+
    ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
                               pad_x = unit(0.15, "in"), pad_y = unit(0.2, "in"),
             style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"),
             line_col = "grey20", text_family = "ArcherPro Book"))+ 
    coord_cartesian(expand = 0)
  return (plt)
}
plot(data,'Observed states') #Figure 10a Page 22
plot(data,'Decoded states') #Figure 10b Page 22
