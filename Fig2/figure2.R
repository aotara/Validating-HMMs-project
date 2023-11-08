

################################################################################
#    Figure 2 Page 8 -  Map for Leith colony
################################################################################

# Load Libraries 
library(ggplot2)
library(tidyverse)
library(marmap)
library(ggsn)
library(dplyr)
# data 
Dcr<-read.csv('Leith_CR_2010.csv') # chick rearing
Din<-read.csv('Leith_IN_2010.csv') # incubation

# replacing 'travelling' with ,not-foraging'
Dcr<-Dcr%>%
  mutate(observed_state=replace(observed_state,observed_state=='Travelling',
                                'Not-foraging'))
Din<-Din%>%
  mutate(observed_state=replace(observed_state,observed_state=='Travelling',
                                'Not-foraging'))
# function to plot map
leithPlot<- function(data,dd,dc){
  # get bathymetry
  leith <- getNOAA.bathy(lon1 = -3.33, lon2 = -2.95, lat1 = 55.945, lat2 =56.155, 
                         resolution = 0.1)
  # set breaks and colours
  breaks <- seq(0, 300, by = 150)
  cols <- c('white','white','grey60','grey60','grey60','grey60','grey60')
  
  autoplot(leith, geom = c('r')) + 
    scale_fill_stepsn(colours = cols, guide = NULL, breaks =breaks)
  # plot
  k <- last_plot()+
    geom_path(data = data, aes(x = LONGITUDE, y = LATITUDE, group = ID,
                               color = `observed_state`), size = 0.7) +
    scale_colour_manual(values = c(Foraging = 'darkorange2',
                                   `Not-foraging` = 'steelblue1'))+
    # site location
    geom_point(aes(x = -3.179, y = 55.99), shape = '\u2605', colour = 'black',
               size = 6) +
    geom_text(aes(x = -3.18, y = 55.975, label = "Leith"), size = 8,
              colour = 'black')+
    # legend
    geom_text(aes(x = -3.22, y = 56.185, label = "Observed State"), size = 8,
              colour = 'black')+
    geom_path(data = dd, aes(x = lng, y = lat), size = 1, colour = 'darkorange2')+
    geom_text(aes(x = -3.25, y = 56.165, label = "Foraging"), size = 7, 
              colour = 'black')+
    geom_path(data = dc, aes(x = lng, y = lat), size = 1,colour = 'steelblue1')+
    geom_text(aes(x = -3.235, y = 56.145, label = "Not-foraging"), size = 7, 
              colour = 'black')+
    # axis labels and plot theme
    xlab("Longitude (\u00B0 W)") + ylab("Latitude (\u00B0 N)") +
    theme(axis.text = element_text(size = 22, colour = "black", face = "bold"),
          axis.title = element_text(size = 25, colour = "black", face = "bold"),
          panel.border = element_rect(colour = "black", fill = NA, size = 2),
          legend.position = 'none')+
    # scalebar and north arrow
    scalebar(x.min = -3.33, x.max = -2.95, y.min = 55.945, y.max = 56.155, 
             st.bottom = T, dist = 5, dist_unit = "km", transform = T,
             model = "WGS84", anchor = c(x = -3.25, y = 55.92), height = 0.03, 
             st.size = 5, st.dist = 0.05, box.fill = c("black", "white",
                                                        "black", "white"))+
    ggspatial::annotation_north_arrow(location = "tl", which_north = "true",
                           pad_x = unit(0.1, "in"), pad_y = unit(0.2, "in"),
          style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"), 
                      line_col = "grey20", text_family = "ArcherPro Book"))+ 
    coord_cartesian(expand = 0)
  return(k)
}
# plot figure 2a
# set data points for legend - incubation
dd <- data.frame(lat = c(56.165, 56.165), lng = c(-3.36, -3.33)) # incubation
dc <- data.frame(lat = c(56.15, 56.15), lng = c(-3.36, -3.33)) # incubation
# plot 
leithPlot(Din,dd,dc) 

# plot figure 2b
# set data points for legend - chick-rearing
dd <- data.frame(lat = c(56.165, 56.165), lng = c(-3.33, -3.30)) # chick-rearing
dc <- data.frame(lat = c(56.145, 56.145), lng = c(-3.33, -3.30)) # chick-rearing
# plot 
leithPlot(Dcr,dd,dc)
