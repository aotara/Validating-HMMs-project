################################################################################
# Load Libraries 
library(ggplot2)
library(tidyverse)
library(marmap)
library(ggsn)
library(dplyr)

################################################################################
# Figure 1 Page 7- Map for study sites
################################################################################

# creating dataframe longitude,latitude of colonies
df1 <- data.frame(lat = c(54.8197, 53.41001, 54.6, 55.38, 57.307165438,
                          56.49629974, 56.3, 56, 56.20), 
                  lng = c(-5.7, -4.51393, -5.5, -1.53850207, -1.98499606,
                          -5.71587596, -2.566667, -3.167832662, -5.6),
                  colony = c('Blue Circle', 'Cemlyn Bay', 'Cockle Island',
                             'Coquet Island', 'Forvie', 'Glas Eileanan',
                             'Isle of May', 'Leith', 'South Shian'),
                  num=c(1,2,3,4,5,6,7,8,9), stringsAsFactors = FALSE)
# legend
df2 <- data.frame(lat = c(59.2,58.8,58.3,57.8,57.3,56.8,56.3,55.8,55.3,54.8), 
                  lng = c(0.8,0.63,0.7,0.8,1.1,0.3,0.8,0.6,0.2,0.68),
                  colony = c('Study site (tern species)','1. Blue Circle (R,S)',
                             '2. Cemlyn Bay (A,C)','3. Cockle Island (A,S)',
                             '4. Coquet Island (A,C,R,S)','5. Forvie (S)',
                             '6. Glas Eileanan (C)','7. Isle of May (A)',
                             '8. Leith (C)','9. South Shian (C)'),
                  stringsAsFactors = FALSE)
# get bathymetry
colony <- getNOAA.bathy(lon1 = -12, lon2 = 3,lat1 = 49.5, lat2 = 59.5, 
                        resolution = 5)
# set breaks and colours
breaks <- seq(-2000, 0, by = 2000)
cols   <- c('white', 'white', 'white', 'grey60', 'grey60')

autoplot(colony, geom = c('r')) +
  scale_fill_stepsn(colours = cols, breaks = breaks, guide = NULL)
# plot
last_plot()+
  geom_text(data = df1, aes(x = lng, y = lat, label = num), size = 7,
            colour = 'blue2',fontface = "bold")+
  geom_text(data = df2,aes(x = lng, y = lat, label = colony), size = 5,
            colour = 'darkblue',fontface = "bold")+
  # axis labels and plot theme
  xlab(expression(paste("Longitude (\u00B0 W)"))) +
  ylab(expression(paste("Latitude (\u00B0 N)"))) +
  theme(axis.text = element_text(size = 20, colour = "black", face = "bold"),
        axis.title = element_text(size = 25,colour = 'black', face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, size = 2))+
  # scalebar and north arrow
  scalebar(x.min = -11, x.max = -7, y.min = 50, y.max = 50.3, st.bottom = T,
           dist = 50, dist_unit = "km", transform = T, model = "WGS84", 
           anchor = c(x = -10, y = 50.3), height = 0.7, st.size = 5, st.dist = 1,
           box.fill = c("black", "white", "black", "white"))+
  ggspatial::annotation_north_arrow( location = "tl", which_north = "true",
                                     pad_x = unit(0.4, "in"), pad_y = unit(0.4, "in"),
                                     style = ggspatial::north_arrow_nautical(fill = c("grey40", "white"),
                                                                             line_col = "grey20", text_family = "ArcherPro Book")+
                                       coord_cartesian(expand = 0)
  )
