# map of SFREC location

library(tidyverse)
library(maps)
library(mapdata)
library(ggmap)


get_map(location = "california", zoom = 6, source = "stamen")
qmplot()

maps::map("states", boundary=FALSE, col="gray", add=TRUE)

maps::map("state", boundary=FALSE, col="gray", add=TRUE)

usa <- map_data('usa')
state <- map_data("state")


ggplot(data=subset(state, grepl("Calif|Oregon", region)), aes(x=long, y=lat, fill=region, group=group)) + 
  geom_polygon(color = "white") + 
  #guides(fill=FALSE)
  theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) + 
  ggtitle('U.S. Map with States') + 
  coord_fixed(1.3)
