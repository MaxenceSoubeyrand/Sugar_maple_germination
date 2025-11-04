rm(list=ls())

library(tidyverse)
library(ggthemes)
library(mapview)
library(viridis)
library(sf)
library(maps)
library(ggspatial)
library(ggnewscale)
library(cowplot)
library(ggpubr)
library(grid)
library(ggplotify)
library(raster)

library(RODBC)

#Site localisation
site <- data.frame(site=c("MUS", "HED", "MON", "KEK", "KIP"), 
                   latitude=c(49.932, 49.243, 48.460, 48.191, 46.740),
                   longitude=c(-78.698, -78.311, -79.418, -79.112, -78.905)) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326)

PdB <- data.frame(site=c("PdB"), 
                   latitude=c(46.032984),
                   longitude=c(-73.184487)) %>% 
  sf::st_as_sf(coords = c("longitude", "latitude"),
               crs = 4326)


#Opening bioclimatic domain
dom_bio <- st_read(dsn = "~/postdoc/germination/github/data/dom_bio/dom_bio.shp",
                   layer = "dom_bio")
dom_bio_map <- subset(dom_bio, dom_bio$DOM_BIO %in% c("3","4", "5", "6"))

lim_bio <- st_bbox(site)

canada <- sf::st_as_sf(map('world', regions="canada", plot = FALSE, fill = TRUE))
world <- sf::st_as_sf(map('world', plot = FALSE, fill = TRUE))


quebec <- st_read(dsn = "~/postdoc/germination/github/data/quebec/quebec.shp",
                  layer = "quebec")

sugar_maple <- st_read(dsn = "~/postdoc/germination/github/data/limite_ERS/ERS_northern_boundary.shp",
                  layer = "ERS_northern_boundary")

draw_key_custom <- function(data, params, size) {
  gp = grid::gpar(
    col = "red",
    fill = "white",
    lwd = 3,
    lty = 1
  )
  grid::grobTree(
    grid::rectGrob(width = 1, height = 1),
    grid::linesGrob(x = c(0.2,0.8), y = 0.5),
    gp = gp
  )
}


main <- ggplot()+ 
  geom_sf(data = canada) +
  geom_sf(data = world) +
  geom_sf(data=dom_bio_map, aes(fill=DOM_BIO)) +
  geom_sf(data = sugar_maple, aes(color=as.character(id)), 
          fill=NA,
          linetype="solid", linewidth=3) +
  geom_sf(data=site, size=5, alpha=1, shape=15) +
  geom_sf(data=PdB, size=5, alpha=1, shape=24, col="red", fill="red") +
  geom_sf_label(data=site, aes(label = site), size=5, nudge_x = 0.6) +
  geom_sf_label(data=PdB, aes(label = site), size=5, nudge_x = 0.6, shape= 24, col="red") +
  scale_fill_manual(
    values = c("#FB65FF", "#C2BF84", "#FFCD5A", "#81ED66", "#63CAE9", "#B3CCFF"),
    breaks = c("1", "2", "3","4", "5", "6"),
    labels = c("Hickory -\nmaple", "Basswood -\nsugar maple", 
               "Yellow birch -\nsugar maple ","Balsam fir -\nyellow birch", 
               "Balsam fir -\nwhite birch", "Spruce -\nmoss")) +
  scale_color_manual(values="red",
                     breaks = c("1"),
                     labels = c("Sugar maple northern\ncontinuous limit"),
                     name="") +
  labs(fill="Bioclimatic domain") +
  coord_sf(xlim = c(lim_bio[1]-0.3, -72), ylim = c(46, lim_bio[4]+0.6), expand=T) +
  theme(panel.background = element_rect(fill = "dodgerblue3"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 15),
        legend.position="bottom",
        legend.box="vertical",
        strip.text = element_text(size = 12),
        plot.title = element_text(size = 16)
  ) + 
  guides(fill=guide_legend(nrow=2,byrow=TRUE, order = 1),
         colour = guide_legend(order = 2)) +
  annotation_north_arrow(location = "bl", which_north = "true", 
                         pad_x = unit(0.05, "in"), pad_y = unit(0.15, "in"),
                         style = north_arrow_fancy_orienteering) +
  annotation_scale(location = "bl", width_hint = 0.35, height = unit(0.1, "cm")) 


countries <- data.frame(lat=c(54, 44.0),
                        long=c(-73, -73.0),
                        country=c("Quebec", "USA"))

emprise <- ggplot(data=world) +
  geom_sf() +
  geom_sf(data = quebec) +
  geom_text(data=countries, aes(x=long, y=lat, label=country), size=4)+
  #???coord_sf(xlim = c(lim_bio[1,1], lim_bio[1,2]+lim_bio[1,2]/7), ylim = c(lim_bio[2,1], lim_bio[2,2]), expand = FALSE) +
  coord_sf(xlim = c(-82, -60), ylim = c(40, 60), expand = FALSE) +
  theme_void() +
  annotate("rect", xmin = lim_bio[1]-0.3, xmax = lim_bio[3]+5, ymin = lim_bio[2]-0.1, ymax = lim_bio[4]+0.6,
           alpha = .6, fill = "grey50", color="black", linewidth=0.5) +
  theme(
    panel.background = element_rect(fill = "dodgerblue3"),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black",
                                fill = NA,
                                size = 1))

map <-
  ggdraw() +
  draw_plot(main) +
  draw_plot(emprise, x=0.588, y=0.68, width=.42, height=.3)


png(filename="~/postdoc/germination/github/figures/map.png",
    unit="in", height=7, width=7, res=500)
map
dev.off()

