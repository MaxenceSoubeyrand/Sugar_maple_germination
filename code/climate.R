#Script for analysis climate data
rm(list=ls())

setwd("~/postdoc/germination/github")

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(readxl)
library(viridis)

#Opening temperature data
temp_moy <- readRDS("data/temperature.rds") %>% 
  #par exemple prenons la moyenne du mois de juin
  filter(month(date) %in% c(5, 6)) %>% 
  group_by(site, year) %>% 
  summarize(
    mean_temperature = mean(temperature),
    temp_lower = mean(temperature) -  sd(temperature),
    temp_upper = mean(temperature) +  sd(temperature)) %>% 
  mutate(year=paste0("20", year))

#Number of frost events
late_frost_event <- readRDS("data/temperature.rds")%>% 
  mutate(year=paste0("20", year)) %>% 
  filter(month(date) %in% c(5, 6),
         temperature<0) %>% 
  mutate(dm=paste0(day, "_", month)) %>% 
  dplyr::select(site, dm, year) %>%
  unique() %>% 
  group_by(site, year) %>% 
  summarize(nb_late_frost_event=n())

#Join temeprature and frost events
clim <- full_join(late_frost_event, temp_moy) %>%
  ungroup() %>% 
  mutate(nb_late_frost_event = if_else(is.na(nb_late_frost_event), 0, nb_late_frost_event))

clim$site <- factor(clim$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Climate data for plotting
clim_plot <- clim %>% 
  rename(`Mean temperature (in °C)`=mean_temperature,
         `# of frost events`=nb_late_frost_event) %>% 
  pivot_longer(cols = `Mean temperature (in °C)`:`# of frost events`,
               names_to = "variable",
               values_to = "values") %>% 
  mutate(temp_lower=case_when(variable=="# of frost events" ~ NA, .default = temp_lower),
         temp_upper=case_when(variable=="# of frost events" ~ NA, .default = temp_upper))

#Rectangle for the plot
rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df$`Sugar maple range` <- factor(rect_df$`Sugar maple range`, levels = c("Within", "Outside (North)"))


#Temperature and number of late frost event figure
#Graphiques climatiques
clim_frost_plot <- ggplot(data=clim_plot, aes(x=site , y= values, ymin=temp_lower, ymax=temp_upper)) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Sugar maple range`), 
            alpha = .3, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  geom_point(size=2) +
  geom_errorbar(width = 0.3)+
  facet_grid(variable~year, scales="free_y") + 
  ylab(NULL) +
  xlab("Sites") +
  ylim(0, NA) +
  theme(strip.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom")


#Data for the late frost event temporality
late_frost_event <- readRDS("~/postdoc/germination/data/temperature/temperature.rds")%>% 
  mutate(year=paste0("20", year)) %>% 
  filter(month(date) %in% c(5, 6),
         temperature<0) %>% 
  mutate(var="Frost events")

year(late_frost_event$date) <- 2005

rect_df2 <- data.frame(ymin = c(0, 2), ymax = c(2, 6),
                       labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df2$`Sugar maple range` <- factor(rect_df2$`Sugar maple range`, levels = c("Within", "Outside (North)"))

late_frost_event$site <- factor(late_frost_event$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

frost_events <- ggplot(late_frost_event, aes(x=date, y=site)) +
  geom_rect(data = rect_df2, aes(ymin = ymin, ymax = ymax, 
                                 xmin = as.Date("2004-05-01"), xmax = as.Date("2006-06-30"), 
                                 fill = `Sugar maple range`), 
            alpha = .3, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  ylab(NULL) +
  xlab(NULL) +
  geom_point(shape=3) + 
  facet_grid(var~year, scales="free_x")+
  coord_cartesian(xlim = c(as.Date("2005-05-01"), as.Date("2005-06-30"))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom")

#Combine the plots
plot_frost_clim <- ggarrange(clim_frost_plot, frost_events,
                             common.legend = TRUE, legend="bottom",
                             nrow=2, ncol=1,
                             heights=c(1.6,1),
                             labels=c("A", "B")) + 
  bgcolor("white")+
  border("white")


ggsave(plot=plot_frost_clim, filename="figures/temp_frost.png", 
       width=8, height=8.1)
