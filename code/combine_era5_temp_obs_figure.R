#Script that compute temperature and precipitation data, and gave figure 2 and S1

rm(list=ls())

library(tidyverse)
theme_set(theme_bw())
library(readxl)
library(ggpubr)
library(corrplot)

#Climate data
temp_obs <- readRDS("~/postdoc/germination/github/data/temperature.rds")
temp_simul <- read_csv("~/postdoc/germination/github/data/temperature_era5.csv")
precip_simul <- read_csv("~/postdoc/germination/github/data/precipitation_era5.csv")

clim_simul <- left_join(temp_simul, precip_simul)

####Mean temperature####
#Mean temperature of min and max temperature by day for observed and simulated temperature
#Simulated data
clim_simul_day <- clim_simul %>%
  mutate(day = as.Date(time),
         temperature=temperature-273.15,
         precip=precip*1000) %>%   
  group_by(site, day) %>%
  summarize(temperature_simul = (max(temperature) + min(temperature)) / 2,
            precip_simul = sum(precip), 
            .groups = "drop") %>% 
  mutate(time=day,
         day = format(time, "%d"),
         month = format(time, "%m"),
         year = format(time, "%Y")) %>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  dplyr::select(site, day, month, year, temperature_simul, precip_simul)

#Observed data
temp_obs_day <- temp_obs %>% 
  group_by(site, day, month, year) %>% 
  summarize(temperature = (max(temperature) + min(temperature)) / 2)%>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  mutate(year=paste0("20", year)) 

#Join observed and simulate
join_temp <- full_join(clim_simul_day,temp_obs_day)%>% 
  mutate(temperature=case_when(temperature < -10 ~ NA,
                               .default=temperature))

#Correcting the biais
bias_model <- lm(temperature ~ temperature_simul, data = join_temp)
summary(bias_model)
summary(bias_model)$r.squared


join_temp <- join_temp %>%
  mutate(temperature_simul = predict(bias_model, newdata = join_temp))

#Fill the gap in the observed serie
mean_temp <- join_temp %>% 
  mutate(temperature=case_when(is.na(temperature) ~ temperature_simul,
                               .default=temperature)) %>% 
  dplyr::select(-temperature_simul)




####Number of late frost events####
obs_hour <- temp_obs %>%
  mutate(datetime = as.POSIXct(paste(year, month, day, hour, minute, second),
                          format="%y %m %d %H %M %S", tz="UTC")) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(site, time) %>%
  summarise(temp_obs = mean(temperature, na.rm = TRUE))

temp_hour <- temp_simul %>%
  mutate(day = as.Date(time),
         temperature_simul=temperature-273.15) %>% 
  mutate(
         day = format(time, "%d"),
         month = format(time, "%m"),
         year = format(time, "%Y")) %>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  dplyr::select(site, time, temperature_simul)

join_temp_frost <- full_join(temp_hour,obs_hour) %>% 
  mutate(temp_obs=case_when(temp_obs < -10 ~ NA,
                               .default=temp_obs)) %>%
  filter(temperature_simul  >= -10 & temperature_simul  <= 10)

bias_model <- lm(temp_obs ~ temperature_simul, data = join_temp_frost)
summary(bias_model)

join_temp_frost <- join_temp_frost %>%
  mutate(temp_simul_corrected = predict(bias_model, newdata = join_temp_frost))

join_temp_frost <- join_temp_frost %>%
  mutate(temp_final = ifelse(is.na(temp_obs), temp_simul_corrected, temp_obs))

nb_late_frost_events <- join_temp_frost %>%
  mutate(date = as.Date(time)) %>%
  group_by(site, date) %>%
  summarise(tmin = min(temp_final), .groups = "drop") %>%
  mutate(year = year(date)) %>%
  group_by(site, year) %>%
  summarise(n_frost = sum(tmin < 0), .groups = "drop")


####For the germination####
clim_germ <- mean_temp %>% 
  filter(month %in% c("05", "06")) %>% 
  group_by(site,year) %>% 
  summarize(temperature=mean(temperature), 
            precipitation=sum(precip_simul)) %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events)

saveRDS(clim_germ, "~/postdoc/germination/github/data/clim_germination.rds")


####For the growth####
clim_growth <- mean_temp %>% 
  filter(month %in% c("05", "06", "07", "08")) %>% 
  group_by(site, year) %>% 
  summarize(temperature=mean(temperature, na.rm=T),
            precipitation=sum(precip_simul, na.rm=T)) %>% 
  mutate(temperature = case_when(
    year == 2006 ~ (temperature[year == 2005] + temperature[year == 2006]) / 2,
    TRUE ~ temperature),
    precipitation = case_when(
      year == 2006 ~ (precipitation[year == 2005] + precipitation[year == 2006]) / 2,
      TRUE ~ precipitation)) %>% 
  filter(year != "2005") %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events)

saveRDS(clim_growth, "~/postdoc/germination/github/data/clim_growth.rds")


####Figure 2####
#####Temperature and Number of late frost event#####
clim_frost <- mean_temp %>% 
  group_by(site, year) %>% 
  summarize(mean_temp = mean(temperature, na.rm = TRUE),
            ecart_type = sd(temperature, na.rm = TRUE),
            temp_lower = mean_temp - ecart_type,
            temp_upper= mean_temp + ecart_type,
            sum_precip=sum(precip_simul, na.rm = TRUE)) %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events) %>% 
  dplyr::select(-ecart_type)

clim_cor <- clim_frost

clim_frost %>% 
  group_by(year) %>% 
  summarize(m=mean(n_frost))


clim_frost$site <- factor(clim_frost$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

clim_frost <- clim_frost %>% 
  rename(`Mean temperature (in °C)`=mean_temp,
         `# of frost events`=n_frost,
         `Precipitation (in mm)`=sum_precip) %>% 
  pivot_longer(cols = c(`Mean temperature (in °C)`,`# of frost events`, `Precipitation (in mm)`),
               names_to = "variable",
               values_to = "values") %>% 
  mutate(temp_lower=case_when(variable=="# of frost events" ~ NA, .default = temp_lower),
         temp_upper=case_when(variable=="# of frost events" ~ NA, .default = temp_upper),
         temp_lower=case_when(variable=="Precipitation (in mm)" ~ NA, .default = temp_lower),
         temp_upper=case_when(variable=="Precipitation (in mm)" ~ NA, .default = temp_upper))


rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("within or", "north of the current sugar maple range")) %>% 
  rename(`Site`=labs)

rect_df$`Site` <- factor(rect_df$`Site`, levels = c("within or", "north of the current sugar maple range"))



#Graphs
#Figure 2.a
temp <- ggplot(data=filter(clim_frost, variable=="Mean temperature (in °C)"), aes(x=site , y= values, ymin=temp_lower, ymax=temp_upper)) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  geom_point(size=2) +
  geom_errorbar(width = 0.3)+
  facet_wrap(~year) + 
  ylab("Temperature (in °C)") +
  xlab(NULL) +
  ylim(0, NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom")

#Figure 2.b
prec <- ggplot(data=filter(clim_frost, variable=="Precipitation (in mm)"), aes(x=site , y= values, ymin=temp_lower, ymax=temp_upper)) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  geom_point(size=2) +
  geom_errorbar(width = 0.3)+
  facet_wrap(~year) + 
  ylab("Precipitation (in mm)") +
  xlab(NULL) +
  ylim(0, NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom")

#Figure 2.c
frost <- ggplot(data=filter(clim_frost, variable=="# of frost events"), aes(x=site , y= values, ymin=temp_lower, ymax=temp_upper)) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE) +
  scale_fill_viridis_d() +
  geom_point(size=2) +
  geom_errorbar(width = 0.3)+
  facet_wrap(~year) + 
  ylab("# of frost events") +
  xlab(NULL) +
  ylim(0, NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(strip.text = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom")

#####Temporality of late frost event#####
late_frost_event <-  join_temp_frost %>%
  mutate(date = as.Date(time)) %>%
  group_by(site, date) %>%
  summarise(tmin = min(temp_final), .groups = "drop") %>%
  filter(tmin<0) %>% 
  mutate(year = year(date)) %>% 
  mutate(var="# frost events")

year(late_frost_event$date) <- 2005

rect_df2 <- data.frame(ymin = c(0, 2), ymax = c(2, 6),
                       labs = c("within or", "north of the current sugar maple range")) %>%
  rename(`Site`=labs)

rect_df2$`Site` <- factor(rect_df2$`Site`, levels = c("within or", "north of the current sugar maple range"))

late_frost_event$site <- factor(late_frost_event$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Figure 2.d
frost_events <- ggplot(late_frost_event, aes(x=date, y=site)) +
  geom_rect(data = rect_df2, aes(ymin = ymin, ymax = ymax, 
                                 xmin = as.Date("2004-05-01"), xmax = as.Date("2006-06-30"), 
                                 fill = `Site`), 
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
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom")

plot_frost_clim <- ggarrange(temp, prec, frost, frost_events,
                             common.legend = TRUE, legend="bottom",
                             nrow = 4, ncol = 1,
                             labels = c("a", "b", "c", "d"),
                             align = "v") +
  bgcolor("white") +
  border("white")

ggsave(plot=plot_frost_clim, filename="~/postdoc/germination/github/figures/temp_frost.png", 
       width=5, height=11)

site <- data.frame(site=c("MUS", "HED", "MON", "KEK", "KIP"), 
                   latitude=c(49.932, 49.243, 48.460, 48.191, 46.740),
                   longitude=c(-78.698, -78.311, -79.418, -79.112, -78.905))

clim_cor <- left_join(site, clim_cor)

clim_cor <- clim_cor[,c(2, 5, 8, 9)]
colnames(clim_cor) <- c("Latitude", "Mean temperature (in °C)", "Precipitation (in mm)", "# of frost events")

M <-  cor(clim_cor)


png(filename="~/postdoc/germination/github/figures/corr_plot.png", width=7, height=7, units="in", res=1000)
corrplot(M, method = 'number',tl.col = "black", cl.cex=1)
dev.off()



