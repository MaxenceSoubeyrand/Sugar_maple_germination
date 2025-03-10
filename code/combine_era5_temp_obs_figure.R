# Script to combine observed temperature data and era5 data 
# to fill the missing data

rm(list=ls())

setwd("~/postdoc/germination/github")

library(tidyverse)
theme_set(theme_bw())
library(readxl)
library(ggpubr)

#Oberserved and era5 climate data
temp_obs <- readRDS("data/temperature.rds")
era5 <- read_csv("~/postdoc/germination/data/era5/temperature.csv")

str(temp_obs)
str(era5)

####Temperature####
#Mean between daily minimal and maximal temperature  
#era5 temperature
era5_day <- era5 %>%
  mutate(day = as.Date(time),
         temperature=temperature-273.15) %>%   
  group_by(site, day) %>%
  summarize(temperature_era5 = (max(temperature) + min(temperature)) / 2, 
            .groups = "drop") %>% 
  mutate(time=day,
         day = format(time, "%d"),
         month = format(time, "%m"),
         year = format(time, "%Y")) %>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  dplyr::select(site, day, month, year, temperature_era5)

#Observed temeprature
temp_obs_day <- temp_obs %>% 
  group_by(site, day, month, year) %>% 
  summarize(temperature = (max(temperature) + min(temperature)) / 2)%>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  mutate(year=paste0("20", year)) 

#Joining observed and era5 temperature
join_temp <- full_join(era5_day,temp_obs_day)%>% 
  mutate(temperature=case_when(temperature < -10 ~ NA,
                               .default=temperature))

#Correcting potential biais that correspond to a mismatch between obserced and era5
bias_model <- lm(temperature ~ temperature_era5, data = join_temp)
summary(bias_model)


join_temp <- join_temp %>%
  mutate(temperature_era5 = predict(bias_model, newdata = join_temp))

mean_temp <- join_temp %>% 
  mutate(temperature=case_when(is.na(temperature) ~ temperature_era5,
                               .default=temperature)) %>% 
  dplyr::select(-temperature_era5)


####Number of late forst event####
#observed
obs_hour <- temp_obs %>%
  mutate(datetime = as.POSIXct(paste(year, month, day, hour, minute, second),
                          format="%y %m %d %H %M %S", tz="UTC")) %>%
  mutate(time = floor_date(datetime, "hour")) %>%
  group_by(site, time) %>%
  summarise(temp_obs = mean(temperature, na.rm = TRUE))

#era5
era5_hour <- era5 %>%
  mutate(day = as.Date(time),
         temperature_era5=temperature-273.15) %>% 
  mutate(
         day = format(time, "%d"),
         month = format(time, "%m"),
         year = format(time, "%Y")) %>%
  filter(month %in% c("05", "06", "07", "08")) %>% 
  dplyr::select(site, time, temperature_era5)

#Joining
join_temp_frost <- full_join(era5_hour,obs_hour) %>% 
  mutate(temp_obs=case_when(temp_obs < -10 ~ NA,
                               .default=temp_obs)) %>%
  filter(temperature_era5  >= -10 & temperature_era5  <= 10)

#Correcting potential biais in temperature
bias_model <- lm(temp_obs ~ temperature_era5, data = join_temp_frost)
summary(bias_model)

join_temp_frost <- join_temp_frost %>%
  mutate(temp_era5_corrected = predict(bias_model, newdata = join_temp_frost))

join_temp_frost <- join_temp_frost %>%
  mutate(temp_final = ifelse(is.na(temp_obs), temp_era5_corrected, temp_obs))

nb_late_frost_events <- join_temp_frost %>%
  mutate(date = as.Date(time)) %>%
  group_by(site, date) %>%
  summarise(tmin = min(temp_final), .groups = "drop") %>%
  mutate(year = year(date)) %>%
  group_by(site, year) %>%
  summarise(n_frost = sum(tmin < 0), .groups = "drop")


####For germination####
#different from height growth because different temporality
clim_germ <- mean_temp %>% 
  filter(month %in% c("05", "06")) %>% 
  group_by(site,year) %>% 
  summarize(temperature=mean(temperature)) %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events)

saveRDS(clim_germ, "data/clim_germination.rds")


####For height growth####
clim_growth <- mean_temp %>% 
  filter(month %in% c("05", "06", "07", "08")) %>% 
  group_by(site, year) %>% 
  summarize(temperature=mean(temperature, na.rm=T)) %>% 
  mutate(temperature = case_when(
    year == 2006 ~ (temperature[year == 2005] + temperature[year == 2006]) / 2,
    TRUE ~ temperature)) %>% 
  filter(year != "2005") %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events)

saveRDS(clim_growth, "data/clim_growth.rds")

####Figure####

#####Temperature and Number of late frost event#####
temp_frost <- mean_temp %>% 
  group_by(site, year) %>% 
  summarize(mean_temp = mean(temperature, na.rm = TRUE),
            ecart_type = sd(temperature, na.rm = TRUE),
            temp_lower = mean_temp - ecart_type,
            temp_upper= mean_temp + ecart_type) %>% 
  mutate(year=as.numeric(year)) %>% 
  left_join(nb_late_frost_events) %>% 
  dplyr::select(-ecart_type)

temp_frost %>% 
  group_by(year) %>% 
  summarize(m=mean(n_frost))


temp_frost$site <- factor(temp_frost$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

temp_frost <- temp_frost %>% 
  rename(`Mean temperature (in °C)`=mean_temp,
         `# of frost events`=n_frost) %>% 
  pivot_longer(cols = c(`Mean temperature (in °C)`,`# of frost events`),
               names_to = "variable",
               values_to = "values") %>% 
  mutate(temp_lower=case_when(variable=="# of frost events" ~ NA, .default = temp_lower),
         temp_upper=case_when(variable=="# of frost events" ~ NA, .default = temp_upper))


rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df$`Sugar maple range` <- factor(rect_df$`Sugar maple range`, levels = c("Within", "Outside (North)"))



#Figure
clim_frost_plot <- ggplot(data=temp_frost, aes(x=site , y= values, ymin=temp_lower, ymax=temp_upper)) +
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

#####Temporality of late frost event#####

late_frost_event <-  join_temp_frost %>%
  mutate(date = as.Date(time)) %>%
  group_by(site, date) %>%
  summarise(tmin = min(temp_final), .groups = "drop") %>%
  filter(tmin<0) %>% 
  mutate(year = year(date)) %>% 
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

plot_frost_clim <- ggarrange(clim_frost_plot, frost_events,
                             common.legend = TRUE, legend="bottom",
                             nrow=2, ncol=1,
                             heights=c(1.6,1),
                             labels=c("A", "B"),
                             align="v") + 
  bgcolor("white")+
  border("white")

ggsave(plot=plot_frost_clim, filename="figures/temp_frost.png", 
       width=8, height=8.1)
