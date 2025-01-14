#Script for germination experiment. 

rm(list=ls())

setwd("~/postdoc/germination/github")

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(readxl)
library(lme4)
library(lmerTest)
library(AICcmodavg)
library(viridis)
library(sjPlot)
library(MuMIn)

#Opening temperature data
temp_moy <- readRDS("data/temperature.rds") %>% 
  #par exemple prenons la moyenne du mois de juin
  filter(month(date) %in% c(5, 6)) %>% 
  group_by(site, year) %>% 
  summarize(temperature=mean(temperature)) %>% 
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

####Germination###
#Open germination data.
germ <- read_excel("data/germination.xlsx") %>% 
  mutate(year=as.factor(year),
         light=as.factor(light))

germ$site <- factor(germ$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

germ <- germ %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(site_block=paste0(site, "_", block))
  

#Random effect selection
mod_re1 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
               (1 | site/block) + (1 | year), data=germ)

mod_re2 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  (1 | site/block) , data=germ)

mod_re3 <- lm(log(germination+1)~ light + temperature + nb_late_frost_event, data=germ)

AIC_re <- data.frame(mod=1:3, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3)))

#Best model is mod_re2

mod_me1 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  light:temperature + temperature:nb_late_frost_event + nb_late_frost_event:light +
                  (1 | site/block), data=germ)

mod_me2 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  light:temperature + temperature:nb_late_frost_event + 
                  (1 | site/block), data=germ)

mod_me3 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  light:temperature +
                  (1 | site/block), data=germ)

mod_me4 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  (1 | site/block), data=germ)

mod_me5 <- lmer(log(germination+1)~ light + temperature + 
                  (1 | site/block), data=germ)

mod_me6 <- lmer(log(germination+1)~ light + 
                  (1 | site/block), data=germ)

mod_me7 <- lmer(log(germination+1)~ 1 +
                  (1 | site/block), data=germ)

liste_mod <- list(
  mod_me1 = mod_me1,
  mod_me2 = mod_me2,
  mod_me3 = mod_me3,
  mod_me4 = mod_me4,
  mod_me5 = mod_me5,
  mod_me6 = mod_me6,
  mod_me7 = mod_me7)

#AIC table
aictab(liste_mod)

#R²
r.squaredGLMM(mod_me6)

#Anova
anova(mod_me6)

#Residuals sum
sum(resid(mod_me6)^2)

#Predictions
nd <- expand.grid(site=c("KIP", "KEK","MON","HED","MUS"),
            light=c("100", "30", "10"), 
            block=c("A", "B", "C"))

pred <- predict(mod_me6, newdata=nd, se.fit = TRUE)

germ_pred <- data.frame(fit=exp(pred$fit)-1,
         se.fit=exp(pred$se.fit)-1)%>%
  bind_cols(nd)


#rectangle for plot design
rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df$`Sugar maple range` <- factor(rect_df$`Sugar maple range`, levels = c("Within", "Outside (North)"))

germ_pred$light <- factor(germ_pred$light, levels = c("10", "30", "100"), 
                          labels = c("10% light level", "30% light level", "100% light level"))

#Median germination
germ_obs <- germ %>% 
  group_by(site, light) %>% 
  summarize(fit=median(germination))

germ_obs$light <- factor(germ_obs$light, levels = c("10", "30", "100"), 
                          labels = c("10% light level", "30% light level", "100% light level"))

plot_light_site <- ggplot()+
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Sugar maple range`), 
            alpha = .3, inherit.aes = FALSE)+
  scale_fill_viridis_d() +
  geom_point(data=germ_pred, aes(x=site, y=fit)) +
  geom_point(data=germ_obs, aes(x=site, y=fit), color="red", shape = 3, size=3) +
  geom_errorbar(data=germ_pred, aes(x=site,ymin=fit-se.fit, ymax=fit+se.fit)) +
  facet_wrap(~light) +
  coord_cartesian(ylim=c(0, 37)) +
  xlab("Sites") + 
  ylab("Germination rate") + 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom")

#Figure 3
ggsave(plot=plot_light_site, filename="figures/germination_site_effect.png", 
       width=8, height=4)
