#Analysis of germination rate in relation to climate

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

#Germination rate data
germ <- read_excel("data/germination.xlsx") %>% 
  mutate(year=as.factor(year),
         light=as.factor(light))
head(germ)

germ$site <- factor(germ$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#temperature data
clim <- readRDS("data/clim_germination.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

####Germination###
germ <- germ %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(site_block=paste0(site, "_", block))

#germaintion rate mean by ligh condition
germ %>% 
  group_by(light) %>% 
  summarize(mean=mean(germination))

#random effect selection
mod_re1 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
               (1 | site/block) + (1 | year), data=germ)

mod_re2 <- lmer(log(germination+1)~ light + temperature + nb_late_frost_event + 
                  (1 | site/block) , data=germ)

mod_re3 <- lm(log(germination+1)~ light + temperature + nb_late_frost_event, data=germ)

AIC_re <- data.frame(mod=1:3, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3)))
AIC_re

#Marginal effect selection

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


aictab(liste_mod)

model_list <- list(mod_me6, mod_me5)

#Model averaging
model_avg <- model.avg(model_list, rank = "AICc", revised.var = TRUE)

sel_mod <- model_avg

summary(sel_mod)

#anova on the model with most effect
anova(mod_me5)

#R²
r2_models <- data.frame(
  model = c("mod_me6", "mod_me5"),
  R2m = c(r.squaredGLMM(mod_me6)[1], r.squaredGLMM(mod_me5)[1]),
  R2c = c(r.squaredGLMM(mod_me6)[2], r.squaredGLMM(mod_me5)[2]),
  weight = c(0.91, 0.09)  # Poids AICc de chaque modèle
)

# Weighted mean
R2m_avg <- sum(r2_models$R2m * r2_models$weight)
R2c_avg <- sum(r2_models$R2c * r2_models$weight)

R2m_avg  # Mean marginal R² 
R2c_avg  # Mean conditionnal R²

nd <- expand.grid(site=c("KIP", "KEK","MON","HED","MUS"),
            light=c("100", "30", "10"), 
            block=c("A", "B", "C")) %>% 
  left_join(germ %>% 
              group_by(site) %>% 
              summarise(temperature=mean(temperature)) %>% 
              select(site,temperature))

pred <- predict(sel_mod, newdata=nd, se.fit = TRUE)

germ_pred <- data.frame(fit=exp(pred$fit)-1,
         se.fit=exp(pred$se.fit)-1)%>%
  bind_cols(nd)

rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df$`Sugar maple range` <- factor(rect_df$`Sugar maple range`, levels = c("Within", "Outside (North)"))

germ_pred$light <- factor(germ_pred$light, levels = c("10", "30", "100"), 
                          labels = c("10% light level", "30% light level", "100% light level"))

germ_pred$site <- factor(germ_pred$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

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
  
ggsave(plot=plot_light_site, filename="figures/germination_site_effect.png", 
       width=8, height=4)
