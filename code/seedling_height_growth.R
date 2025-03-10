rm(list=ls())

setwd("~/postdoc/germination/github")

library(tidyverse)
theme_set(theme_bw())
library(readxl)
library(lme4)
library(lmerTest)
library(MuMIn)
library(ggpubr)
library(AICcmodavg)

#Les données de croissance
height <- read_excel("data/height.xlsx")

head(height)

#Growth data
growth <- height %>% 
  pivot_longer(cols=c(`HT 2005`, `HT 2006`, `HT 2007`), 
               names_to="year", values_to = "height_note") %>% 
  mutate(height=ifelse(grepl("[a-zA-Z]", height_note), NA, height_note),
         height=as.numeric(height),
         year=gsub("HT ", "", year)) %>% 
  arrange(site, block, light ,seedling_ID, year) %>%
  group_by(site, seedling_ID, block, light) %>% 
  mutate(growth=height-lag(height),
         growth=case_when(year=="2006"~growth/2,
                          year=="2007"~growth)) %>% 
  na.omit()


growth$site <- factor(growth$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Climate data
clim <- readRDS("data/clim_growth.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

#Join climate and height data
growth <- growth %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(seedling_ID=paste(site, seedling_ID, light, block, sep="_")) %>% 
  mutate(site_block=paste0(site, "_", block))

#Random effect model selection
mod_re1 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                (1 | site/block) + (1 | seedling_ID) + (1 | year), data=growth)

mod_re2 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                  (1 | site/block) + (1 | seedling_ID), data=growth)

mod_re3 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                  (1 | site/block) , data=growth)

mod_re4 <- lm(log(growth+1)~ height + light + temperature + nb_late_frost_event, data=growth)

AIC_re <- data.frame(mod=1:4, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3), AIC(mod_re4)))
AIC_re

#Marginal effect model selection
mod_me1 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                  height:light + height:temperature + light:temperature +
                  (1 | site/block) + (1 | seedling_ID) + (1 | year), data=growth)

mod_me2 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                  height:light + height:temperature + 
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me3 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event + 
                  height:light + 
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me4 <- lmer(log(growth+1)~ height + light + temperature + nb_late_frost_event +
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me5 <- lmer(log(growth+1)~ height + light + temperature + 
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me6 <- lmer(log(growth+1)~ height + light + 
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me7 <- lmer(log(growth+1)~ height + 
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

mod_me8 <- lmer(log(growth+1)~ 1 +
                  (1 | site/block) + (1 | seedling_ID)+ (1 | year), data=growth)

liste_mod <- list(
  mod_me1 = mod_me1,
  mod_me2 = mod_me2,
  mod_me3 = mod_me3,
  mod_me4 = mod_me4,
  mod_me5 = mod_me5,
  mod_me6 = mod_me6,
  mod_me7 = mod_me7,
  mod_me8 = mod_me8)



aictab(liste_mod)

model_list <- list(mod_me6, mod_me5, mod_me4)

#Model averaging
model_avg <- model.avg(model_list, rank = "AICc", revised.var = TRUE)

sel_mod <- model_avg

summary(sel_mod)


#anova on the model with most effect
anova(mod_me4)

#R²
r2_models <- data.frame(
  model = c("mod_me6", "mod_me5", "mod_me4"),
  R2m = c(r.squaredGLMM(mod_me6)[1], r.squaredGLMM(mod_me5)[1], r.squaredGLMM(mod_me4)[1]),
  R2c = c(r.squaredGLMM(mod_me6)[2], r.squaredGLMM(mod_me5)[2], r.squaredGLMM(mod_me4)[2]),
  weight = c(0.53, 0.29, 0.18)  # Poids AICc de chaque modèle
)

# R² weighted mean
R2m_avg <- sum(r2_models$R2m * r2_models$weight)
R2c_avg <- sum(r2_models$R2c * r2_models$weight)

R2m_avg  #Marginal mean R²
R2c_avg  #Conditionnal mean R²

#Figures
#Temperature
new_data <- expand.grid(height=mean(growth$height), 
            nb_late_frost_event= mean(growth$nb_late_frost_event),
            light=c("SH-1", "SH-2", "SH-3"),
            temperature=seq(min(growth$temperature), max(growth$temperature), length.out = 100))

predicted_values <- predict(sel_mod, newdata = new_data, re.form = NA, se=T)

d <- bind_cols(new_data, pred=predicted_values$fit,
               pred_min=predicted_values$fit-predicted_values$se.fit, 
               pred_max=predicted_values$fit+predicted_values$se.fit) %>% 
  rename(`Light level`=light)


d$`Light level` <- case_when(
  d$`Light level` == "SH-1" ~ "100%",
  d$`Light level` == "SH-2" ~ "30%",
  d$`Light level` == "SH-3" ~ "10%")
d$`Light level` <- factor(d$`Light level`, levels = c("10%", "30%", "100%"))

temp <- ggplot(d, aes(x=temperature,
              y= exp(pred)+1,
              ymin=exp(pred_min)+1,
              ymax=exp(pred_max)+1,
              color=`Light level`, fill=`Light level`)) + 
  geom_ribbon(alpha=0.5, color=NA) +
  geom_line(linewidth=1) +
  ylab("Height growth (cm.year⁻¹)") + 
  xlab("Temperature (°C)")+ 
  ylim(0,100)+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12))


#Height
new_data <- expand.grid(height=seq(min(growth$height), max(growth$height), length.out = 100), 
                        nb_late_frost_event= mean(growth$nb_late_frost_event),
                        light=c("SH-1", "SH-2", "SH-3"),
                        temperature=mean(growth$temperature))

predicted_values <- predict(sel_mod, newdata = new_data, re.form = NA, se=T)

d <- bind_cols(new_data, pred=predicted_values$fit,
               pred_min=predicted_values$fit-predicted_values$se.fit, 
               pred_max=predicted_values$fit+predicted_values$se.fit)%>% 
  rename(`Light level`=light)

d$`Light level` <- case_when(
  d$`Light level` == "SH-1" ~ "100%",
  d$`Light level` == "SH-2" ~ "30%",
  d$`Light level` == "SH-3" ~ "10%")
d$`Light level` <- factor(d$`Light level`, levels = c("10%", "30%", "100%"))

height <- ggplot(d, aes(x=height,
              y= exp(pred)+1,
              ymin=exp(pred_min)+1,
              ymax=exp(pred_max)+1,
              color=`Light level`, fill=`Light level`)) + 
  geom_ribbon(alpha=0.5, color=NA) +
  geom_line(linewidth=1) +
  ylab("Height growth (cm.year⁻¹)") + 
  xlab("Seedling height (in cm)")+ 
  ylim(0,100)+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12))


#Number of frost day
new_data <- expand.grid(height=mean(growth$height), 
                        nb_late_frost_event= seq(min(growth$nb_late_frost_event), max(growth$nb_late_frost_event), length.out = 100),
                        light=c("SH-1", "SH-2", "SH-3"),
                        temperature=mean(growth$temperature))

predicted_values <- predict(sel_mod, newdata = new_data, re.form = NA, se=T)

d <- bind_cols(new_data, pred=predicted_values$fit,
               pred_min=predicted_values$fit-predicted_values$se.fit, 
               pred_max=predicted_values$fit+predicted_values$se.fit)%>% 
  rename(`Light level`=light)

d$`Light level` <- case_when(
  d$`Light level` == "SH-1" ~ "100%",
  d$`Light level` == "SH-2" ~ "30%",
  d$`Light level` == "SH-3" ~ "10%")
d$`Light level` <- factor(d$`Light level`, levels = c("10%", "30%", "100%"))

frost_day <- ggplot(d, aes(x=nb_late_frost_event,
              y= exp(pred)+1,
              ymin=exp(pred_min)+1,
              ymax=exp(pred_max)+1,
              color=`Light level`, fill=`Light level`)) + 
  geom_ribbon(alpha=0.5, color=NA) +
  geom_line(linewidth=1) +
  ylab("Height growth (cm.year⁻¹)") + 
  xlab("# of frost day")+ 
  ylim(0,100)+ 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=9),
        axis.title=element_text(size=10),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12))

plot <- ggarrange(height, temp, frost_day, 
          common.legend = TRUE, legend="bottom", 
          labels="AUTO", nrow=1) + 
  bgcolor("white")+
  border("white")

ggsave(plot=plot, filename="figures/growth_height_effect.png", 
       width=8, height=3)



#Height growth differrences between sites
growth$fitted <-  unname(exp(predict(model_avg, newdata = growth, type = "response")) - 1)

site_effect <- growth %>% 
  group_by(site, light) %>% 
  summarize(mean = mean(fitted),
            sd = sd(fitted)) %>% 
  rename(`Light level` = light)

site_effect$site <- factor(site_effect$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))
site_effect$`Light level` <- case_when(
  site_effect$`Light level` == "SH-1" ~ "100% light level",
  site_effect$`Light level` == "SH-2" ~ "30% light level",
  site_effect$`Light level` == "SH-3" ~ "10% light level")
site_effect$`Light level` <- factor(site_effect$`Light level`, levels = c("10% light level", "30% light level", "100% light level"))

rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("Within", "Outside (North)")) %>% 
  rename(`Sugar maple range`=labs)

rect_df$`Sugar maple range` <- factor(rect_df$`Sugar maple range`, levels = c("Within", "Outside (North)"))

growth_obs <- growth %>% 
  group_by(site, light) %>% 
  summarize(med_growth=median(growth)) %>% 
  rename(`Light level`=light)

growth_obs$site <- factor(growth_obs$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))
growth_obs$`Light level` <- case_when(
  growth_obs$`Light level` == "SH-1" ~ "100% light level",
  growth_obs$`Light level` == "SH-2" ~ "30% light level",
  growth_obs$`Light level` == "SH-3" ~ "10% light level")
growth_obs$`Light level` <- factor(growth_obs$`Light level`, levels = c("10% light level", "30% light level", "100% light level"))

site_plot <- ggplot()+
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Sugar maple range`), 
            alpha = .3, inherit.aes = FALSE)+
  scale_fill_viridis_d() +
  geom_point(data=site_effect, aes(x=site, y= mean)) +
  geom_errorbar(data=site_effect, aes(x=site, ymin=mean-sd, ymax=mean+sd)) +
  geom_point(data=growth_obs, aes(x=site, y= med_growth), color="red", shape=3, size=3) +
  facet_wrap(~`Light level`) +
  ylab("Height growth (cm.year⁻¹)") + 
  xlab("Sites") + 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.position = "bottom")

site_plot

ggsave(plot=site_plot, filename="figures/growth_site_effect.png", 
       width=8, height=4)

