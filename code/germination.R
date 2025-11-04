#Script taht model the germination rate of sugar maple. 

rm(list=ls())

Sys.setlocale("LC_TIME", "en_US.UTF-8")

library(tidyverse)
theme_set(theme_bw())
library(ggpubr)
library(readxl)
library(lme4)
library(lmerTest)
library(car)
library(AICcmodavg)
library(viridis)
library(sjPlot)
library(MuMIn)
library(ggeffects)

germ <- read_excel("~/postdoc/germination/github/data/germination.xlsx") %>% 
  mutate(year=as.factor(year),
         light=as.factor(light))
head(germ)

germ$site <- factor(germ$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Climate
clim <- readRDS("~/postdoc/germination/github/data/clim_germination.rds") %>% 
  mutate(year=as.character(year)) %>% 
  rename(nb_late_frost_event=n_frost)

#Join climate and germination
germ <- germ %>% 
  left_join(clim, by = join_by(site, year)) %>% 
  na.omit() %>% 
  mutate(site_block=paste0(site, "_", block))

#Mean germiantion rate accroding to light
germ %>%
  group_by(light) %>%
  summarize(
    mean = mean(germination, na.rm = TRUE),
    sd   = sd(germination, na.rm = TRUE),
    n    = n(),
    se   = sd / sqrt(n),  # erreur standard
    ci_lower = mean - qt(0.975, df = n - 1) * se,
    ci_upper = mean + qt(0.975, df = n - 1) * se
  )

germ$site <- relevel(as.factor(germ$site), ref = "KIP")

germ$light <- factor(germ$light, levels = c("10", "30", "100"), 
                          labels = c("Dense-canopy", "Intermediate", "Open-canopy"))


#Random effect selection model
mod_re1 <- glmer(cbind(germination,100-germination) ~ light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) + 
                   (1 | site/block) + (1 | year), family = binomial, data=germ)

mod_re2 <- glmer(cbind(germination,100-germination)~ light  + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                   (1 | site/block), family = binomial, data=germ)

mod_re3 <- glm(cbind(germination,100-germination)~ light + scale(temperature) + 
                 scale(nb_late_frost_event) + scale(precipitation), family = binomial, data=germ)

AIC_re <- data.frame(mod=1:3, AIC=c(AIC(mod_re1), AIC(mod_re2), AIC(mod_re3)))
AIC_re

#Fixed effect selection
mod_me1 <- glmer(cbind(germination,100-germination)~ light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                  light:scale(temperature) + light:scale(precipitation) + light:scale(nb_late_frost_event)+
                 (1 | site/block), 
                 family = binomial, data=germ)

mod_me2 <- glmer(cbind(germination,100-germination)~ light +scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                  light:scale(temperature) + light:scale(precipitation) + 
                  (1 | site/block), family = binomial, data=germ)

mod_me3 <- glmer(cbind(germination,100-germination)~ light + scale(temperature)+ 
                  scale(nb_late_frost_event) + scale(precipitation)+
                  light:scale(temperature) +
                  (1 | site/block), family = binomial, data=germ)

mod_me4 <- glmer(cbind(germination,100-germination)~ light + scale(temperature) + 
                   scale(nb_late_frost_event) + scale(precipitation) +
                  (1 | site/block), family = binomial, data=germ)

mod_me5 <- glmer(cbind(germination,100-germination)~ light + scale(temperature) + 
                   scale(precipitation) +
                  (1 | site/block), family = binomial, data=germ)

mod_me6 <- glmer(cbind(germination,100-germination)~ light + scale(temperature) + 
                  (1 | site/block), family = binomial, data=germ)

mod_me7 <- glmer(cbind(germination,100-germination)~ light +
                  (1 | site/block), family = binomial, data=germ)

mod_me8 <- glmer(cbind(germination,100-germination)~ 1 +
                   (1 | site/block), family = binomial, data=germ)

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

sel_mod <- mod_me1

summary(sel_mod)
car::Anova(sel_mod, type=2)
anova(sel_mod)
#Marginal and conditional R²
r2 <- r.squaredGLMM(sel_mod)
r2

#Prediction on the model by site (Figure 4)
nd <- expand.grid(site=c("KIP", "KEK","MON","HED","MUS"),
            light=c("Open-canopy", "Intermediate", "Dense-canopy"), 
            block=c("A", "B", "C")) %>% 
  left_join(germ %>% 
              group_by(site) %>% 
              summarise(temperature=mean(temperature),
                        precipitation=mean(precipitation),
                        nb_late_frost_event=mean(nb_late_frost_event)) %>% 
              dplyr::select(site,temperature, precipitation, nb_late_frost_event))

pred <- predict(sel_mod, newdata=nd, se.fit = TRUE)

germ_pred <- data.frame(
  fit = plogis(pred$fit)*100,           
  lower = plogis(pred$fit - 1.96*pred$se.fit)*100,  
  upper = plogis(pred$fit + 1.96*pred$se.fit)*100
) %>%
  bind_cols(nd) %>%
  group_by(site, light) %>%
  summarise(
    fit = mean(fit),                  
    lower = mean(lower),                    
    upper = mean(upper)                      
  ) %>%
  ungroup()

rect_df <- data.frame(xmin = c(0, 2), xmax = c(2, 6),
                      labs = c("within or", "north of the current sugar maple range")) %>% 
  rename(`Site`=labs)

rect_df$`Site` <- factor(rect_df$`Site`, levels = c("within or", "north of the current sugar maple range"))

germ_pred$site <- factor(germ_pred$site, levels = c("KIP", "KEK", "MON", "HED", "MUS"))

#Get observed median
germ_obs <- germ %>% 
  group_by(site, light) %>% 
  summarize(fit=median(germination))

germ_pred$light <- factor(germ_pred$light, levels = c("Dense-canopy", "Intermediate", "Open-canopy"))


plot_light_site <- ggplot()+
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, 
                                ymin = -Inf, ymax = Inf, 
                                fill = `Site`), 
            alpha = .3, inherit.aes = FALSE)+
  scale_fill_viridis_d() +
  geom_point(data=germ_pred, aes(x=site, y=fit)) +
  geom_errorbar(data=germ_pred, aes(x=site,ymin=lower, ymax=upper)) +
  geom_point(data=germ_obs, aes(x=site, y=fit), color="red", shape = 3, size=3) +
  facet_wrap(~light) +
  xlab("Sites") + 
  ylab("Germination rate") + 
  theme(strip.text.x = element_text(size = 12),
        axis.text=element_text(size=12),
        axis.title=element_text(size=14),
        legend.title = element_text(size = 12),
        legend.text = element_text(size=12),
        legend.position = "bottom")
  
ggsave(plot=plot_light_site, filename="~/postdoc/germination/github/figures/germination_site_effect.png", 
       width=8, height=4)






#Prédiction for each climate effect with interaction (Figure 3)
pred_temp  <- ggpredict(sel_mod, terms = c("temperature [all]", "light"),
                        back_transform = TRUE)

pred_prec  <- ggpredict(sel_mod, terms = c("precipitation [all]", "light"))

pred_frost <- ggpredict(sel_mod, terms = c("nb_late_frost_event [all]", "light"))

make_plot <- function(pred, xlab) {
  ggplot(pred, aes(x = x, y = predicted, color = group)) +
    geom_line(linewidth = 1) +
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = group),
                alpha = 0.2, color = NA) +
    labs(x = xlab, y = "Predicted germination probability", color = "Light treatment", fill = "Light treatment") +
    theme_minimal()+ 
    theme(strip.text.x = element_text(size = 14),
          axis.text=element_text(size=14),
          axis.title=element_text(size=15),
          legend.title = element_text(size = 16),
          legend.text = element_text(size=14),
          legend.position = "bottom")
}

p1 <- make_plot(pred_temp,  "Temperature")
p2 <- make_plot(pred_prec,  "Precipitation")
p3 <- make_plot(pred_frost, "# late frost events")

#Combine plots
plot_germ_effect <- ggarrange(p1, p2, p3, nrow=1, ncol=3, 
                              common.legend = T, legend="bottom",
                              labels=c("a","b","c"),
                              label.x = c(-0.02),  # décalage spécifique pour chaque
                              label.y = c(1.02)) + 
  bgcolor("white")+
  border("white")
# combinaison en une seule figure

ggsave(plot=plot_germ_effect, filename="~/postdoc/germination/github/figures/germination_interaction_effect.png", 
       width=10, height=4)

