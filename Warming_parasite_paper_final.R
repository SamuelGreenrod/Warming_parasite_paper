## Warming alters life-history traits and competition in a microbial parasite community
## Samuel T.E Greenrod, Daniel Cazares Lopez, Serena Johnson, Tobias E Hector, Emily J Stevens, R Craig MacLean, and Kayla C King

# This R script contains the code used to generate figures and carry out statistical tests for the paper titled above.

# Dependencies

require("car")
require("ggplot2")
require("RColorBrewer")
require("patchwork")
require("tidyverse") 
require("zoo") 
require("broom")
require("ggbeeswarm")
require("nlme")
require("emmeans")
require("ggsignif")
require("FSA")
require("pgirmess")
require("DescTools")
require("ggpmisc")
require("devtools")
require("ggradar")
require("scales")
require(colorblindr)
require("ggtext")



#Figure S1 - Optical density is correlated with colony-forming units
CFU_OD_regression <- read_csv("CFU_OD_regression.csv")

model <- CFU_OD_regression %>%
  lm(sqrt(CFU_ml) ~ OD, data = .)

autoplot(model)
summary(model)

model <- CFU_OD_regression %>%
  lm(sqrt(CFU_ml) ~ OD, data = .)

predslm <- predict(model,type = "response", se.fit=T)

CFU_OD_regression <- CFU_OD_regression %>%
  mutate(predictions = predslm$fit) %>%
  mutate(se_fit = predslm$se.fit) %>%
  mutate(conf_int = se_fit * 1.96) %>%
  mutate(upper_fit = predictions + conf_int) %>%
  mutate(lower_fit = predictions - conf_int)

OD_CFU <- ggplot(CFU_OD_regression, aes(OD, sqrt(CFU_ml))) +
  geom_ribbon( aes(ymin = lower_fit, ymax = upper_fit), alpha = .15) +
  geom_line( aes(y = predictions), size = 0.8, col = "darkred")+
  geom_point()+
  labs(x="Optical density (OD595)", y="Colony-forming units (CFU/ml) (square-root)", colour = "Time point")+
  theme_bw()+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=12,face="bold"))

ggsave("CFU_OD.tiff")



## Figure 1 - High temperatures can restrict phage infectivity

# Figure 1A - Bacterial growth in presence of phage (OD)

BactOD <- read.csv("OD_readings.csv",fileEncoding="UTF-8-BOM")

BactOD$OD <- as.numeric(BactOD$OD)
BactOD$Temp <- as.factor(BactOD$Temp)
BactOD$Phage <- factor(BactOD$Phage, levels = c("14_1", "3-phage_cocktail", "LUZ19" ,"LUZ19_14-1" ,  "No_phage" ,   "PEV2" , "PEV2_14-1" , "PEV2_LUZ19"))

nophage <- BactOD %>%
  filter(Phage=="No_phage")

singlephage <- BactOD %>%
  filter(Phage %in% c("No_phage", "PEV2", "LUZ19", "14_1"))
    
# Single phage

singlephage$Phage <- factor(singlephage$Phage, levels=c("No_phage", "PEV2", "LUZ19","14_1"))
singlephage$Temp <- factor(singlephage$Temp, levels=c("37","40","42"), labels=c("37\u00B0C","40\u00B0C","42\u00B0C"))

bacterialgrowth_wphage <- ggplot(singlephage, aes(hour, OD, colour = as.factor(Phage))) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Phage), stat = "summary", fun = mean, size=0.8) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1,size=0.8)+
  theme_bw()+
  ylab("Bacterial density (OD595)") +
  xlab("Time (hours)")+
  labs(colour="Treatment")+
  scale_colour_manual(values = c("#fde725", "#a0da39", "#277f8e","#440154"),labels=c("No phage",expression(phi*PEV2),expression(phi*LUZ19),expression(phi*'14-1')))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12),legend.text.align = 0)+
  theme(axis.title=element_text(size=12,face="bold"))+
  theme(axis.title.x = element_blank())+
  scale_x_continuous(breaks = 0:5)+
  facet_grid(~Temp)+
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_rect(fill="lightgrey", colour="black",size=0.3))

## OD stats for single phage

# Prepare the data
OD_stats <- singlephage %>%
  group_by(Temp) %>%
  filter(hour == "5") %>%
  nest()

# Run model

OD_models <- OD_stats %>%
  mutate(model = map(data, ~lme(OD ~ Phage, random= ~1|Batch, data= ., method = "REML")))

# Diagnostic plots showing unequal variances
T37 <- OD_models$model[[1]]
plot(T37) 

T40 <- OD_models$model[[2]]
plot(T40)

T42 <- OD_models$model[[3]]
plot(T42)


# Revised model accounting for unequal variances including pairwise comparisons

revised_OD_models <- OD_stats %>%
  mutate(model = map(data, ~lme(OD ~ Phage, random= ~1|Batch, data= ., method = "REML", weights = varIdent(form = ~1|Phage)))) %>%
  mutate(model_summary = map(model, summary)) %>%
  mutate(model_pairwise = map(model, ~ {
    emmeans(.x, pairwise ~ Phage)
    }))

T37_new <- revised_OD_models$model[[1]]
plot(T37_new) 

T40_new <- revised_OD_models$model[[2]]
plot(T40_new)

T42_new <- revised_OD_models$model[[3]]
plot(T42_new)

revised_OD_models$model_summary
revised_OD_models$model_pairwise



# Figure 1B - Plaque assay results

pfu <- read.csv("Plaque_assay_data.csv",fileEncoding="UTF-8-BOM")

pfu$Rel_pfu_ml <- as.numeric(pfu$Rel_pfu_ml)
pfu$Time_hour <- as.numeric(pfu$Time_hour)

pfu$Temp <- as.factor(pfu$Temp)
pfu$Phage <- factor(pfu$Phage, levels=c("PEV2","LUZ19","phage_14_1"), labels=c("\U03D5PEV2","\U03D5LUZ19","\U03D5 14-1"))


phagegrowth_temp_pfu <- ggplot(pfu, aes(Time_hour, Rel_pfu_ml, colour = as.factor(Temp))) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Temp), stat = "summary", fun = mean, size = 0.8) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1, size = 0.8)+
  scale_y_continuous(trans='log10', limits = c(0.01,10^5)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  ylab("Relative phage density") +
  xlab("Time (hours)")+
  labs(colour="Temperature (°C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=12,face="bold"))+
  theme(axis.title.x = element_blank())+
  scale_x_continuous(breaks = 0:5)+
  facet_grid(~Phage)+
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_rect(fill="lightgrey", colour="black",size=0.3))


# Make models
pfu_grouped_data <- pfu %>%
  group_by(Temp, Phage) %>%
  filter(Time %in% c("T0", "T5")) %>%
  nest()
  
pfu_models <- pfu_grouped_data %>%
  mutate(model = map(data, ~ lme(pfu_ml ~ Time, random= ~1|Rep, data= ., method = "REML", weights = varIdent(form = ~1|Time)))) %>%
  mutate(model_summaries = map(model, summary))

# Diagnostic plots
length(pfu_models$model)

plots <- lapply(1:9, function(i) plot(pfu_models$model[[i]], main = paste(pfu_models$Phage[[i]], pfu_models$Temp[[i]])))
cowplot::plot_grid(plotlist = plots)


# Generate model summaries
pfu_models$model_summaries





## Figure 1C - phage densities based on DNA concentrations with qPCR

qPCR_data <- read.csv("qPCR_data.csv",fileEncoding="UTF-8-BOM")

qPCR_data$Fold_change_DNAconc_to_zero <- as.numeric(qPCR_data$Fold_change_DNAconc_to_zero)
qPCR_data$Time <- as.numeric(qPCR_data$Time)
qPCR_data$Temp <- as.factor(qPCR_data$Temp)

qPCR_data_single <- qPCR_data %>%
  filter(Treatment %in% c('P','L','phage14'))

options(scipen=999) # To make the Y-axis correct numbers

qPCR_data_single$Phage <- factor(qPCR_data_single$Phage, levels=c("PEV2","LUZ19","phage14"), labels=c("\U03D5PEV2","\U03D5LUZ19","\U03D5 14-1"))

phagegrowth_wtemp <- ggplot(qPCR_data_single, aes(Time, Fold_change_DNAconc_to_zero, colour = as.factor(Temp))) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Temp), stat = "summary", fun = mean, size = 0.8) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1, size = 0.8)+
  scale_y_continuous(trans='log10', limits = c(0.1,10^4)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  ylab("Relative\n phage DNA concentration") +
  xlab("Time (hours)")+
  labs(colour="Temperature (°C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=12,face="bold"))+
  scale_x_continuous(breaks = 0:5)+
  facet_grid(~Phage)+
  theme(strip.text = element_text(face="bold", size=9),
        strip.background = element_rect(fill="lightgrey", colour="black",size=0.3))


## qPCR stats for single phage over time

# Make models
qpcr_grouped_data <- qPCR_data_single %>%
  group_by(Temp, Treatment) %>%
  filter(Time %in% c("0", "5")) %>%
  nest()

qpcr_models <- qpcr_grouped_data %>%
  mutate(model = map(data, ~ lme(log(DNA_conc) ~ Time, random= ~1|Rep, data= ., method = "REML"))) %>%
  mutate(model_summaries = map(model, summary))

# Diagnostic plots
length(qpcr_models$model)

plots <- lapply(1:9, function(i) plot(qpcr_models$model[[i]], main = paste(qpcr_models$Treatment[[i]], qpcr_models$Temp[[i]])))
cowplot::plot_grid(plotlist = plots)


# Generate model summaries
qpcr_models$model_summaries


# Plot
bacterialgrowth_wphage / phagegrowth_temp_pfu / phagegrowth_wtemp + plot_layout(guides="collect")  & theme(legend.position = 'bottom')

ggsave("Fig_1.tiff", width = 9, height = 9)



## Figure S2. Temperature has a phage-specific impact on phage lytic activity

CFU_wphage <- read.csv("CFU_sodiumcitrate.csv",fileEncoding="UTF-8-BOM")

CFU_wphage$Phage <- factor(CFU_wphage$Phage, levels=c("No_phage","pev2","luz19","14_1"), labels=c("No phage","PEV2","LUZ19","14-1"))
CFU_wphage$Temp <- as.factor(CFU_wphage$Temp)


# Set 0s to 1 so data can be plotted on log scale
CFU_wphage <- CFU_wphage %>%
  dplyr::mutate(standard_CFU_ml = ifelse(CFU_ml == 0, 1, CFU_ml))

ggplot(CFU_wphage, aes(hour, standard_CFU_ml, colour = as.factor(Temp))) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Temp), stat = "summary", fun = mean, size=0.8) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1,size=0.8)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x), 
                labels = trans_format("log10", scales::math_format(10^.x)))+
  annotation_logticks(sides="l")+
  theme_bw()+
  ylab("Bacterial density (CFU/ml)") +
  xlab("Time (hours)")+
  labs(colour="Temperature (C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12),legend.text.align = 0)+
  theme(axis.title=element_text(size=11,face="bold"))+
  facet_grid(~Phage)

ggsave("Fig_S2.tiff")


## Figure S3. Bacterial growth in absence of phage is similar between temperatures

# Figure S3A

BactOD <- read_csv("OD_readings.csv",fileEncoding="UTF-8-BOM")

BactOD$OD <- as.numeric(BactOD$OD)
BactOD$Temp <- as.factor(BactOD$Temp)
BactOD$Phage <- factor(BactOD$Phage, levels = c("14_1", "3-phage_cocktail", "LUZ19" ,"LUZ19_14-1" ,  "No_phage" ,   "PEV2" , "PEV2_14-1" , "PEV2_LUZ19"))

nophage <- BactOD %>%
  filter(Phage=="No_phage")

# Model
nophage_model <- lme(OD ~ Temp + time_point + Temp:time_point, random= ~1|Batch, data= nophage)

summary(model)
anova(model)

nophage_OD <- ggplot(nophage, aes(hour, OD, colour = as.factor(Temp))) +
  geom_jitter(width = 0.2) +
  geom_line(aes(group = Temp), stat = "summary", fun = mean, size = 0.8) +
  geom_errorbar(stat = "summary", fun.data = function(x) {
    data.frame(ymin = mean(x) - sd(x), ymax = mean(x) + sd(x))
  }, width = 0.1, linewidth = 0.8)+
  theme_bw()+
  ylab("Bacterial density (OD595)") +
  xlab("Time (hours)")+
  labs(colour="Temperature (C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  theme(legend.title = element_text(colour="black", size=12,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=12,face="bold"))



# Figure S3B

CFU_nophage <- read_csv("T0_T5_CFU.csv")

nophage_CFU <- ggplot(data=CFU_nophage, aes(x = as.factor(hour), y= CFU_ml, fill=as.factor(Temp))) + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1,lwd=1, alpha = 0.7,show.legend = F,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=1,show.legend = F)+
  ylab("Bacterial density (CFU/ml)") +
  xlab("Time point") +
  #ylim((10^8),(10^10))+
  scale_y_continuous(trans='log10', limits = c(10^7,10^9)) +
  annotation_logticks(sides="l")+
  scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  theme_bw()+
  geom_smooth(method = "lm")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=12,face="bold"))+
  theme(legend.position="none")

nophage_OD + nophage_CFU+ plot_layout(guides="collect")


# No phage CFU stats

CFU_nophage_T5 <- CFU_nophage %>%
  filter(hour =="5")

model <- lme(CFU_ml ~ Temp, random= ~1|Biol_rep, data= CFU_nophage_T5, weights = varIdent(form = ~1|Temp))
plot(model)
summary(model)
anova(model)
emmeans(model, pairwise ~ Temp)

# Get number of bacterial doublings
group_mean <- CFU_nophage %>%
  dplyr::group_by(hour, Temp) %>%
  dplyr::summarise(mean_value = mean(CFU_ml))

group_mean %>%
  dplyr::group_by(Temp) %>%
  dplyr::summarise(Bacterial_doublings = log2(mean_value[hour == 5] / mean_value[hour == 0]))




### Figure 2 - Phage life-history traits across temperatures

### Figure S4 - Survivorship

decayrate_assay <- read.csv("Decay_rate_120523.csv",fileEncoding="UTF-8-BOM")

decayrate_assay$PFU_ML <- as.numeric(decayrate_assay$PFU_ML)
decayrate_assay$Time <- as.numeric(decayrate_assay$Time)
decayrate_assay$Temp <- as.factor(decayrate_assay$Temp)

decayrate_assay <- decayrate_assay %>%
  filter(Time<25)

decayrate_assay$Phage <- factor(decayrate_assay$Phage, levels=c("PEV2","LUZ19","phage14_1"), labels=c("PEV2","LUZ19","14-1"))


## Survivorship stats

# For survivorship and attachment stats, the 95% confidence intervals of time to 50% decay/adsorption was predicted from models

decay_rate_data <- decayrate_assay %>%
  mutate(rounded_pfu = round(PFU_ML)) %>%
  group_by(Phage) %>%
  nest()

decay_rate_models <- decay_rate_data %>%
  mutate(model = map(data, ~ glm(rounded_pfu ~ Time + Temp + Time:Temp, data= ., family="poisson")))%>%
  mutate(model_summaries = map(model, summary))

# Diagnostic plots
length(decay_rate_models$model)

plots <- lapply(1:3, function(i) plot(decay_rate_models$model[[i]], main = paste(decay_rate_models$Phage[[i]])))
cowplot::plot_grid(plotlist = plots)

decay_rate_models$model_summaries

# Model predictions between 1 and 100 hours

new_x <- seq(0, 1000, 0.1)

predict_for_temp <- function(model, temp, new_x) {
  predict(model, 
          newdata = data.frame(Time = new_x, Temp = factor(rep(temp, length(new_x)))),
          type = "response", se.fit = TRUE)
}

decay_rate_models_predict <- decay_rate_models %>%
  mutate(T37 = map(model, ~ predict_for_temp(.x, "37", new_x)),
         T40 = map(model, ~ predict_for_temp(.x, "40", new_x)),
         T42 = map(model, ~ predict_for_temp(.x, "42", new_x))) %>%
  select(Phage, data, T37, T40, T42) %>%
  pivot_longer(
    cols = starts_with("T"),
    names_to = "Temp",
    values_to = "predictions"
  )


# Extract fitted values and calculate upper and lower confidence interval bounds

model_data <- decay_rate_models_predict %>%
  mutate(
    fitted = map(predictions, "fit"),
    se = map(predictions, "se.fit")
    ) %>%
  unnest(cols = c(fitted, se)) %>%
  mutate(conf_int = se * 1.96,
         upper_fit = fitted + conf_int,
         lower_fit = fitted - conf_int,
         ) %>%
  dplyr::group_by(Phage, Temp) %>%
  dplyr::mutate(Time = new_x) %>%
  ungroup() %>%
  dplyr::group_by(Phage) %>%
  dplyr::mutate(half_dens = map(data, ~max(.x$PFU_ML/2))) %>% # Include the value when phage density has decreased by half
  ungroup() %>%
  select(Phage, Temp, fitted, se, upper_fit, lower_fit, half_dens, Time)

model_data$Temp <- factor(model_data$Temp, levels=c("T37","T40", "T42"), labels=c("37","40","42"))


# Calculate the time predicted for phage densities to decrease by 50%
# Note values are not provided for 14-1 37C as decay is expected to take longer than 1000 hours

model_data %>%
  group_by(Phage, Temp) %>%
  filter(upper_fit > half_dens & lower_fit < half_dens) %>%
  mutate(time_to_50 = paste(min(Time), "-", max(Time))) %>%
  select(Phage, Temp, time_to_50) %>%
  distinct()


# Plot of phage decay

model_for_plotting <- model_data %>%
  filter(Time <= 24)

ggplot(decayrate_assay, aes(Time, PFU_ML, group = Temp)) +
  geom_jitter(aes(shape = Temp, fill = Temp),width = 1, size=2.2,col = "black") +
  geom_line(data = model_for_plotting, aes(y = fitted, color = Temp), size = 0.8,show.legend = FALSE) +
  geom_ribbon(data = model_for_plotting, aes(x = Time, ymin = lower_fit, ymax = upper_fit, fill = Temp), alpha = 0.2, inherit.aes = FALSE) +
  xlim(-2,26)+
  scale_y_continuous(trans='log10', limits = c(10^4,10^6)) +
  annotation_logticks(sides="l")+
  #ylim(0,500000)+
  ylab("Phage density (PFU/ml)") +
  xlab("Time (hours)")+
  labs(fill="Temperature (\u00B0C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"), guide="none")+
  scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  scale_shape_manual(values = c(23,21,24), guide="none")+
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=11,face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape=c(23,21,24)))) +
  facet_grid(~ Phage)

ggsave("Fig_S4_update.tiff")


### Figure S5 - Host attachment

adsorption_assay <- read.csv("Adsorption_assay_020523.csv",fileEncoding="UTF-8-BOM")

adsorption_assay <- adsorption_assay %>%
  mutate(across(c(PFU_ML, PFU_ML_round,Rel_PFU_ml,Time), as.numeric)) %>%
  mutate(Temp = as.factor(Temp))

adsorption_assay$Phage <- factor(adsorption_assay$Phage, levels=c("PEV2","LUZ19","phage14_1"), labels=c("PEV2","LUZ19","14-1"))

adsorption_assay_filtered <- adsorption_assay %>%
  filter((Phage %in% c("PEV2", "LUZ19") & Time < 25) |
           (Phage == "14-1" & Time < 10))


## Attachment stats


adsorption_rate_models <- adsorption_assay_filtered %>%
  group_by(Phage) %>%
  nest() %>%
  mutate(model = map2(data, Phage, ~ {
    if (.y %in% c("PEV2", "LUZ19")) {
      glm(PFU_ML_round ~ Time + Temp + Time:Temp, data = ., family = "poisson")
    } else {
      lm(PFU_ML_round ~ Time + Temp + Time:Temp, data = .)
    }
  })) %>%
  mutate(model_summaries = map(model, summary))



# Diagnostic plots
length(adsorption_rate_models_P_L$model)

plots <- lapply(1:3, function(i) plot(adsorption_rate_models_P_L$model[[i]], main = paste(adsorption_rate_models$Phage[[i]])))
cowplot::plot_grid(plotlist = plots)

adsorption_rate_models$model_summaries

# Model predictions between 1 and 100 hours

new_x <- seq(0, 50, 0.1)

predict_for_temp <- function(model, temp, new_x) {
  predict(model, 
          newdata = data.frame(Time = new_x, Temp = factor(rep(temp, length(new_x)))),
          type = "response", se.fit = TRUE)
}


# PEV2 and LUZ19 GLMs
adsorption_rate_models_predict <- adsorption_rate_models %>%
  mutate(T37 = map(model, ~ predict_for_temp(.x, "37", new_x)),
         T40 = map(model, ~ predict_for_temp(.x, "40", new_x)),
         T42 = map(model, ~ predict_for_temp(.x, "42", new_x))) %>%
  select(Phage, data, T37, T40, T42) %>%
  pivot_longer(
    cols = starts_with("T"),
    names_to = "Temp",
    values_to = "predictions"
  )


# Extract fitted values and calculate upper and lower confidence interval bounds

model_data <- adsorption_rate_models_predict %>%
  mutate(
    fitted = map(predictions, "fit"),
    se = map(predictions, "se.fit")
  ) %>%
  unnest(cols = c(fitted, se)) %>%
  mutate(conf_int = se * 1.96,
         upper_fit = fitted + conf_int,
         lower_fit = fitted - conf_int,
  ) %>%
  dplyr::group_by(Phage, Temp) %>%
  dplyr::mutate(Time = new_x) %>%
  ungroup() %>%
  dplyr::group_by(Phage) %>%
  dplyr::mutate(half_dens = map(data, ~max(.x$PFU_ML_round/2))) %>% # Include the value when phage density has decreased by half
  ungroup() %>%
  select(Phage, Temp, fitted, se, upper_fit, lower_fit, half_dens, Time)

model_data$Temp <- factor(model_data$Temp, levels=c("T37","T40", "T42"), labels=c("37","40","42"))


# Calculate the time predicted for phage densities to decrease by 50% (50% attachment)
# Note values are not provided for 14-1 37C as decay is expected to take longer than 1000 hours

model_data %>%
  group_by(Phage, Temp) %>%
  filter(upper_fit > half_dens & lower_fit < half_dens) %>%
  mutate(time_to_50 = paste(min(Time), "-", max(Time))) %>%
  select(Phage, Temp, time_to_50) %>%
  distinct()



# Plot of phage adsorption

model_for_plotting <- model_data %>%
  filter(Time <= 25 & fitted > 0)

data_for_plotting <- adsorption_assay %>%
  filter(Time < 25)

ggplot(data_for_plotting, aes(Time, PFU_ML, group = Temp)) +
  geom_jitter(aes(shape = Temp, fill = Temp),width = 1, size=2.2,col = "black") +
  geom_line(data = model_for_plotting, aes(y = fitted, color = Temp), size = 0.8,show.legend = FALSE) +
  geom_ribbon(data = model_for_plotting, aes(x = Time, ymin = lower_fit, ymax = upper_fit, fill = Temp), alpha = 0.2, inherit.aes = FALSE) +
  xlim(-2,25)+
  #scale_y_continuous(trans='log10', limits = c(10^3,10^5)) +
  #annotation_logticks(sides="l")+
  #ylim(0,500000)+
  ylab("Phage density (PFU/ml)") +
  xlab("Time (hours)")+
  labs(fill="Temperature (\u00B0C)")+
  scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"), guide="none")+
  scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  scale_shape_manual(values = c(23,21,24), guide="none")+
  theme_bw() +
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=12))+
  theme(axis.title=element_text(size=11,face="bold"))+
  guides(fill = guide_legend(override.aes = list(shape=c(23,21,24)))) +
  facet_grid(~ Phage)

ggsave("Fig_S5_update.tiff")



## Figure 2A - Radar charts showing life-history traits

radar <- read_csv("Radar_charts_correct_new.csv")

radar_rename <- radar %>%
  select(Phage, Temperature, Decay_prop, Recip_adsorption_time_prop, log_T5_T0_PFU_prop_norm, log_T5_T0_CFU_prop_norm) %>%
  rename("Stability" = Decay_prop,
         "Host_attachment" = Recip_adsorption_time_prop,
         "Reproductive_capacity" = log_T5_T0_PFU_prop_norm,
         "Virulence" = log_T5_T0_CFU_prop_norm)

radar_rename <- radar_rename %>%
  mutate(across(c(Stability, Host_attachment,Reproductive_capacity,Virulence), as.numeric))


# Make radar plots

radar_rename %>%
  dplyr::group_by(Phage) %>%
  dplyr::group_map(~ {

    p <- ggradar(., 
                values.radar = c(0, 0.5, 1),
                background.circle.colour = "grey90",
                axis.line.colour = "gray60",
                plot.extent.x.sf =1.2,
                gridline.min.colour = "gray60",
                gridline.mid.colour = "gray60",
                gridline.max.colour = "gray60",
                gridline.min.linetype = 1,
                gridline.max.linetype = 1,
                axis.label.size = 4, 
                grid.label.size = 3)+
      scale_color_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))
    
    # Save the plot
    ggsave(filename = paste0("plot_", unique(.y$Phage), ".tiff"), plot = p, width = 8, height = 5)
  })



### Figure 2B - Correlation plot of virulence and population growth at T2

radar$Temperature <- as.factor(radar$Temperature)
radar$Temperature <- factor(radar$Temperature, levels = c("37","40","42"), labels = c("37°C","40°C","42°C"))
radar$Phage <- factor(radar$Phage, levels=c("PEV2", "LUZ19", "phage14"), labels= c(expression('\u03d5PEV2'), expression("\u03d5LUZ19"), expression("\u03d514-1")))


ggplot(radar,aes(x=Log_average_ratio_T2_and_T0_CFU,y=Log_average_ratio_T2_and_T0_PFU, group=Phage))+ 
  geom_point(aes(shape = Phage, fill = Temperature), col="black",size=5)+
  geom_errorbar(aes(ymin=Log_average_ratio_T2_and_T0_PFU-Log_average_ratio_T2_and_T0_PFU_SE, ymax=Log_average_ratio_T2_and_T0_PFU+Log_average_ratio_T2_and_T0_PFU_SE), width=.01,col="darkgrey",
                position=position_dodge(.01)) +
  geom_errorbar(aes(xmin=Log_average_ratio_T2_and_T0_CFU-Log_average_ratio_T2_and_T0_CFU_SE, xmax=Log_average_ratio_T2_and_T0_CFU + Log_average_ratio_T2_and_T0_CFU_SE), height=.03,col="darkgrey",
                position=position_dodge(.01)) +
  ylab("Phage population growth (phage doublings after 2h)")+
  xlab("Phage virulence (bacterial doublings after 2h)") +
  scale_x_reverse()+
  theme_bw()+
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+
  scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  scale_shape_manual(values = c(23,21,24))+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=9))+
  guides(fill = guide_legend(override.aes = list(shape=22)))

ggsave("Fig_2B.tiff")


## Virulence stats - calculating standard errors of the difference between means for each temperature at T5

# T5 virulence stats

CFU_wphage <- read_csv("CFU_sodiumcitrate.csv")

CFU_wphage$Phage <- factor(CFU_wphage$Phage, levels=c("No_phage","pev2","luz19","14_1"), labels=c("No phage","PEV2","LUZ19","14-1"))
CFU_wphage$Temp <- as.factor(CFU_wphage$Temp)

CFU_group <- CFU_wphage %>%
  filter(hour < 1 | hour == 5, Phage != "No phage") %>%
  group_by(Phage, Temp, hour) %>%
  summarise(mean = mean(CFU_ml, na.rm = TRUE),
            se = sd(CFU_ml, na.rm = TRUE) / sqrt(n()))
  
CFU_group %>%
  pivot_wider(names_from = hour, values_from = c(mean, se)) %>%
  mutate(mean_ratio = log2(mean_5/mean_0),
         se_ratio = sqrt(((se_0^2)/(mean_0^2)) + ((se_5^2)/(mean_5^2)))) %>%
  select(Phage, Temp, mean_ratio, se_ratio)


# Virulence stats after 2h

CFU_group <- CFU_wphage %>%
  filter(hour < 1 | hour == 2, Phage != "No phage") %>%
  group_by(Phage, Temp, hour) %>%
  summarise(mean = mean(CFU_ml, na.rm = TRUE),
            se = sd(CFU_ml, na.rm = TRUE) / sqrt(n()))

CFU_group %>%
  pivot_wider(names_from = hour, values_from = c(mean, se)) %>%
  mutate(mean_ratio = log2(mean_2/mean_0),
         se_ratio = sqrt(((se_0^2)/(mean_0^2)) + ((se_2^2)/(mean_2^2)))) %>%
  select(Phage, Temp, mean_ratio, se_ratio)



## Population growth stats - calculating standard errors of the difference between means for each temperature at T5

# T5 population growth

pfu <- read.csv("Plaque_assay_data.csv",fileEncoding="UTF-8-BOM")

pfu$Rel_pfu_ml <- as.numeric(pfu$Rel_pfu_ml)
pfu$Time_hour <- as.numeric(pfu$Time_hour)

pfu$Temp <- as.factor(pfu$Temp)
pfu$Phage <- factor(pfu$Phage, levels=c("PEV2","LUZ19","phage_14_1"), labels=c("\U03D5PEV2","\U03D5LUZ19","\U03D5 14-1"))

PFU_group <- pfu %>%
  filter(Time_hour < 1 | Time_hour == 5, Phage != "No phage") %>%
  group_by(Phage, Temp, Time_hour) %>%
  summarise(mean = mean(pfu_ml, na.rm = TRUE),
            se = sd(pfu_ml, na.rm = TRUE) / sqrt(n()))

PFU_group %>%
  pivot_wider(names_from = Time_hour, values_from = c(mean, se)) %>%
  mutate(mean_ratio = log2(mean_5/mean_0),
         se_ratio = sqrt(((se_0^2)/(mean_0^2)) + ((se_5^2)/(mean_5^2)))) %>%
  select(Phage, Temp, mean_ratio, se_ratio)



# Population growth stats after 2h

PFU_group <- pfu %>%
  filter(Time_hour < 1 | Time_hour == 2, Phage != "No phage") %>%
  group_by(Phage, Temp, Time_hour) %>%
  summarise(mean = mean(pfu_ml, na.rm = TRUE),
            se = sd(pfu_ml, na.rm = TRUE) / sqrt(n()))

PFU_group %>%
  pivot_wider(names_from = Time_hour, values_from = c(mean, se)) %>%
  mutate(mean_ratio = log2(mean_2/mean_0),
         se_ratio = sqrt(((se_0^2)/(mean_0^2)) + ((se_2^2)/(mean_2^2)))) %>%
  select(Phage, Temp, mean_ratio, se_ratio)



### Figure 3 - Competition between phages across temperatures

# Figure 3A - Phage competition at 37C and 40C

qPCR_data <- read_csv("qPCR_data.csv")

qPCR_data$Fold_change_DNAconc_to_zero <- as.numeric(qPCR_data$Fold_change_DNAconc_to_zero)
qPCR_data$Time <- as.numeric(qPCR_data$Time)
qPCR_data$Temp <- as.factor(qPCR_data$Temp)


qPCR_data_group <- qPCR_data %>%
  group_by(Phage) %>%
  filter(Time == "5" & Temp != 42) %>%
  nest()

qpcr_models <- qPCR_data_group %>%
  mutate(model = map(data, ~ lme(DNA_conc ~ Competitor + Temp + Competitor:Temp, random= ~1|Batch, data= ., method = "REML", weights = varIdent(form = ~1|Competitor)))) %>%
  mutate(model_summaries = map(model, summary))%>%
  mutate(model_pairwise = map(model, ~ {
    emmeans(.x, pairwise ~ Competitor:Temp)
  }))

# Diagnostic plots
length(qpcr_models$model)

plots <- lapply(1:3, function(i) plot(qpcr_models$model[[i]], main = paste(qpcr_models$Treatment[[i]], qpcr_models$Temp[[i]])))
cowplot::plot_grid(plotlist = plots)

# Generate model summaries
qpcr_models$model_summaries
qpcr_models$model_pairwise


# Plots - separate plots needed for each phage as different competitor order and different significance indicators

# Filter dataset
PEV2_qPCR <- qPCR_data %>%
  filter(Phage=="PEV2" & Time == "5" & Temp != 42)

LUZ19_qPCR <- qPCR_data %>%
  filter(Phage=="LUZ19" & Time == "5" & Temp != 42)

phage14_qPCR <- qPCR_data %>%
  filter(Phage=="phage14" & Time == "5" & Temp != 42)

# Order competitors

PEV2_qPCR$Competitor <- factor(PEV2_qPCR$Competitor, levels=c("No_competitor","LUZ19", "phage14","3P"), labels=c("No competitor","LUZ19","14-1","3-phage"))
LUZ19_qPCR$Competitor <- factor(LUZ19_qPCR$Competitor, levels=c("No_competitor","PEV2", "phage14","3P"), labels=c("No competitor","PEV2","14-1","3-phage"))
phage14_qPCR$Competitor <- factor(phage14_qPCR$Competitor, levels=c("No_competitor","PEV2", "LUZ19","3P"), labels=c("No competitor","PEV2","LUZ19","3-phage"))

#Plots

pev2_plot <- ggplot(data=PEV2_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.5,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=3, dodge.width = .75, color="black",alpha=.7,show.legend = F)+
  ylab("Fold increase in phage DNA (T0-T5)") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(100,22000)) +
  annotation_logticks(sides="l")+
  geom_signif(y_position = 2.8, xmin = 2.8, 
              xmax = 3.17, annotation = "*",
              tip_length = 0.05,col = "firebrick")+
  geom_signif(y_position = 2.8, xmin = 0.8, 
              xmax = 1.8, annotation = "*",
              tip_length = 0.05)+
  theme_bw()+
  labs(fill="Temperature (\u00B0C)")+
  scale_fill_manual(values = c("#ffeda0","#feb24c"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=9))


luz19_plot <-ggplot(data=LUZ19_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=3, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Relative phage DNA concentration") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(100,22000)) +
  annotation_logticks(sides="l")+
  geom_signif(y_position = 2.8, xmin = 0.8, 
              xmax = 1.8, annotation = "***",
              tip_length = 0.05)+
  geom_signif(y_position = 3, xmin = 0.8, 
              xmax = 2.8, annotation = "*",
              tip_length = 0.05)+
  geom_signif(y_position = 3.2, xmin = 0.8, 
              xmax = 3.8, annotation = "**",
              tip_length = 0.05)+
  theme_bw()+
  labs(fill="Temperature (\u00B0C)")+
  scale_fill_manual(values = c("#ffeda0","#feb24c"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) + 
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=10,face="bold"),legend.text = element_text(size=9))


phage14_plot <- ggplot(data=phage14_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=3, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Phage growth (fold increase) across 5 hour period") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(100,22000)) +
  annotation_logticks(sides="l")+
  geom_signif(y_position = c(3.3,3.15), xmin = c(1.77,3.78), 
              xmax = c(2.18,4.2), annotation = c("***","***"),
              tip_length = 0.02, col = "firebrick")+
  geom_signif(y_position = c(4,4.2), xmin = c(0.8,0.8), 
              xmax = c(1.77,3.75), annotation = c("**","**"),
              tip_length = 0.02)+
  theme_bw()+
  labs(fill="Temperature (\u00B0C)")+
  scale_fill_manual(values = c("#ffeda0","#feb24c"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) + 
  theme(axis.title=element_text(size=12,face="bold")) + 
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=9))

pev2_plot + luz19_plot + phage14_plot + plot_layout(guides="collect")& theme(legend.position="bottom")

cvd_grid(pev2_plot) # Check for colour blind accessibility

ggsave("Fig_3A.tiff")



# Figure 3B - Summary plot of phage competitiveness

average_res <- read.csv("Phage_impact_on_phage_replication_average.csv",fileEncoding="UTF-8-BOM")

average_res$Temp <- as.factor(average_res$Temp)
average_res$Temp <- factor(average_res$Temp, levels = c("37C","40C","42C"), labels = c("37°C","40°C","42°C"))
average_res$Phage <- factor(average_res$Phage, levels=c("PEV2", "LUZ19", "14_1"), labels= c(expression('\u03d5PEV2'), expression("\u03d5LUZ19"), expression("\u03d514-1")))

ggplot(average_res,aes(x=Average_restricting_opposite,y=Average_restricted, group = Phage))+ 
  geom_point(aes(shape = Phage, fill = Temp), col="black",size=5)+
  geom_path(aes(x = Average_restricting_opposite, y = Average_restricted, group = Phage), 
            arrow = arrow(length = unit(0.35, "cm")), col="black")+
  labs(y="Resistance to competitors\n(Growth with competitors relative to no competition)",x="Suppression of competitors\n(Reduction in competitor growth relative to no competition)", fill="Temperature")+
  theme_bw()+
  ylim(0, 1.25)+
  xlim(-0.1,1)+
  scale_fill_manual(values = c("#ffeda0","#feb24c","#fc4e2a"))+
  scale_shape_manual(values = c(23,21,24))+
  geom_hline(yintercept=1, linetype="dashed", color = "black")+
  geom_vline(xintercept=0, linetype="dashed", color = "black")+
  theme(axis.title=element_text(size=11,face="bold"))+
  theme(legend.title = element_text(colour="black", size=11,face="bold"),legend.text = element_text(size=9))+
  guides(fill = guide_legend(override.aes = list(shape=22)))

ggsave("Fig_3C.tiff")



## Figure S6 - Phage competition at 37C after 2h

qPCR_data <- read_csv("qPCR_data.csv")

qPCR_data$Fold_change_DNAconc_to_zero <- as.numeric(qPCR_data$Fold_change_DNAconc_to_zero)
qPCR_data$Time <- as.numeric(qPCR_data$Time)
qPCR_data$Temp <- as.factor(qPCR_data$Temp)

qPCR_data_group <- qPCR_data %>%
  group_by(Phage) %>%
  filter(Time == "2" & Temp == 37) %>%
  nest()

qpcr_models <- qPCR_data_group %>%
  mutate(model = map(data, ~ lme(DNA_conc ~ Competitor, random= ~1|Batch, data= ., method = "REML", weights = varIdent(form = ~1|Competitor)))) %>%
  mutate(model_summaries = map(model, summary))%>%
  mutate(model_pairwise = map(model, ~ {
    emmeans(.x, pairwise ~ Competitor)
  }))

# Diagnostic plots
length(qpcr_models$model)

plots <- lapply(1:3, function(i) plot(qpcr_models$model[[i]], main = paste(qpcr_models$Treatment[[i]], qpcr_models$Temp[[i]])))
cowplot::plot_grid(plotlist = plots)

# Generate model summaries
qpcr_models$model_summaries
qpcr_models$model_pairwise


# Plots - separate plots needed for each phage as different competitor order and different significance indicators

# Filter dataset
PEV2_qPCR <- qPCR_data %>%
  filter(Phage=="PEV2" & Time == "2" & Temp == 37)

LUZ19_qPCR <- qPCR_data %>%
  filter(Phage=="LUZ19" & Time == "2" & Temp == 37)

phage14_qPCR <- qPCR_data %>%
  filter(Phage=="phage14" & Time == "2" & Temp == 37)

# Order competitors

PEV2_qPCR$Competitor <- factor(PEV2_qPCR$Competitor, levels=c("No_competitor","LUZ19", "phage14","3P"), labels=c("No competitor","LUZ19","14-1","3-phage"))
LUZ19_qPCR$Competitor <- factor(LUZ19_qPCR$Competitor, levels=c("No_competitor","PEV2", "phage14","3P"), labels=c("No competitor","PEV2","14-1","3-phage"))
phage14_qPCR$Competitor <- factor(phage14_qPCR$Competitor, levels=c("No_competitor","PEV2", "LUZ19","3P"), labels=c("No competitor","PEV2","LUZ19","3-phage"))

#Plots

pev2_plot <- ggplot(data=PEV2_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(10,1000)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  scale_fill_manual(values = "#ffeda0") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(axis.title.x=element_blank()) + 
  theme(legend.position = "none")


luz19_plot <- ggplot(data=LUZ19_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(10,1000)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  scale_fill_manual(values = "#ffeda0") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) + 
  theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.position = "none")


phage14_plot <- ggplot(data=phage14_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(10,1000)) +
  annotation_logticks(sides="l")+
  geom_signif(y_position = c(1.8), xmin = c(1), 
              xmax = c(4), annotation = c("*"),textsize=5,
              tip_length = 0.03)+
  theme_bw()+
  scale_fill_manual(values = "#ffeda0") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) +
  theme(axis.title.x=element_blank()) + 
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.position = "none")


pev2_plot + luz19_plot + phage14_plot

ggsave("Fig_S6.tiff")



## Figure S7 - Phage competition at 42C

qPCR_data <- read_csv("qPCR_data.csv")

qPCR_data$Fold_change_DNAconc_to_zero <- as.numeric(qPCR_data$Fold_change_DNAconc_to_zero)
qPCR_data$Time <- as.numeric(qPCR_data$Time)
qPCR_data$Temp <- as.factor(qPCR_data$Temp)

qPCR_data_group <- qPCR_data %>%
  group_by(Phage) %>%
  filter(Time == "5" & Temp == 42) %>%
  nest()

qpcr_models <- qPCR_data_group %>%
  mutate(model = map(data, ~ lme(DNA_conc ~ Competitor, random= ~1|Batch, data= ., method = "REML", weights = varIdent(form = ~1|Competitor)))) %>%
  mutate(model_summaries = map(model, summary))%>%
  mutate(model_pairwise = map(model, ~ {
    emmeans(.x, pairwise ~ Competitor)
  }))

# Diagnostic plots
length(qpcr_models$model)

plots <- lapply(1:3, function(i) plot(qpcr_models$model[[i]], main = paste(qpcr_models$Treatment[[i]], qpcr_models$Temp[[i]])))
cowplot::plot_grid(plotlist = plots)

# Generate model summaries
qpcr_models$model_summaries
qpcr_models$model_pairwise


# Plots - separate plots needed for each phage as different competitor order and different significance indicators

# Filter dataset
PEV2_qPCR <- qPCR_data %>%
  filter(Phage=="PEV2" & Time == "5" & Temp == 42)

LUZ19_qPCR <- qPCR_data %>%
  filter(Phage=="LUZ19" & Time == "5" & Temp == 42)

phage14_qPCR <- qPCR_data %>%
  filter(Phage=="phage14" & Time == "5" & Temp == 42)

# Order competitors

PEV2_qPCR$Competitor <- factor(PEV2_qPCR$Competitor, levels=c("No_competitor","LUZ19", "phage14","3P"), labels=c("No competitor","LUZ19","14-1","3-phage"))
LUZ19_qPCR$Competitor <- factor(LUZ19_qPCR$Competitor, levels=c("No_competitor","PEV2", "phage14","3P"), labels=c("No competitor","PEV2","14-1","3-phage"))
phage14_qPCR$Competitor <- factor(phage14_qPCR$Competitor, levels=c("No_competitor","PEV2", "LUZ19","3P"), labels=c("No competitor","PEV2","LUZ19","3-phage"))

#Plots

pev2_plot <- ggplot(data=PEV2_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(0.1,1000)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  scale_fill_manual(values = "#fc4e2a") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.position = "none")


luz19_plot <- ggplot(data=LUZ19_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(0.1,1000)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  scale_fill_manual(values = "#fc4e2a") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) + 
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.position = "none")


phage14_plot <- ggplot(data=phage14_qPCR, aes(Competitor, Fold_change_DNAconc_to_zero, fill = Temp))  + 
  geom_boxplot(width=0.5,notch = FALSE,  outlier.size = -1, color="black",lwd=1, alpha = 0.7,show.legend = F, varwidth=FALSE,position = position_dodge(width = .75))+
  ggbeeswarm::geom_quasirandom(shape = 21,size=2, dodge.width = .75, color="black",alpha=.5,show.legend = F)+
  ylab("Fold change in DNA copies relative to T0") +
  xlab("Phage competitor") +
  scale_y_continuous(trans='log10', limits = c(0.1,1000)) +
  annotation_logticks(sides="l")+
  theme_bw()+
  scale_fill_manual(values = "#fc4e2a") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.y=element_blank()) + 
  theme(axis.title=element_text(size=10,face="bold")) + 
  theme(legend.position = "none")


pev2_plot + luz19_plot + phage14_plot

ggsave("Fig_S7.tiff")

