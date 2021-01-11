########## Rice Project - Beta Regressions on the Rice Herbivore Consumption by Predators ##########
########## Author: Gen-Chang Hsu ###########
library(ggplot2)
library(emmeans)
library(betareg)
library(multcompView)
library(multcomp)
library(lmtest)
library(glmmTMB)
library(visreg)
library(car)
library(afex)
library(reshape2)
library(dplyr)
library(ggpubr)

### Define a ggplot theme
thm <- theme(axis.text.x = element_text(size = 12, color = "black"),
             axis.text.y = element_text(size = 12, color = "black"),
             axis.title.x = element_text(size = 15, margin = margin(t = 12)),
             axis.title.y = element_text(size = 15, margin = margin(r = 10, l = 0)),
             plot.title = element_text(hjust = 0.5, size = 18),
             plot.margin = unit(c(0.1, 0, 0.1, 0.1),"null"),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = NA),
             panel.border = element_rect(color = "black", fill = NA),
             strip.text = element_text(size = 15),
             legend.background = element_rect(color = "transparent"),
             strip.background = element_rect(color = "black"),
             legend.key = element_rect(fill = "transparent"),
             legend.key.height = unit(0.65, "cm"),
             legend.title.align = 0.5)


### Load the data
Summary <- read.csv("Output/Data_clean/Mix_mod_out_summary.csv", header = T)
Summary$Farm.ID <- ordered(Summary$Farm.ID, levels = c("MO1", "MO2", "MO3", "MC1", "MC2", "MC3","LO1", "LO2", "LO3", "LC1", "LC2", "LC3", "SO1", "SC1"))
Summary$Farmtype <- ordered(Summary$Farmtype, levels = c("Or", "Cv"))
Summary$Landscape <- ordered(Summary$Landscape, levels = c("Mount", "Land", "Sea"))
Summary$Stage <- ordered(Summary$Stage, levels = c("Seedling", "Tillering", "Flowering", "Ripening"))
Summary$Source <- ordered(Summary$Source, levels = c("Rice.herb", "Tour.herb", "Detritivore"))
Prop.Rice.herb <- subset(Summary, Source == "Rice.herb")

Source.Rel.Abd <- read.csv("Output/Data_clean/Diet_prop_by_farm.csv", header = T, row.names = 1)
Rel.Abd_Rice.herb <- subset(Source.Rel.Abd, Trophic == "Rice.herb")
Rel.Abd_Rice.herb$Farm.ID <- gsub("-", "", Rel.Abd_Rice.herb$Farm.ID)

index <- sapply(1:nrow(Rel.Abd_Rice.herb), function(x){
  which(Rel.Abd_Rice.herb[x, ]$Farm.ID == Prop.Rice.herb$Farm.ID &
        as.character(Rel.Abd_Rice.herb[x, ]$Stage) == Prop.Rice.herb$Stage)
})

empty <- which(1:41 %in% unlist(index) == F)
index2 <- which(index !=0)
Prop.Rice.herb$Rel.Abd <- NA
Prop.Rice.herb$Rel.Abd[-empty] <- Rel.Abd_Rice.herb$Rel.Abd[index2]

# Assign proportion of 1.0 to 0.999
Prop.Rice.herb$Rel.Abd[which(Prop.Rice.herb$Rel.Abd == 1)] <- 0.999


### Beta regression model 1: Farm type * crop stage
beta.out <- betareg(Mean ~ Farmtype * Stage, data = Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",], type = "ML")
summary(beta.out)
lrtest(beta.out)

# Anova table
Anova(beta.out, test.statistic = "F")

# Model diagnostics
plot(scale(beta.out$fitted.values), scale(beta.out$residuals))
summary(lm(scale(beta.out$residuals)~scale(beta.out$fitted.values)))
qqnorm(scale(beta.out$residuals))

# Estimated marginal means (EMMs) with multiple comparisons
refgrid <- ref_grid(beta.out)
summary(refgrid)
test(refgrid)
confint(refgrid, adjust = "tukey")                    # Make multiple comparisons among all grids
confint(refgrid, adjust = "tukey", by = "Farmtype")   # Make multiple comparisons between "Stage" split by "Farmtype"
confint(refgrid, adjust = "tukey", by = "Stage")      # Make multiple comparisons between "Farmtype" split by "Stage"

EMM <- emmeans(beta.out, c("Farmtype", "Stage"), data = Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",])
CI.adj <- cld(EMM, Letters = letters)

# Plot the Tukey-adjusted C.I.
P1 <- ggplot(data = as.data.frame(CI.adj), aes(x = Stage, y = emmean, fill = Farmtype)) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(shape = 21, size = 5, position = position_dodge(width = 0.5)) +
  scale_fill_manual(values = c("black", "white"), label = c("Organic", "Conventional"), name = NULL) +
  ylab(paste("Proportion of rice herbivores in \n predators' diet", " (Mean ", "\u00B1", " 95% CI)", sep = "")) +
  xlab("Crop stage") +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1.1)) +
  thm +
  theme(legend.position = c(0.8, 0.3),
        legend.background = element_rect(color = "transparent")) +
  geom_text(data = as.data.frame(CI.adj), aes(x = Stage, y = asymp.UCL+0.05, label = gsub(" ", "", .group, fixed = TRUE)), position = position_dodge(width = 0.5))

P1
ggsave("Output/Figures/Beta_reg_farmtype_cropstage.tiff", width = 6, height = 5, dpi = 600)


### Beta regression model 2: relative abundance of rice herbivores
Beta.all <- betareg(Mean ~ Rel.Abd, data = Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",], type = "ML")
Beta.Or <- betareg(Mean ~ Rel.Abd, data = subset(Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",], Farmtype == "Or"), type = "ML")
Beta.Cv <- betareg(Mean ~ Rel.Abd, data = subset(Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",], Farmtype == "Cv"), type = "ML")

summary(Beta.all)
summary(Beta.Or)
summary(Beta.Cv)

Beta.all.fit <- predict(Beta.all, data.frame(Rel.Abd = seq(0.01, 0.99, 0.05)))
Beta.Or.fit <- predict(Beta.Or, data.frame(Rel.Abd = seq(0.01, 0.99, 0.05)))
Beta.Cv.fit <- predict(Beta.Cv, data.frame(Rel.Abd = seq(0.01, 0.99, 0.05)))
fit.dat <- data.frame(Beta.all.fit, Beta.Or.fit, Beta.Cv.fit, x = seq(0.01, 0.99, 0.05))
fit.dat2 <- melt(fit.dat, id.vars = "x")

# Plot the regression lines
P2 <- ggplot() +
  geom_point(data = Prop.Rice.herb[Prop.Rice.herb$Stage != "Seedling",], aes(x = Rel.Abd, y = Mean, shape = Farmtype), size = 2) +
  geom_line(data = filter(fit.dat2, variable != "Beta.all.fit"), aes(x = x, y = value, group = variable, linetype = variable), size = 1) +
  geom_line(data = filter(fit.dat2, variable == "Beta.all.fit"), aes(x = x, y = value), size = 1, linetype = 4) +
  ylim(0, 1) +
  xlim(0, 1) +
  thm +
  theme(legend.background = element_rect(color = "transparent", fill = "transparent"),
        legend.position = c(0.225, 0.9),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.size = unit(0.4, "line"),
        legend.text = element_text(size = 12)) +
  xlab("Relative abundance of rice herbivores") +
  ylab("Proportion of rice herbivores in \n predators' diet") +
  scale_shape_manual(values = c(16, 21), label = c("Organic", "Conventional"), name = "") +
  scale_linetype_manual(values = c(1, 5), label = c("Organic", "Conventional"), name = "") +
  guides(linetype = guide_legend(override.aes = list(linetype = c(1, 2))))

P2
ggsave("Output/Figures/Beta_reg_rel_abd.tiff", width = 6, height = 5, dpi = 600)

# # Arrange the two plots
# ggarrange(P1, P2, labels = c("(a)", "(b)"), label.x = 0.2)
# ggsave("Output/Figures/Beta_reg.tiff", width = 10, height = 4.5, dpi = 600)

