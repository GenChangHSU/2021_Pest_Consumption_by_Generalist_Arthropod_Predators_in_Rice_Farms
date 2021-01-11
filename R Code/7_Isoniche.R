########## Rice Project - Isotopic Niche of Predators ##########
########## Author: Gen-Chang Hsu ###########
library(ggplot2)
library(vegan)
library(dplyr)
library(ggpubr)

### Load the data
SID <- read.csv("Output/Data_clean/SID.csv", header = T)
SID$Farmtype <- ordered(SID$Farmtype, levels = c("Or", "Cv"))
SID$Landscape <- ordered(SID$Landscape, levels = c("Mount", "Land", "Sea"))
SID$Stage <- ordered(SID$Stage, levels = c("Seedling", "Tillering", "Flowering", "Ripening"))


### PERMANOVA to test the difference in isotopic signatures
Predator.SID <- subset(SID, Family == "ARA"|
                            Family == "COC"|
                            Family == "TET")
PER.dat <- Predator.SID[, c("d13C", "d15N", "Farmtype", "Stage")]
PER.dat2 <- subset(PER.dat, Stage != "Seedling")

# Euclidean distance
Dis.mat <- vegdist(PER.dat2[, 1:2], method = "eu")
set.seed(111)
PER.out <- adonis(Dis.mat ~ Farmtype*Stage, data = PER.dat2[, 3:4], perm = 9999)


### Test for homogeneous multivariate dispersions
dps.out.Farmtype <- betadisper(Dis.mat, group = PER.dat2$Farmtype, type = "centroid")
dps.out.Stage <- betadisper(Dis.mat, group = PER.dat2$Stage, type = "centroid")
permutest(dps.out.Farmtype, pairwise = F)
permutest(dps.out.Stage, pairwise = F)


### Mean and SD of distance-to-centroid in organic and conventional farms
Dist_to_centroid <- data.frame(Distance = dps.out.Farmtype$distances, Farmtype = dps.out.Farmtype$group)
Dist_to_centroid %>% group_by(Farmtype) %>% summarise(Mean = mean(Distance), SD = sd(Distance), n = n())
t.test(Distance~Farmtype, data = Dist_to_centroid)


### Isotopic niche of predators
# 1. Farm type + Crop stage
ggplot(Predator.SID[Predator.SID$Stage != "Seedling", ], aes(x = d13C, y = d15N, color = Stage, shape = Farmtype)) +
  geom_point(size = 1.25, alpha = 0.5, stroke = 2) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.715, 0.775),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.box.just = "center",
        legend.justification = c(0, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), y = expression(paste(delta^{15}, "N (\u2030)", sep = "")))+
  stat_ellipse(aes(lty = Farmtype), type = "t", level = 0.5, size = 0.8) +
  scale_linetype_manual(values = c(1, 2), labels = c("Organic", "Conventional")) +
  scale_color_manual(values = c("red2", "blue", "green3", "black"), label = c("Tillering", "Flowering", "Ripening")) +
  scale_shape_manual(values = c(16, 6), labels = c("Organic", "Conventional")) +
  guides(color = guide_legend(title = "Crop stage", override.aes = list(shape = NA)),
         linetype = guide_legend(title = "Farm type"),
         shape = guide_legend(title = "Farm type")) +
  scale_y_continuous(limits = c(5, 13), breaks = seq(6, 12, 2), labels = c("6.0", "8.0", "10.0", "12.0"))

# 2.Farm type only
P_farmtype <- ggplot(Predator.SID[Predator.SID$Stage != "Seedling", ], aes(x = d13C, y = d15N, shape = Farmtype)) +
  geom_point(size = 2.5, alpha = 0.5, stroke = 1) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.55, 0.9),
        legend.spacing.x = unit(0.2, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(1.5, "cm"),
        legend.key.height = unit(0.65, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.box.just = "center",
        legend.justification = c(0, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), y = expression(paste(delta^{15}, "N (\u2030)", sep = "")))+
  stat_ellipse(aes(lty = Farmtype), type = "t", level = 0.5, size = 0.8) +
  scale_linetype_manual(values = c(1, 2), labels = c("Organic", "Conventional")) +
  scale_shape_manual(values = c(16, 1), labels = c("Organic", "Conventional")) +
  guides(linetype = guide_legend(title = ""),
         shape = guide_legend(title = "")) +
  scale_y_continuous(limits = c(5, 13), breaks = seq(6, 12, 2), labels = c("6.0", "8.0", "10.0", "12.0"))

# 3. Crop stage only
P_cropstage <- ggplot(Predator.SID[Predator.SID$Stage != "Seedling", ], aes(x = d13C, y = d15N, color = Stage)) +
  geom_point(size = 2, alpha = 0.5, stroke = 1) +
  theme(axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = 10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.title = element_text(hjust = 0.5, size = 18),
        plot.margin = rep(unit(0.05,"null"), 4),
        panel.background = element_rect(fill = "transparent"),
        plot.background = element_rect(colour = "transparent"),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        legend.position = c(0.6, 0.86),
        legend.spacing.x = unit(0.1, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        legend.key.width = unit(1.3, "cm"),
        legend.key.size = unit(1, "line"),
        legend.key.height = unit(0.65, "cm"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.box.just = "center",
        legend.justification = c(0, 0.5),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), y = expression(paste(delta^{15}, "N (\u2030)", sep = "")))+
  stat_ellipse(type = "t", level = 0.5, size = 0.8) +
  scale_color_manual(values = c("#d95f02", "#7570b3", "#1b9e77", "black"), label = c("Tillering", "Flowering", "Ripening"), name = "") +
  scale_y_continuous(limits = c(5, 13), breaks = seq(6, 12, 2), labels = c("6.0", "8.0", "10.0", "12.0"))



ggarrange(P_farmtype, P_cropstage, labels = c("(a)", "(b)"),
          label.x = 0.1)

ggsave("Isospace2.tiff", width = 10, height = 4.5, dpi = 600)

