########## Rice Project - Stable Isotope Mixing Model ##########
########## Author: Gen-Chang Hsu ##########
library(ggplot2)
library(plyr)
library(dplyr)
library(MixSIAR)
library(xlsx)
library(vegan)

### Define a ggplot theme
thm <- theme(axis.text.x = element_text(size = 12, color = "black"),
             axis.text.y = element_text(size = 12, color = "black"),
             axis.title.x = element_text(size = 15),
             axis.title.y = element_text(size = 15),
             plot.title = element_text(hjust = 0.5, size = 18),
             plot.margin = rep(unit(0.1,"null"),4),
             panel.grid.major.x = element_blank(),
             panel.grid.minor.x = element_blank(),
             panel.grid.major.y = element_blank(),
             panel.grid.minor.y = element_blank(),
             panel.background = element_rect(fill = NA),
             panel.border = element_rect(color = "black", fill = NA),
             strip.text = element_text(size = 15),
             legend.background = element_rect(color = "black"),
             strip.background = element_rect(color = "black"),
             legend.key = element_rect(fill = "transparent"))

### Load the data
SID <- read.csv("Output/Data_clean/SID.csv", header = T)
SID$C_conc <- SID$C.amount_ug/SID$Total.amount_mg/1000
SID$N_conc <- SID$N.amount_ug/SID$Total.amount_mg/1000

### Subset to different trophic guilds
# Predator
Predator.dat <- subset(SID, Family == "ARA"|
                            Family == "COC"|
                            Family == "TET")

# Rice herbivore
Rice.herb.dat <- subset(SID, Family == "DEL"|
                             Family == "CIC"|
                             Family == "PEN"|
                             Family == "ALY"|
                             Family == "LYG"|
                             Family == "PYR"|
                             Family == "HES"|
                             Family == "PYL")

# Tourist herbivore
Tour.herb.dat <- subset(SID, Family == "ACR"|
                             Family == "CHR")

# Detritivore
Detritivore.dat <- subset(SID, Family == "CHI"|
                               Family == "SCI"|
                               Family == "MUS"|
                               Family == "EPH"|
                               Family == "EMP"|
                               Family == "STR"|
                               Family == "CHL"|
                               Family == "TER")


# ### Check the multivariate differences between farm types for each prey source
# adonis(Rice.herb.dat[, c(3, 5)]~Rice.herb.dat$Farm.ID, method = "euclidean", permutations = 999)
# adonis(Tour.herb.dat[, c(3, 5)]~Tour.herb.dat$Farm.ID, method = "euclidean", permutations = 999)
# adonis(Detritivore.dat[, c(3, 5)]~Detritivore.dat$Farm.ID, method = "euclidean", permutations = 999)

### Predator sample sizes
Predator.dat %>% group_by(Farmtype, Stage) %>% summarise(N = n())


### Mixture data
Mixture <- data.frame(d13C = Predator.dat$d13C,
                      d15N = Predator.dat$d15N,
                      Farm = factor(Predator.dat$Farm.ID, levels = c("MC1", "MC2", "MC3", "MO1", "MO2", "MO3", "LC1", "LC2", "LC3", "LO1", "LO2", "LO3", "SO1", "SC1"), ordered = T),
                      Stage = factor(Predator.dat$Stage, levels = c("Seedling", "Tillering", "Flowering", "Ripening"), ordered = T))

write.csv(Mixture, "Output/Data_clean/Mixture.csv", row.names = FALSE)

# Mixture data setup
mix <- load_mix_data(filename = "Output/Data_clean/Mixture.csv",
                     iso_names = c("d13C","d15N"),
                     factors = c("Farm", "Stage"),
                     fac_random = c(F, F),
                     fac_nested = c(F, F),
                     cont_effects = NULL)

### Source data
Source <- data.frame(Source = factor(c(rep("Rice.herb", nrow(Rice.herb.dat)), rep("Tour.herb", nrow(Tour.herb.dat)), rep("Detritivore", nrow(Detritivore.dat))), levels = c("Rice.herb", "Tour.herb", "Detritivore"), ordered = T),
                     d13C = c(Rice.herb.dat$d13C, Tour.herb.dat$d13C, Detritivore.dat$d13C),
                     d15N = c(Rice.herb.dat$d15N, Tour.herb.dat$d15N, Detritivore.dat$d15N),
                     Concd13C = c(Rice.herb.dat$C_conc, Tour.herb.dat$C_conc, Detritivore.dat$C_conc),
                     Concd15N = c(Rice.herb.dat$N_conc, Tour.herb.dat$N_conc, Detritivore.dat$N_conc),
                     Farm = factor(c(as.character(Rice.herb.dat$Farm.ID), as.character(Tour.herb.dat$Farm.ID), as.character(Detritivore.dat$Farm.ID)), levels = c("MC1", "MC2", "MC3", "MO1", "MO2", "MO3", "LC1", "LC2", "LC3", "LO1", "LO2", "LO3", "SO1", "SC1"), ordered = T))

write.csv(Source, "Output/Data_clean/Source.csv", row.names = FALSE)

# Source data setup
source <- load_source_data(filename = "Output/Data_clean/Source.csv",
                           source_factors = NULL,
                           conc_dep = TRUE,
                           data_type = "raw",
                           mix)

### Discrimination factors: Diet-Dependent Discrimination Factor (DDDF) from Caut et al. (2009)
DF_C.fun <- function(x){-0.113*x - 1.916}
DF_N.fun <- function(x){-0.311*x + 4.065}

Discrimination <- ddply(.data = Source, .variables = "Source", summarize, Meand13C = mean(DF_C.fun(d13C)), SDd13C = sd(DF_C.fun(d13C)),
                                                                          Meand15N = mean(DF_N.fun(d15N)), SDd15N = sd(DF_N.fun(d15N)))
write.csv(Discrimination, "Output/Data_clean/Discrimination.csv", row.names = F)
discr <- load_discr_data(filename = "Output/Data_clean/Discrimination.csv", mix)


### Prey source biplot (TDF-corrected)
d13C_correct <- sapply(1:nrow(Source), function(x){Source[x, 2] + Discrimination[which(Discrimination$Source == Source[x, 1]), 2]})
d15N_correct <- sapply(1:nrow(Source), function(x){Source[x, 3] + Discrimination[which(Discrimination$Source == Source[x, 1]), 4]})
Source_correct <- data.frame(Source = Source$Source, d13C_correct, d15N_correct)
Source_correct <- ddply(Source_correct, "Source", summarise, Mean_d13C = mean(d13C_correct), Se_d13C = sd(d13C_correct)/sqrt(length(d13C_correct)), Mean_d15N = mean(d15N_correct), Se_d15N = sd(d15N_correct)/sqrt(length(d15N_correct)))

# Stable isotope signatures of primary producer
SID2017 <- read.xlsx("Data_raw/SID2017.xlsx", header = T, sheetIndex = 1)
SID2017.Rice <- subset(SID2017, Species == "Os")
SID2017.Rice$d_13C <- as.numeric(as.character(SID2017.Rice$d_13C))
SID2017.Rice$d_15N <- as.numeric(as.character(SID2017.Rice$d_15N))

Mean.d13C.rice <- mean(SID2017.Rice$d_13C)
SE.d13C.rice <- sd(SID2017.Rice$d_13C)/sqrt(length(SID2017.Rice$d_13C))
Mean.d15N.rice <- mean(SID2017.Rice$d_15N)
SE.d15N.rice <- sd(SID2017.Rice$d_15N)/sqrt(length(SID2017.Rice$d_15N))
Rice.all <- as.data.frame(cbind(Mean.d13C.rice, SE.d13C.rice, Mean.d15N.rice, SE.d15N.rice))

ggplot() +
  geom_point(data = Source_correct, aes(x = Mean_d13C, y = Mean_d15N, color = Source, shape = Source), size = 4) +
  geom_errorbar(data = Source_correct, aes(x = Mean_d13C, ymin = Mean_d15N-Se_d15N, ymax = Mean_d15N+Se_d15N, color = Source), width = 0.1, size = 1.2) +
  geom_errorbarh(data = Source_correct, aes(y = Mean_d15N, xmin = Mean_d13C-Se_d13C, xmax = Mean_d13C+Se_d13C, color = Source), height = 0.1, size = 1.2) +
  geom_point(data = Rice.all, aes(x = Mean.d13C.rice, y = Mean.d15N.rice), size = 1) +
  geom_errorbar(data = Rice.all, aes(x = Mean.d13C.rice, ymin = Mean.d15N.rice-SE.d15N.rice, ymax = Mean.d15N.rice+SE.d15N.rice), width = 0.1, size = 1.2) +
  geom_errorbarh(data = Rice.all, aes(y = Mean.d15N.rice, xmin = Mean.d13C.rice-SE.d13C.rice, xmax = Mean.d13C.rice+SE.d13C.rice), height = 0.1, size = 1.2) +
  ylim(5, 11) +
  xlim(-31, -21) +
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
        legend.position = c(0, 1.05),
        legend.spacing.x = unit(0.25, "cm"),
        legend.key.width = unit(0.7, "cm"),
        legend.key.size = unit(1.2, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.box.just = "center",
        legend.justification = c(-0.2, 1.2),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  labs(x = expression(paste(delta^{13}, "C (\u2030)", sep = "")), y = expression(paste(delta^{15}, "N (\u2030)", sep = ""))) +
  scale_color_manual(values=c("#00BA38", "#619CFF", "#993300"), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore"), name = "") +
  scale_shape_manual(values=c(15, 16, 17), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore"), name = "") +
  annotate(geom = "text", x = -28.3, y = 5.5, label = "Rice plant", size = 4.5)

ggsave("Output/Figures/Biplot.tiff", width = 6, height = 5, dpi = 600)


### Write JAGS file
dr <- getwd()
setwd(dir = paste0(dr, "/Output/data_clean"))
model_filename <- "JAGS.txt"
resid_err <- T
process_err <- T
write_JAGS_model(model_filename, resid_err, process_err, mix, source)

### Run JAGS file
jags <- run_model(run = "short", mix, source, discr, model_filename, alpha.prior = 1, resid_err, process_err)
setwd(dr)

### Evaluate JAGS output
dr <- getwd()
setwd(dir = paste0(dr, "/Output/data_clean"))
options(max.print = 1000000)
output_JAGS(jags, mix, source, output_options = list(summary_save = T,
                                                        summary_name = "Mix_mod_out",
                                                        sup_post = F,
                                                        plot_post_save_pdf = F,
                                                        plot_post_name = "Posterior_density",
                                                        sup_pairs = F,
                                                        plot_pairs_save_pdf = F,
                                                        plot_pairs_name = "Pairs_plot",
                                                        sup_xy = T,
                                                        plot_xy_save_pdf = F,
                                                        plot_xy_name = "xy_plot",
                                                        gelman = F,
                                                        heidel = F,
                                                        geweke = F,
                                                        diag_save = F,
                                                        diag_name = "Diagnostics",
                                                        indiv_effect = F,
                                                        plot_post_save_png = F,
                                                        plot_pairs_save_png = F,
                                                        plot_xy_save_png = F))


setwd(dr)

### Read model output summary data
Sum.stat <- read.table("Output/Data_clean/Mix_mod_out.txt", header = F, fill = TRUE)
Sum.stat <- Sum.stat[-c(1, 3, 4), ]
colnames(Sum.stat) <- c("ID", "Mean", "SD", "2.5%", "5%", "25%", "50%", "75%", "95%", "97.5%")
row.names(Sum.stat) <- NULL
Sum.stat <- Sum.stat[-1, ]

temp1 <- strsplit(as.character(Sum.stat$ID), "[.]")
temp2 <- ldply(temp1, .fun = function(x) c(x[2], x[3], x[4]))
temp2$V3 <- replace(temp2$V3, temp2$V3 == "Rice", "Rice.herb") %>%
            replace(temp2$V3 == "Tour", values = "Tour.herb")
temp2$Farmtype <- substr(temp2$V1, start = 2, stop = 2) %>%
  replace(substr(temp2$V1, start = 2, stop = 2) == "O", "Or") %>%
  replace(substr(temp2$V1, start = 2, stop = 2) == "C", "Cv")
temp2$Landscape <- substr(temp2$V1, start = 1, stop = 1) %>%
  replace(substr(temp2$V1, start = 1, stop = 1) == "M", "Mount") %>%
  replace(substr(temp2$V1, start = 1, stop = 1) == "L", "Land") %>%
  replace(substr(temp2$V1, start = 1, stop = 1) == "S", "Sea")

names(temp2)[1:3] <- c("Farm.ID", "Stage", "Source")
Summary <- cbind(temp2[, c("Farm.ID", "Landscape", "Farmtype", "Stage", "Source")], Sum.stat[, -1])
Summary$Farm.ID <- ordered(Summary$Farm.ID, levels = c("MO1", "MO2", "MO3", "MC1", "MC2", "MC3","LO1", "LO2", "LO3", "LC1", "LC2", "LC3", "SO1", "SC1"))
Summary$Farmtype <- ordered(Summary$Farmtype, levels = c("Or", "Cv"))
Summary$Landscape <- ordered(Summary$Landscape, levels = c("Mount", "Land", "Sea"))
Summary$Stage <- ordered(Summary$Stage, levels = c("Seedling", "Tillering", "Flowering", "Ripening"))
Summary$Source <- ordered(Summary$Source, levels = c("Rice.herb", "Tour.herb", "Detritivore"))
Summary <- arrange(Summary, Farm.ID, Stage, Source)

write.csv(Summary, "Output/Data_clean/Mix_mod_out_summary.csv", row.names = F)


### Figures for mixing model results
### Load the data
Summary <- read.csv('Output/Data_clean/Mix_mod_out_summary.csv', header = T)
Summary$Farm.ID <- ordered(Summary$Farm.ID, levels = c("MO1", "MO2", "MO3", "MC1", "MC2", "MC3","LO1", "LO2", "LO3", "LC1", "LC2", "LC3", "SO1", "SC1"))
Summary$Farmtype <- ordered(Summary$Farmtype, levels = c("Or", "Cv"))
Summary$Landscape <- ordered(Summary$Landscape, levels = c("Mount", "Land", "Sea"))
Summary$Stage <- ordered(Summary$Stage, levels = c("Seedling", "Tillering", "Flowering", "Ripening"))
Summary$Source <- ordered(Summary$Source, levels = c("Rice.herb", "Tour.herb", "Detritivore"))

### Dietary proportions of predators by farm
ggplot(Summary, aes(x=Stage, y=Mean, group=Source, color=Source, shape=Source)) +
  geom_line(position = position_dodge(0.1), size = 1.2) +
  geom_point(aes(x = Stage, y = Mean), position = position_dodge(0.1), size = 3)+
  facet_wrap(~Farm.ID, nrow = 3, ncol = 6) +
  geom_errorbar(aes(ymin = ifelse((Mean-SD)>0, Mean-SD, 0), ymax = ifelse((Mean+SD)<1, Mean+SD, 1)), width = 0.2, position = position_dodge(0.1), size = 1) +
  ylim(0, 1) +
  labs(title = "Bayesian Posterior Estimates") +
  xlab("Stage") +
  ylab(paste("Proportion", "(mean", "\u00B1", "SD)", spe = "")) +
  theme(axis.text.x = element_text(size = 8, color = "black", angle = 45, hjust = 0.75),
        axis.text.y = element_text(size = 8, color = "black"),
        axis.title.x = element_text(size = 15, margin = margin(t = -10)),
        axis.title.y = element_text(size = 15, margin = margin(r = 8)),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "cm"),
        panel.background = element_rect(fill = NA),
        panel.border = element_rect(colour = "black", fill = NA, size = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "grey"),
        panel.grid.minor.y = element_blank(),
        strip.background = element_rect(color = "black"),
        legend.position = c(0.65, 0.15),
        legend.spacing.x = unit(1, "cm"),
        legend.key.width = unit(0.8, "cm"),
        legend.key.size = unit(1.2, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 10),
        legend.direction = "horizontal",
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "black")) +
  scale_color_manual(values=c("#00BA38", "#619CFF", "#993300"), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore"))+
  scale_shape_manual(values=c(15, 16, 17), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore")) +
  guides(colour = guide_legend(title.position = "top"))+
  guides(shape = guide_legend(title.position = "top"))

# ggsave("Output/Figures/Diet_prop_by_farm.tiff", width = 11, height = 6, dpi = 600)

### Organic vs Conventional
Summary2 <- ddply(Summary, c("Farmtype", "Stage", "Source"), summarise, mean = mean(Mean), Sd = sd(Mean), n = length(Mean), Se = Sd/sqrt(n))
Summary2$Farmtype2 <- as.character(Summary2$Farmtype) %>%
  replace(as.character(Summary2$Farmtype) == "Or", "Organic") %>%
  replace(as.character(Summary2$Farmtype) == "Cv", "Conventional")
Summary2$Farmtype2 <- ordered(Summary2$Farmtype2, levels = c("Organic", "Conventional"))
write.csv(Summary2, "Output/Data_clean/Mix_mod_out_summary2.csv")

Labs <- data.frame(Farmtype2 = c("Organic", "Conventional"),
                   x = c(2, 2),
                   y = c(1, 1),
                   Source = c("Rice.herb", "Rice.herb"))

ggplot(Summary2[Summary2$Stage != "Seedling", ] , aes(x = Stage, y = mean, color = Source, shape = Source, group = Source)) +
  geom_point(size = 3, position = position_dodge(0.1)) +
  geom_line(size = 1, position = position_dodge(0.1)) +
  facet_grid(~Farmtype2) +
  geom_errorbar(aes(ymin = ifelse((mean-Se)>0, mean-Se, 0), ymax = ifelse((mean+Se)<1, mean+Se, 1)), width = 0.2, position = position_dodge(0.1), size = 0.6) +
  xlab("Crop stage") +
  ylab(paste("Proportion of prey sources in \n predators' diet", " (Mean ", "\u00B1", " SE)", sep = "")) +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10, color = "black"),
        axis.text.y = element_text(size = 10, color = "black"),
        axis.title.x = element_text(size = 12.5, margin = margin(t = 12)),
        axis.title.y = element_text(size = 12.5, margin = margin(r = 10)),
        plot.margin = unit(c(0.5, 0.2, 0.2, 0.2), "cm"),
        panel.background = element_rect(fill = NA),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        strip.background = element_rect(fill = NA, color = "transparent"),
        strip.text.x = element_text(size = 10),
        legend.position = c(0.5, 1.15),
        legend.direction = "horizontal",
        legend.spacing.x = unit(0.2, "cm"),
        legend.key.width = unit(1, "cm"),
        legend.key.size = unit(0, "line"),
        legend.key = element_blank(),
        legend.text = element_text(size = 8.5),
        legend.title = element_text(size = 10),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "transparent", size = 0.5, linetype = "solid", colour = "transparent")) +
  scale_color_manual(values=c("#00BA38", "#619CFF", "#993300"), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore"), name = "") +
  scale_shape_manual(values=c(15, 16, 17), labels = c("Rice herbivore", "Tourist herbivore", "Detritivore"), name = "") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf)

ggsave("Output/Figures/Diet_prop_Or_Cv.tiff", width = 6, height = 3.5, dpi = 600)



