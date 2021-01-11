########## Rice Project - Guild Assignment ##########
########## Author: Gen-Chang Hsu ##########
library(ggplot2)
library(plyr)

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

### Exclude predators
SID2 <- subset(SID, Family != "TET"&
                    Family != "ARA"&
                    Family != "COC"&
                    Family != "STA"&
                    Family != "DOL"&
                    Family != "CAR"&
                    Family != "GEO"&
                    Family != "ICH", drop=TRUE)
SID2$Family <- as.factor(as.character(SID2$Family))

### Kmeans clustering
set.seed(1)
kmeans <- kmeans(SID2[, c("d13C", "d15N")] , centers = 3)
SID2$kmeans <- as.factor(kmeans$cluster)

### SI biplot of prey sources
ggplot(data = SID2, aes(x = d13C, y = d15N, color = kmeans)) +
  geom_point() +
  geom_text(aes(label = as.character(Family)), show.legend = F) +
  xlab(expression(paste(delta^'13'*C, sep = ""))) +
  ylab(expression(paste(delta^'15'*N, sep = ""))) +
  scale_color_discrete(name = "Group") +
  thm +
  theme(legend.position=c(1, 1),
        legend.justification=c(1.2, 1.2))

# ggsave("Output/Figures/Source_biplot.tiff", width = 7, height = 5, dpi = 600)

### Numbers of capsules in each cluster for each prey family
Group.dat <- tapply(SID2$Family, INDEX = SID2$kmeans, table)
Group.dat <- unlist(Group.dat)
Group.dat <- data.frame(Group = substr(names(Group.dat), 1, 1),
                        Family = substr(names(Group.dat), 3, 5),
                        N = Group.dat)
Prop.dat <- ddply(Group.dat, "Family", transform, Prop = N/sum(N))

row.names(Group.dat) <- NULL
Group.dat$X <- paste(c(rep("", 18), rep(" ", 18), rep("  ", 18)), Group.dat$Family, sep = "")
Group.dat$X <- factor(Group.dat$X, levels = Group.dat[order(Group.dat$N, decreasing = T), 'X'])
Group.dat$X <- ordered(Group.dat$X)

ggplot(data = Group.dat, aes(x = X, y = N)) +
  geom_col() +
  facet_wrap(~Group, scales = "free_x") +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 70)) +
  xlab("Family") +
  ylab("Number of Capsules") +
  thm +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        plot.margin = rep(unit(0.1, "inch"),4))

# ggsave("Output/Figures/Capsules.tiff", width = 8, height = 4, dpi = 600)

### Proportions of the numbers of capsules in each cluster for each prey family
ggplot(data = Prop.dat, aes(x = Family, y = Prop, fill = Group)) +
  geom_bar(stat = "identity") +
  xlab("Family") +
  ylab("Proportion") +
  scale_y_continuous(expand = c(0, 0), limit = c(0, 1)) +
  thm +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        plot.margin = rep(unit(0.1, "inch"),4))

# ggsave("Output/Figures/Proportion.tiff", width = 6, height = 4, dpi = 600)






