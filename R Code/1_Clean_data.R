########## Rice Project - Data Cleaning ##########
########## Author: Gen-Chang Hsu ##########
library(xlsx)
library(dplyr)

### Load and clean the data
SID.1st <- read.xlsx("Stable Isotope Data - 1st Survey.xls", sheetIndex = 1, startRow = 7, header = T)
SID.1st <- SID.1st[, c(1, 2, 3, 5, 6, 14)]
SID.1st <- SID.1st[complete.cases(SID.1st), ]
names(SID.1st) <- c("Sample.ID", "d13C", "C.amount_ug", "d15N", "N.amount_ug", "Total.amount_mg")
SID.1st$Sample.ID <- substr(SID.1st$Sample.ID, start = 1, stop = 14)
SID.1st$Stage <- rep("Seedling", nrow(SID.1st))
SID.1st$Family <- substr(SID.1st$Sample.ID, start = 9, stop = 11)
SID.1st$Farmtype <- substr(SID.1st$Sample.ID, start = 13, stop = 13) %>%
  replace(substr(SID.1st$Sample.ID, start = 13, stop = 13) == "O", "Or") %>%
  replace(substr(SID.1st$Sample.ID, start = 13, stop = 13) == "C", "Cv")
SID.1st$Landscape <- substr(SID.1st$Sample.ID, start = 12, stop = 12) %>%
  replace(substr(SID.1st$Sample.ID, start = 12, stop = 12) == "M", "Mount") %>%
  replace(substr(SID.1st$Sample.ID, start = 12, stop = 12) == "L", "Land") %>%
  replace(substr(SID.1st$Sample.ID, start = 12, stop = 12) == "S", "Sea")
SID.1st$Date <- substr(SID.1st$Sample.ID, start = 1, stop = 8)
SID.1st$Farm.ID <- substr(SID.1st$Sample.ID, start = 12, stop = 14)

SID.2nd <- read.xlsx("Stable Isotope Data - 2nd Survey.xls", sheetIndex = 1, startRow = 7, header = T)
SID.2nd <- SID.2nd[, c(1, 2, 3, 5, 6, 14)]
SID.2nd <- SID.2nd[complete.cases(SID.2nd), ]
names(SID.2nd) <- c("Sample.ID", "d13C", "C.amount_ug", "d15N", "N.amount_ug", "Total.amount_mg")
SID.2nd$Sample.ID <- substr(SID.2nd$Sample.ID, start = 1, stop = 14)
SID.2nd$Stage <- rep("Tillering", nrow(SID.2nd))
SID.2nd$Family <- substr(SID.2nd$Sample.ID, start = 9, stop = 11)
SID.2nd$Farmtype <- substr(SID.2nd$Sample.ID, start = 13, stop = 13) %>%
  replace(substr(SID.2nd$Sample.ID, start = 13, stop = 13) == "O", "Or") %>%
  replace(substr(SID.2nd$Sample.ID, start = 13, stop = 13) == "C", "Cv")
SID.2nd$Landscape <- substr(SID.2nd$Sample.ID, start = 12, stop = 12) %>%
  replace(substr(SID.2nd$Sample.ID, start = 12, stop = 12) == "M", "Mount") %>%
  replace(substr(SID.2nd$Sample.ID, start = 12, stop = 12) == "L", "Land") %>%
  replace(substr(SID.2nd$Sample.ID, start = 12, stop = 12) == "S", "Sea")
SID.2nd$Date <- substr(SID.2nd$Sample.ID, start = 1, stop = 8)
SID.2nd$Farm.ID <- substr(SID.2nd$Sample.ID, start = 12, stop = 14)

SID.3rd <- read.xlsx("Stable Isotope Data - 3rd Survey.xls", sheetIndex = 1, startRow = 7, header = T)
SID.3rd <- SID.3rd[, c(1, 2, 3, 5, 6, 14)]
SID.3rd <- SID.3rd[complete.cases(SID.3rd), ]
names(SID.3rd) <- c("Sample.ID", "d13C", "C.amount_ug", "d15N", "N.amount_ug", "Total.amount_mg")
SID.3rd$Sample.ID <- substr(SID.3rd$Sample.ID, start = 1, stop = 14)
SID.3rd$Stage <- rep("Flowering", nrow(SID.3rd))
SID.3rd$Family <- substr(SID.3rd$Sample.ID, start = 9, stop = 11)
SID.3rd$Farmtype <- substr(SID.3rd$Sample.ID, start = 13, stop = 13) %>%
  replace(substr(SID.3rd$Sample.ID, start = 13, stop = 13) == "O", "Or") %>%
  replace(substr(SID.3rd$Sample.ID, start = 13, stop = 13) == "C", "Cv")
SID.3rd$Landscape <- substr(SID.3rd$Sample.ID, start = 12, stop = 12) %>%
  replace(substr(SID.3rd$Sample.ID, start = 12, stop = 12) == "M", "Mount") %>%
  replace(substr(SID.3rd$Sample.ID, start = 12, stop = 12) == "L", "Land") %>%
  replace(substr(SID.3rd$Sample.ID, start = 12, stop = 12) == "S", "Sea")
SID.3rd$Date <- substr(SID.3rd$Sample.ID, start = 1, stop = 8)
SID.3rd$Farm.ID <- substr(SID.3rd$Sample.ID, start = 12, stop = 14)

SID.4th <- read.xlsx("Stable Isotope Data - 4th Survey.xls", sheetIndex = 1, startRow = 7, header = T)
SID.4th <- SID.4th[, c(1, 2, 3, 5, 6, 14)]
SID.4th <- SID.4th[complete.cases(SID.4th), ]
names(SID.4th) <- c("Sample.ID", "d13C", "C.amount_ug", "d15N", "N.amount_ug", "Total.amount_mg")
SID.4th$Sample.ID <- substr(SID.4th$Sample.ID, start = 1, stop = 14)
SID.4th$Stage <- rep("Ripening", nrow(SID.4th))
SID.4th$Family <- substr(SID.4th$Sample.ID, start = 9, stop = 11)
SID.4th$Farmtype <- substr(SID.4th$Sample.ID, start = 13, stop = 13) %>%
  replace(substr(SID.4th$Sample.ID, start = 13, stop = 13) == "O", "Or") %>%
  replace(substr(SID.4th$Sample.ID, start = 13, stop = 13) == "C", "Cv")
SID.4th$Landscape <- substr(SID.4th$Sample.ID, start = 12, stop = 12) %>%
  replace(substr(SID.4th$Sample.ID, start = 12, stop = 12) == "M", "Mount") %>%
  replace(substr(SID.4th$Sample.ID, start = 12, stop = 12) == "L", "Land") %>%
  replace(substr(SID.4th$Sample.ID, start = 12, stop = 12) == "S", "Sea")
SID.4th$Date <- substr(SID.4th$Sample.ID, start = 1, stop = 8)
SID.4th$Farm.ID <- substr(SID.4th$Sample.ID, start = 12, stop = 14)

SID <- rbind(SID.1st, SID.2nd, SID.3rd, SID.4th)
SID <- SID[, c(1, 8, 2:6, 9, 10, 7, 12, 11)]
write.csv(SID, "SID.csv", row.names = F)




