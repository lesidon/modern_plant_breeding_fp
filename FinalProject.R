######### PREPARATION OF THE WORKING INTERFACE IN R ######################################
### I. Set working directory
#On RStudio: tab 'Session'-> Set Working Directory -> Choose Directory.
#Choose the directory containing the R script.

### II. Installation R packages needed for the analysis on RStudio:
#Click on the 'Packages' tab in the bottom-right window of R Studio interface->'Install Packages'
#Comment #1: R package installation requires a connection to internet
#Comment #2: Once packages have been installed, NO need to re-install them again 
# when you close and open again RStudio.

### III. Initialisation of the working space:clean workspace, clean memory 
# To erase all graphs
graphics.off()
# To erase objects from the working space - Clean up of the memory
rm(list = ls())
# garbage collector. To use to free memory
gc() 

################################################################################

## Loading of the R packages needed for the analysis.
library(lme4) 
library(ggplot2)
library(tidyverse)
library(lmerTest)
library(pracma)
library(multcomp)
library(lattice)
library(car)    # Levene test
library(agricolae)   # Newman-Keuls & Tukey test
library(ggplot2)
library(gridExtra)
library(dplyr)
library(MASS)#BoxCox

#Set working directory

maindir <- here::here()
setwd(maindir)

###################
# Data import in R
###################
#########

ayn16 <- read.table("ayn.byd0.16.csv", sep = ",", header = TRUE,
                    dec = ".", stringsAsFactors = TRUE)

ayn17 <- read.table("ayn.byd0.17.csv", sep = ",", header = TRUE, 
                    dec = ".", stringsAsFactors = TRUE)

ayn18 <- read.table("ayn.byd0.18.csv", sep = ",", header = TRUE, 
                    dec = ".", stringsAsFactors = TRUE)

ayn19 <- read.table("ayn.byd0.19.csv", sep = ",", header = TRUE, 
                    dec = ".", stringsAsFactors = TRUE)

ayn20 <- read.table("ayn.byd0.20.csv", sep = ",", header = TRUE, 
                    dec = ".", stringsAsFactors = TRUE)

#merge datasets

df <- c("ayn.byd0.16.csv", "ayn.byd0.17.csv", "ayn.byd0.18.csv", 
        "ayn.byd0.19.csv", "ayn.byd0.20.csv")

#create empty datafram to store data
merged_data <- data.frame()

#loop datasets from df
for (file in df) {
  data <- read.csv(file)
  merged_data <- bind_rows(merged_data,data)
}

#write csv for merged data
write.csv(merged_data, "merged_data.csv", row.names = FALSE)

str(merged_data)

merged_data$yield <- as.numeric(merged_data$yield)
merged_data$byd <- as.numeric(merged_data$byd)
merged_data$trt <- as.factor(merged_data$trt)
merged_data$location <- as.factor(merged_data$location)
merged_data$planting_date <- as.factor(merged_data$planting_date)
