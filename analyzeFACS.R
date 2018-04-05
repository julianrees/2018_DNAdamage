#---- Header content ----
library(ggplot2)
library(BiocInstaller)
library(ggcyto)
library(flowCore)
library(outliers)

# source("https://bioconductor.org/biocLite.R")
# biocLite()

setwd('./')
#---- Data import ----
# Set 1 data, Control, High, Low, Set 2 data, same order
rdata <- list()
rdata[[1]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S1 CTRL.csv")[,3])
rdata[[2]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S1 HD.csv")[,3])
rdata[[3]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S1 LD.csv")[,3])
rdata[[4]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S2 CTRL.csv")[,3])
rdata[[5]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S2 HD.csv")[,3])
rdata[[6]] <- c(read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S2 LD.csv")[,3])

#---- Data processing ----
logdata <- list()
logmedians <- list()
for (i in seq(length(rdata))){
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 4)
  logmedians[[i]] <- median(logdata[[i]])
}

# make the average of the medians of the control data sets
ctrl_abs_mean <- mean(c(logmedians[[1]], logmedians[[4]]))
controls <- as.list(c(1,1,1,4,4,4))

# adjust the data sets to the absolute mean using the median absolute deviation
normdata <- list()
logmads <- list()
normfactors <- list()
temp_normdata <- list()
for (i in seq(length(logdata))){
  logmads[[i]] <- mad(logdata[[i]], constant = ctrl_abs_mean)
  temp_normdata[[i]] <- logdata[[i]] - logmads[[controls[[i]]]]
  normfactors[[i]] <- median(temp_normdata[[controls[[i]]]])
  normdata[[i]] <- temp_normdata[[i]] / normfactors[[i]]
}



#---- Plotting the data ----

dfs <- list()
for (i in seq(length(normdata))){
  dfs[[i]] <- data.frame(fl = normdata[[i]])
}

ggplot(dfs[[1]], aes(x = fl)) + 
  geom_bar(aes(fill = 'Control')) + 
  geom_bar(data = dfs[[2]], aes(x = fl, fill = 'High Dose')) + 
  geom_bar(data = dfs[[3]], aes(x = fl, fill = 'Low Dose'))

