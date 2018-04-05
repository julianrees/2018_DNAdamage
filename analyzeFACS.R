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
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 6)
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
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())


dfs <- list()
for (i in seq(length(normdata))){
  dfs[[i]] <- data.frame(fl = normdata[[i]], set = i)
}

ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = 'Control')) + 
  geom_rug(aes(x = fl, y = 0), position = position_jitter(height = 0))

ggplot(dfs[[1]], aes(fl)) + 
  geom_histogram(aes(fill = 'Control'), binwidth = 0.002)

ggplot(dfs[[1]], aes(fl)) + 
  geom_freqpoly(aes(fill = 'Control'), binwidth = 0.005)

# set the plotting options - alp is transparency, bw is the bandwidth multiplier
alp = 0.2
bw = 0.75
ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = 'Control'), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = 'High Dose'), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = 'Low Dose'), alpha = alp,  adjust = bw)# +
#facet_grid(~ set)

