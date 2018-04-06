#---- Header content ----
library(ggplot2)
library(BiocInstaller)
library(ggcyto)
library(flowCore)
library(outliers)

# source("https://bioconductor.org/biocLite.R")
# biocLite()

setwd('./')
#---- Data import from the listed folders ----
folders <- c("Data/4h_gH2aX_BT549/", "Data/24h_gH2aX_BT549/")

temp_rdata <- list()
rdata <- list()
dataset_name <- list()
k = 0
for (i in seq(length(folders))){
  files <- dir(folders[i])
  for (j in seq(length(files))){
    rdata[[j+k]] <- c(read.csv(paste(folders[i], files[j], sep = ''))[,3])
    dataset_name[[j+k]] <- files[j]
  }
  k = k + j
  # NEED TO CHECK THIS
}


#---- Data processing ----
logdata <- list()
logmedians <- list()
for (i in seq(length(rdata))){
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 6)
  logmedians[[i]] <- median(logdata[[i]])
}

# make the average of the medians of the control data sets


# NEED TO FIRST SPLIT BY CELL LINE, THEN BY TIME POINT FOR THE CONTROLS
# COMPARE TREATMENTS TO CONTROLS, ALSO HAVE DIFFERENT TARGETS

controls <- c(1,1,1,4,4,4,7,7,7,10,10,10,13,13,13)
# timepoint 1 is 4 hrs, 2 is 24 hrs, 3 is 4_24, and 4 is 24_24
timepoints <- c('4 hrs','24 hrs','4 hrs + 24 hrs','24 hrs + 24 hrs')
timepoint <- c(1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
# antibody 1 is gH2aX, 2 is ATF2
antibodies <- c('gH2aX','ATF2')
antibody <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
# cellline 1 is BT549, 2 is SKBR3
celllines <- c('BT549','SKBR3')
cellline <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)


# select which controls to include in MAD via logical bit
includecontrols <- c(1,1,1,0,0,0,1,1,1,1,1,1,0,0,0)
includecontrols <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)


normdata <- list()
logmads <- list()
maddata <- list()
ctrl_abs_mean <- array()
normfactors <- array()

# adjusts the datasets to the control absolute means for each time point. may need to do this differently
for (j in seq(max(timepoint))){
  ctrl_abs_mean[j] <- mean(as.numeric(logmedians[as.numeric(unique(controls[which(
    timepoint == j & includecontrols == 1)]))]))
  for (i in seq(length(logdata))){
    logmads[[i]] <- mad(logdata[[i]], constant = ctrl_abs_mean[j])
    maddata[[i]] <- logdata[[i]] - logmads[[controls[[i]]]]
    normfactors[i] <- median(maddata[[controls[[i]]]])
    normdata[[i]] <- maddata[[i]] / normfactors[i]
  }
}

ctrl_abs_mean

#---- Make the log, MAD and normalized dataframes for plotting ----
r_dfs <- list()
log_dfs <- list()
mad_dfs <- list()
dfs <- list()
for (i in seq(length(maddata))){
  r_dfs[[i]] <- data.frame(fl = rdata[[i]], set = dataset_name[[i]])
  log_dfs[[i]] <- data.frame(fl = logdata[[i]], set = dataset_name[[i]])
  mad_dfs[[i]] <- data.frame(fl = maddata[[i]], set = dataset_name[[i]])
  dfs[[i]] <- data.frame(fl = normdata[[i]], set = dataset_name[[i]])
}


#---- generate geometric means of data from corresponding sets ----



#---- Plotting the data ----
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())
alp = 0.2 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

# compare the control groups in log space and not normalized
ggplot(log_dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = log_dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_vline(xintercept = ctrl_abs_mean, aes(color))


ggplot(mad_dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = mad_dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_vline(xintercept = unique(normfactors))


ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw)

ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = 'Control')) + 
  geom_rug(aes(x = fl, y = 0), position = position_jitter(height = 0)) + 
  geom_histogram(aes(fill = 'Black'), binwidth = 0.002)

ggplot(dfs[[1]], aes(fl)) + 
  geom_histogram(aes(fill = 'Control'), binwidth = 0.002)

ggplot(dfs[[1]], aes(fl)) + 
  geom_freqpoly(aes(fill = 'Control'), binwidth = 0.005)

ggplot(dfs[[13]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = set), alpha = alp,  adjust = bw)# +
  facet_grid(~ set)
  
  ggplot(dfs[[10]], aes(fl)) + 
    geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
    geom_density(data = dfs[[11]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
    geom_density(data = dfs[[12]], aes(fl, fill = set), alpha = alp,  adjust = bw)# +
  facet_grid(~ set)

ggplot(r_dfs[[13]], aes(fl)) + 
    geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
    geom_density(data = r_dfs[[14]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
    geom_density(data = r_dfs[[15]], aes(fl, fill = set), alpha = alp,  adjust = bw)# +
  facet_grid(~ set)
  
                       