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
folders <- c("Data/pATF2/4h_BT549_ATF2/",
             "Data/pATF2/4h_24h_BT549_ATF2/",
             "Data/pATF2/4h_SKBR3_ATF2/",
             "Data/pATF2/4h_24h_SKBR3_ATF2/")
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
dataset_name

controls <- c(1,1,1,4,4,4,7,7,7,10,10,10,13,13,13,16,16,16,19,19,19,
              22,22,22,25,25,25,28,28,28,31,31,31,34,34,34)
# timepoint 1 is 4 hrs, 2 is 24 hrs, 3 is 4_24, and 4 is 24_24
timepoints <- c('4 hrs','4 hrs + 24 hrs', '24 hrs','24 hrs + 24 hrs')
timepoint <- c(1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,
               1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2)
# antibody 1 is gH2aX, 2 is ATF2
antibodies <- c('gH2aX','ATF2')
antibody <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,
              2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
# cellline 1 is BT549, 2 is SKBR3
celllines <- c('BT549','SKBR3')
cellline <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
              2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)
# dose 1 is control, 2 is high dose, 3 is low dose
doses <- c('Control','High','Low')
dose <- c(1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,1,2,3,
          1,2,3,1,2,3,1,2,3,1,2,3,1,2,3)
set <- c(1,1,1,2,2,2,3,3,3,1,1,1,2,2,2,3,3,3,1,1,1,
         2,2,2,3,3,3,1,1,1,2,2,2,3,3,3)

# select which controls to include in MAD via logical bit
includecontrols <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,
                     1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)

normdata <- list()
logmads <- list()
maddata <- list()
ctrl_abs_mean <- matrix(nrow = max(timepoint), ncol = max(cellline))
normfactors <- array()

# adjusts the datasets to the control absolute means for each time point. may need to do this differently
for (j in seq(max(timepoint))){
  for (k in seq(max(cellline))){
    ctrl_abs_mean[j,k] <- mean(as.numeric(logmedians[as.numeric(unique(controls[which(
      timepoint == j & includecontrols == 1 & antibody == 2 & cellline == k)]))]))
  }
  for (i in seq(length(logdata))){
    logmads[[i]] <- mad(logdata[[i]], constant = ctrl_abs_mean[j,k])
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
  r_dfs[[i]] <- data.frame(fl = rdata[[i]], set = dataset_name[[i]], 
                           antibody = antibodies[antibody[i]], 
                           cell = celllines[cellline[i]],
                           timepoint = timepoints[timepoint[i]],
                           dose = doses[dose[i]],
                           sets = set[i])
  log_dfs[[i]] <- data.frame(fl = logdata[[i]], set = dataset_name[[i]], 
                             antibody = antibodies[antibody[i]], 
                             cell = celllines[cellline[i]],
                             timepoint = timepoints[timepoint[i]],
                             dose = doses[dose[i]],
                             sets = set[i])
  mad_dfs[[i]] <- data.frame(fl = maddata[[i]], set = dataset_name[[i]], 
                             antibody = antibodies[antibody[i]], 
                             cell = celllines[cellline[i]],
                             timepoint = timepoints[timepoint[i]],
                             dose = doses[dose[i]],
                             sets = set[i])
  dfs[[i]] <- data.frame(fl = normdata[[i]], set = dataset_name[[i]], 
                         antibody = antibodies[antibody[i]], 
                         cell = celllines[cellline[i]],
                         timepoint = timepoints[timepoint[i]],
                         dose = doses[dose[i]],
                         sets = as.factor(set[i]))
}


#---- generate geometric means of data from corresponding sets ----



#---- Plotting the data ----
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())
alp = 0.4 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

# compare the control groups in log space and not normalized
ggplot(log_dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = log_dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[16]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[19]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[22]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[25]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[28]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = log_dfs[[31]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = log_dfs[[34]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = as.vector(ctrl_abs_mean), aes(color)) + 
  facet_grid(cell ~ timepoint) + 
  ggsave(filename = 'Figures/controls_logspace.pdf',
                 width = 8.5, height = 5.5, units = "in")




ggplot(mad_dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = mad_dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[16]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[19]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[22]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[25]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[28]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = mad_dfs[[31]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = mad_dfs[[34]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = unique(normfactors), aes(color)) + 
  facet_grid(cell ~ timepoint) + 
  ggsave(filename = 'Figures/maddif_logspace.pdf',
         width = 8.5, height = 5.5, units = "in")


ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = set), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = set), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = set), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(cell ~ timepoint) + 
  ggsave(filename = 'Figures/controls_normalized.pdf',
         width = 8.5, height = 5.5, units = "in")

alp = 0.4
ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = sets), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[5]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[6]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[8]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[9]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[11]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[12]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[17]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[18]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[20]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[21]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[23]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[24]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[26]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[27]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[29]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[30]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[32]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[33]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[35]], aes(fl, fill = sets), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[36]], aes(fl, fill = sets), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(dose~cell) +
  ggsave(filename = 'Figures/ATF2_cell_by_dose2.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = dose), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[5]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[6]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[8]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[9]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[11]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[12]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[17]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[18]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[20]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[21]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[23]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[24]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[26]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[27]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[29]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[30]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[32]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[33]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[35]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[36]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(sets~cell) + 
  ggsave(filename = 'Figures/ATF2_cell_by_set.pdf',
         width = 8.5, height = 5.5, units = "in")


ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = dose), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[5]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[6]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[8]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[9]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[11]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[12]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[17]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[18]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[20]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[21]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[23]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[24]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[26]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[27]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[29]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[30]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[32]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[33]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[35]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[36]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(sets~cell) + 
  ggsave(filename = 'Figures/ATF2_cell_by_set.pdf',
         width = 8.5, height = 5.5, units = "in")


ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = dose), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[5]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[6]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[8]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[9]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[11]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[12]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[17]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[18]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[20]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[21]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[23]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[24]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[26]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[27]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[29]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[30]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[32]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[33]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[35]], aes(fl, fill = dose), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[36]], aes(fl, fill = dose), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(timepoint~cell) + 
  ggsave(filename = 'Figures/ATF2_cell_by_timepoint.pdf',
         width = 8.5, height = 5.5, units = "in")

alp = 0.6
ggplot(dfs[[1]], aes(fl)) + 
  geom_density(aes(fill = cell), alpha = alp, adjust = bw) + 
  geom_density(data = dfs[[2]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[3]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[4]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[5]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[6]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[7]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[8]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[9]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[10]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[11]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[12]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[13]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[14]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[15]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[16]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[17]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[18]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[19]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[20]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[21]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[22]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[23]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[24]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[25]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[26]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[27]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[28]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[29]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[30]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[31]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[32]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[33]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[34]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_density(data = dfs[[35]], aes(fl, fill = cell), alpha = alp,  adjust = bw) + 
  geom_density(data = dfs[[36]], aes(fl, fill = cell), alpha = alp,  adjust = bw) +
  geom_vline(xintercept = 1) + 
  facet_grid(timepoint~dose) + 
  ggsave(filename = 'Figures/ATF2_dose_by_timepoint.pdf',
         width = 8.5, height = 5.5, units = "in")
