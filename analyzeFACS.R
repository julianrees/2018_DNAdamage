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
             "Data/pATF2/4h_24h_SKBR3_ATF2/",
             "Data/pATF2/24h_BT549_ATF2/",
             "Data/pATF2/24h_24h_BT549_ATF2/",
             "Data/pATF2/24h_SKBR3_ATF2/",
             "Data/pATF2/24h_24h_SKBR3_ATF2/")

# read the directory contents in, populate raw data in the rdata, generate names in dataset_names
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
}

# split the filenames at the underscores, read the categorical info from the text strings into the following arrays
timepoint <- array(dim = length(rdata))
cellline <- array(dim = length(rdata))
antibody <- array(dim = length(rdata))
replicate <- array(dim = length(rdata))
dose <- array(dim = length(rdata))
control <- array(dim = length(rdata))

for (i in seq(length(rdata))){
  if (length(unlist(strsplit(dataset_name[[i]], split = "_"))) == 5){
    timepoint[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[1]
    cellline[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[2]
    antibody[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[3]
    replicate[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[4]
    dose[i] <- unlist(strsplit(unlist(strsplit(dataset_name[[i]], split = c("_")))[5], "[.]"))
  }
  else {
    timepoint[i] <- paste(unlist(strsplit(dataset_name[[i]], split = "_"))[1], 
                      unlist(strsplit(dataset_name[[i]], split = "_"))[2])
    cellline[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[3]
    antibody[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[4]
    replicate[i] <- unlist(strsplit(dataset_name[[i]], split = "_"))[5]
    dose[i] <- unlist(strsplit(unlist(strsplit(dataset_name[[i]], split = c("_")))[6], "[.]"))
  }
}

# create array for the control ids, assuming control is the first one in each replicate
for (i in seq(length(rdata))){
  if (dose[i] == "CTRL"){
    j = i
  }
  control[i] <- j
}

# generate the array of "include" bits for the control groups
includecontrol <- array(dim = length(rdata), 1)


#---- Data processing ----
logdata <- list()
logmedians <- list()
for (i in seq(length(rdata))){
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 6)
  logmedians[[i]] <- median(logdata[[i]])
}

normdata <- list()
logmads <- list()
maddata <- list()
ctrl_abs_mean <- matrix(nrow = length(unique(timepoint)), ncol = length(unique(cellline)))
normfactors <- array()

# adjusts the datasets to the control absolute means for each time point. may need to do this differently
for (j in seq(NROW(ctrl_abs_mean))){
  for (k in seq(NCOL(ctrl_abs_mean))){
    ctrl_abs_mean[j,k] <- mean(as.numeric(logmedians[as.numeric(unique(control[which(
        timepoint == unique(timepoint)[j] & 
        includecontrol == 1 & 
        antibody == unique(antibody[2]) & 
        cellline == unique(cellline)[k])]))]))
    print(j)
  }
  for (i in seq(length(logdata))){
    logmads[[i]] <- mad(logdata[[i]], constant = ctrl_abs_mean[j,k])
    maddata[[i]] <- logdata[[i]] - logmads[[control[[i]]]]
    normfactors[i] <- median(maddata[[control[[i]]]])
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
                           antibody = antibody[i], 
                           cellline = cellline[i],
                           timepoint = timepoint[i],
                           dose = dose[i],
                           replicate = replicate[i])
  log_dfs[[i]] <- data.frame(fl = logdata[[i]], set = dataset_name[[i]], 
                             antibody = antibody[i], 
                             cellline = cellline[i],
                             timepoint = timepoint[i],
                             dose = dose[i],
                             replicate = replicate[i])
  mad_dfs[[i]] <- data.frame(fl = maddata[[i]], set = dataset_name[[i]], 
                             antibody = antibody[i], 
                             cellline = cellline[i],
                             timepoint = timepoint[i],
                             dose = dose[i],
                             replicate = replicate[i])
  dfs[[i]] <- data.frame(fl = normdata[[i]], set = dataset_name[[i]], 
                         antibody = antibody[i], 
                         cellline = cellline[i],
                         timepoint = timepoint[i],
                         dose = dose[i],
                         replicate = replicate[i])
}

# generate "total" dataframes for each df
tdf <- dfs[[1]]
tr_df <- r_dfs[[1]]
tlog_df <- log_dfs[[1]]
tmad_df <- mad_dfs[[1]]

for (i in seq(2,length(dfs))){
  tdf <- rbind(tdf, dfs[[i]])
  # tr_df <- rbind(tr_df, r_dfs[[i]])
  # tlog_df <- rbind(tlog_df, log_dfs[[i]])
  # tmad_df <- rbind(tmad_df, mad_dfs[[i]])
}




#---- Plotting the data ----
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

alp = 0.5 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

pfill = 'cellline' # select the fill 
pdf = dfs[[1]]
  
ggplot(data = pdf, aes_string("fl", fill = pfill), alpha = alp, adjust = bw) + 
  geom_density()


tdf <- dfs[[1]]
for (i in seq(2,length(dfs))){
  tdf <- rbind(tdf, dfs[[i]])
}

ggplot(tdf, aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) + 
  geom_vline(xintercept = 1)

alp = 0.8
ggplot(tdf, aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) + 
  geom_vline(xintercept = 1)







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

pdf = dfs
pfill = 'cell'
palpha = 0.5


ggplot(data = dfs[[1]], aes_string("fl", fill = pfill), alpha = alp, adjust = bw) + 
  geom_density()
