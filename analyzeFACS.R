#---- Header content ----
library(ggplot2)
library(BiocInstaller)
library(ggcyto)
library(flowCore)
library(outliers)
library(fitdistrplus)
library(clusterSim)
library(dplyr)
library(multcomp)

# source("https://bioconductor.org/biocLite.R")
# biocLite()

setwd('./')
#---- Data import from the listed folders ----
folders <- dir("Data/041218_case2/")
for (i in seq(length(folders))){
  folders[i] <- paste("Data/041218_case2/", folders[i], '/', sep = "")
}
# folders <- c("Data/pATF2/4h_BT549_ATF2/",
#              "Data/pATF2/4h_24h_BT549_ATF2/",
#              "Data/pATF2/4h_SKBR3_ATF2/",
#              "Data/pATF2/4h_24h_SKBR3_ATF2/",
#              "Data/pATF2/24h_BT549_ATF2/",
#              "Data/pATF2/24h_24h_BT549_ATF2/",
#              "Data/pATF2/24h_SKBR3_ATF2/",
#              "Data/pATF2/24h_24h_SKBR3_ATF2/")

# read the directory contents in, populate raw data into rdata, generate names in dataset_names
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

#---- Data processing ----
# take the logs of the data and the medians of the log data
logdata <- list()
logmedians <- list()
for (i in seq(length(rdata))){
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 6)
  logmedians[[i]] <- median(logdata[[i]])
}

# normalize to 1 (or standardize) the data by taking 1 + (x - median(control set for x))/mad(control set for x)
normdata <- list()
logmads <- list()
maddata <- list()
normfactors <- array()

for (i in seq(length(logdata))){
  logmads[[i]] <- mad(logdata[[i]])
  normdata[[i]] <- 1 + (logdata[[i]] - logmedians[[control[[i]]]])/logmads[[control[[i]]]]
  normfactors[i] <- median(normdata[[control[[i]]]])
}

#---- Make the log and normalized dataframes ----
r_dfs <- list()
log_dfs <- list()
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

for (i in seq(2,length(dfs))){
  tdf <- rbind(tdf, dfs[[i]])
  #tr_df <- rbind(tr_df, r_dfs[[i]])
  #tlog_df <- rbind(tlog_df, log_dfs[[i]])
}

#---- Statistical description and analysis ----
indexer <- which(tdf$dose == "LD" &
                   tdf$cellline == "SKBR3" &
                   tdf$timepoint == "4h 24h")

subdata <- tdf$fl[indexer]

# distribution analysis
plotdist(subdata, histo = TRUE, demp = TRUE)
descdist(subdata, boot = 100)

# T-test
t.test(x = tdf$fl[which(tdf$cellline == "BT549" &
                          tdf$timepoint == "4h" &
                          tdf$dose == "CTRL")], 
       y = tdf$fl[which(tdf$cellline == "BT549" &
                          tdf$timepoint == "4h" &
                          tdf$dose == "LD")], 
       conf.level = 0.95)


# ANOVA
fit <- aov(fl ~ timepoint, data = tdf[which(tdf$dose == 'CTRL' & tdf$cellline == "BT549"),])
summary(fit)
summary(glht(fit, linfct=mcp(timepoint="Tukey")))


#---- Plotting the data ----
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

alp = 0.5 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

pfill = 'cellline' # select the fill 
pdf = dfs[[1]]
  
# ggplot(data = pdf, aes_string("fl", fill = pfill), alpha = alp, adjust = bw) + 
#   geom_density()
# 
# 
# ggplot(tlog_df[which(tlog_df$dose == "CTRL"),], aes(fl, fill = cellline, by = replicate)) +
#   geom_density(alpha = alp,  adjust = bw) + 
#   facet_grid(timepoint ~ dose)



# plot the grouped and individual control groups
# ggplot(tlog_df[which(tlog_df$dose == 'CTRL'), ], aes(fl, fill = cellline)) +
#   geom_density(alpha = alp,  adjust = bw) + 
#   facet_grid(timepoint ~ antibody) + 
#   geom_vline(xintercept = as.vector(ctrl_abs_mean))
# 
# ggplot(tlog_df[which(tlog_df$dose == 'CTRL'), ], aes(fl, fill = cellline, by = replicate)) +
#   geom_density(alpha = alp,  adjust = bw) + 
#   facet_grid(timepoint ~ dose) + 
#   geom_vline(xintercept = as.vector(ctrl_abs_mean)) #+ 
#   # ggsave(filename = 'Figures/controls_logspace.pdf',
#   #        width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$dose == 'CTRL'), ], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ antibody) + 
  geom_vline(xintercept = 1) + 
  ggsave(filename = 'Figures/controls_normalized.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$dose == 'CTRL'), ], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ antibody) + 
  geom_vline(xintercept = 1) + 
  ggsave(filename = 'Figures/controls_normalized_grouped.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "H2aX"),], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  geom_vline(xintercept = 1) +
  ggsave(filename = 'Figures/H2aX_dose_by_timepoint.pdf',
         width = 8.5, height = 5.5, units = "in")
  
ggplot(tdf[which(tdf$antibody == "ATF2"),], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  geom_vline(xintercept = 1) + 
  ggsave(filename = 'Figures/ATF2_dose_by_timepoint.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "H2aX"),], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  geom_vline(xintercept = 1) +
  ggsave(filename = 'Figures/H2aX_dose_by_timepoint_averaged.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "ATF2"),], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) + 
  facet_grid(timepoint ~ dose) +                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                     
  geom_vline(xintercept = 1) + 
  ggsave(filename = 'Figures/ATF2_dose_by_timepoint_averaged.pdf',
         width = 8.5, height = 5.5, units = "in")


ggplot(tdf, aes(x = dose, y = fl, fill = cellline, by = replicate)) + 
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = 1) +
  scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = 'Figures/boxplot_doses_timepoint_by_antibody.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf, aes(x = dose, y = fl, fill = cellline)) + 
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = 1) +
  scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = 'Figures/boxplot_doses_timepoint_by_antibody_averaged.pdf',
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf, aes(x = cellline, y = fl, fill = dose)) + 
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = 1) +
  # scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = 'Figures/boxplot_celllines_timepoint_by_antibody_averaged.pdf',
         width = 8.5, height = 5.5, units = "in")


ggplot(tdf, aes(x = cellline, y = fl, fill = dose)) + 
  geom_violin(adjust = bw) + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = 1) +
  #scale_x_discrete(limits = rev(levels(tdf$dose))[1:2]) + 
  ggsave(filename = 'Figures/violin_celllines_timepoint_by_antibody_averaged.pdf',
         width = 8.5, height = 5.5, units = "in")
