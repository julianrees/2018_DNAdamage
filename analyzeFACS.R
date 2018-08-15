#---- Header content ----
library(ggplot2)
library(BiocInstaller)
library(ggcyto)
library(flowCore)
library(outliers)
library(fitdistrplus)
library(clusterSim)
library(dplyr)
library(plyr)
library(multcomp)

# source("https://bioconductor.org/biocLite.R")
# biocLite()

# custom function definitions
abcellplot <- function(Vab, Vcell) {
  ggplot(tlog_df[which(tlog_df$antibody == Vab & tlog_df$cellline == Vcell), ], aes(fl, fill = replicate, by = replicate)) +
    geom_density(alpha = alp,  adjust = bw) +
    facet_grid(timepoint ~ dose) +
    geom_vline(aes(xintercept = mean), log_ctrl_means[which(log_ctrl_means$cellline == Vcell &
                                                              log_ctrl_means$antibody == Vab),]) +
    geom_vline(aes(xintercept = mean+sd), log_ctrl_means[which(log_ctrl_means$cellline == Vcell &
                                                                 log_ctrl_means$antibody == Vab),],
               linetype = 2) +
    geom_vline(aes(xintercept = mean-sd), log_ctrl_means[which(log_ctrl_means$cellline == Vcell &
                                                                 log_ctrl_means$antibody == Vab),],
               linetype = 2)
}











#for dir in ./*; do (cd "$dir" && bulk_rename _S _B csv); done

setwd('./')
#---- Data import from the listed folders ----
folders <- dir("Data/merged_expts//")
for (i in seq(length(folders))){
  folders[i] <- paste("Data/merged_expts//", folders[i], '/', sep = "")
}

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
# take the logs of the data and the medians and means of the log data
logdata <- list()
logmedians <- list()
logmeans <- list()
for (i in seq(length(rdata))){
  logdata[[i]] <- signif(log(rdata[[i]]), digits = 6)
  #logmedians[[i]] <- median(logdata[[i]])
  logmeans[[i]] <- mean(logdata[[i]])
}


# make the raw and log dataframes
r_dfs <- list()
log_dfs <- list()

# make lists of individual data frams for the raw and log data
for (i in seq(length(logdata))){
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
}



# generate "total" dataframes for the raw and log data
tr_df <- r_dfs[[1]]
tlog_df <- log_dfs[[1]]

for (i in seq(2,length(log_dfs))){
  #tr_df <- rbind(tr_df, r_dfs[[i]])
  tlog_df <- rbind(tlog_df, log_dfs[[i]])

}


# Make a new factor column for the different experiments

expts <- array(dim = nrow(tlog_df))

sel = which(tlog_df$replicate == 'S1'
            | tlog_df$replicate == 'S2'
            | tlog_df$replicate == 'S3')
expts[sel] <- 'S'

sel = which(tlog_df$replicate == 'T1'
            | tlog_df$replicate == 'T2'
            | tlog_df$replicate == 'T3'
            | tlog_df$replicate == 'T4')
expts[sel] <- 'T'

sel = which(tlog_df$replicate == 'U1'
            | tlog_df$replicate == 'U2'
            | tlog_df$replicate == 'U3'
            | tlog_df$replicate == 'U4')
expts[sel] <- 'U'

tlog_df <- cbind(tlog_df, experiment = as.factor(expts))
rm(expts)

# make a data frame of run details
runs <- ddply(tlog_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
              log.mean = round(mean(fl), 3),
              sd = round(sd(fl), 3))
runlength <- array()
for (i in seq(1,length(logdata))){
  runlength[i] <- length(logdata[[i]])
}

runscounts <- cbind(runs, length = runlength)

#---- Start looking at statistics and sets to toss ----

# calculate the mean and SD of all distributions - which are outside of normal shape?
mean(runs$sd)
sd(runs$sd)

ggplot(runs, aes(x = antibody, y = sd)) +
  geom_boxplot()

ggplot(runs, aes(sd)) +
  geom_histogram() +
  facet_grid(cellline~antibody)

ggplot(runs, aes(sd)) +
  geom_density() +
  facet_grid(cellline~antibody)


#---- SET REMOVAL ----

nrow(tlog_df[which(tlog_df$antibody == 'H2aX' &
                     tlog_df$cellline == 'HCC' &
                     tlog_df$timepoint == '24h 24h' &
                     tlog_df$dose == 'LD' &
                     tlog_df$replicate == 'S1'),])


#build a dataframe of sets to remove, either based on cell count or set SD
removals <- runs[1,1:5]
removals <- rbind(removals, c('H2aX','SKBR3','4h','U3','CTRL'))
removals <- rbind(removals, c('ATF2','HCC','24h','T3','HD'))
removals <- rbind(removals, c('ATF2','BT549','4h','S2','HD'))
removals <- rbind(removals, c('ATF2','BT549','4h','S2','LD'))
removals <- rbind(removals, c('ATF2','BT549','24h 24h','S3','CTRL'))
# less than 10,000 cells
removals <- rbind(removals, c('ATF2','HCC','24h','S1','HD'))
removals <- rbind(removals, c('ATF2','HCC','24h','S1','LD'))
removals <- rbind(removals, c('ATF2','HCC','24h','S2','LD'))
removals <- rbind(removals, c('H2aX','HCC','24h 24h','S1','CTRL'))
removals <- rbind(removals, c('H2aX','HCC','24h 24h','S2','CTRL'))

removals <- removals[-1,]


# remove sets from tlog_df, rebuild runs & log_ctrl_means

for (r in seq(nrow(removals))){
tlog_df <- tlog_df[-which(tlog_df$antibody == as.character(removals$antibody[r]) &
        tlog_df$cellline == as.character(removals$cellline[r]) &
        tlog_df$timepoint == as.character(removals$timepoint[r]) &
        tlog_df$replicate == as.character(removals$replicate[r]) &
        tlog_df$dose == as.character(removals$dose[r])),]
}

runs <- ddply(tlog_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
              log.mean = round(mean(fl), 3),
              sd = round(sd(fl), 3))

log_ctrl_means <- ddply(runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        sd = round(sd(log.mean), 3))

# ask for a summary of the means and standard deviations within each condition
# for each experiment
log_fl_means <- ddply(tlog_df, .(antibody, cellline, timepoint, dose, experiment), summarize,
                      mean = round(mean(fl), 3),
                      sd = round(sd(fl), 3))

# pull the means and SDs by cellline and antibody of the controls for each timepoint
log_mean_mean <- ddply(runs, .(antibody, cellline, timepoint, dose, experiment), summarize,
                       mean = round(mean(log.mean), 3),
                       sd = round(sd(log.mean), 3))

# pull the means and SDs by cellline and antibody of the controls across all timepoints
log_ctrl_means <- ddply(runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        sd = round(sd(log.mean), 3))

# runs_norms <- merge(runs, log_ctrl_means, by = c('timepoint',
#                                                  'antibody',
#                                                  'cellline'))
# colnames(runs_norms)[8:9] <- c('ctrl_mean','ctrl_sd')


runs[which(runs$antibody == pAb & runs$cellline == pcell &
             runs$dose == 'CTRL' &
           runs$timepoint == '24h'),]


which(log_fl_means$mean[which(log_fl_means$antibody == 'ATF2')] >
        atf2_mean + atf2_sd)

which(log_fl_means$mean[which(log_fl_means$antibody == 'H2aX')] >
        h2ax_mean + h2ax_sd)

# inspect the data for systematic shifts
# plotting variables
pAb = 'H2aX'
pcell = 'HCC'


abcellplot(pAb, pcell)

# SKBR3 + ATF2 + 4h + Experiment T has a systematic shift. Correcting it, keeping variance within set.
# BT549 + ATF2 + 4h + Experiment S has a systematic shift. Will correct it too.
# HCC + ATF2 + 24h 24h / 4h / 4h 24h + Experiment S has a systematic shift.
fixAb = 'ATF2'
fixCell = 'HCC'
fixTime = '4h 24h'
fixExp = 'S'

goodmean <- mean(runs$log.mean[which(runs$cellline == fixCell &
                                 runs$antibody == fixAb &
                                 runs$timepoint == fixTime &
                                 runs$experiment != fixExp)])
setmean <- mean(runs$log.mean[which(runs$cellline == fixCell &
                                      runs$antibody == fixAb &
                                      runs$timepoint == fixTime &
                                      runs$experiment == fixExp)])

# now need to fix tlog_df, runs, log_mean_means, log_ctrl_means,
fixFactor = goodmean - setmean

tlog_df$fl[which(tlog_df$antibody == fixAb &
                   tlog_df$cellline == fixCell &
                   tlog_df$timepoint == fixTime &
                   tlog_df$experiment == fixExp)] <-
  tlog_df$fl[which(tlog_df$antibody == fixAb &
                     tlog_df$cellline == fixCell &
                     tlog_df$timepoint == fixTime &
                     tlog_df$experiment == fixExp)] + fixFactor

runs <- ddply(tlog_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
      log.mean = round(mean(fl), 3),
      sd = round(sd(fl), 3))

log_ctrl_means <- ddply(runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        sd = round(sd(log.mean), 3))

# pull the means and SDs by cellline and antibody of the controls for each timepoint
log_mean_mean <- ddply(runs, .(timepoint, antibody, cellline, dose, experiment), summarize,
                       mean = round(mean(log.mean), 3),
                       sd = round(sd(log.mean), 3))

# check to make sure it worked by plotting

abcellplot(fixAb, fixCell)






mean(tlog_df$fl[which(tlog_df$antibody == 'ATF2' & tlog_df$dose == 'CTRL' &
                tlog_df$experiment == 'U')])








runs_norms <- merge(runs, log_mean_mean[which(log_mean_mean$dose == "CTRL"),], by = c('timepoint',
                                                                                      'antibody'))
colnames(runs_norms)[6] <- 'dose'
colnames(runs_norms)[9:10] <- c('normfactor','normfacSD')


normdata <- list()
for (i in seq(length(logdata))){
  normdata[[i]] <- 1 + logdata[[i]] - runs_norms$normfactor[[i]]
}


# make lists of individual data frams for the normalized data
for (i in seq(length(logdata))){
  dfs[[i]] <- data.frame(fl = normdata[[i]], set = dataset_name[[i]],
                         antibody = antibody[i],
                         cellline = cellline[i],
                         timepoint = timepoint[i],
                         dose = dose[i],
                         replicate = replicate[i],
                         experiment = expts2[i])
}







tdf <- dfs[[1]]
for (i in seq(2,length(dfs))){
  tdf <- rbind(tdf, dfs[[i]])
  #tr_df <- rbind(tr_df, r_dfs[[i]])
  #tlog_df <- rbind(tlog_df, log_dfs[[i]])
}





# normalize to 1 (or standardize) the data by taking 1 + (x - median(control set for x))/mad(control set for x)

# normfactors <- array()
#
# for (i in seq(length(logdata))){
#   logmads[[i]] <- mad(logdata[[i]])
#   normdata[[i]] <- 1 + (logdata[[i]] - logmedians[[control[[i]]]])/logmads[[control[[i]]]]
#   normfactors[i] <- median(normdata[[control[[i]]]])
# }







# need to remove sets which don't make sense here, before making means of them
# calculate the mean and standard deviation of the means of alllllll data sets
# for a given antibody
h2ax_mean <- mean(runs$log.mean[which(runs$antibody == 'H2aX')])
atf2_mean <- mean(runs$log.mean[which(runs$antibody == 'ATF2')])
# h2ax_median <- median(runs$log.mean[which(runs$antibody == 'H2aX')])
# atf2_median <- median(runs$log.mean[which(runs$antibody == 'ATF2')])
h2ax_sd <- sd(runs$log.mean[which(runs$antibody == 'H2aX')])
atf2_sd <- sd(runs$log.mean[which(runs$antibody == 'ATF2')])

# repeat above, but only for control groups
h2ax_ctrl_mean <- mean(runs$log.mean[which(runs$antibody == 'H2aX' & runs$dose == 'CTRL')])
atf2_ctrl_mean <- mean(runs$log.mean[which(runs$antibody == 'ATF2' & runs$dose == 'CTRL')])
# h2ax_ctrl_median <- median(runs$log.mean[which(runs$antibody == 'H2aX' & runs$dose == 'CTRL')])
# atf2_ctrl_median <- median(runs$log.mean[which(runs$antibody == 'ATF2' & runs$dose == 'CTRL')])
h2ax_ctrl_sd <- sd(runs$log.mean[which(runs$antibody == 'H2aX' & runs$dose == 'CTRL')])
atf2_ctrl_sd <- sd(runs$log.mean[which(runs$antibody == 'ATF2' & runs$dose == 'CTRL')])

dixon.test(runs$log.mean[which(runs$antibody == 'H2aX' & runs$dose == 'CTRL')])
dixon.test(runs$log.mean[which(runs$antibody == 'ATF2' & runs$dose == 'CTRL')])






#---- Statistical description and analysis ----

# define the sample size n
n = 5000

Dose = "CTRL"
Cell = "SKBR3"
Time = "24h 24h"
Ab = "H2aX"


indexer <- which(tlog_df$dose == Dose &
                   tlog_df$cellline == Cell &
                   tlog_df$timepoint == Time &
                   tlog_df$antibody == Ab &
                   tlog_df$replicate == 'S1')

ref_indexer <- which(tlog_df$dose == Dose &
                   tlog_df$cellline == Cell &
                   tlog_df$timepoint == Time &
                   tlog_df$antibody == Ab
 #                  & tlog_df$replicate == 'S1'
                   )

fitdata <- sample_n(tlog_df[indexer,], n)
refdata <- sample_n(tlog_df[ref_indexer,], n)

disdata <- fitdist(fitdata$fl, 'norm')
refdisdata <- fitdist(refdata$fl, 'norm')
plot(disdata)
plot(refdisdata)


# Kolmogorov-Smirnov test for normality (D closer to 0 is more normal / the same distribution)
ks.test(unique(disdata$data), 'pnorm')
ks.test(unique(disdata$data), unique(refdisdata$data))



#densityplot(fitdata$fl)

#scores(disdata$data)



#plot(ecdf(disdata$data))






kruskal.test(fl ~ replicate, data = fitdata)



#fit <- aov(fl ~ timepoint, data = tdf[which(tdf$dose == 'CTRL' & tdf$cellline == "BT549"),])

fit <- aov(fl ~ replicate, data = fitdata)

summary(fit)
summary(glht(fit, linfct=mcp(replicate = "Tukey")))



summary(tr_df$replicate[indexer])

subdata <- tdf$fl[indexer]
subdata2 <- tdf$fl[indexer2]

# distribution analysis
qqnorm(subdata); qqline(subdata)
#qqplot(subdata, 'pnorm')
plotdist(subdata, histo = TRUE, demp = TRUE)
descdist(ja3$data, boot = 100)

ja <- fitdist(subdata, "norm")
ja2 <- fitdist(subdata2, "norm")
ja3 <- fitdist(tr_df$fl[indexer], 'lnorm')
denscomp(ja3)
cdfcomp(ja3)

plot(ja3)

gofstat(ja3, fitnames = "norm")


# Kolmogorov-Smirnov test for normality (D closer to 0 is more normal)
plot(ecdf(disdata$data))
ks.test(unique(disdata$data), 'pnorm')
ks.test(unique(ja3$data), unique(subdata2))


# Wilcoxon Rank Sum Test
wilcox.test(subdata, subdata2)
kruskal.test(fl ~ dose, data = tdf[indexer,])


str(subdata)

# T-test

t.test(x = tdf$fl[which(tdf$dose == "HD" &
                    tdf$cellline == "HCC" &
                    tdf$timepoint == "24h" &
                    tdf$antibody == 'H2aX')],
       y = tdf$fl[which(tdf$dose == "LD" &
                          tdf$cellline == "HCC" &
                          tdf$timepoint == "24h" &
                          tdf$antibody == 'H2aX')],
       conf.level = 0.95)

# ANOVA



# ==== PLOT THE DATA ====
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

alp = 0.5 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

pfill = 'replicate' # select the fill
pdf = dfs[[1]]

figprefix <- 'Figures/merged_expts/'

# ---- Plot the raw data ----
ggplot(tr_df[which(tr_df$dose == 'CTRL' & tr_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggsave(filename = paste(figprefix, 'raw_controls_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$dose == 'CTRL' & tr_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  xlim(4e4,3e5) +
  ggtitle("raw controls H2aX") +
  ggsave(filename = paste(figprefix, 'raw_controls_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'raw_BT549_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'raw_BT549_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'raw_SKBR3_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'raw_SKBR3_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'raw_HCC_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  xlim(4e4,3e5)
  ggsave(filename = paste(figprefix, 'raw_HCC_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")



ggplot(tr_df[which(tr_df$dose == 'LD' & tr_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline)

ggplot(tr_df[which(tr_df$dose == 'LD' & tr_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  xlim(4e4,1.7e5)



ggplot(tr_df[which(tr_df$dose == 'HD' & tr_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline)

ggplot(tr_df[which(tr_df$dose == 'HD' & tr_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  xlim(4e4,1.7e5)


# ---- Plot the log data ----

ggplot(tlog_df[which(tlog_df$dose == 'CTRL' & tlog_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggsave(filename = paste(figprefix, 'log_controls_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")



ggplot(tlog_df[which(tlog_df$dose == 'CTRL' & tlog_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggtitle('log congrols H2aX') +
  ggsave(filename = paste(figprefix, 'log_controls_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'log_BT549_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'log_BT549_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'log_SKBR3_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'log_SKBR3_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'log_HCC_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
ggsave(filename = paste(figprefix, 'log_HCC_H2aX.pdf', sep = ""),
       width = 8.5, height = 5.5, units = "in")







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


# ---- Plot the normalized data ----

ggplot(tdf[which(tdf$dose == 'CTRL'), ], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ antibody) +
  geom_vline(xintercept = 1) +
  ggsave(filename = paste(figprefix, 'controls_normalized.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$dose == 'CTRL'), ], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ antibody) +
  geom_vline(xintercept = 1) +
  ggsave(filename = paste(figprefix,'controls_normalized_grouped.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "H2aX"),], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  geom_vline(xintercept = 1) +
  #xlim(-8,10) +
  ggsave(filename = paste(figprefix,'H2aX_dose_by_timepoint.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "ATF2"),], aes(fl, fill = cellline, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  geom_vline(xintercept = 1) +
  ggsave(filename = paste(figprefix,'ATF2_dose_by_timepoint.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "H2aX"),], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  geom_vline(xintercept = 1) +
  ggsave(filename = paste(figprefix,'H2aX_dose_by_timepoint_averaged.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf[which(tdf$antibody == "ATF2"),], aes(fl, fill = cellline)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  geom_vline(xintercept = 1) +
  ggsave(filename = paste(figprefix,'ATF2_dose_by_timepoint_averaged.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tdf, aes(x = dose, y = fl, fill = cellline, by = replicate, color = replicate)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody~timepoint) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = paste(figprefix, 'boxplot_doses_timepoint_by_antibody.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

# ggplot(tlog_df[which(tlog_df$dose == "CTRL"),], aes(x = dose, y = fl, fill = cellline, by = replicate, color = replicate)) +
#   geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
#   facet_grid(antibody~timepoint) +
#   ggsave(filename = paste(figprefix, 'boxplot_doses_timepoint_by_antibody_unnormalized_controls.pdf', sep = ""),
#          width = 8.5, height = 5.5, units = "in")

ggplot(tdf, aes(x = dose, y = fl, fill = cellline)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody~timepoint) +
  geom_hline(yintercept = 1) +
  scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = paste(figprefix, 'boxplot_doses_timepoint_by_antibody_averaged.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tdf, aes(x = cellline, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody~timepoint) +
  geom_hline(yintercept = 1) +
  # scale_x_discrete(limits = levels(tdf$dose)[c(1,3,2)]) +
  ggsave(filename = paste(figprefix, 'boxplot_celllines_timepoint_by_antibody_averaged.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tdf, aes(x = cellline, y = fl, fill = dose)) +
  geom_violin(adjust = bw) +
  facet_grid(antibody~timepoint) +
  geom_hline(yintercept = 1) +
  #scale_x_discrete(limits = rev(levels(tdf$dose))[1:2]) +
  ggsave(filename = paste(figprefix, 'violin_celllines_timepoint_by_antibody_averaged.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")
