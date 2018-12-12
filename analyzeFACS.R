#---- Header content ----
## @knitr header
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
library(plotly)
library(lme)
library(nlme)
library(FSA)

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
               linetype = 2) + 
    ggtitle(paste(Vab, Vcell))
}

tabcellplot <- function(Vab, Vcell) {
  ggplot(tdf[which(tdf$antibody == Vab & tdf$cellline == Vcell), ], aes(fl, fill = replicate, by = replicate)) +
    geom_density(alpha = alp,  adjust = bw) +
    facet_grid(timepoint ~ dose) +
    geom_vline(aes(xintercept = mean), ctrl_means[which(ctrl_means$cellline == Vcell &
                                                              ctrl_means$antibody == Vab),]) +
    geom_vline(aes(xintercept = mean+sd), ctrl_means[which(ctrl_means$cellline == Vcell &
                                                                 ctrl_means$antibody == Vab),],
               linetype = 2) +
    geom_vline(aes(xintercept = mean-sd), ctrl_means[which(ctrl_means$cellline == Vcell &
                                                                 ctrl_means$antibody == Vab),],
               linetype = 2) + 
    ggtitle(paste(Vab, Vcell))
}

abcellbox <- function(Vab, Vcell) {
  ggplot(trunc_log_data[which(trunc_log_data$antibody == Vab & trunc_log_data$cellline == Vcell), ], aes(fl, fill = replicate, by = replicate)) +
    geom_density(alpha = alp,  adjust = bw) +
    facet_grid(timepoint ~ dose) +
    geom_vline(aes(xintercept = mean), trunc_log_ctrl_means[which(trunc_log_ctrl_means$cellline == Vcell &
                                                          trunc_log_ctrl_means$antibody == Vab),]) +
    geom_vline(aes(xintercept = mean+sd), trunc_log_ctrl_means[which(trunc_log_ctrl_means$cellline == Vcell &
                                                             trunc_log_ctrl_means$antibody == Vab),],
               linetype = 2) +
    geom_vline(aes(xintercept = mean-sd), trunc_log_ctrl_means[which(trunc_log_ctrl_means$cellline == Vcell &
                                                             trunc_log_ctrl_means$antibody == Vab),],
               linetype = 2) + 
    ggtitle(paste(Vab, Vcell))
}




# plotting preferences
theme_set(theme_bw())
theme_update(plot.title = element_text(hjust = 0.5),
             panel.grid.major = element_blank(),
             panel.grid.minor = element_blank())

alp = 0.5 # set transparency
bw = 0.5 # multiplier for bandwidth relative to defualt (SD)

pfill = 'replicate' # select the fill

# command for bulk renaming within a directory
#for dir in ./*; do (cd "$dir" && bulk_rename _S _B csv); done

setwd('./')

#---- Data import from the listed folders ----
folders <- dir("Data/merged_expts/")
for (i in seq(length(folders))){
  folders[i] <- paste("Data/merged_expts/", folders[i], '/', sep = "")
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

#---- Data cleanup ----
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

# make lists of individual data frames for the raw and log data
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

log_ctrl_means <- ddply(runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        sd = round(sd(log.mean), 3))

# make a static table of original runs, including cell counts
runlength <- array()
for (i in seq(1,length(logdata))){
  runlength[i] <- length(logdata[[i]])
}
runscounts <- cbind(runs[,-7:-8], cellcount = runlength)
runs <- merge(runs, runscounts, by = c('antibody',
                                       'cellline',
                                       'timepoint',
                                       'replicate',
                                       'dose',
                                       'experiment'))

#---- Data processing ----
# calculate the mean and SD of all distributions - which are outside of normal shape?
mean(runs$sd)
runs$sd
## @knitr load_preprocessed
load('Workspaces/180827_pre_processing.RData')

#---

ggplot(runs, aes(x = antibody, y = sd)) +
  geom_boxplot(aes(fill = cellline))

ggplot(runs, aes(x = antibody, y = sd)) +
  geom_jitter(aes(color = cellline)) +
  facet_grid(~cellline)


#---- SET REMOVAL ----

# bulk removal of sets, based on larger standard deviations among sets
setremovals <- matrix(nrow = 0, ncol = 4)
setremovals <- rbind(setremovals, c('H2aX','SKBR3','24h 24h','T'))
setremovals <- rbind(setremovals, c('H2aX','SKBR3','24h','T'))
setremovals <- rbind(setremovals, c('H2aX','SKBR3','4h 24h','T'))
setremovals <- rbind(setremovals, c('H2aX','SKBR3','4h','U'))

setremovals <- rbind(setremovals, c('ATF2','SKBR3','24h','T'))
setremovals <- rbind(setremovals, c('ATF2','SKBR3','4h 24h','T'))

setremovals <- rbind(setremovals, c('H2aX','BT549','24h 24h','T'))
setremovals <- rbind(setremovals, c('H2aX','BT549','4h','T')) # SD among control means much higher than T

for (r in seq(nrow(setremovals))){
  tlog_df <- tlog_df[-which(tlog_df$antibody == setremovals[r,1] &
                              tlog_df$cellline == setremovals[r,2] &
                              tlog_df$timepoint == setremovals[r,3] &
                              tlog_df$experiment == setremovals[r,4]),]
}

runs <- ddply(tlog_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
              log.mean = round(mean(fl), 3),
              sd = round(sd(fl), 3))
runs <- merge(runs, runscounts, by = c('antibody',
                                       'cellline',
                                       'timepoint',
                                       'replicate',
                                       'dose',
                                       'experiment'))

pcell = 'SKBR3'
pAb = 'ATF2'

ggplot(tlog_df[which(tlog_df$antibody == pAb),], aes(x = dose, y = fl, fill = dose, by = replicate)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(cellline ~ timepoint)









# build a table of runs to remove, based on e.g. cell count, SD
removals <- runs[1,1:5]
removals <- rbind(removals, c('ATF2','HCC','24h','T3','HD')) # SD over .5
removals <- rbind(removals, c('ATF2','BT549','4h','S2','HD')) # SD over .5
removals <- rbind(removals, c('ATF2','BT549','4h','S2','LD')) # SD over .5
#removals <- rbind(removals, c('H2aX','SKBR3','4h 24h','U3','HD')) # not a normal distrubtion
#removals <- rbind(removals, c('H2aX','BT549','4h','U4','CTRL')) # very low cell count
removals <- removals[-1,]

cellnumber = 2000
removals <- rbind(removals, runs[which(runs$cellcount < cellnumber),1:5])
dropped <- nrow(runs[which(runs$cellcount < cellnumber),1:5])

for (r in seq(nrow(removals))){
  if (removals$replicate[r] != 'T4' & removals$replicate[r] != 'U4' & removals$dose[r] == 'CTRL') {
    removals <- rbind(removals, cbind(removals[r,1:4], dose = levels(removals$dose)[2]))
    removals <- rbind(removals, cbind(removals[r,1:4], dose = levels(removals$dose)[3]))
  }
}

removals <- rbind(removals, runs[which(runs$replicate == 'U4' | runs$replicate == 'T4'), 1:5])
removals <- unique(removals)


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
              log.median = round(median(fl),3),
              sd = round(sd(fl), 3))

runs <- merge(runs, runscounts, by = c('antibody',
                                       'cellline',
                                       'timepoint',
                                       'replicate',
                                       'dose',
                                       'experiment'))

log_ctrl_means <- ddply(runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        median = round(mean(log.median),3),
                        sd = round(sd(log.mean), 3))

trunc_log_data <- ddply(tlog_df, .(antibody, cellline, timepoint, dose, replicate), 
                        function(x) head(x, cellnumber))

trunc_runs <- ddply(trunc_log_data, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
                    log.mean = round(mean(fl), 3),
                    log.median = round(median(fl),3),
                    sd = round(sd(fl), 3))
trunc_log_ctrl_means <- ddply(trunc_runs[which(runs$dose == 'CTRL'),],
                              .(antibody, cellline, timepoint), summarize,
                              mean = round(mean(log.mean), 3),
                              median = round(mean(log.median),3),
                              sd = round(sd(log.mean), 3))
trunc_log_means <- ddply(trunc_runs[which(trunc_runs$dose != 'CTRL'),], 
                         .(antibody, cellline), summarize, 
                         mean = mean(log.mean), 
                         sd = sd(log.mean))

plot(runs$log.mean[order(runs$log.mean)], trunc_runs$log.mean[order(trunc_runs$log.mean)])
plot(runs$log.mean[order(runs$log.mean)], trunc_runs$log.mean[order(trunc_runs$log.mean)] - 
       runs$log.mean[order(runs$log.mean)])

plot(runs$log.median[order(runs$log.median)], trunc_runs$log.median[order(trunc_runs$log.median)])
plot(runs$log.median[order(runs$log.median)], trunc_runs$log.median[order(trunc_runs$log.median)] - 
       runs$log.median[order(runs$log.median)])


plot(runs$sd[order(runs$sd)], trunc_runs$sd[order(trunc_runs$sd)])
plot(runs$sd[order(runs$sd)], trunc_runs$sd[order(trunc_runs$sd)] - runs$sd[order(runs$sd)])



sd(norm_runs$log.mean[which(norm_runs$antibody == 'H2aX' & norm_runs$cellline == '' & norm_runs$dose != 'CTRL')])

pAb = 'H2aX'
pcell = 'BT549'
conf = 1.5
# abcellbox(pAb, pcell)

ggplot(trunc_log_data[which(trunc_log_data$cellline == pcell),], aes(x = replicate, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) + 
  geom_hline(aes(yintercept = mean), trunc_log_means[which(trunc_log_means$cellline == pcell),]) + 
  geom_hline(aes(yintercept = mean+sd*conf), trunc_log_means[which(trunc_log_means$cellline == pcell),], 
             linetype = 2) + 
  geom_hline(aes(yintercept = mean-sd*conf), trunc_log_means[which(trunc_log_means$cellline == pcell),], 
             linetype = 2) 



trunc_norm_data <- merge(trunc_log_data, ddply(trunc_runs[which(trunc_runs$dose == 'CTRL'),],
                                           .(antibody, cellline, timepoint, replicate), summarize,
                                           mean = round(mean(log.mean), 3),
                                           median = round(mean(log.median),3),
                                           sd = round(sd(log.mean), 3)), by = c('antibody','cellline','timepoint','replicate'))
trunc_norm_data$fl <- (trunc_norm_data$fl-trunc_norm_data$median)+1

trunc_norm_runs <- ddply(trunc_norm_data, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
                    log.mean = round(mean(fl), 3),
                    log.median = round(median(fl),3),
                    sd = round(sd(fl), 3))
trunc_norm_ctrl_means <- ddply(trunc_norm_runs[which(runs$dose == 'CTRL'),],
                              .(antibody, cellline, timepoint), summarize,
                              mean = round(mean(log.mean), 3),
                              median = round(mean(log.median),3),
                              sd = round(sd(log.mean), 3))
trunc_norm_means <- ddply(trunc_norm_runs[which(trunc_norm_runs$dose != 'CTRL'),], 
                         .(antibody, cellline), summarize, 
                         mean = mean(log.mean), 
                         median = mean(log.median),
                         sd = sd(log.mean))

pAb = 'H2aX'
pcell = 'BT549'
conf = 1.5
# abcellbox(pAb, pcell)

ggplot(trunc_norm_data[which(trunc_norm_data$cellline == pcell),], aes(x = replicate, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) + 
  geom_hline(aes(yintercept = median), trunc_norm_means[which(trunc_norm_means$cellline == pcell),]) + 
  geom_hline(aes(yintercept = median+sd*conf), trunc_norm_means[which(trunc_norm_means$cellline == pcell),], 
             linetype = 2) + 
  geom_hline(aes(yintercept = median-sd*conf), trunc_norm_means[which(trunc_norm_means$cellline == pcell),], 
             linetype = 2)

trunc_norm_runs <- merge(trunc_norm_runs, trunc_norm_means, 
                         by = c('antibody','cellline'))
colnames(trunc_norm_runs) <- c(colnames(trunc_norm_runs)[1:6],'norm.mean','norm.median','norm.sd','global.mean','global.median','global.sd')

trunc_norm_runs[which((trunc_norm_runs$norm.mean < (trunc_norm_runs$global.median - trunc_norm_runs$global.sd*conf) |
  trunc_norm_runs$norm.median > (trunc_norm_runs$global.median + trunc_norm_runs$global.sd*conf))),1:5]

removals2 <- trunc_norm_runs[which((trunc_norm_runs$norm.median < (trunc_norm_runs$global.median - trunc_norm_runs$global.sd*conf) |
                                      trunc_norm_runs$norm.median > (trunc_norm_runs$global.median + trunc_norm_runs$global.sd*conf))),1:5]

for (r in seq(nrow(removals2))){
  trunc_norm_data <- trunc_norm_data[-which(trunc_norm_data$antibody == as.character(removals2$antibody[r]) &
                              trunc_norm_data$cellline == as.character(removals2$cellline[r]) &
                              trunc_norm_data$timepoint == as.character(removals2$timepoint[r]) &
                              trunc_norm_data$replicate == as.character(removals2$replicate[r]) &
                              trunc_norm_data$dose == as.character(removals2$dose[r])),]
}

trunc_norm_runs2 <- ddply(trunc_norm_data, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
                         log.mean = round(mean(fl), 3),
                         log.median = round(median(fl),3),
                         sd = round(sd(fl), 3))

trunc_norm_stats <- ddply(trunc_norm_runs2, .(antibody, cellline, timepoint, dose), summarize, 
                          stat.mean = mean(log.mean), 
                          stat.median = mean(log.median), 
                          stat.meansd = sd(log.mean), 
                          stat.mediansd = sd(log.median))


ggplot(trunc_norm_runs[which(trunc_norm_runs$dose != 'CTRL'),], aes(x = cellline, y = norm.mean)) + 
  geom_boxplot(aes(fill = antibody))






pcell = 'HCC'
pAb = 'ATF2'

ggplot(trunc_norm_data[which(trunc_norm_data$cellline == pcell),], aes(x = replicate, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) + 
  geom_hline(aes(yintercept = median), trunc_norm_means[which(trunc_norm_means$cellline == pcell),]) + 
  geom_hline(aes(yintercept = median+sd*conf), trunc_norm_means[which(trunc_norm_means$cellline == pcell),], 
             linetype = 2) + 
  geom_hline(aes(yintercept = median-sd*conf), trunc_norm_means[which(trunc_norm_means$cellline == pcell),], 
              linetype = 2) #+ 
  # ggsave(filename = paste(figprefix, 'normalized_boxplot_byset_HCC.pdf', sep = ""),
  #        width = 8.5, height = 5.5, units = "in")


ggplot(trunc_norm_data[which(trunc_norm_data$cellline == pcell),], aes(x = dose, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) + 
  geom_hline(yintercept = 1)


ggplot(trunc_norm_data[which(trunc_norm_data$antibody == pAb),], aes(x = dose, y = fl, fill = dose)) +
  geom_violin(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(cellline ~ timepoint) + 
  geom_hline(yintercept = 1) + 
  xlab('Dose') + 
  ylab('Normalized log(fluorescence)') + 
  ggtitle('ATF2 Dose Response') +
  ggsave(filename = paste(figprefix, 'ATF2_combined_violin.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(trunc_norm_data[which(trunc_norm_data$antibody == pAb),], aes(x = dose, y = fl, fill = dose)) +
  geom_violin(outlier.color = NULL, position = "dodge") +
  facet_grid(cellline ~ timepoint) + 
  geom_hline(yintercept = 1) 

ggplot(trunc_norm_stats[which(trunc_norm_stats$antibody == pAb),], 
       aes(x = dose, y = stat.median)) + 
  geom_col() + 
  geom_errorbar(aes(ymin = stat.median-stat.mediansd, ymax = stat.median+stat.mediansd)) +
  facet_grid(cellline~timepoint) + 
  geom_point(data = trunc_norm_runs2[which(trunc_norm_runs2$antibody == pAb & trunc_norm_runs2$dose != 'CTRL'),],
             aes(y = log.median, color = experiment, shape = experiment, group = experiment), position = 'dodge', size = 3)

ggplot(trunc_norm_runs2[which(trunc_norm_runs2$dose != 'CTRL'),], aes(x = log.mean, y = log.median)) + 
  geom_point(aes(color = dose), size = 3) + 
  facet_wrap(~cellline) + 
  geom_abline(slope = 1, intercept = 0)



write.csv(trunc_norm_runs2, file = 'FinalRuns_181030.csv', sep = ',', quote = FALSE)
write.csv(ddply(trunc_norm_runs2, .(antibody, timepoint, cellline, dose), summarize, 
      ave_median = mean(log.median), 
      sd_median = sd(log.median)), 
      file = 'FinalRuns_Summary_181030.csv', sep = ',', quote = FALSE)


# inspect the data for systematic shifts
# plotting variables
pAb = 'H2aX'
pcell = 'HCC'

abcellplot(pAb, pcell)


ggplot(runs[which(runs$cellline == 'HCC'),], aes(x = replicate, y = cellcount)) + 
  geom_col(aes(fill = dose), position = 'dodge') + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = cellnumber)

runs[which(runs$cellcount < cellnumber),1:5]

# 
# 
# trunc_log_data <- ddply(tlog_df, .(antibody, cellline, timepoint, dose, replicate), 
#           function(x) head(x, cellnumber))
# 
# trunc_runs <- ddply(trunc_log_data, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
#               log.mean = round(mean(fl), 3),
#               log.median = round(median(fl),3),
#               sd = round(sd(fl), 3))
# trunc_log_ctrl_means <- ddply(trunc_runs[which(runs$dose == 'CTRL'),],
#                         .(antibody, cellline, timepoint), summarize,
#                         mean = round(mean(log.mean), 3),
#                         median = round(mean(log.median),3),
#                         sd = round(sd(log.mean), 3))
# 
# 
# 
# plot(runs$log.mean[order(runs$log.mean)], trunc_runs$log.mean[order(trunc_runs$log.mean)])
# plot(runs$log.median[order(runs$log.median)], trunc_runs$log.median[order(trunc_runs$log.median)])
# plot(runs$sd[order(runs$sd)], trunc_runs$sd[order(trunc_runs$sd)])
# 
# 
# 
# 
# 
# 
# mean(tlog_df$fl[which(tlog_df$antibody == 'ATF2' & tlog_df$dose == 'CTRL' &
#                         tlog_df$experiment == 'U')])
# 
# 
# 
# 
# 
# 
# norm_runs <- ddply(trunc_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
#                            log.mean = round(mean(fl), 3),
#                            sd = round(sd(fl), 3))
# ctrl_means <- ddply(norm_runs[which(runs$dose == 'CTRL'),],
#                         .(antibody, cellline, timepoint), summarize,
#                         mean = round(mean(log.mean), 3),
#                         sd = round(sd(log.mean), 3))


pAb = 'ATF2'
pcell = 'SKBR3'
tabcellplot(pAb, pcell)




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

# this part is totally not working right now, need to evaluate HOW to do stats testing
# define the sample size n
n = 4000

Dose = "HD"
Cell = "HCC"
Time = "24h"
Ab = "H2aX"


indexer <- which(trunc_norm_data$cellline == Cell &
                   trunc_norm_data$timepoint == Time &
                   trunc_norm_data$antibody == Ab)

ref_indexer <- which(trunc_norm_data$dose == Dose &
                       trunc_norm_data$cellline == Cell &
                       trunc_norm_data$timepoint == Time &
                       trunc_norm_data$antibody == Ab
                     #                  & trunc_norm_data$replicate == 'S1'
)

fitdata <- trunc_norm_data[indexer,]
refdata <- sample_n(trunc_norm_data[ref_indexer,], n)

disdata <- fitdist(fitdata$fl, 'norm')
refdisdata <- fitdist(refdata$fl, 'norm')
plot(disdata)
plot(refdisdata)

# Kolmogorov-Smirnov test for normality (D closer to 0 is more normal / the same distribution)
ks.test(unique(refdisdata$data), 'pnorm')
ks.test(unique(disdata$data), unique(refdisdata$data))

#densityplot(fitdata$fl)
#scores(disdata$data)
#plot(ecdf(disdata$data))
kruskal.test(fl ~ dose, data = fitdata)
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
kruskal.test(fl ~ dose, data = trunc_norm_data[indexer,])

str(subdata)

# T-test

t.test(x = trunc_norm_data$fl[which(trunc_norm_data$dose == "HD" &
                                      trunc_norm_data$cellline == "HCC" &
                                      trunc_norm_data$timepoint == "24h" &
                                      trunc_norm_data$antibody == 'H2aX')],
       y = trunc_norm_data$fl[which(trunc_norm_data$dose == "LD" &
                                      trunc_norm_data$cellline == "HCC" &
                                      trunc_norm_data$timepoint == "24h" &
                                      trunc_norm_data$antibody == 'H2aX')],
       conf.level = 0.95)

# ANOVA

ano_hcc_h2ax <- aov(data = trunc_norm_data[which(trunc_norm_data$antibody == 'ATF2' 
                                                 & trunc_norm_data$cellline == 'BT549' 
                                           & trunc_norm_data$timepoint == '24h 24h'),], 
                    fl ~ dose)

summary(glht(ano_hcc_h2ax, linfct=mcp(dose = "Dunnett")))

anova(ano_hcc_h2ax)

ano_hcc_h2ax <- lm(data = trunc_norm_runs2[which(trunc_norm_runs2$antibody == 'ATF2' 
                                                 #& trunc_norm_runs2$cellline == 'HCC' 
                                                 #& trunc_norm_runs2$dose == 'HD'
                                                 & trunc_norm_runs2$timepoint == '24h 24h'
                                                 ),], 
                    log.median ~ dose + cellline + dose:cellline)

#summary(ano_hcc_h2ax)
summary(glht(ano_hcc_h2ax, linfct=mcp(dose = "Tukey")))
Anova(ano_hcc_h2ax, type = 'II')
summary(ano_hcc_h2ax)




# ==== PLOT THE DATA ====
## @knitr dataload 
load("/Volumes/Seagate Backup Plus Drive/Projects/DNA Damage/FACS/R/Workspaces/180821_4000cells.RData")


# ---- Plot the raw data ----
figprefix <- 'Figures/Deepa/'

ggplot(tr_df[which(tr_df$dose == 'CTRL' & tr_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggsave(filename = paste(figprefix, 'controls_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$dose == 'CTRL' & tr_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  xlim(4e4,3e5) +
  ggtitle("raw controls H2aX") +
  ggsave(filename = paste(figprefix, 'controls_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'BT549_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'BT549_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'SKBR3_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'SKBR3_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'ATF2' & tr_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'HCC_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tr_df[which(tr_df$antibody == 'H2aX' & tr_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  xlim(4e4,3e5)
ggsave(filename = paste(figprefix, 'HCC_H2aX.pdf', sep = ""),
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

figprefix <- 'Figures/merged_expts/log/'

ggplot(tlog_df[which(tlog_df$dose == 'CTRL' & tlog_df$antibody == 'ATF2'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggsave(filename = paste(figprefix, 'controls_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")



ggplot(tlog_df[which(tlog_df$dose == 'CTRL' & tlog_df$antibody == 'H2aX'), ], aes(fl, fill = replicate, by = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ cellline) +
  ggtitle('log congrols H2aX') +
  ggsave(filename = paste(figprefix, 'controls_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'BT549_ATF2_bimodal.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'BT549'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'BT549_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")


ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'SKBR3_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'SKBR3'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'SKBR3_H2aX.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'ATF2' & tlog_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'HCC_ATF2.pdf', sep = ""),
         width = 8.5, height = 5.5, units = "in")

ggplot(tlog_df[which(tlog_df$antibody == 'H2aX' & tlog_df$cellline == 'HCC'), ], aes(fl, fill = replicate)) +
  geom_density(alpha = alp,  adjust = bw) +
  facet_grid(timepoint ~ dose) +
  ggsave(filename = paste(figprefix, 'HCC_H2aX.pdf', sep = ""),
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
  ggtitle('H2aX') + 
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


ggplot(tdf, aes(x = dose, y = fl, fill = cellline, by = experiment, color = experiment)) +
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




ggplot(trunc_runs[which(trunc_runs$antibody == 'H2aX' & trunc_runs$cellline == 'HCC'),], aes(x = log.mean)) + 
  geom_histogram(aes(fill = dose))



mean(trunc_runs$log.mean[which(trunc_runs$antibody == 'H2aX' & trunc_runs$cellline == 'HCC' & trunc_runs$dose != 'CTRL')])


