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

#---- Start looking at statistics and sets to toss ----

# calculate the mean and SD of all distributions - which are outside of normal shape?
mean(runs$sd)
runs$sd
## @knitr deviation_plots
ggplot(runs, aes(x = antibody, y = sd)) +
  geom_boxplot(aes(fill = cellline))



ggplot(runs, aes(x = antibody, y = sd)) +
  geom_jitter(aes(color = cellline)) +
  facet_grid(~cellline)

##
# plot_ly(runscounts, x = ~antibody, y = ~sd,
#         type = 'box',
#         color = ~cellline)


ggplot(runs, aes(sd)) +
  geom_histogram() +
  facet_grid(cellline~antibody)

ggplot(runs, aes(sd)) +
  geom_density() +
  facet_grid(cellline~antibody)

#---- SET REMOVAL ----

# bulk removal of sets 
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



# build a table of runs to remove, based on e.g. cell count, SD
removals <- matrix(nrow = 0, ncol = 5)
removals <- rbind(removals, c('ATF2','HCC','24h','T3','HD')) # SD over .5
removals <- rbind(removals, c('ATF2','BT549','4h','S2','HD')) # SD over .5
removals <- rbind(removals, c('ATF2','BT549','4h','S2','LD')) # SD over .5
removals <- rbind(removals, c('H2aX','SKBR3','4h 24h','U3','HD')) # not a normal distrubtion
removals <- rbind(removals, c('H2aX','BT549','4h','U4','CTRL')) # very low cell count


cellnumber = 4000




remv <- data.frame(lapply(runs[which(runs$cellcount < cellnumber),1:5], as.character), stringsAsFactors = FALSE)
rbind(removals, remv)


removals <- rbind(removals, c('H2aX','BT549','24h 24h','U2','CTRL'))
removals <- rbind(removals, c('H2aX','BT549','24h 24h','U2','LD'))
removals <- rbind(removals, c('H2aX','BT549','24h 24h','U2','HD'))
removals <- rbind(removals, c('H2aX','BT549','24h 24h','U3','CTRL')) # could possibly replace U3 ctrl with U4
removals <- rbind(removals, c('H2aX','BT549','24h 24h','U3','LD'))
removals <- rbind(removals, c('H2aX','BT549','24h 24h','U3','HD')) 

removals <- rbind(removals, c('ATF2','SKBR3','24h','S1','CTRL'))
removals <- rbind(removals, c('ATF2','SKBR3','24h','S1','LD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h','S1','HD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h','S2','CTRL'))
removals <- rbind(removals, c('ATF2','SKBR3','24h','S2','LD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h','S2','HD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T2','CTRL'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T2','LD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T2','HD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T3','CTRL'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T3','LD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T3','HD'))
removals <- rbind(removals, c('ATF2','SKBR3','24h 24h','T4','CTRL'))

removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S2','CTRL'))
removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S2','LD'))
removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S2','HD'))
removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S3','CTRL'))
removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S3','LD'))
removals <- rbind(removals, c('H2aX','SKBR3','24h 24h','S3','HD'))

removals <- rbind(removals, c('ATF2','HCC','24h 24h','T1','CTRL'))
removals <- rbind(removals, c('ATF2','HCC','24h 24h','T1','LD'))
removals <- rbind(removals, c('ATF2','HCC','24h 24h','T1','HD'))
removals <- rbind(removals, c('ATF2','HCC','24h','S2','CTRL'))
removals <- rbind(removals, c('ATF2','HCC','24h','S2','LD'))
removals <- rbind(removals, c('ATF2','HCC','24h','T3','CTRL'))
removals <- rbind(removals, c('ATF2','HCC','24h','T3','LD'))

removals <- rbind(removals, c('H2aX','HCC','4h','T1','CTRL'))
removals <- rbind(removals, c('H2aX','HCC','4h','T1','LD'))
removals <- rbind(removals, c('H2aX','HCC','4h','T1','HD'))

removals <- rbind(removals, c('H2aX','HCC','4h 24h','T3','CTRL'))
removals <- rbind(removals, c('H2aX','HCC','4h 24h','T3','LD'))
removals <- rbind(removals, c('H2aX','HCC','4h 24h','T3','HD'))
#removals <- rbind(removals, c('ATF2','HCC','24h','S1','LD'))
#removals <- rbind(removals, c('ATF2','HCC','24h','S2','LD'))
#removals <- rbind(removals, c('H2aX','HCC','24h 24h','S1','CTRL'))
#removals <- rbind(removals, c('H2aX','HCC','24h 24h','S2','CTRL'))

# remove sets from tlog_df, rebuild runs & log_ctrl_means
for (r in seq(nrow(removals))){
  tlog_df <- tlog_df[-which(tlog_df$antibody == removals[r,1] &
                              tlog_df$cellline == removals[r,2] &
                              tlog_df$timepoint == removals[r,3] &
                              tlog_df$replicate == removals[r,4] &
                              tlog_df$dose == removals[r,5]),]
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

# inspect the data for systematic shifts
# plotting variables
pAb = 'H2aX'
pcell = 'SKBR3'

abcellplot(pAb, pcell)


ggplot(runs[which(runs$cellline == 'HCC'),], aes(x = replicate, y = cellcount)) + 
  geom_col(aes(fill = dose), position = 'dodge') + 
  facet_grid(antibody~timepoint) + 
  geom_hline(yintercept = cellnumber)

runs[which(runs$cellcount < cellnumber),1:5]



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



plot(runs$log.mean[order(runs$log.mean)], trunc_runs$log.mean[order(trunc_runs$log.mean)])
plot(runs$log.median[order(runs$log.median)], trunc_runs$log.median[order(trunc_runs$log.median)])
plot(runs$sd[order(runs$sd)], trunc_runs$sd[order(trunc_runs$sd)])






mean(tlog_df$fl[which(tlog_df$antibody == 'ATF2' & tlog_df$dose == 'CTRL' &
                        tlog_df$experiment == 'U')])





trunc_norms <- merge(trunc_log_data, ddply(trunc_runs[which(trunc_runs$dose == 'CTRL'),],
                           .(antibody, cellline, timepoint, replicate), summarize,
                           mean = round(mean(log.mean), 3),
                           median = round(mean(log.median),3),
                           sd = round(sd(log.mean), 3)), by = c('antibody','cellline','timepoint','replicate'))

trunc_df <- trunc_norms
str(tdf)
trunc_df$fl <- (trunc_norms$fl-trunc_norms$mean)+1

norm_runs <- ddply(trunc_df, .(antibody, cellline, timepoint, replicate, dose, experiment), summarize,
                           log.mean = round(mean(fl), 3),
                           sd = round(sd(fl), 3))
ctrl_means <- ddply(norm_runs[which(runs$dose == 'CTRL'),],
                        .(antibody, cellline, timepoint), summarize,
                        mean = round(mean(log.mean), 3),
                        sd = round(sd(log.mean), 3))


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
Cell = "SKBR3"
Time = "4h 24h"
Ab = "H2aX"


indexer <- which(tlog_df$dose == Dose &
                   tlog_df$cellline == Cell &
                   tlog_df$timepoint == Time &
                   tlog_df$antibody == Ab &
                   tlog_df$replicate == 'U3')

ref_indexer <- which(tlog_df$dose == Dose &
                       tlog_df$cellline == Cell &
                       tlog_df$timepoint == Time &
                       tlog_df$antibody == Ab
                     #                  & tlog_df$replicate == 'S1'
)

fitdata <- sample_n(tlog_df[indexer,], n, replace = TRUE)
refdata <- sample_n(tlog_df[ref_indexer,], n)

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
## @knitr dataload 
load("/Volumes/Seagate Backup Plus Drive/Projects/DNA Damage/FACS/R/Workspaces/180821_4000cells.RData")


# ---- Plot the raw data ----
figprefix <- 'Figures/merged_expts/raw/'

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
  ggsave(filename = paste(figprefix, 'BT549_ATF2.pdf', sep = ""),
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


sd(norm_runs$log.mean[which(norm_runs$antibody == 'ATF2' & norm_runs$cellline == 'BT549' & norm_runs$dose != 'CTRL')])


ggplot(trunc_df[which(trunc_df$cellline == "HCC"),], aes(x = dose, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) +
  geom_hline(yintercept = 1) + 
  geom_hline(yintercept = 1.5 * 0.2185)


ggplot(trunc_df[which(trunc_df$cellline == "BT549"),], aes(x = replicate, y = fl, fill = dose)) +
  geom_boxplot(notch = TRUE, notchwidth = 0.25, outlier.color = NULL, position = "dodge") +
  facet_grid(antibody ~ timepoint) + 
  geom_hline(yintercept = 1.5 * 0.194 + 0.936) + 
  geom_hline(yintercept = 0.936 - 1.5 * 0.194)
  
  


ggplot(trunc_runs[which(trunc_runs$antibody == 'H2aX' & trunc_runs$cellline == 'HCC'),], aes(x = log.mean)) + 
  geom_histogram(aes(fill = dose))



mean(trunc_runs$log.mean[which(trunc_runs$antibody == 'H2aX' & trunc_runs$cellline == 'HCC' & trunc_runs$dose != 'CTRL')])


