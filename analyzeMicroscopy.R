#---- Header content ----
## @knitr header
library(ggplot2)
#library(BiocInstaller)
#library(ggcyto)
#library(flowCore)
library(outliers)
library(reshape2)
library(fitdistrplus)
#library(clusterSim)
library(plyr)
library(dplyr)
library(multcomp)



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

se <- function(x) sqrt(var(x)/length(x))


setwd('./')

#---- Data import and cleaning from the listed folders ----
datafold <- 'confocal_data/csv_data/'
csvfiles <- dir(datafold)

xlsheets <- array()
rdata <- list()
filenames <- array()
k = 0

pb = txtProgressBar(min = 1, max = length(csvfiles), style = 3)
j = 1
for (i in csvfiles){
  filenames[j] <- i
  rdata[[i]] <- data.frame(read.csv(paste(datafold, i, sep = ''), sep = ',', header = FALSE))
  j = j+1
}
  
timepoint <- array(dim = length(rdata))
cellline <- array(dim = length(rdata))
antibody <- array(dim = length(rdata))
parameter <- array(dim = length(rdata))

for (i in seq(length(filenames))){
  filelength <- length(unlist(strsplit(filenames[[i]], split = "_")))
   timepoint[i] <- unlist(strsplit(unlist(strsplit(filenames[i], split = "_"))[filelength], split = "\\."))[1]
   cellline[i] <- unlist(strsplit(filenames[i], split = "_"))[1]
   antibody[i] <- unlist(strsplit(filenames[i], split = "_"))[3]
   parameter[i] <- paste(unlist(strsplit(filenames[i], split = "_"))[4:(filelength-2)], collapse = '_')
}


binned_data <- cbind(rdata[[1]], 
                 Antibody = antibody[1], 
                 Cellline = cellline[1],
                 Timepoint = timepoint[1], 
                 Parameter = parameter[1])
levels(binned_data$Antibody) <- unique(antibody)
levels(binned_data$Cellline) <- unique(cellline)
levels(binned_data$Timepoint) <- unique(timepoint)
levels(binned_data$Parameter) <- unique(parameter)

binned_data4 <- cbind(binned_data[0,1:7], binned_data[0,1:2],binned_data[0,7:11])
unbinned_data <- binned_data[0,-1:-4]

pb = txtProgressBar(min = 1, max = length(csvfiles), style = 3)
for (i in 2:length(rdata)){
  if (ncol(rdata[[i]]) == 7){
    binned_data <- rbind(binned_data, cbind(rdata[[i]], 
                               Antibody = antibody[i], 
                               Cellline = cellline[i],
                               Timepoint = timepoint[i], 
                               Parameter = parameter[i]))
  } else if (ncol(rdata[[i]]) == 3){
    unbinned_data <- rbind(unbinned_data, cbind(rdata[[i]], 
                                            Antibody = antibody[i], 
                                            Cellline = cellline[i],
                                            Timepoint = timepoint[i], 
                                            Parameter = parameter[i]))
  } else if (ncol(rdata[[i]]) == 9){
    binned_data4 <- rbind(binned_data4, cbind(rdata[[i]], 
                                                Antibody = antibody[i], 
                                                Cellline = cellline[i],
                                                Timepoint = timepoint[i], 
                                                Parameter = parameter[i]))
  } else {
    print(i)
  }
  setTxtProgressBar(pb, i)
}

colnames(binned_data) <- c('Dose','V1','V2','V3','SE1','SE2','SE3',colnames(binned_data)[8:11])
colnames(binned_data4) <- c('Dose','V1','V2','V3','V4','SE1','SE2','SE3','SE4',colnames(binned_data4)[10:13])
colnames(unbinned_data) <- c('Dose','value','SE',colnames(unbinned_data)[4:7])


binned_data$Dose <- factor(binned_data$Dose, 
                           levels = levels(binned_data$Dose)[c(1,3,2)])
binned_data4$Dose <- factor(binned_data4$Dose, 
                           levels = levels(binned_data4$Dose)[c(1,3,2)])
unbinned_data$Dose <- factor(unbinned_data$Dose, 
                           levels = levels(unbinned_data$Dose)[c(1,3,2)])


#### START MAKING FIGURES ####


unbinpars <- c(3,13,8,18)

plotdata <- dplyr::filter(unbinned_data, 
                           Antibody == unique(antibody)[2] &
                           Parameter == unique(unbinned_data$Parameter)[18] & 
                           Cellline != 'HCCold')


pos <- position_dodge(width = 0.9)
ggplot(plotdata, aes(x = Dose, y = value, by = Dose, ymin = value - SE, ymax = value + SE)) + 
  geom_col(aes(fill = Dose), position = pos) + 
  geom_errorbar(width = 0.2, position = pos) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  facet_grid(Cellline~Timepoint) + 
  ggtitle(paste(plotdata$Antibody, plotdata$Parameter, sep = ' ')) + 
  ggsave(filename = paste('confocal_figures/barplot',
                          plotdata$Antibody[1],
                          plotdata$Parameter[1],
                          'byCellTime.pdf',
                          sep = '_'),
         width = 8.5, height = 5.5, units = "in")




plotdata <- dplyr::filter(binned_data, Antibody == unique(antibody)[3] &
                            Parameter == unique(binned_data$Parameter)[3] & 
                            Cellline != 'HCCnew')

mplotdata <- melt(plotdata[,-5:-7], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))
mplotdata <- cbind(mplotdata, 
                   SEvar = melt(plotdata[,-2:-4], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))[,6],
                   SEval = melt(plotdata[,-2:-4], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))[,7])
cplotdata <- ddply(mplotdata, Dose + Antibody + Timepoint + Cellline + Parameter ~ variable)


pos <- position_dodge(width = 0.9)
ggplot(cplotdata, aes(x = Dose, y = value, by = variable, ymin = value - SEval, ymax = value + SEval)) + 
  geom_col(aes(fill = variable), position = pos) + 
  geom_errorbar(width = 0.2, position = pos) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  facet_grid(Cellline~Timepoint) + 
  ggtitle(paste(plotdata$Antibody, plotdata$Parameter, sep = ' ')) + 
  ggsave(filename = paste('confocal_figures/barplot',
                          plotdata$Antibody[1],
                          plotdata$Parameter[1],
                          'bin_byCellTime.pdf',
                          sep = '_'),
         width = 8.5, height = 5.5, units = "in")



bin4pars <- c(1)

plotdata <- dplyr::filter(binned_data4, Antibody == unique(antibody)[2] &
                            Parameter == unique(binned_data4$Parameter)[1] & 
                            Cellline != 'HCCold')

mplotdata <- melt(plotdata[,-6:-9], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))
mplotdata <- cbind(mplotdata, 
                   SEvar = melt(plotdata[,-2:-5], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))[,6],
                   SEval = melt(plotdata[,-2:-5], id=c('Dose','Antibody','Timepoint','Parameter','Cellline'))[,7])
cplotdata <- ddply(mplotdata, Dose + Antibody + Timepoint + Cellline + Parameter ~ variable)


pos <- position_dodge(width = 0.9)
ggplot(cplotdata, aes(x = Dose, y = value, by = variable, ymin = value - SEval, ymax = value + SEval)) + 
  geom_col(aes(fill = variable), position = pos) + 
  geom_errorbar(width = 0.2, position = pos) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  facet_grid(Cellline~Timepoint) + 
  ggtitle(paste(plotdata$Antibody, plotdata$Parameter, sep = ' ')) + 
  ggsave(filename = paste('confocal_figures/barplot',
                          plotdata$Antibody[1],
                          plotdata$Parameter[1],
                          'dose_byCellTime.pdf',
                          sep = '_'),
         width = 8.5, height = 5.5, units = "in")


pos <- position_dodge(width = 0.9)
ggplot(cplotdata, aes(x = variable, y = value, by = Dose, ymin = value - SEval, ymax = value + SEval)) + 
  geom_col(aes(fill = Dose), position = pos) + 
  geom_errorbar(width = 0.2, position = pos) +
  theme(axis.text.x = element_text(angle = 20, hjust = 1)) + 
  facet_grid(Cellline~Timepoint) + 
  ggtitle(paste(plotdata$Antibody, plotdata$Parameter, sep = ' ')) + 
  ggsave(filename = paste('confocal_figures/barplot',
                          plotdata$Antibody[1],
                          plotdata$Parameter[1],
                          'bin_byCellTime.pdf',
                          sep = '_'),
         width = 8.5, height = 5.5, units = "in")
