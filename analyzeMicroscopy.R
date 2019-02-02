#---- Header content ----
## @knitr header
library(ggplot2)
#library(BiocInstaller)
#library(ggcyto)
#library(flowCore)
library(outliers)
library(fitdistrplus)
#library(clusterSim)
library(plyr)
library(dplyr)
library(multcomp)
library(plotly)



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

unbinned_data <- binned_data[0,-1:-4]

pb = txtProgressBar(min = 1, max = length(csvfiles), style = 3)
for (i in 2:length(rdata)){
  if (ncol(rdata[[i]]) == 7){
    binned_data <- rbind(binned_data, cbind(rdata[[i]], 
                               Antibody = antibody[i], 
                               Cellline = cellline[i],
                               Timepoint = timepoint[i], 
                               Parameter = parameter[i]))
  }
  if (ncol(rdata[[i]]) == 3){
    unbinned_data <- rbind(unbinned_data, cbind(rdata[[i]], 
                                            Antibody = antibody[i], 
                                            Cellline = cellline[i],
                                            Timepoint = timepoint[i], 
                                            Parameter = parameter[i]))
  }
  setTxtProgressBar(pb, i)
}

colnames(binned_data) <- c('Dose','V1','V2','V3','SE1','SE2','SE3',colnames(binned_data)[8:11])
colnames(unbinned_data) <- c('Dose','value','SE',colnames(unbinned_data)[4:7])

unique(parameter)

subdata <- dplyr::filter(unbinned_data, Cellline == unique(cellline[1]))

ggplot(subdata, aes(x = Parameter, y = value)) + 
  geom_point() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
