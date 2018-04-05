#---- Header content ----
library(ggplot2)
library(BiocInstaller)
library(ggcyto)
library(flowCore)
library(outliers)

source("https://bioconductor.org/biocLite.R")
biocLite()

setwd('./')
#---- Data import ----
data_4h_BT549_gH2aX_S1_CTRL <- read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S1 CTRL.csv")[,3]
data_4h_BT549_gH2aX_S1_HD <- read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S1 HD.csv")[,3]
# data_4h_24h_BT549_gH2aX_S1_LD <- read.csv("Data/4h_24h BT549  gH2aX/4h_24h BT549  gH2aX_S1 LD.csv")[,3]
data_4h_BT549_gH2aX_S2_CTRL <- read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S2 CTRL.csv")[,3]
data_4h_BT549_gH2aX_S2_HD <- read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S2 HD.csv")[,3]
# data_4h_24h_BT549_gH2aX_S2_LD <- read.csv("Data/4h_24h BT549  gH2aX/4h_24h BT549  gH2aX_S2 LD.csv")[,3]
#data_4h_BT549_gH2aX_S3_CTRL <- read.csv("Data/4h_gH2aX_BT549/4h_gH2aX_BT549_S3 CTRL.csv")[,3]
# data_4h_24h_BT549_gH2aX_S3_HD <- read.csv("Data/4h_24h BT549  gH2aX/4h_24h BT549  gH2aX_S3 HD.csv")[,3]
# data_4h_24h_BT549_gH2aX_S3_LD <- read.csv("Data/4h_24h BT549  gH2aX/4h_24h BT549  gH2aX_S3 LD.csv")[,3]

S1_CTRL <- data.frame(data_4h_BT549_gH2aX_S1_CTRL)
colnames(S1_CTRL) <- "Alexa488"

S1_CTRL$Alexa488 <- signif(log(S1_CTRL$Alexa488), digits = 4)

ggplot(data = S1_CTRL, aes(x =Alexa488)) + 
  geom_bar()

median(S1_CTRL$Alexa488)


#---- S2 Control ----
S2_CTRL <- data.frame(data_4h_BT549_gH2aX_S2_CTRL)
colnames(S2_CTRL) <- "Alexa488"

S2_CTRL$Alexa488 <- signif(log(S2_CTRL$Alexa488), digits = 4)

ggplot(data = S2_CTRL, aes(x =Alexa488)) + 
  geom_bar()

median(S2_CTRL$Alexa488)

#---- Master control ----
# make the average of the medians of the control data sets
ctrl_abs_mean <- mean(c(median(S1_CTRL$Alexa488), median(S2_CTRL$Alexa488)))

# adjust the control data sets to the absolute mean using the median absolute deviation
norm_S1_CTRL <- data.frame(norm = S1_CTRL$Alexa488 - mad(S1_CTRL$Alexa488, constant = ctrl_abs_mean))
norm_S2_CTRL <- data.frame(norm = S2_CTRL$Alexa488 - mad(S2_CTRL$Alexa488, constant = ctrl_abs_mean))

# write the normalization divisors for each set 
S1_norm <- median(norm_S1_CTRL$norm)
S2_norm <- median(norm_S2_CTRL$norm)

# normalize the control sets 
norm_S1_CTRL$norm <- norm_S1_CTRL$norm / S1_norm
norm_S2_CTRL$norm <- norm_S2_CTRL$norm / S2_norm
                           

ggplot(data = norm_S1_CTRL, aes(x = norm)) + 
  geom_bar()
ggplot(data = norm_S2_CTRL, aes(x = norm)) + 
  geom_bar()

#---- S1 High dose ----
S1_HD <- data.frame(data_4h_BT549_gH2aX_S1_HD)
colnames(S1_HD) <- "Alexa488"

S1_HD$Alexa488 <- signif(log(S1_HD$Alexa488), digits = 4)

ggplot(data = S1_HD, aes(x =Alexa488)) + 
  geom_bar()

# adjust the data set using the S1 control MAD and normalize using the S1 control factor
norm_S1_HD <- data.frame(norm = S1_HD$Alexa488 - mad(S1_CTRL$Alexa488, constant = ctrl_abs_mean)) 
norm_S1_HD$norm <- norm_S1_HD$norm / S1_norm

ggplot(data = norm_S1_HD, aes(x = norm)) + 
  geom_bar()
median(norm_S1_HD$norm)


#---- S2 High dose ----
S2_HD <- data.frame(data_4h_BT549_gH2aX_S2_HD)
colnames(S2_HD) <- "Alexa488"

S2_HD$Alexa488 <- signif(log(S2_HD$Alexa488), digits = 4)

ggplot(data = S2_HD, aes(x =Alexa488)) + 
  geom_bar()


# adjust the data set using the S1 control MAD and normalize using the S1 control factor
norm_S2_HD <- data.frame(norm = S2_HD$Alexa488 - mad(S2_CTRL$Alexa488, constant = ctrl_abs_mean)) 
norm_S2_HD$norm <- norm_S2_HD$norm / S2_norm

ggplot(data = norm_S2_HD, aes(x = norm)) + 
  geom_bar()
median(norm_S2_HD$norm)

ggplot(data = norm_S1_CTRL, aes(x = norm)) + 
  geom_bar(aes(fill = 'red')) + 
  geom_bar(data = norm_S1_HD, aes(x = norm, fill = 'green'))

S1_norm
S2_norm
