SUBJID <- 303
session <- "fMRI" # "beh" or "fMRI"
########################################################################################

library(plyr)
library(tidyverse)


if (session == "fMRI"){
  path <- paste("data/fMRI/", SUBJID, "ordch_fmri.txt", sep = "")
}else{
  path <- paste("data/behavioral/", SUBJID, "ordch_beh.txt", sep = "")
}

data <- read.table(path, header = TRUE)
data <- data[data$SUBID != "SUBID",] 

if (session != "fMRI"){
  data <- data[data$BLOCK %in% seq(from = 2, to = 72, by = 2),] #include only the test blocks in the beh cond
}

data$BLOCK <- as.numeric(data$BLOCK)

#calculate blockwise accuracy
blockwise <- ddply(data, .(BLOCK), summarize, ACC = mean(as.numeric(ACC)))

#plot blockwise accuracy
ggplot(blockwise, aes(x = BLOCK, y = ACC, group = 1))+
  geom_point()+
  geom_line()+
  theme_bw()

#return overall accuracy
print(paste("Overall accuracy is ", round(mean(as.numeric(data$ACC))*100 ,3), " %.", sep = ""))


