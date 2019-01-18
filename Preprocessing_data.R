# Preprocess data for Bayesian model fit
# Jan-17-2019
# Set up ----------------------------------------------------------------------
rm(list=ls()) # Clears everything from environment

# Below: function to load necessary packages and install if necessary
init <- function(need) {
  ip <- .packages(all.available = T)
  if (any((need %in% ip) == F)) {
    install.packages(need[!(need %in% ip)])
  }
  ok <- sapply(1:length(need), function(p) require(need[[p]], 
                                                   character.only = T))
}

init(c("plyr","reshape2","car"))

# Import files -----------------------------------------------------------------


cidrs18 <- read.csv("cidrs18_fulldatascreened.csv",header=TRUE)
head(cidrs18)

cidrs18 <- cidrs18[cidrs18$ExclReason=="valid",]

cidrs18binnedfreq <- melt(cidrs18, id.vars = c("SubjID","Condition","StimType"), measure.vars = c("RecConf"))
cidrs18binnedfreq2 <- dcast(cidrs18binnedfreq, Condition + SubjID + StimType ~ value, length)

cidrs18Item<-cidrs18binnedfreq2[cidrs18binnedfreq2$Condition==1,]

cidrs18Item_info <- list(
  numberParticipants = length(unique(cidrs18Item$SubjID)),
  testData_signalitems = as.matrix(cidrs18Item[cidrs18Item$StimType=="1",as.character(1:6)]),
  testData_noiseitems = as.matrix(cidrs18Item[cidrs18Item$StimType=="2",as.character(1:6)])
)

numberParticipants <- length(unique(cidrs18Item$SubjID))

cidrs18Context<-cidrs18binnedfreq2[cidrs18binnedfreq2$Condition==2,]

cidrs18Context_info <- list(
  numberParticipants = length(unique(cidrs18Context$SubjID)),
  testData_signalitems = as.matrix(cidrs18Context[cidrs18Context$StimType=="1",as.character(1:6)]),
  testData_noiseitems = as.matrix(cidrs18Context[cidrs18Context$StimType=="2",as.character(1:6)])
)

numberParticipants <- length(unique(cidrs18Context$SubjID))