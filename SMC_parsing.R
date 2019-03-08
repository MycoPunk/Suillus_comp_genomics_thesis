
#data was generated using GPCRHMM
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#load packages 
library(tidyr)
library(reshape2)

#read in the input files
SMCs_DF<- read.csv("SMCs_Suillus_vs_OtherECM.csv")
#shrink DF
SMC_DF_small<- SMCs_DF[,5:21]

#shrink DF to the important stuff
SMCs_DF_Suillus<- SMCs_DF[SMCs_DF$group == "S",]
SMCs_DF_Suillus2<- SMCs_DF_Suillus[,5:20]
SMCs_DF_Suillus3<- data.frame(round(colMeans(SMCs_DF_Suillus2),0))
names<- rownames(SMCs_DF_Suillus4)
SMCs_DF_Suillus4<- cbind(SMCs_DF_Suillus3, names, rep(("a_Suillus"), length(SMCs_DF_Suillus3)))
colnames(SMCs_DF_Suillus4)<- c("counts", "names", "group")



SMCs_DF_Other<- SMCs_DF[SMCs_DF$group == "O",]
SMCs_DF_Other2<- SMCs_DF_Other[,5:20]
SMCs_DF_Other3<- data.frame(round(colMeans(SMCs_DF_Other2),0))
names<- rownames(SMCs_DF_Other4)
SMCs_DF_Other4<- cbind(SMCs_DF_Other3, names, rep(("b_Other"), length(SMCs_DF_Other3)))
colnames(SMCs_DF_Other4)<- c("counts", "names", "group")


#bind them 
SMCs_DF<- rbind(SMCs_DF_Suillus4, SMCs_DF_Other4)

#make repeating rows by SMC counts
counts_gather2 <- as.data.frame(lapply(SMCs_DF, rep, SMCs_DF$counts))

counts_gather_table<- table(counts_gather2$group, counts_gather2$names)

groups<- c("Suillus", "Other")
#spine plot
par(las = 2)
spineplot(counts_gather_table, main = " ", col=c("#5677A1", "#405952", "#9B9A79", "#FED393","#F79552", "#F0502B", "#A82B0E"), border = NA, cex.axis = 0.66, ylab = "", xlab = "", yaxlabels = NA, xaxlabels = groups)
#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , fill=c( "#5677A1", "#405952", "#9B9A79", "#FED393","#F79552", "#F0502B", "#A82B0E"),
       ncol = 1, cex = .4, lwd = 0, box.lwd = 0, box.lty =0, box.col =0, xjust = 1)



dim(counts_gather_table)



#or do we want it equalized?
barplot(t(prop.table(counts_gather_table,1)), col=c("#4F5C69", "#838558", "#A49C4C","#E4D9AC", "#F9F1D2", "#FBE898", "#FFCC6E", "#E6A871"), 
        border = NA, cex.axis = 0.66, names = c("Suillus", "Other ECM"))




