
#data was generated using fungal antiSMASH 
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
names<- rownames(SMCs_DF_Suillus3)
SMCs_DF_Suillus4<- cbind(SMCs_DF_Suillus3, names, rep(("a_Suillus"), length(SMCs_DF_Suillus3)))
colnames(SMCs_DF_Suillus4)<- c("counts", "names", "group")


SMCs_DF_Other<- SMCs_DF[SMCs_DF$group == "O",]
SMCs_DF_Other2<- SMCs_DF_Other[,5:20]
SMCs_DF_Other3<- data.frame(round(colMeans(SMCs_DF_Other2),0))
names<- rownames(SMCs_DF_Other3)
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
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
spineplot(counts_gather_table, main = " ", col=c("#5677A1", "#405952", "#9B9A79", "#FED393","#F79552", "#F0502B", "#A82B0E"), border = NA, cex.axis = 0.66, ylab = "", xlab = "", yaxlabels = NA, xaxlabels = groups)
#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , fill=c( "#5677A1", "#405952", "#9B9A79", "#FED393","#F79552", "#F0502B", "#A82B0E"),
       ncol = 1, cex = .4, lwd = 0, box.lwd = 0, box.lty =0, box.col =0, xjust = 1)

#get n and N
counts_gather_table
N<- rowSums(counts_gather_table)

#make grouped boxplots
#example

#seperate Suillus and Other from counts_gather2

#names unique(counts_gather2$names)
counts_gather2_Suillus<- counts_gather2[counts_gather2$group == "a_Suillus",]
counts_gather2_Other<- counts_gather2[counts_gather2$group == "b_Other",]






#make groped box plot
#example
library(ggplot2)

#first put into long form
library(tidyr)
library(dplyr)
#add repeat col (ours is only 1)
SMC_DF_small_no_zeros$id<- seq(1:nrow(SMC_DF_small_no_zeros))
SMC_DF_small_no_zeros_long <- SMC_DF_small_no_zeros %>% group_by(group) %>%
  gather(data = SMC_DF_small_no_zeros, id, cf_fatty_acid, cf_putative, nrps, other, t1pks, terpene, indole)

#check class
class(SMC_DF_small_no_zeros[1,1])

#make col names match
colnames(SMC_DF_small_no_zeros_long)<- c("Group", "SMC", "Count")


#make grouped box plot
p<- ggplot(SMC_DF_small_no_zeros_long, aes(x=SMC, y=Count, fill=Group)) + 
  geom_boxplot()
p+scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()

#add points?
p + geom_point(position=position_jitterdodge(),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()



help(position_jitterdodge)

p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()


#OR, do we want to make this horizontal with the colors to match the SMC categories?
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal() + coord_flip()

#how to set colors by shade? 





#####stats start here 
#####ANOVA 
#shrink to relevent content 
SMC_DF_small_no_zeros<- SMC_DF_small[, c(1:6, 9, 17)]

#test for variance heterogenety 
var.test((as.numeric(SMC_DF_small_no_zeros$cf_fatty_acid) [SMC_DF_small_no_zeros$group == "S"]), (as.numeric(SMC_DF_small_no_zeros$cf_fatty_acid) [SMC_DF_small_no_zeros$group == "O"]))
#not equal variance 

library(outliers)
cochran.test(as.vector(SMC_DF_small_no_zeros[,1]) ~ SMC_DF_small_no_zeros$group, inlying = FALSE)
help(cochran.test)


#run cochrans in a loop
for (i in SMC_DF_small_no_zeros) {
  print(var.test(as.formula(paste(i,'~group')),data))
}


#cochran.test(object, data, inlying = FALSE)
help(cochran.test)








#set factors
options(contrasts=c("contr.sum","contr.poly"))
SMC_DF_small_no_zeros_long$SMC <- factor(SMC_DF_small_no_zeros_long$SMC,levels=c("cf_putative", "nrps", "other", "t1pks", "terpene", "indole", "cf_fatty_acid"))
names(SMC_DF_small_no_zeros_long) <- c("Group", "SMC", "Count")
SMC_DF_small_no_zeros_long$Group <- factor(SMC_DF_small_no_zeros_long$Group, levels=c("S", "O"))


aov_results<- aov(SMC_DF_small_no_zeros_long$Count ~ SMC_DF_small_no_zeros_long$SMC * SMC_DF_small_no_zeros_long$Group, data = SMC_DF_small_no_zeros_long)
summary(aov_results)

TukeyHSD(aov_results)


