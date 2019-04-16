
#data was generated using GPCRHMM
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#load packages 
library(tidyr)
library(reshape2)

#read in the input files
SMCs_DF<- read.csv("SMCs_Suillus_vs_OtherECM.csv")


#shrink DF
SMC_DF_small<- SMCs_DF[,5:21]

View(SMC_DF_small)
#SMC_DF_small_no_zeros<- SMC_DF_small[, c(1:6, 9, 17)]

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
#wants to inherit zeros in the following- get rid of them by switching types (I dono)
counts_gather3<- as.matrix(counts_gather2)
counts_gather4<- as.data.frame(counts_gather3)

#call table
counts_gather_table<- table(counts_gather4$group, counts_gather4$names)

groups<- c("Suillus", "Other")
#spine plot
par(las = 1)
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
spineplot(counts_gather_table, main = " ", 
          col=c("#5676A1", "#405952","#969574", "#FDD191", "#885053", "#B2675E", "#442B47"), 
          border = NA, 
          #cex.axis = 0.66, 
          #ylab = "", #xlab = "", 
          yaxlabels = NA, 
          xaxlabels = groups)

#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , 
       fill=c("#5676A1", "#405952","#969574", "#FDD191", "#885053", "#B2675E", "#442B47"),
       ncol = 1, 
       cex = .4, 
       lwd = 0, 
       box.lwd = 0, 
       box.lty =0, 
       box.col =0, 
       xjust = 1)


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

#add points
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal()


#horizontal
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666")) + theme_minimal() + coord_flip()


#####stats start here 
#####ANOVA 

#inspect data. 
#test for variance heterogenety 
var.test((as.numeric(SMC_DF_small_no_zeros$cf_fatty_acid) [SMC_DF_small_no_zeros$group == "S"]), (as.numeric(SMC_DF_small_no_zeros$cf_fatty_acid) [SMC_DF_small_no_zeros$group == "O"]))
#not equal variance 

View(SMC_DF_small_no_zeros)

#build model
genome_size_model <- lm(all_df$genome_size ~ all_df$group)

#test normality and var
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(genome_size_model, lambda = seq(-4, 1, 1/10)) #what does boxcox suggest? log transformation works. 

#log transform to improve normality?
all_df$genome_size_log<- log(all_df$genome_size)

#see if the ransformation helped
genome_size_model2 <- lm(all_df$genome_size_log ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#not a lot, try -2 ==(1/x^2) - inverse square transformation 

all_df$genome_size_inv_suqare<- all_df$genome_size^(-2)

genome_size_model3 <- lm(all_df$genome_size_inv_suqare ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(genome_size_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#that looks a little better 

#Run anova
summary(aov(genome_size_model3))







#set factors
options(contrasts=c("contr.sum","contr.poly"))
SMC_DF_small_no_zeros_long$SMC <- factor(SMC_DF_small_no_zeros_long$SMC,levels=c("cf_putative", "nrps", "other", "t1pks", "terpene", "indole", "cf_fatty_acid"))
names(SMC_DF_small_no_zeros_long) <- c("Group", "SMC", "Count")
SMC_DF_small_no_zeros_long$Group <- factor(SMC_DF_small_no_zeros_long$Group, levels=c("S", "O"))


aov_results<- aov(SMC_DF_small_no_zeros_long$Count ~ SMC_DF_small_no_zeros_long$SMC * SMC_DF_small_no_zeros_long$Group, data = SMC_DF_small_no_zeros_long)
summary(aov_results)

TukeyHSD(aov_results)

#type 2 ANOVA
library("car")
#make model
model.for.t2<-lm(SMC_DF_small_no_zeros_long$Count ~ SMC_DF_small_no_zeros_long$SMC * SMC_DF_small_no_zeros_long$Group)

model.for.t3<- Anova(model.for.t2, type=2)

summary(SMC_DF_small_no_zeros_long)

#run TukeyHSD
library("cfcdae")

#define interaction term
SMC_DF_small_no_zeros_long$Group_SMC_interaction <- as.factor(interaction(SMC_DF_small_no_zeros_long$SMC, SMC_DF_small_no_zeros_long$Group))

#re-run model with interaction term
model.for.t4<- lm(SMC_DF_small_no_zeros_long$Count ~ SMC_DF_small_no_zeros_long$Group_SMC_interaction)

#take a look - lines in common mean not significantly different 
sidelines(pairwise(model.for.t4, SMC_DF_small_no_zeros_long$Group_SMC_interaction,confidence = 0.95, type = "hsd"))

#full result with p-vals:
pairwise(model.for.t4, SMC_DF_small_no_zeros_long$Group_SMC_interaction,confidence = 0.95, type = "hsd")

#full result with p-vals:
pairs_out<- pairwise(model.for.t4, SMC_DF_small_no_zeros_long$Group_SMC_interaction,confidence = 0.95, type = "regwr")

sidelines(pairs_out)
#terpenes and Other are dignificantly different 

#get p-value
p_vals<-with(SMC_DF_small_no_zeros_long, pairwise.t.test(SMC_DF_small_no_zeros_long$Count,SMC_DF_small_no_zeros_long$Group_SMC_interaction,
                                     p.adjust.method="holm"))

