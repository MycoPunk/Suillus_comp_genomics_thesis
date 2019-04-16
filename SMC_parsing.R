
#data was generated using fungal antiSMASH
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#load packages 
library(tidyr)
library(reshape2)
library(dplyr)

#read in the input files
SMCs_DF<- read.csv("SMCs_Suillus_vs_OtherECM.csv")


#shrink DF
SMC_DF_small<- SMCs_DF[,5:21]

SMC_DF_small_no_zeros<- SMC_DF_small[, c(1:6, 9, 17)]

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

#full result with hsd -- note you can't get p-vals with hsd though with cfcdae. 
pairwise(model.for.t4, SMC_DF_small_no_zeros_long$Group_SMC_interaction,confidence = 0.95, type = "hsd")

#run model for multiple t tests so that you can get p-vals
pairs_out<- pairwise(model.for.t4, SMC_DF_small_no_zeros_long$Group_SMC_interaction,confidence = 0.95, type = "regwr")

sidelines(pairs_out)
#terpenes and Other are dignificantly different 

#get p-value using mulitple t tests and holm adjustment for multiple comparisons. 
p_vals<-with(SMC_DF_small_no_zeros_long, pairwise.t.test(SMC_DF_small_no_zeros_long$Count,SMC_DF_small_no_zeros_long$Group_SMC_interaction,
                                     p.adjust.method="holm"))





###### Within Suillus by host comparison ######

#parse dataframe by host association
#red
red_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "R",]
#white
white_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "W",]
#larch
larch_SMC_df<- SMCs_DF_Suillus[SMCs_DF_Suillus$host == "L",]

#format dataframs to relevent info
#red
red_SMC_df2<- red_SMC_df[,5:20]
red_SMC_df3<- data.frame(round(colMeans(red_SMC_df2),0))
names<- rownames(red_SMC_df3)
red_SMC_df4<- cbind(red_SMC_df3, names, rep(("a_Red"), length(red_SMC_df3)))
colnames(red_SMC_df4)<- c("counts", "names", "group")

#white
white_SMC_df2<- white_SMC_df[,5:20]
white_SMC_df3<- data.frame(round(colMeans(white_SMC_df2),0))
names<- rownames(white_SMC_df3)
white_SMC_df4<- cbind(white_SMC_df3, names, rep(("b_White"), length(white_SMC_df3)))
colnames(white_SMC_df4)<- c("counts", "names", "group")

#larch
larch_SMC_df2<- larch_SMC_df[,5:20]
larch_SMC_df3<- data.frame(round(colMeans(larch_SMC_df2),0))
names<- rownames(larch_SMC_df3)
larch_SMC_df4<- cbind(larch_SMC_df3, names, rep(("c_Larch"), length(larch_SMC_df3)))
colnames(larch_SMC_df4)<- c("counts", "names", "group")

#bind them 
SMCs_Suillus_DF<- rbind(red_SMC_df4, white_SMC_df4, larch_SMC_df4)

#make repeating rows by SMC counts
counts_gather2 <- as.data.frame(lapply(SMCs_Suillus_DF, rep, SMCs_Suillus_DF$counts))
#wants to inherit zeros in the following- get rid of them by switching types (I dono)
counts_gather3<- as.matrix(counts_gather2)
counts_gather4<- as.data.frame(counts_gather3)

#call table
counts_gather_table<- table(counts_gather4$group, counts_gather4$names)

groups<- c("Red", "White", "Larch")
#spine plot - note - I removed hte light green color for "indole" as there are not indols in the Suillus group 
par(las = 1)
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
spineplot(counts_gather_table, main = " ", 
          col=c("#5676A1", "#405952", "#FDD191", "#885053", "#B2675E", "#442B47"), 
          border = NA, 
          #cex.axis = 0.66, 
          #ylab = "", #xlab = "", 
          yaxlabels = NA, 
          xaxlabels = groups)

#legend for spineplot
legend("topleft", legend = colnames(counts_gather_table) , 
       fill=c("#5676A1", "#405952", "#FDD191", "#885053", "#B2675E", "#442B47"),
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
N






####make grouped box plot
#first put into long form
SMCs_DF_Suillus2<- SMCs_DF_Suillus[c(SMCs_DF_Suillus$host == "W" | SMCs_DF_Suillus$host == "R" | SMCs_DF_Suillus$host == "L"),]

SMC_DF_small_no_zeros_Suillus<- SMCs_DF_Suillus2[, c(5:10, 22)]

#add repeat col (ours is only 1)
SMC_DF_small_no_zeros_Suillus$id<- seq(1:nrow(SMC_DF_small_no_zeros_Suillus))
SMC_DF_small_no_zeros_Suillus <- SMC_DF_small_no_zeros_Suillus %>% group_by(host) %>%
  gather(data = SMC_DF_small_no_zeros_Suillus, id, cf_fatty_acid, cf_putative, nrps, other, t1pks, terpene)

#check class
class(SMC_DF_small_no_zeros[1,1])

#make col names match
colnames(SMC_DF_small_no_zeros_Suillus)<- c("Group", "SMC", "Count")


#set order of boxes
SMC_DF_small_no_zeros_Suillus$Group <- factor(SMC_DF_small_no_zeros_Suillus$Group,
                       levels = c('L','W', 'R'),ordered = TRUE)

###make grouped box plot
p<- ggplot(SMC_DF_small_no_zeros_Suillus, aes(x=SMC, y=Count, fill=Group)) + 
  geom_boxplot()
p+scale_fill_manual(values=c("#B1913A", "#273666", "#FCD090")) + theme_minimal()

#horizontal
p + geom_point(position=position_jitterdodge(jitter.width = .13),alpha=0.3) +
  theme_bw(base_size = 3) + scale_fill_manual(values=c("#B1913A", "#273666", "#FCD090")) + theme_minimal() + coord_flip()


####### Stats #######
#type 2 ANOVA
library("car")
#make model
model.for.t2_Suillus_only<-lm(SMC_DF_small_no_zeros_Suillus$Count ~ SMC_DF_small_no_zeros_Suillus$SMC * SMC_DF_small_no_zeros_Suillus$Group)

model.for.t3<- Anova(model.for.t2_Suillus_only, type=2)

summary(SMC_DF_small_no_zeros_long)

#no significant difference between groups. 
