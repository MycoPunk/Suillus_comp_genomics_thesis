


#load libraries
library("seqinr")
library("data.table")
library("stringr")
library("car")

options(stringsAsFactors = FALSE)
setwd("~/Desktop/Project_Suillus_comp_genomics/R")

#read it back in 
SSPs_coded_within_Suillus<- read.csv("SSPs_coded_within_Suillus.csv")

#split the two data categories for the three groups and get stats for each 
Red_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="Red",]
White_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="White",]
Larch_df<- SSPs_coded_within_Suillus[SSPs_coded_within_Suillus$group=="Larch",]

#means
Red_genome_size_mean<- mean(Red_df$genome_size)
White_genome_size_mean<- mean(White_df$genome_size)
Larch_genome_size_mean<- mean(Larch_df$genome_size)

#SD and SE
std <- function(x) sd(x)/sqrt(length(x))

Red_genome_size_SD<- sd(Red_df$genome_size)
Red_genome_size_SE<- std(Red_df$genome_size)  

White_genome_size_SD<- sd(White_df$genome_size)
White_genome_size_SE<- std(White_df$genome_size)  

Larch_genome_size_SD<- sd(Larch_df$genome_size)
Larch_genome_size_SE<- std(Larch_df$genome_size)  

all_df<- rbind(Red_df, White_df, Larch_df)
#test normality 
shapiro.test(all_df$genome_size)
#not normal, need to transform

#build model
genome_size_model <- lm(all_df$genome_size ~ all_df$group)

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


####proteins
#means
Red_prot_mean<- mean(Red_df$n_proteins_from_gene_cat)
White_prot_mean<- mean(White_df$n_proteins_from_gene_cat)
Larch_prot_mean<- mean(Larch_df$n_proteins_from_gene_cat)

#SD and SE
Red_prot_SD<- sd(Red_df$n_proteins_from_gene_cat)
Red_prot_SE<- std(Red_df$n_proteins_from_gene_cat)  

White_prot_SD<- sd(White_df$n_proteins_from_gene_cat)
White_prot_SE<- std(White_df$n_proteins_from_gene_cat)  

Larch_prot_SD<- sd(Larch_df$n_proteins_from_gene_cat)
Larch_prot_SE<- std(Larch_df$n_proteins_from_gene_cat)  

#test normality 
shapiro.test(all_df$n_proteins_from_gene_cat)
#normal

#build model
prot_model <- lm(all_df$n_proteins_from_gene_cat ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(prot_model) #what does boxcox suggest? log transformation works. 
powerTransform(prot_model) #how about powerTransform? #log should be find for this.

#log transform to improve normality?
all_df$n_proteins_from_gene_cat_log<- log(all_df$n_proteins_from_gene_cat)

#see if the transformation helped
prot_model2 <- lm(all_df$n_proteins_from_gene_cat_log ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#well, that didn't help...

#try the boc recomended transformation
prot_model3 <- lm(all_df$n_proteins_from_gene_cat^-.5 ~ all_df$group)
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(prot_model3) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#worse! go back. do not pass go. do not colelct $100. Use original data. 

#run anova
summary(aov(prot_model))
#not significant


####SSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)
White_ssp_mean<- mean(White_df$n_SSPs_signalP_TMHMM_lt_300aa)
Larch_ssp_mean<- mean(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)

#SD and SE
Red_ssp_SD<- sd(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)
Red_ssp_SE<- std(Red_df$n_SSPs_signalP_TMHMM_lt_300aa)  

White_ssp_SD<- sd(White_df$n_SSPs_signalP_TMHMM_lt_300aa)
White_ssp_SE<- std(White_df$n_SSPs_signalP_TMHMM_lt_300aa)  

Larch_ssp_SD<- sd(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)
Larch_ssp_SE<- std(Larch_df$n_SSPs_signalP_TMHMM_lt_300aa)  

#build model
ssp_model <- lm(all_df$n_SSPs_signalP_TMHMM_lt_300aa ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(ssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(ssp_model) #how about powerTransform? # ^-2 is what's recommended 

all_df$n_SSPs_signalP_TMHMM_lt_300aa_inv_square<- all_df$n_SSPs_signalP_TMHMM_lt_300aa^(-2)

ssp_model2 <- lm(all_df$n_SSPs_signalP_TMHMM_lt_300aa_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(ssp_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(ssp_model2))
#not significant



####SSSPs
#means
Red_ssp_mean<- mean(Red_df$n_SSSPs)
White_ssp_mean<- mean(White_df$n_SSSPs)
Larch_ssp_mean<- mean(Larch_df$n_SSSPs)

#SD and SE
Red_sssp_SD<- sd(Red_df$n_SSSPs)
Red_sssp_SE<- std(Red_df$n_SSSPs)  

White_sssp_SD<- sd(White_df$n_SSSPs)
White_sssp_SE<- std(White_df$n_SSSPs)  

Larch_sssp_SD<- sd(Larch_df$n_SSSPs)
Larch_sssp_SE<- std(Larch_df$n_SSSPs)  

#build model
sssp_model <- lm(all_df$n_SSSPs ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sssp_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(sssp_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(sssp_model) #how about powerTransform? # ^-1 is what's recommended 

#transform
all_df$n_SSSPs_inv_one<- all_df$n_SSSPs^(-1)

#make new model
sssp_model2 <- lm(all_df$n_SSSPs_inv_one ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(sssp_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(sssp_model2))
#significant!

#what's different?
sssp_aov<- aov(sssp_model2)
TukeyHSD(sssp_aov)
#red and larch are sig. different 


####Effectors
#means
Red_effector_mean<- mean(Red_df$n_effectors_from_EffectorP)
White_effector_mean<- mean(White_df$n_effectors_from_EffectorP)
Larch_effector_mean<- mean(Larch_df$n_effectors_from_EffectorP)

#SD and SE
Red_effector_SD<- sd(Red_df$n_effectors_from_EffectorP)
Red_effector_SE<- std(Red_df$n_effectors_from_EffectorP)  

White_effector_SD<- sd(White_df$n_effectors_from_EffectorP)
White_effector_SE<- std(White_df$n_effectors_from_EffectorP)  

Larch_effector_SD<- sd(Larch_df$n_effectors_from_EffectorP)
Larch_effector_SE<- std(Larch_df$n_effectors_from_EffectorP)  

#build model
effector_model <- lm(all_df$n_effectors_from_EffectorP ~ all_df$group)

#vis normality and variance
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
boxCox(effector_model, lambda = seq(-6,2,2)) #what does boxcox suggest?  
powerTransform(effector_model) #how about powerTransform? # ^-2 is what's recommended 

#transform
all_df$n_effectors_from_EffectorP_inv_square<- all_df$n_effectors_from_EffectorP^(-2)

#make new model
effector_model2 <- lm(all_df$n_effectors_from_EffectorP_inv_square ~ all_df$group)

#see if the transformation helped
par(mfrow = c(2,2)) #prints up to 4 plots on a 2 x 2 page
plot(effector_model2) #data is non normal with non constant variance
par(mfrow = c(1,1)) #back to one image per page
#looks better

#run anova
summary(aov(effector_model2))
#not significant



