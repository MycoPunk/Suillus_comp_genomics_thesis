#run TOPCONS on GPCR lists, in this example we're running it on the crystal validated GPCR's. 

#http://single.topcons.net///index.php
#import "all topologies" outputflle, save as a .csv file
setwd("~/Desktop/")

#kill those factors
options(stringsAsFactors = FALSE)
#import files 
TOPCONS_validated_GPCRs<- read.table("Crystal_validated_GPCR_topologies.txt", fill = TRUE, header = FALSE)
#format TOPCONS files
validated_GPCRs.names<- TOPCONS_validated_GPCRs[grep(">", TOPCONS_validated_GPCRs$V1),1]
validated_GPCRs.seq<- TOPCONS_validated_GPCRs[ - grep(">", TOPCONS_validated_GPCRs$V1),1]
TOPCONS_validated_GPCRs<-cbind(validated_GPCRs.names, validated_GPCRs.seq)
TOPCONS_validated_GPCRs.2<- data.frame(gsub(">", "", TOPCONS_validated_GPCRs))
colnames(TOPCONS_validated_GPCRs.2)<- c("name", "seq")

#names
View(validated_GPCRs.names)
#positions
View(TOPCONS_validated_GPCRs.2)



#import validated GPCR fasta (format like Thega1_pos_fastas_from_Phobius)
#View(Thega1_pos_fastas_from_Phobius)
validated_GPCRs_pos_fastas_from_Phobius<- read.fasta("Crystal_validated_GPCRs.txt", as.string = TRUE, seqtype = "AA")

#format the fasta files to match against
validated_GPCRs_pos_fastas_from_Phobius_M<- as.matrix(validated_GPCRs_pos_fastas_from_Phobius)
validated_GPCRs_pos_fastas_from_Phobius_df<- as.data.frame(validated_GPCRs_pos_fastas_from_Phobius_M)

#rename fastas 
seq.df<- data.frame(validated_GPCRs_pos_fastas_from_Phobius_df)
map.key<- data.frame(TOPCONS_validated_GPCRs.2[,2])

#get indices of "M" regions 
map.list.list<- apply(map.key, 2, function (x) str_locate_all(pattern ='M', x))
#get starts only
map.list.list<- map.list.list[[1]]
map.df<- mapply('[', map.list.list, TRUE, 1)

#itterate over all 
validated_GPCRs_cat_TM_list<- unlist(Map(function(x, y) paste(substring(x, unlist(y), 
                                                                        unlist(y)), collapse=""), seq.df[[1]], map.df))
##re-attach sequence names and format back into fasta files. 
#get list for fasta headers
names<- TOPCONS_validated_GPCRs.2$name

#rename rows
validated_GPCRs_cat_TM_list<- data.frame(validated_GPCRs_cat_TM_list)
rownames(validated_GPCRs_cat_TM_list)<- TOPCONS_validated_GPCRs.2$name

#add the carrot back into the row names
rownames(validated_GPCRs_cat_TM_list) <- paste0(">",  rownames(validated_GPCRs_cat_TM_list))

#format and print this to a file 
write.table(validated_GPCRs_cat_TM_list, file = "validated_GPCRs_cat_TM_list.fasta", col.names=FALSE, quote = FALSE, sep = "\n")

#quick QC check - make sure the pro length is the same as the respective M lists. 
str_count(map.key$TOPCONS_validated_GPCRs.2...2., "M")
nchar(validated_GPCRs_cat_TM_list[1,1])
nchar(validated_GPCRs_cat_TM_list[2,1])
#looks good. 

#align the outputfiles using clustal OMEGA, and build using RAxML-HPC2 w/ boot strapped over 10000 generations. 




