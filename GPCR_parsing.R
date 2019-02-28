########
#in MSI directory
kennedyp/llofgren/COMP/GPCRHMM/gpcrhmm

#remove stop codons
sed 's/*//g' <species>aa.fasta > <species>_wo_stops.fasta 

#form command line
#gpcrhmm.pl Suibr0_without_stops.fasta > Suibr0test.fasta
#note- due to licensing issues with the software on MSI, I ended up running it on the wed server: https://gpcrhmm.sbc.su.se 

#copy interactive output into excel, 
#then:
#1) delete the first 4 rows
#2) concatenate the header from "Sequence identifier" to "Sequence_identifier" 
#3) replace "Too short sequence" "Tooshortsequence, -,"
#4) text to columns separating on spaces and "," 
#5) save as "<speciesID>_GPCRHMM_output.csv
########

#now, in R
#load libraries
library("seqinr")

#data was generated using GPCRHMM
setwd("~/Desktop/Project_Suillus_comp_genomics/R")
#read in the input files
Suivar1_GPCRHMM<- read.csv("Suivar1_GPCRHMM_output.csv")
Suitom1_GPCRHMM<- read.csv("Suitom1_GPCRHMM_output.csv")
Suisub1_GPCRHMM<- read.csv("Suisub1_GPCRHMM_output.csv")
Suisu1_GPCRHMM<- read.csv("Suisu1_GPCRHMM_output.csv")
Suipla1_GPCRHMM<- read.csv("Suipla1_GPCRHMM_output.csv")
Suipic1_GPCRHMM<- read.csv("Suipic1_GPCRHMM_output.csv")
Suipal1_GPCRHMM<- read.csv("Suipal1_GPCRHMM_output.csv")
Suiocc1_GPCRHMM<- read.csv("Suiocc1_GPCRHMM_output.csv")
Suilu4_GPCRHMM<- read.csv("Suilu4_GPCRHMM_output.csv")
Suilak1_GPCRHMM<- read.csv("Suilak1_GPCRHMM_output.csv")
Suihi1_GPCRHMM<- read.csv("Suihi1_GPCRHMM_output.csv")
#suigr1 here
Suidec1_GPCRHMM<- read.csv("Suidec1_GPCRHMM_output.csv")
Suicot1_GPCRHMM<- read.csv("Suicot1_GPCRHMM_output.csv")
Suicli1_GPCRHMM<- read.csv("Suicli1_GPCRHMM_output.csv")
Suibr2_GPCRHMM<- read.csv("Suibr2_GPCRHMM_output.csv")
Suibov1_GPCRHMM<- read.csv("Suibov1_GPCRHMM_output.csv")
Suiamp1_GPCRHMM<- read.csv("Suiamp1_GPCRHMM_output.csv")
Suiame1_GPCRHMM<- read.csv("Suiame1_GPCRHMM_output.csv")

#shrink files to contain only positive hits
Suivar1_GPCRHMM_pos<- Suivar1_GPCRHMM[Suivar1_GPCRHMM$pred =='GPCR',]
Suitom1_GPCRHMM_pos<- Suitom1_GPCRHMM[Suitom1_GPCRHMM$pred =='GPCR',]
Suisub1_GPCRHMM_pos<- Suisub1_GPCRHMM[Suisub1_GPCRHMM$pred =='GPCR',]
Suisu1_GPCRHMM_pos<- Suisu1_GPCRHMM[Suisu1_GPCRHMM$pred =='GPCR',]
Suipla1_GPCRHMM_pos<- Suipla1_GPCRHMM[Suipla1_GPCRHMM$pred =='GPCR',]
Suipic1_GPCRHMM_pos<- Suipic1_GPCRHMM[Suipic1_GPCRHMM$pred =='GPCR',]
Suipal1_GPCRHMM_pos<- Suipal1_GPCRHMM[Suipal1_GPCRHMM$pred =='GPCR',]
Suiocc1_GPCRHMM_pos<- Suiocc1_GPCRHMM[Suiocc1_GPCRHMM$pred =='GPCR',]
Suilu4_GPCRHMM_pos<- Suilu4_GPCRHMM[Suilu4_GPCRHMM$pred =='GPCR',]
Suilak1_GPCRHMM_pos<- Suilak1_GPCRHMM[Suilak1_GPCRHMM$pred =='GPCR',]
Suihi1_GPCRHMM_pos<- Suihi1_GPCRHMM[Suihi1_GPCRHMM$pred =='GPCR',]
#guigr1 here
Suidec1_GPCRHMM_pos<- Suidec1_GPCRHMM[Suidec1_GPCRHMM$pred =='GPCR',]
Suicot1_GPCRHMM_pos<- Suicot1_GPCRHMM[Suicot1_GPCRHMM$pred =='GPCR',]
Suicli1_GPCRHMM_pos<- Suicli1_GPCRHMM[Suicli1_GPCRHMM$pred =='GPCR',]
Suibr2_GPCRHMM_pos<- Suibr2_GPCRHMM[Suibr2_GPCRHMM$pred =='GPCR',]
Suibov1_GPCRHMM_pos<- Suibov1_GPCRHMM[Suibov1_GPCRHMM$pred =='GPCR',]
Suiamp1_GPCRHMM_pos<- Suiamp1_GPCRHMM[Suiamp1_GPCRHMM$pred =='GPCR',]
Suiame1_GPCRHMM_pos<- Suiame1_GPCRHMM[Suiame1_GPCRHMM$pred =='GPCR',]



#read in the whole proteome fastas
Suivar1_in<- seqinr::read.fasta(file = "Suivar1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suitom1_in<- seqinr::read.fasta(file = "Suitom1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisub1_in<- seqinr::read.fasta(file = "Suisub1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisu1_in<- seqinr::read.fasta(file = "Suisu1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipla1_in<- seqinr::read.fasta(file = "Suipla1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipic1_in<- seqinr::read.fasta(file = "Suipic1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipal1_in<- seqinr::read.fasta(file = "Suipal1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiocc1_in<- seqinr::read.fasta(file = "Suiocc1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilu4_in<- seqinr::read.fasta(file = "Suilu4_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilak1_in<- seqinr::read.fasta(file = "Suilak1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suihi1_in<- seqinr::read.fasta(file = "Suihi1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
#Suigr1_in<- seqinr::read.fasta(file = "Suigr1_wo_stops.fasta", 
 #                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suidec1_in<- seqinr::read.fasta(file = "Suidec1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicot1_in<- seqinr::read.fasta(file = "Suicot1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicli1_in<- seqinr::read.fasta(file = "Suicli1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibr2_in<- seqinr::read.fasta(file = "Suibr2_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibov1_in<- seqinr::read.fasta(file = "Suibov1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiamp1_in<- seqinr::read.fasta(file = "Suiamp1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiame1_in<- seqinr::read.fasta(file = "Suiame1_wo_stops.fasta", 
                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)



#get fastas of only positive hits for GPCR analysis 
Suivar1_pos_fastas<- Suivar1_in[c(which(names(Suivar1_in) %in% Suivar1_GPCRHMM_pos$Sequence_identifier))]
Suitom1_pos_fastas<- Suitom1_in[c(which(names(Suitom1_in) %in% Suitom1_GPCRHMM_pos$Sequence_identifier))]
Suisub1_pos_fastas<- Suisub1_in[c(which(names(Suisub1_in) %in% Suisub1_GPCRHMM_pos$Sequence_identifier))]
Suisu1_pos_fastas<- Suisu1_in[c(which(names(Suisu1_in) %in% Suisu1_GPCRHMM_pos$Sequence_identifier))]
Suipla1_pos_fastas<- Suipla1_in[c(which(names(Suipla1_in) %in% Suipla1_GPCRHMM_pos$Sequence_identifier))]
Suipic1_pos_fastas<- Suipic1_in[c(which(names(Suipic1_in) %in% Suipic1_GPCRHMM_pos$Sequence_identifier))]
Suipal1_pos_fastas<- Suipal1_in[c(which(names(Suipal1_in) %in% Suipal1_GPCRHMM_pos$Sequence_identifier))]
Suiocc1_pos_fastas<- Suiocc1_in[c(which(names(Suiocc1_in) %in% Suiocc1_GPCRHMM_pos$Sequence_identifier))]
Suilu4_pos_fastas<- Suilu4_in[c(which(names(Suilu4_in) %in% Suilu4_GPCRHMM_pos$Sequence_identifier))]
Suilak1_pos_fastas<- Suilak1_in[c(which(names(Suilak1_in) %in% Suilak1_GPCRHMM_pos$Sequence_identifier))]
Suihi1_pos_fastas<- Suihi1_in[c(which(names(Suihi1_in) %in% Suihi1_GPCRHMM_pos$Sequence_identifier))]
#Suigr1_pos_fastas<- Suigr1_in[c(which(names(Suigr1_in) %in% Suigr1_GPCRHMM_pos$Sequence_identifier))]
Suidec1_pos_fastas<- Suidec1_in[c(which(names(Suidec1_in) %in% Suidec1_GPCRHMM_pos$Sequence_identifier))]
Suicot1_pos_fastas<- Suicot1_in[c(which(names(Suicot1_in) %in% Suicot1_GPCRHMM_pos$Sequence_identifier))]
Suicli1_pos_fastas<- Suicli1_in[c(which(names(Suicli1_in) %in% Suicli1_GPCRHMM_pos$Sequence_identifier))]
Suibr2_pos_fastas<- Suibr2_in[c(which(names(Suibr2_in) %in% Suibr2_GPCRHMM_pos$Sequence_identifier))]
Suibov1_pos_fastas<- Suibov1_in[c(which(names(Suibov1_in) %in% Suibov1_GPCRHMM_pos$Sequence_identifier))]
Suiamp1_pos_fastas<- Suiamp1_in[c(which(names(Suiamp1_in) %in% Suiamp1_GPCRHMM_pos$Sequence_identifier))]
Suiame1_pos_fastas<- Suiame1_in[c(which(names(Suiame1_in) %in% Suiame1_GPCRHMM_pos$Sequence_identifier))]




#print the fastas to files
write.fasta(Suivar1_pos_fastas, names = names(Suivar1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suivar1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suitom1_pos_fastas, names = names(Suitom1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suitom1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suisub1_pos_fastas, names = names(Suisub1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisub1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suisu1_pos_fastas, names = names(Suisu1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisu1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suipla1_pos_fastas, names = names(Suipla1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipla1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suipic1_pos_fastas, names = names(Suipic1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipic1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suipal1_pos_fastas, names = names(Suipal1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipal1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suiocc1_pos_fastas, names = names(Suiocc1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiocc1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suilu4_pos_fastas, names = names(Suilu4_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilu4_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suilak1_pos_fastas, names = names(Suilak1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilak1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suihi1_pos_fastas, names = names(Suihi1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suihi1_GPCRs_from_GPCRHMM.fasta")
#write.fasta(Suigr1_pos_fastas, names = names(Suigr1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suigr1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suidec1_pos_fastas, names = names(Suidec1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suidec1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suicot1_pos_fastas, names = names(Suicot1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicot1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suicli1_pos_fastas, names = names(Suicli1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicli1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suibr2_pos_fastas, names = names(Suibr2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibr2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suibov1_pos_fastas, names = names(Suibov1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibov1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suiamp1_pos_fastas, names = names(Suiamp1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiamp1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suiame1_pos_fastas, names = names(Suiame1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiame1_GPCRs_from_GPCRHMM.fasta")


#basic numbers on how many GPCR's precited per genome with GPCRHMM
Suivar1<- nrow(Suivar1_GPCRHMM_pos)
Suitom1<- nrow(Suitom1_GPCRHMM_pos)
Suisub1<- nrow(Suisub1_GPCRHMM_pos)
Suisu1<- nrow(Suisu1_GPCRHMM_pos)
Suipla1<- nrow(Suipla1_GPCRHMM_pos)
Suipic1<- nrow(Suipic1_GPCRHMM_pos)
Suipal1<- nrow(Suipal1_GPCRHMM_pos)
Suiocc1<- nrow(Suiocc1_GPCRHMM_pos)
Suilu4<- nrow(Suilu4_GPCRHMM_pos)
Suilak1<- nrow(Suilak1_GPCRHMM_pos)
Suihi1<- nrow(Suihi1_GPCRHMM_pos)
#Suigr1<- nrow(Suigr1_GPCRHMM_pos)
Suidec1<- nrow(Suidec1_GPCRHMM_pos)
Suicot1<- nrow(Suicot1_GPCRHMM_pos)
Suicli1<- nrow(Suicli1_GPCRHMM_pos)
Suibr2<- nrow(Suibr2_GPCRHMM_pos)
Suibov1<- nrow(Suibov1_GPCRHMM_pos)
Suiamp1<- nrow(Suiamp1_GPCRHMM_pos)
Suiame1<- nrow(Suiame1_GPCRHMM_pos)

totals<- data.frame(cbind(Suivar1, 
                          Suitom1, 
                          Suisub1, 
                          Suisu1, 
                          Suipla1, 
                          Suipic1, 
                          Suipal1, 
                          Suiocc1, 
                          Suilu4, 
                          Suilak1, 
                          Suihi1, 
                          Suidec1, 
                          Suicot1, 
                          Suicli1, 
                          Suibr2, 
                          Suibov1, 
                          Suiamp1, 
                          Suiame1),
                     row.names = "#GPCR's_w_GPCRHMM")



 

#######
#run TMHMM analysis in bash on MSI then re-import the output files
# On MSI:
# module load tmhmm
# 
# #copy these the first time
# #cp /panfs/roc/msisoft/tmhmm/2.0c/lib/TMHMM2.0.model /home/kennedyp/llofgren/COMP/TMHMM
# #cp /panfs/roc/msisoft/tmhmm/2.0c/lib/TMHMM2.0.options /home/kennedyp/llofgren/COMP/TMHMM
# 
# #Then you can create a variable holding the path to the program:
# DECODE="/panfs/roc/msisoft/tmhmm/2.0c/bin/decodeanhmm"
# FORMAT="/panfs/roc/msisoft/tmhmm/2.0c/bin/tmhmmformat.pl"
# 
# #run tmhmm in the format: 
# #cat my.fasta | ${DECODE} -f /my/destination_path/TMHMM2.0.options
# #-modelfile /my/destination_path/TMHMM2.0.model -N1 -PrintNumbers
# #-PrintScore -PrintStat | perl ${FORMAT} > formatted_output.txt
# 
# #like so:
# #cat Suivar1_GPCRs_from_GPCRHMM.fasta | ${DECODE} -f /home/kennedyp/llofgren/COMP/TMHMM/TMHMM2.0.options -modelfile /home/kennedyp/llofgren/COMP/TMHMM/TMHMM2.0.model -N1 -PrintNumbers -PrintScore -PrintStat | perl ${FORMAT}> Suivar1_TMHMM_output.txt
#######

#get only the subset of seqs that have 7 TMH's

#read in files 
Suivar1_TMHMM<- data.frame(read.csv("Suivar1_TMHMM_output.txt", col.names = "header"))
Suitom1_TMHMM<- data.frame(read.csv("Suitom1_TMHMM_output.txt", col.names = "header"))
Suisub1_TMHMM<- data.frame(read.csv("Suisub1_TMHMM_output.txt", col.names = "header"))
Suisu1_TMHMM<- data.frame(read.csv("Suisu1_TMHMM_output.txt", col.names = "header"))
Suipla1_TMHMM<- data.frame(read.csv("Suipla1_TMHMM_output.txt", col.names = "header"))
Suipic1_TMHMM<- data.frame(read.csv("Suipic1_TMHMM_output.txt", col.names = "header"))
Suipal1_TMHMM<- data.frame(read.csv("Suipal1_TMHMM_output.txt", col.names = "header"))
Suiocc1_TMHMM<- data.frame(read.csv("Suiocc1_TMHMM_output.txt", col.names = "header"))
Suilu4_TMHMM<- data.frame(read.csv("Suilu4_TMHMM_output.txt", col.names = "header"))
Suilak1_TMHMM<- data.frame(read.csv("Suilak1_TMHMM_output.txt", col.names = "header"))
Suihi1_TMHMM<- data.frame(read.csv("Suihi1_TMHMM_output.txt", col.names = "header"))
#Suigr1_TMHMM<- data.frame(read.csv("Suigr1_TMHMM_output.txt", col.names = "header"))
Suidec1_TMHMM<- data.frame(read.csv("Suidec1_TMHMM_output.txt", col.names = "header"))
Suicot1_TMHMM<- data.frame(read.csv("Suicot1_TMHMM_output.txt", col.names = "header"))
Suicli1_TMHMM<- data.frame(read.csv("Suicli1_TMHMM_output.txt", col.names = "header"))
Suibr2_TMHMM<- data.frame(read.csv("Suibr2_TMHMM_output.txt", col.names = "header"))
Suibov1_TMHMM<- data.frame(read.csv("Suibov1_TMHMM_output.txt", col.names = "header"))
Suiamp1_TMHMM<- data.frame(read.csv("Suiamp1_TMHMM_output.txt", col.names = "header"))
Suiame1_TMHMM<- data.frame(read.csv("Suiame1_TMHMM_output.txt", col.names = "header"))


#isolate hashed lines with relevent info
Suivar1_TMHMM_sm<- data.frame(Suivar1_TMHMM[ grep("Number of predicted TMHs:  7", Suivar1_TMHMM$header),])
Suitom1_TMHMM_sm<- data.frame(Suitom1_TMHMM[ grep("Number of predicted TMHs:  7", Suitom1_TMHMM$header),])
Suisub1_TMHMM_sm<- data.frame(Suisub1_TMHMM[ grep("Number of predicted TMHs:  7", Suisub1_TMHMM$header),])
Suisu1_TMHMM_sm<- data.frame(Suisu1_TMHMM[ grep("Number of predicted TMHs:  7", Suisu1_TMHMM$header),])
Suipla1_TMHMM_sm<- data.frame(Suipla1_TMHMM[ grep("Number of predicted TMHs:  7", Suipla1_TMHMM$header),])
Suipic1_TMHMM_sm<- data.frame(Suipic1_TMHMM[ grep("Number of predicted TMHs:  7", Suipic1_TMHMM$header),])
Suipal1_TMHMM_sm<- data.frame(Suipal1_TMHMM[ grep("Number of predicted TMHs:  7", Suipal1_TMHMM$header),])
Suiocc1_TMHMM_sm<- data.frame(Suiocc1_TMHMM[ grep("Number of predicted TMHs:  7", Suiocc1_TMHMM$header),])
Suilu4_TMHMM_sm<- data.frame(Suilu4_TMHMM[ grep("Number of predicted TMHs:  7", Suilu4_TMHMM$header),])
Suilak1_TMHMM_sm<- data.frame(Suilak1_TMHMM[ grep("Number of predicted TMHs:  7", Suilak1_TMHMM$header),])
Suihi1_TMHMM_sm<- data.frame(Suihi1_TMHMM[ grep("Number of predicted TMHs:  7", Suihi1_TMHMM$header),])
#Suigr1_TMHMM_sm<- data.frame(Suigr1_TMHMM[ grep("Number of predicted TMHs:  7", Suigr1_TMHMM$header),])
Suidec1_TMHMM_sm<- data.frame(Suidec1_TMHMM[ grep("Number of predicted TMHs:  7", Suidec1_TMHMM$header),])
Suicot1_TMHMM_sm<- data.frame(Suicot1_TMHMM[ grep("Number of predicted TMHs:  7", Suicot1_TMHMM$header),])
Suicli1_TMHMM_sm<- data.frame(Suicli1_TMHMM[ grep("Number of predicted TMHs:  7", Suicli1_TMHMM$header),])
Suibr2_TMHMM_sm<- data.frame(Suibr2_TMHMM[ grep("Number of predicted TMHs:  7", Suibr2_TMHMM$header),])
Suibov1_TMHMM_sm<- data.frame(Suibov1_TMHMM[ grep("Number of predicted TMHs:  7", Suibov1_TMHMM$header),])
Suiamp1_TMHMM_sm<- data.frame(Suiamp1_TMHMM[ grep("Number of predicted TMHs:  7", Suiamp1_TMHMM$header),])
Suiame1_TMHMM_sm<- data.frame(Suiame1_TMHMM[ grep("Number of predicted TMHs:  7", Suiame1_TMHMM$header),])


#get summary numbers for each of these and attach them to the previous output file. 
Suivar1<- nrow(Suivar1_TMHMM_sm)
Suitom1<- nrow(Suitom1_TMHMM_sm)
Suisub1<- nrow(Suisub1_TMHMM_sm)
Suisu1<- nrow(Suisu1_TMHMM_sm)
Suipla1<- nrow(Suipla1_TMHMM_sm)
Suipic1<- nrow(Suipic1_TMHMM_sm)
Suipal1<- nrow(Suipal1_TMHMM_sm)
Suiocc1<- nrow(Suiocc1_TMHMM_sm)
Suilu4<- nrow(Suilu4_TMHMM_sm)
Suilak1<- nrow(Suilak1_TMHMM_sm)
Suihi1<- nrow(Suihi1_TMHMM_sm)
#Suigr1<- nrow(Suigr1_TMHMM_sm)
Suidec1<- nrow(Suidec1_TMHMM_sm)
Suicot1<- nrow(Suicot1_TMHMM_sm)
Suicli1<- nrow(Suicli1_TMHMM_sm)
Suibr2<- nrow(Suibr2_TMHMM_sm)
Suibov1<- nrow(Suibov1_TMHMM_sm)
Suiamp1<- nrow(Suiamp1_TMHMM_sm)
Suiame1<- nrow(Suiame1_TMHMM_sm)


totals.2<- data.frame(cbind(Suivar1, 
                            Suitom1, 
                            Suisub1, 
                            Suisu1, 
                            Suipla1, 
                            Suipic1, 
                            Suipal1, 
                            Suiocc1, 
                            Suilu4, 
                            Suilak1, 
                            Suihi1, 
                            Suidec1, 
                            Suicot1, 
                            Suicli1, 
                            Suibr2, 
                            Suibov1, 
                            Suiamp1, 
                            Suiame1),
                      row.names = "#GPCR's_w_GPCRHMM_andTMHMM")
  
#write out totals to a file
results<- rbind(totals, totals.2[1,])

totals.2[1,]
write.csv(results, quote = FALSE, file = "GPCR_totals.csv")
