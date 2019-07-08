########
#in MSI directory
kennedyp/llofgren/COMP/GPCRHMM/gpcrhmm

#remove stop codons
#sed 's/*//g' <species>aa.fasta > <species>_wo_stops.fasta 

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
library("plyr")
library("waffle")
library("vioplot")
library("tidyr")
library("dplyr")


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
Suigr1_GPCRHMM<- read.csv("Suigr1_GPCRHMM_output.csv")
Suidec1_GPCRHMM<- read.csv("Suidec1_GPCRHMM_output.csv")
Suicot1_GPCRHMM<- read.csv("Suicot1_GPCRHMM_output.csv")
Suicli1_GPCRHMM<- read.csv("Suicli1_GPCRHMM_output.csv")
Suibr2_GPCRHMM<- read.csv("Suibr2_GPCRHMM_output.csv")
Suibov1_GPCRHMM<- read.csv("Suibov1_GPCRHMM_output.csv")
Suiamp1_GPCRHMM<- read.csv("Suiamp1_GPCRHMM_output.csv")
Suiame1_GPCRHMM<- read.csv("Suiame1_GPCRHMM_output.csv")
Suifus1_GPCRHMM<- read.csv("Suifus1_GPCRHMM_output.csv")


#non_Suillus set
Rhivul1_GPCRHMM<- read.csv("Rhivul1_GPCRHMM_output.csv")
Rhitru1_GPCRHMM<- read.csv("Rhitru1_GPCRHMM_output.csv")
Amamu1_GPCRHMM<- read.csv("Amamu1_GPCRHMM_output.csv")
Hebcy2_GPCRHMM<- read.csv("Hebcy2_GPCRHMM_output.csv")
Lacbi2_GPCRHMM<- read.csv("Lacbi2_GPCRHMM_output.csv")
Paxin1_GPCRHMM<- read.csv("Paxin1_GPCRHMM_output.csv")
Pilcr1_GPCRHMM<- read.csv("Pilcr1_GPCRHMM_output.csv")
Pismi1_GPCRHMM<- read.csv("Pismi1_GPCRHMM_output.csv")
Sclci1_GPCRHMM<- read.csv("Sclci1_GPCRHMM_output.csv")
Cananz1_GPCRHMM<- read.csv("Cananz1_GPCRHMM_output.csv")
Gyrli1_GPCRHMM<- read.csv("Gyrli1_GPCRHMM_output.csv")
Gaumor1_GPCRHMM<- read.csv("Gaumor1_GPCRHMM_output.csv")
Hydru2_GPCRHMM<- read.csv("Hydru2_GPCRHMM_output.csv")
Hyssto1_GPCRHMM<- read.csv("Hyssto1_GPCRHMM_output.csv")
Lacam2_GPCRHMM<- read.csv("Lacam2_GPCRHMM_output.csv")
Pisti1_GPCRHMM<- read.csv("Pisti1_GPCRHMM_output.csv")
Rusbre1_GPCRHMM<- read.csv("Rusbre1_GPCRHMM_output.csv")
Ruscom1_GPCRHMM<- read.csv("Ruscom1_GPCRHMM_output.csv")
Rhisa1_GPCRHMM<- read.csv("Rhisa1_GPCRHMM_output.csv")
Rhives1_GPCRHMM<- read.csv("Rhives1_GPCRHMM_output.csv")
Rhivi1_GPCRHMM<- read.csv("Rhivi1_GPCRHMM_output.csv")
Theter1_GPCRHMM<- read.csv("Theter1_GPCRHMM_output.csv")
Thega1_GPCRHMM<- read.csv("Thega1_GPCRHMM_output.csv")


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
Suigr1_GPCRHMM_pos<- Suigr1_GPCRHMM[Suigr1_GPCRHMM$pred =='GPCR',]
Suidec1_GPCRHMM_pos<- Suidec1_GPCRHMM[Suidec1_GPCRHMM$pred =='GPCR',]
Suicot1_GPCRHMM_pos<- Suicot1_GPCRHMM[Suicot1_GPCRHMM$pred =='GPCR',]
Suicli1_GPCRHMM_pos<- Suicli1_GPCRHMM[Suicli1_GPCRHMM$pred =='GPCR',]
Suibr2_GPCRHMM_pos<- Suibr2_GPCRHMM[Suibr2_GPCRHMM$pred =='GPCR',]
Suibov1_GPCRHMM_pos<- Suibov1_GPCRHMM[Suibov1_GPCRHMM$pred =='GPCR',]
Suiamp1_GPCRHMM_pos<- Suiamp1_GPCRHMM[Suiamp1_GPCRHMM$pred =='GPCR',]
Suiame1_GPCRHMM_pos<- Suiame1_GPCRHMM[Suiame1_GPCRHMM$pred =='GPCR',]
Suifus1_GPCRHMM_pos<- Suifus1_GPCRHMM[Suifus1_GPCRHMM$pred =='GPCR',]
#non-Suillus set
Rhivul1_GPCRHMM_pos<- Rhivul1_GPCRHMM[Rhivul1_GPCRHMM$pred =='GPCR',]
Rhitru1_GPCRHMM_pos<- Rhitru1_GPCRHMM[Rhitru1_GPCRHMM$pred =='GPCR',]
Amamu1_GPCRHMM_pos<- Amamu1_GPCRHMM[Amamu1_GPCRHMM$pred =='GPCR',]
Hebcy2_GPCRHMM_pos<- Hebcy2_GPCRHMM[Hebcy2_GPCRHMM$pred =='GPCR',]
Lacbi2_GPCRHMM_pos<- Lacbi2_GPCRHMM[Lacbi2_GPCRHMM$pred =='GPCR',]
Paxin1_GPCRHMM_pos<- Paxin1_GPCRHMM[Paxin1_GPCRHMM$pred =='GPCR',]
Pilcr1_GPCRHMM_pos<- Pilcr1_GPCRHMM[Pilcr1_GPCRHMM$pred =='GPCR',]
Pismi1_GPCRHMM_pos<- Pismi1_GPCRHMM[Pismi1_GPCRHMM$pred =='GPCR',]
Sclci1_GPCRHMM_pos<- Sclci1_GPCRHMM[Sclci1_GPCRHMM$pred =='GPCR',]
Cananz1_GPCRHMM_pos<- Cananz1_GPCRHMM[Cananz1_GPCRHMM$pred =='GPCR',]
Gyrli1_GPCRHMM_pos<- Gyrli1_GPCRHMM[Gyrli1_GPCRHMM$pred =='GPCR',]
Gaumor1_GPCRHMM_pos<- Gaumor1_GPCRHMM[Gaumor1_GPCRHMM$pred =='GPCR',]
Hydru2_GPCRHMM_pos<- Hydru2_GPCRHMM[Hydru2_GPCRHMM$pred =='GPCR',]
Hyssto1_GPCRHMM_pos<- Hyssto1_GPCRHMM[Hyssto1_GPCRHMM$pred =='GPCR',]
Lacam2_GPCRHMM_pos<- Lacam2_GPCRHMM[Lacam2_GPCRHMM$pred =='GPCR',]
Pisti1_GPCRHMM_pos<- Pisti1_GPCRHMM[Pisti1_GPCRHMM$pred =='GPCR',]
Rusbre1_GPCRHMM_pos<- Rusbre1_GPCRHMM[Rusbre1_GPCRHMM$pred =='GPCR',]
Ruscom1_GPCRHMM_pos<- Ruscom1_GPCRHMM[Ruscom1_GPCRHMM$pred =='GPCR',]
Rhisa1_GPCRHMM_pos<- Rhisa1_GPCRHMM[Rhisa1_GPCRHMM$pred =='GPCR',]
Rhives1_GPCRHMM_pos<- Rhives1_GPCRHMM[Rhives1_GPCRHMM$pred =='GPCR',]
Rhivi1_GPCRHMM_pos<- Rhivi1_GPCRHMM[Rhivi1_GPCRHMM$pred =='GPCR',]
Theter1_GPCRHMM_pos<- Theter1_GPCRHMM[Theter1_GPCRHMM$pred =='GPCR',]
Thega1_GPCRHMM_pos<- Thega1_GPCRHMM[Thega1_GPCRHMM$pred =='GPCR',]




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
Suigr1_in<- seqinr::read.fasta(file = "Suigr1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
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
Suifus1_in<- seqinr::read.fasta(file = "Suifus1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
#non-Suillus set
Rhivul1_in<- seqinr::read.fasta(file = "Rhivul1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhitru1_in<- seqinr::read.fasta(file = "Rhitru1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Amamu1_in<- seqinr::read.fasta(file = "Amamu1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hebcy2_in<- seqinr::read.fasta(file = "Hebcy2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacbi2_in<- seqinr::read.fasta(file = "Lacbi2_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Paxin1_in<- seqinr::read.fasta(file = "Paxin1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pilcr1_in<- seqinr::read.fasta(file = "Pilcr1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pismi1_in<- seqinr::read.fasta(file = "Pismi1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Sclci1_in<- seqinr::read.fasta(file = "Sclci1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Cananz1_in<- seqinr::read.fasta(file = "Cananz1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gyrli1_in<- seqinr::read.fasta(file = "Gyrli1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Gaumor1_in<- seqinr::read.fasta(file = "Gaumor1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hydru2_in<- seqinr::read.fasta(file = "Hydru2_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Hyssto1_in<- seqinr::read.fasta(file = "Hyssto1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Lacam2_in<- seqinr::read.fasta(file = "Lacam2_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Pisti1_in<- seqinr::read.fasta(file = "Pisti1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rusbre1_in<- seqinr::read.fasta(file = "Rusbre1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Ruscom1_in<- seqinr::read.fasta(file = "Ruscom1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhisa1_in<- seqinr::read.fasta(file = "Rhisa1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhives1_in<- seqinr::read.fasta(file = "Rhives1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Rhivi1_in<- seqinr::read.fasta(file = "Rhivi1_wo_stops.fasta", 
                                seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Theter1_in<- seqinr::read.fasta(file = "Theter1_wo_stops.fasta", 
                               seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Thega1_in<- seqinr::read.fasta(file = "Thega1_wo_stops.fasta", 
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
Suigr1_pos_fastas<- Suigr1_in[c(which(names(Suigr1_in) %in% Suigr1_GPCRHMM_pos$Sequence_identifier))]
Suidec1_pos_fastas<- Suidec1_in[c(which(names(Suidec1_in) %in% Suidec1_GPCRHMM_pos$Sequence_identifier))]
Suicot1_pos_fastas<- Suicot1_in[c(which(names(Suicot1_in) %in% Suicot1_GPCRHMM_pos$Sequence_identifier))]
Suicli1_pos_fastas<- Suicli1_in[c(which(names(Suicli1_in) %in% Suicli1_GPCRHMM_pos$Sequence_identifier))]
Suibr2_pos_fastas<- Suibr2_in[c(which(names(Suibr2_in) %in% Suibr2_GPCRHMM_pos$Sequence_identifier))]
Suibov1_pos_fastas<- Suibov1_in[c(which(names(Suibov1_in) %in% Suibov1_GPCRHMM_pos$Sequence_identifier))]
Suiamp1_pos_fastas<- Suiamp1_in[c(which(names(Suiamp1_in) %in% Suiamp1_GPCRHMM_pos$Sequence_identifier))]
Suiame1_pos_fastas<- Suiame1_in[c(which(names(Suiame1_in) %in% Suiame1_GPCRHMM_pos$Sequence_identifier))]
Suifus1_pos_fastas<- Suifus1_in[c(which(names(Suifus1_in) %in% Suifus1_GPCRHMM_pos$Sequence_identifier))]
#non-Suillus set
Rhivul1_pos_fastas<- Rhivul1_in[c(which(names(Rhivul1_in) %in% Rhivul1_GPCRHMM_pos$Sequence_identifier))]
Rhitru1_pos_fastas<- Rhitru1_in[c(which(names(Rhitru1_in) %in% Rhitru1_GPCRHMM_pos$Sequence_identifier))]
Amamu1_pos_fastas<- Amamu1_in[c(which(names(Amamu1_in) %in% Amamu1_GPCRHMM_pos$Sequence_identifier))]
Hebcy2_pos_fastas<- Hebcy2_in[c(which(names(Hebcy2_in) %in% Hebcy2_GPCRHMM_pos$Sequence_identifier))]
Lacbi2_pos_fastas<- Lacbi2_in[c(which(names(Lacbi2_in) %in% Lacbi2_GPCRHMM_pos$Sequence_identifier))]
Paxin1_pos_fastas<- Paxin1_in[c(which(names(Paxin1_in) %in% Paxin1_GPCRHMM_pos$Sequence_identifier))]
Pilcr1_pos_fastas<- Pilcr1_in[c(which(names(Pilcr1_in) %in% Pilcr1_GPCRHMM_pos$Sequence_identifier))]
Pismi1_pos_fastas<- Pismi1_in[c(which(names(Pismi1_in) %in% Pismi1_GPCRHMM_pos$Sequence_identifier))]
Sclci1_pos_fastas<- Sclci1_in[c(which(names(Sclci1_in) %in% Sclci1_GPCRHMM_pos$Sequence_identifier))]
Cananz1_pos_fastas<- Cananz1_in[c(which(names(Cananz1_in) %in% Cananz1_GPCRHMM_pos$Sequence_identifier))]
Gyrli1_pos_fastas<- Gyrli1_in[c(which(names(Gyrli1_in) %in% Gyrli1_GPCRHMM_pos$Sequence_identifier))]
Gaumor1_pos_fastas<- Gaumor1_in[c(which(names(Gaumor1_in) %in% Gaumor1_GPCRHMM_pos$Sequence_identifier))]
Hydru2_pos_fastas<- Hydru2_in[c(which(names(Hydru2_in) %in% Hydru2_GPCRHMM_pos$Sequence_identifier))]
Hyssto1_pos_fastas<- Hyssto1_in[c(which(names(Hyssto1_in) %in% Hyssto1_GPCRHMM_pos$Sequence_identifier))]
Lacam2_pos_fastas<- Lacam2_in[c(which(names(Lacam2_in) %in% Lacam2_GPCRHMM_pos$Sequence_identifier))]
Pisti1_pos_fastas<- Pisti1_in[c(which(names(Pisti1_in) %in% Pisti1_GPCRHMM_pos$Sequence_identifier))]
Rusbre1_pos_fastas<- Rusbre1_in[c(which(names(Rusbre1_in) %in% Rusbre1_GPCRHMM_pos$Sequence_identifier))]
Ruscom1_pos_fastas<- Ruscom1_in[c(which(names(Ruscom1_in) %in% Ruscom1_GPCRHMM_pos$Sequence_identifier))]
Rhisa1_pos_fastas<- Rhisa1_in[c(which(names(Rhisa1_in) %in% Rhisa1_GPCRHMM_pos$Sequence_identifier))]
Rhives1_pos_fastas<- Rhives1_in[c(which(names(Rhives1_in) %in% Rhives1_GPCRHMM_pos$Sequence_identifier))]
Rhivi1_pos_fastas<- Rhivi1_in[c(which(names(Rhivi1_in) %in% Rhivi1_GPCRHMM_pos$Sequence_identifier))]
Theter1_pos_fastas<- Theter1_in[c(which(names(Theter1_in) %in% Theter1_GPCRHMM_pos$Sequence_identifier))]
Thega1_pos_fastas<- Thega1_in[c(which(names(Thega1_in) %in% Thega1_GPCRHMM_pos$Sequence_identifier))]


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
write.fasta(Suigr1_pos_fastas, names = names(Suigr1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suigr1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suidec1_pos_fastas, names = names(Suidec1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suidec1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suicot1_pos_fastas, names = names(Suicot1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicot1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suicli1_pos_fastas, names = names(Suicli1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicli1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suibr2_pos_fastas, names = names(Suibr2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibr2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suibov1_pos_fastas, names = names(Suibov1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibov1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suiamp1_pos_fastas, names = names(Suiamp1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiamp1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suiame1_pos_fastas, names = names(Suiame1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiame1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Suifus1_pos_fastas, names = names(Suifus1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suifus1_GPCRs_from_GPCRHMM.fasta")

#non_Suillus set
write.fasta(Rhivul1_pos_fastas, names = names(Rhivul1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivul1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Rhitru1_pos_fastas, names = names(Rhitru1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhitru1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Amamu1_pos_fastas, names = names(Amamu1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Amamu1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Hebcy2_pos_fastas, names = names(Hebcy2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hebcy2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Lacbi2_pos_fastas, names = names(Lacbi2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacbi2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Paxin1_pos_fastas, names = names(Paxin1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Paxin1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Pilcr1_pos_fastas, names = names(Pilcr1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pilcr1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Pismi1_pos_fastas, names = names(Pismi1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pismi1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Sclci1_pos_fastas, names = names(Sclci1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Sclci1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Cananz1_pos_fastas, names = names(Cananz1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Cananz1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Gyrli1_pos_fastas, names = names(Gyrli1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gyrli1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Gaumor1_pos_fastas, names = names(Gaumor1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gaumor1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Hydru2_pos_fastas, names = names(Hydru2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hydru2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Hyssto1_pos_fastas, names = names(Hyssto1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hyssto1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Lacam2_pos_fastas, names = names(Lacam2_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacam2_GPCRs_from_GPCRHMM.fasta")
write.fasta(Pisti1_pos_fastas, names = names(Pisti1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pisti1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Rusbre1_pos_fastas, names = names(Rusbre1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rusbre1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Ruscom1_pos_fastas, names = names(Ruscom1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Ruscom1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Rhisa1_pos_fastas, names = names(Rhisa1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhisa1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Rhives1_pos_fastas, names = names(Rhives1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhives1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Rhivi1_pos_fastas, names = names(Rhivi1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivi1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Theter1_pos_fastas, names = names(Theter1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Theter1_GPCRs_from_GPCRHMM.fasta")
write.fasta(Thega1_pos_fastas, names = names(Thega1_pos_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Thega1_GPCRs_from_GPCRHMM.fasta")




##you are here in updated species list parcing analysis...




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
Suigr1_TMHMM<- data.frame(read.csv("Suigr1_TMHMM_output.txt", col.names = "header"))
Suidec1_TMHMM<- data.frame(read.csv("Suidec1_TMHMM_output.txt", col.names = "header"))
Suicot1_TMHMM<- data.frame(read.csv("Suicot1_TMHMM_output.txt", col.names = "header"))
Suicli1_TMHMM<- data.frame(read.csv("Suicli1_TMHMM_output.txt", col.names = "header"))
Suibr2_TMHMM<- data.frame(read.csv("Suibr2_TMHMM_output.txt", col.names = "header"))
Suibov1_TMHMM<- data.frame(read.csv("Suibov1_TMHMM_output.txt", col.names = "header"))
Suiamp1_TMHMM<- data.frame(read.csv("Suiamp1_TMHMM_output.txt", col.names = "header"))
Suiame1_TMHMM<- data.frame(read.csv("Suiame1_TMHMM_output.txt", col.names = "header"))
Suifus1_TMHMM<- data.frame(read.csv("Suifus1_TMHMM_output.txt", col.names = "header"))

#non_Suillus set
Rhivul1_TMHMM<- data.frame(read.csv("Rhivul1_TMHMM_output.txt", col.names = "header"))
Rhitru1_TMHMM<- data.frame(read.csv("Rhitru1_TMHMM_output.txt", col.names = "header"))
Amamu1_TMHMM<- data.frame(read.csv("Amamu1_TMHMM_output.txt", col.names = "header"))
Hebcy2_TMHMM<- data.frame(read.csv("Hebcy2_TMHMM_output.txt", col.names = "header"))
Lacbi2_TMHMM<- data.frame(read.csv("Lacbi2_TMHMM_output.txt", col.names = "header"))
Paxin1_TMHMM<- data.frame(read.csv("Paxin1_TMHMM_output.txt", col.names = "header"))
Pilcr1_TMHMM<- data.frame(read.csv("Pilcr1_TMHMM_output.txt", col.names = "header"))
Pismi1_TMHMM<- data.frame(read.csv("Pismi1_TMHMM_output.txt", col.names = "header"))
Sclci1_TMHMM<- data.frame(read.csv("Sclci1_TMHMM_output.txt", col.names = "header"))
Cananz1_TMHMM<- data.frame(read.csv("Cananz1_TMHMM_output.txt", col.names = "header"))
Gyrli1_TMHMM<- data.frame(read.csv("Gyrli1_TMHMM_output.txt", col.names = "header"))
Gaumor1_TMHMM<- data.frame(read.csv("Gaumor1_TMHMM_output.txt", col.names = "header"))
Hydru2_TMHMM<- data.frame(read.csv("Hydru2_TMHMM_output.txt", col.names = "header"))
Hyssto1_TMHMM<- data.frame(read.csv("Hyssto1_TMHMM_output.txt", col.names = "header"))
Lacam2_TMHMM<- data.frame(read.csv("Lacam2_TMHMM_output.txt", col.names = "header"))
Pisti1_TMHMM<- data.frame(read.csv("Pisti1_TMHMM_output.txt", col.names = "header"))
Rusbre1_TMHMM<- data.frame(read.csv("Rusbre1_TMHMM_output.txt", col.names = "header"))
Ruscom1_TMHMM<- data.frame(read.csv("Ruscom1_TMHMM_output.txt", col.names = "header"))
Rhisa1_TMHMM<- data.frame(read.csv("Rhisa1_TMHMM_output.txt", col.names = "header"))
Rhives1_TMHMM<- data.frame(read.csv("Rhives1_TMHMM_output.txt", col.names = "header"))
Rhivi1_TMHMM<- data.frame(read.csv("Rhivi1_TMHMM_output.txt", col.names = "header"))
Theter1_TMHMM<- data.frame(read.csv("Theter1_TMHMM_output.txt", col.names = "header"))
Thega1_TMHMM<- data.frame(read.csv("Thega1_TMHMM_output.txt", col.names = "header"))


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
Suigr1_TMHMM_sm<- data.frame(Suigr1_TMHMM[ grep("Number of predicted TMHs:  7", Suigr1_TMHMM$header),])
Suidec1_TMHMM_sm<- data.frame(Suidec1_TMHMM[ grep("Number of predicted TMHs:  7", Suidec1_TMHMM$header),])
Suicot1_TMHMM_sm<- data.frame(Suicot1_TMHMM[ grep("Number of predicted TMHs:  7", Suicot1_TMHMM$header),])
Suicli1_TMHMM_sm<- data.frame(Suicli1_TMHMM[ grep("Number of predicted TMHs:  7", Suicli1_TMHMM$header),])
Suibr2_TMHMM_sm<- data.frame(Suibr2_TMHMM[ grep("Number of predicted TMHs:  7", Suibr2_TMHMM$header),])
Suibov1_TMHMM_sm<- data.frame(Suibov1_TMHMM[ grep("Number of predicted TMHs:  7", Suibov1_TMHMM$header),])
Suiamp1_TMHMM_sm<- data.frame(Suiamp1_TMHMM[ grep("Number of predicted TMHs:  7", Suiamp1_TMHMM$header),])
Suiame1_TMHMM_sm<- data.frame(Suiame1_TMHMM[ grep("Number of predicted TMHs:  7", Suiame1_TMHMM$header),])
Suifus1_TMHMM_sm<- data.frame(Suifus1_TMHMM[ grep("Number of predicted TMHs:  7", Suifus1_TMHMM$header),])
#non_Suillus set
Rhivul1_TMHMM_sm<- data.frame(Rhivul1_TMHMM[ grep("Number of predicted TMHs:  7", Rhivul1_TMHMM$header),])
Rhitru1_TMHMM_sm<- data.frame(Rhitru1_TMHMM[ grep("Number of predicted TMHs:  7", Rhitru1_TMHMM$header),])
Amamu1_TMHMM_sm<- data.frame(Amamu1_TMHMM[ grep("Number of predicted TMHs:  7", Amamu1_TMHMM$header),])
Hebcy2_TMHMM_sm<- data.frame(Hebcy2_TMHMM[ grep("Number of predicted TMHs:  7", Hebcy2_TMHMM$header),])
Lacbi2_TMHMM_sm<- data.frame(Lacbi2_TMHMM[ grep("Number of predicted TMHs:  7", Lacbi2_TMHMM$header),])
Paxin1_TMHMM_sm<- data.frame(Paxin1_TMHMM[ grep("Number of predicted TMHs:  7", Paxin1_TMHMM$header),])
Pilcr1_TMHMM_sm<- data.frame(Pilcr1_TMHMM[ grep("Number of predicted TMHs:  7", Pilcr1_TMHMM$header),])
Pismi1_TMHMM_sm<- data.frame(Pismi1_TMHMM[ grep("Number of predicted TMHs:  7", Pismi1_TMHMM$header),])
Sclci1_TMHMM_sm<- data.frame(Sclci1_TMHMM[ grep("Number of predicted TMHs:  7", Sclci1_TMHMM$header),])
Cananz1_TMHMM_sm<- data.frame(Cananz1_TMHMM[ grep("Number of predicted TMHs:  7", Cananz1_TMHMM$header),])
Gyrli1_TMHMM_sm<- data.frame(Gyrli1_TMHMM[ grep("Number of predicted TMHs:  7", Gyrli1_TMHMM$header),])
Gaumor1_TMHMM_sm<- data.frame(Gaumor1_TMHMM[ grep("Number of predicted TMHs:  7", Gaumor1_TMHMM$header),])
Hydru2_TMHMM_sm<- data.frame(Hydru2_TMHMM[ grep("Number of predicted TMHs:  7", Hydru2_TMHMM$header),])
Hyssto1_TMHMM_sm<- data.frame(Hyssto1_TMHMM[ grep("Number of predicted TMHs:  7", Hyssto1_TMHMM$header),])
Lacam2_TMHMM_sm<- data.frame(Lacam2_TMHMM[ grep("Number of predicted TMHs:  7", Lacam2_TMHMM$header),])
Pisti1_TMHMM_sm<- data.frame(Pisti1_TMHMM[ grep("Number of predicted TMHs:  7", Pisti1_TMHMM$header),])
Rusbre1_TMHMM_sm<- data.frame(Rusbre1_TMHMM[ grep("Number of predicted TMHs:  7", Rusbre1_TMHMM$header),])
Ruscom1_TMHMM_sm<- data.frame(Ruscom1_TMHMM[ grep("Number of predicted TMHs:  7", Ruscom1_TMHMM$header),])
Rhisa1_TMHMM_sm<- data.frame(Rhisa1_TMHMM[ grep("Number of predicted TMHs:  7", Rhisa1_TMHMM$header),])
Rhives1_TMHMM_sm<- data.frame(Rhives1_TMHMM[ grep("Number of predicted TMHs:  7", Rhives1_TMHMM$header),])
Rhivi1_TMHMM_sm<- data.frame(Rhivi1_TMHMM[ grep("Number of predicted TMHs:  7", Rhivi1_TMHMM$header),])
Theter1_TMHMM_sm<- data.frame(Theter1_TMHMM[ grep("Number of predicted TMHs:  7", Theter1_TMHMM$header),])
Thega1_TMHMM_sm<- data.frame(Thega1_TMHMM[ grep("Number of predicted TMHs:  7", Thega1_TMHMM$header),])



#######
##Run phobius on the GPCRHMM results 
##available as a module at MSI
#module load phobius

##in the form: phobius -short Inputfile > outputfile
##phobius -short Suivar1_GPCRs_from_GPCRHMM.fasta > phobius_test.csv - does not run
##to run phobius
##perl /panfs/roc/msisoft/phobius/1.01/phobius.pl my.fasta

#perl /panfs/roc/msisoft/phobius/1.01/phobius.pl -short Sclci1_GPCRs_from_GPCRHMM.fasta > Phobius_Sclci1_output.txt

#to get a list of options you can run:
#perl /panfs/roc/msisoft/phobius/1.01/phobius.pl -h
#######



#read in phobius output files 
Suivar1_Phobius<- data.frame(read.csv("Phobius_Suivar1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suitom1_Phobius<- data.frame(read.csv("Phobius_Suitom1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suisub1_Phobius<- data.frame(read.csv("Phobius_Suisub1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suisu1_Phobius<- data.frame(read.csv("Phobius_Suisu1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suipla1_Phobius<- data.frame(read.csv("Phobius_Suipla1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suipic1_Phobius<- data.frame(read.csv("Phobius_Suipic1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suipal1_Phobius<- data.frame(read.csv("Phobius_Suipal1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suiocc1_Phobius<- data.frame(read.csv("Phobius_Suiocc1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suilu4_Phobius<- data.frame(read.csv("Phobius_Suilu4_output.txt", sep = "", header = TRUE, fill = TRUE))
Suilak1_Phobius<- data.frame(read.csv("Phobius_Suilak1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suihi1_Phobius<- data.frame(read.csv("Phobius_Suihi1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suigr1_Phobius<- data.frame(read.csv("Phobius_Suigr1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suidec1_Phobius<- data.frame(read.csv("Phobius_Suidec1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suicot1_Phobius<- data.frame(read.csv("Phobius_Suicot1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suicli1_Phobius<- data.frame(read.csv("Phobius_Suicli1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suibr2_Phobius<- data.frame(read.csv("Phobius_Suibr2_output.txt", sep = "", header = TRUE, fill = TRUE))
Suibov1_Phobius<- data.frame(read.csv("Phobius_Suibov1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suiamp1_Phobius<- data.frame(read.csv("Phobius_Suiamp1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suiame1_Phobius<- data.frame(read.csv("Phobius_Suiame1_output.txt", sep = "", header = TRUE, fill = TRUE))
Suifus1_Phobius<- data.frame(read.csv("Phobius_Suifus1_output.txt", sep = "", header = TRUE, fill = TRUE))

#non_Suillus set
Rhivul1_Phobius<- data.frame(read.csv("Phobius_Rhivul1_output.txt", sep = "", header = TRUE, fill = TRUE))
Rhitru1_Phobius<- data.frame(read.csv("Phobius_Rhitru1_output.txt", sep = "", header = TRUE, fill = TRUE))
Amamu1_Phobius<- data.frame(read.csv("Phobius_Amamu1_output.txt", sep = "", header = TRUE, fill = TRUE))
Hebcy2_Phobius<- data.frame(read.csv("Phobius_Hebcy2_output.txt", sep = "", header = TRUE, fill = TRUE))
Lacbi2_Phobius<- data.frame(read.csv("Phobius_Lacbi2_output.txt", sep = "", header = TRUE, fill = TRUE))
Paxin1_Phobius<- data.frame(read.csv("Phobius_Paxin1_output.txt", sep = "", header = TRUE, fill = TRUE))
Pilcr1_Phobius<- data.frame(read.csv("Phobius_Pilcr1_output.txt", sep = "", header = TRUE, fill = TRUE))
Pismi1_Phobius<- data.frame(read.csv("Phobius_Pismi1_output.txt", sep = "", header = TRUE, fill = TRUE))
Sclci1_Phobius<- data.frame(read.csv("Phobius_Sclci1_output.txt", sep = "", header = TRUE, fill = TRUE))
Cananz1_Phobius<- data.frame(read.csv("Phobius_Cananz1_output.txt", sep = "", header = TRUE, fill = TRUE))
Gyrli1_Phobius<- data.frame(read.csv("Phobius_Gyrli1_output.txt", sep = "", header = TRUE, fill = TRUE))
Gaumor1_Phobius<- data.frame(read.csv("Phobius_Gaumor1_output.txt", sep = "", header = TRUE, fill = TRUE))
Hydru2_Phobius<- data.frame(read.csv("Phobius_Hydru2_output.txt", sep = "", header = TRUE, fill = TRUE))
Hyssto1_Phobius<- data.frame(read.csv("Phobius_Hyssto1_output.txt", sep = "", header = TRUE, fill = TRUE))
Lacam2_Phobius<- data.frame(read.csv("Phobius_Lacam2_output.txt", sep = "", header = TRUE, fill = TRUE))
Pisti1_Phobius<- data.frame(read.csv("Phobius_Pisti1_output.txt", sep = "", header = TRUE, fill = TRUE))
Rusbre1_Phobius<- data.frame(read.csv("Phobius_Rusbre1_output.txt", sep = "", header = TRUE, fill = TRUE))
Ruscom1_Phobius<- data.frame(read.csv("Phobius_Ruscom1_output.txt", sep = "", header = TRUE, fill = TRUE))
Rhisa1_Phobius<- data.frame(read.csv("Phobius_Rhisa1_output.txt", sep = "", header = TRUE, fill = TRUE))
Rhives1_Phobius<- data.frame(read.csv("Phobius_Rhives1_output.txt", sep = "", header = TRUE, fill = TRUE))
Rhivi1_Phobius<- data.frame(read.csv("Phobius_Rhivi1_output.txt", sep = "", header = TRUE, fill = TRUE))
Theter1_Phobius<- data.frame(read.csv("Phobius_Theter1_output.txt", sep = "", header = TRUE, fill = TRUE))
Thega1_Phobius<- data.frame(read.csv("Phobius_Thega1_output.txt", sep = "", header = TRUE, fill = TRUE))


#get supset of each dataframe where TM == 7 and == 0
Suivar1_Phobius_sm1<- data.frame(Suivar1_Phobius[ grep("7", Suivar1_Phobius$TM),])
Suivar1_Phobius_sm<- data.frame(Suivar1_Phobius_sm1[ grep("0", Suivar1_Phobius_sm1$SP),])

Suitom1_Phobius_sm1<- data.frame(Suitom1_Phobius[ grep("7", Suitom1_Phobius$TM),])
Suitom1_Phobius_sm<- data.frame(Suitom1_Phobius_sm1[ grep("0", Suitom1_Phobius_sm1$SP),])

Suisub1_Phobius_sm1<- data.frame(Suisub1_Phobius[ grep("7", Suisub1_Phobius$TM),])
Suisub1_Phobius_sm<- data.frame(Suisub1_Phobius_sm1[ grep("0", Suisub1_Phobius_sm1$SP),])

Suisu1_Phobius_sm1<- data.frame(Suisu1_Phobius[ grep("7", Suisu1_Phobius$TM),])
Suisu1_Phobius_sm<- data.frame(Suisu1_Phobius_sm1[ grep("0", Suisu1_Phobius_sm1$SP),])

Suipla1_Phobius_sm1<- data.frame(Suipla1_Phobius[ grep("7", Suipla1_Phobius$TM),])
Suipla1_Phobius_sm<- data.frame(Suipla1_Phobius_sm1[ grep("0", Suipla1_Phobius_sm1$SP),])

Suipic1_Phobius_sm1<- data.frame(Suipic1_Phobius[ grep("7", Suipic1_Phobius$TM),])
Suipic1_Phobius_sm<- data.frame(Suipic1_Phobius_sm1[ grep("0", Suipic1_Phobius_sm1$SP),])

Suipal1_Phobius_sm1<- data.frame(Suipal1_Phobius[ grep("7", Suipal1_Phobius$TM),])
Suipal1_Phobius_sm<- data.frame(Suipal1_Phobius_sm1[ grep("0", Suipal1_Phobius_sm1$SP),])

Suiocc1_Phobius_sm1<- data.frame(Suiocc1_Phobius[ grep("7", Suiocc1_Phobius$TM),])
Suiocc1_Phobius_sm<- data.frame(Suiocc1_Phobius_sm1[ grep("0", Suiocc1_Phobius_sm1$SP),])

Suilu4_Phobius_sm1<- data.frame(Suilu4_Phobius[ grep("7", Suilu4_Phobius$TM),])
Suilu4_Phobius_sm<- data.frame(Suilu4_Phobius_sm1[ grep("0", Suilu4_Phobius_sm1$SP),])

Suilak1_Phobius_sm1<- data.frame(Suilak1_Phobius[ grep("7", Suilak1_Phobius$TM),])
Suilak1_Phobius_sm<- data.frame(Suilak1_Phobius_sm1[ grep("0", Suilak1_Phobius_sm1$SP),])

Suihi1_Phobius_sm1<- data.frame(Suihi1_Phobius[ grep("7", Suihi1_Phobius$TM),])
Suihi1_Phobius_sm<- data.frame(Suihi1_Phobius_sm1[ grep("0", Suihi1_Phobius_sm1$SP),])

Suigr1_Phobius_sm1<- data.frame(Suigr1_Phobius[ grep("7", Suigr1_Phobius$TM),])
Suigr1_Phobius_sm<- data.frame(Suigr1_Phobius_sm1[ grep("0", Suigr1_Phobius_sm1$SP),])

Suidec1_Phobius_sm1<- data.frame(Suidec1_Phobius[ grep("7", Suidec1_Phobius$TM),])
Suidec1_Phobius_sm<- data.frame(Suidec1_Phobius_sm1[ grep("0", Suidec1_Phobius_sm1$SP),])

Suicot1_Phobius_sm1<- data.frame(Suicot1_Phobius[ grep("7", Suicot1_Phobius$TM),])
Suicot1_Phobius_sm<- data.frame(Suicot1_Phobius_sm1[ grep("0", Suicot1_Phobius_sm1$SP),])

Suicli1_Phobius_sm1<- data.frame(Suicli1_Phobius[ grep("7", Suicli1_Phobius$TM),])
Suicli1_Phobius_sm<- data.frame(Suicli1_Phobius_sm1[ grep("0", Suicli1_Phobius_sm1$SP),])

Suibr2_Phobius_sm1<- data.frame(Suibr2_Phobius[ grep("7", Suibr2_Phobius$TM),])
Suibr2_Phobius_sm<- data.frame(Suibr2_Phobius_sm1[ grep("0", Suibr2_Phobius_sm1$SP),])

Suibov1_Phobius_sm1<- data.frame(Suibov1_Phobius[ grep("7", Suibov1_Phobius$TM),])
Suibov1_Phobius_sm<- data.frame(Suibov1_Phobius_sm1[ grep("0", Suibov1_Phobius_sm1$SP),])

Suiamp1_Phobius_sm1<- data.frame(Suiamp1_Phobius[ grep("7", Suiamp1_Phobius$TM),])
Suiamp1_Phobius_sm<- data.frame(Suiamp1_Phobius_sm1[ grep("0", Suiamp1_Phobius_sm1$SP),])

Suiame1_Phobius_sm1<- data.frame(Suiame1_Phobius[ grep("7", Suiame1_Phobius$TM),])
Suiame1_Phobius_sm<- data.frame(Suiame1_Phobius_sm1[ grep("0", Suiame1_Phobius_sm1$SP),])

Suifus1_Phobius_sm1<- data.frame(Suifus1_Phobius[ grep("7", Suifus1_Phobius$TM),])
Suifus1_Phobius_sm<- data.frame(Suifus1_Phobius_sm1[ grep("0", Suifus1_Phobius_sm1$SP),])


#non_Suillus set
Rhivul1_Phobius_sm1<- data.frame(Rhivul1_Phobius[ grep("7", Rhivul1_Phobius$TM),])
Rhivul1_Phobius_sm<- data.frame(Rhivul1_Phobius_sm1[ grep("0", Rhivul1_Phobius_sm1$SP),])

Rhitru1_Phobius_sm1<- data.frame(Rhitru1_Phobius[ grep("7", Rhitru1_Phobius$TM),])
Rhitru1_Phobius_sm<- data.frame(Rhitru1_Phobius_sm1[ grep("0", Rhitru1_Phobius_sm1$SP),])

Amamu1_Phobius_sm1<- data.frame(Amamu1_Phobius[ grep("7", Amamu1_Phobius$TM),])
Amamu1_Phobius_sm<- data.frame(Amamu1_Phobius_sm1[ grep("0", Amamu1_Phobius_sm1$SP),])

Hebcy2_Phobius_sm1<- data.frame(Hebcy2_Phobius[ grep("7", Hebcy2_Phobius$TM),])
Hebcy2_Phobius_sm<- data.frame(Hebcy2_Phobius_sm1[ grep("0", Hebcy2_Phobius_sm1$SP),])

Lacbi2_Phobius_sm1<- data.frame(Lacbi2_Phobius[ grep("7", Lacbi2_Phobius$TM),])
Lacbi2_Phobius_sm<- data.frame(Lacbi2_Phobius_sm1[ grep("0", Lacbi2_Phobius_sm1$SP),])

Paxin1_Phobius_sm1<- data.frame(Paxin1_Phobius[ grep("7", Paxin1_Phobius$TM),])
Paxin1_Phobius_sm<- data.frame(Paxin1_Phobius_sm1[ grep("0", Paxin1_Phobius_sm1$SP),])

Pilcr1_Phobius_sm1<- data.frame(Pilcr1_Phobius[ grep("7", Pilcr1_Phobius$TM),])
Pilcr1_Phobius_sm<- data.frame(Pilcr1_Phobius_sm1[ grep("0", Pilcr1_Phobius_sm1$SP),])

Pismi1_Phobius_sm1<- data.frame(Pismi1_Phobius[ grep("7", Pismi1_Phobius$TM),])
Pismi1_Phobius_sm<- data.frame(Pismi1_Phobius_sm1[ grep("0", Pismi1_Phobius_sm1$SP),])

Sclci1_Phobius_sm1<- data.frame(Sclci1_Phobius[ grep("7", Sclci1_Phobius$TM),])
Sclci1_Phobius_sm<- data.frame(Sclci1_Phobius_sm1[ grep("0", Sclci1_Phobius_sm1$SP),])

Cananz1_Phobius_sm1<- data.frame(Cananz1_Phobius[ grep("7", Cananz1_Phobius$TM),])
Cananz1_Phobius_sm<- data.frame(Cananz1_Phobius_sm1[ grep("0", Cananz1_Phobius_sm1$SP),])

Gyrli1_Phobius_sm1<- data.frame(Gyrli1_Phobius[ grep("7", Gyrli1_Phobius$TM),])
Gyrli1_Phobius_sm<- data.frame(Gyrli1_Phobius_sm1[ grep("0", Gyrli1_Phobius_sm1$SP),])

Gaumor1_Phobius_sm1<- data.frame(Gaumor1_Phobius[ grep("7", Gaumor1_Phobius$TM),])
Gaumor1_Phobius_sm<- data.frame(Gaumor1_Phobius_sm1[ grep("0", Gaumor1_Phobius_sm1$SP),])

Hydru2_Phobius_sm1<- data.frame(Hydru2_Phobius[ grep("7", Hydru2_Phobius$TM),])
Hydru2_Phobius_sm<- data.frame(Hydru2_Phobius_sm1[ grep("0", Hydru2_Phobius_sm1$SP),])

Hyssto1_Phobius_sm1<- data.frame(Hyssto1_Phobius[ grep("7", Hyssto1_Phobius$TM),])
Hyssto1_Phobius_sm<- data.frame(Hyssto1_Phobius_sm1[ grep("0", Hyssto1_Phobius_sm1$SP),])

Lacam2_Phobius_sm1<- data.frame(Lacam2_Phobius[ grep("7", Lacam2_Phobius$TM),])
Lacam2_Phobius_sm<- data.frame(Lacam2_Phobius_sm1[ grep("0", Lacam2_Phobius_sm1$SP),])

Pisti1_Phobius_sm1<- data.frame(Pisti1_Phobius[ grep("7", Pisti1_Phobius$TM),])
Pisti1_Phobius_sm<- data.frame(Pisti1_Phobius_sm1[ grep("0", Pisti1_Phobius_sm1$SP),])

Rusbre1_Phobius_sm1<- data.frame(Rusbre1_Phobius[ grep("7", Rusbre1_Phobius$TM),])
Rusbre1_Phobius_sm<- data.frame(Rusbre1_Phobius_sm1[ grep("0", Rusbre1_Phobius_sm1$SP),])

Ruscom1_Phobius_sm1<- data.frame(Ruscom1_Phobius[ grep("7", Ruscom1_Phobius$TM),])
Ruscom1_Phobius_sm<- data.frame(Ruscom1_Phobius_sm1[ grep("0", Ruscom1_Phobius_sm1$SP),])

Rhisa1_Phobius_sm1<- data.frame(Rhisa1_Phobius[ grep("7", Rhisa1_Phobius$TM),])
Rhisa1_Phobius_sm<- data.frame(Rhisa1_Phobius_sm1[ grep("0", Rhisa1_Phobius_sm1$SP),])

Rhives1_Phobius_sm1<- data.frame(Rhives1_Phobius[ grep("7", Rhives1_Phobius$TM),])
Rhives1_Phobius_sm<- data.frame(Rhives1_Phobius_sm1[ grep("0", Rhives1_Phobius_sm1$SP),])

Rhivi1_Phobius_sm1<- data.frame(Rhivi1_Phobius[ grep("7", Rhivi1_Phobius$TM),])
Rhivi1_Phobius_sm<- data.frame(Rhivi1_Phobius_sm1[ grep("0", Rhivi1_Phobius_sm1$SP),])

Theter1_Phobius_sm1<- data.frame(Theter1_Phobius[ grep("7", Theter1_Phobius$TM),])
Theter1_Phobius_sm<- data.frame(Theter1_Phobius_sm1[ grep("0", Theter1_Phobius_sm1$SP),])

Thega1_Phobius_sm1<- data.frame(Thega1_Phobius[ grep("7", Thega1_Phobius$TM),])
Thega1_Phobius_sm<- data.frame(Thega1_Phobius_sm1[ grep("0", Thega1_Phobius_sm1$SP),])

#get summary numbers for how many GPCR's were predicted with GPCRHMM 
Suivar1<- length(Suivar1_pos_fastas)
Suitom1<- length(Suitom1_pos_fastas)
Suisub1<- length(Suisub1_pos_fastas)
Suisu1<- length(Suisu1_pos_fastas)
Suipla1<- length(Suipla1_pos_fastas)
Suipic1<- length(Suipic1_pos_fastas)
Suipal1<- length(Suipal1_pos_fastas)
Suiocc1<- length(Suiocc1_pos_fastas)
Suilu4<- length(Suilu4_pos_fastas)
Suilak1<- length(Suilak1_pos_fastas)
Suihi1<- length(Suihi1_pos_fastas)
Suigr1<- length(Suigr1_pos_fastas)
Suidec1<- length(Suidec1_pos_fastas)
Suicot1<- length(Suicot1_pos_fastas)
Suicli1<- length(Suicli1_pos_fastas)
Suibr2<- length(Suibr2_pos_fastas)
Suibov1<- length(Suibov1_pos_fastas)
Suiamp1<- length(Suiamp1_pos_fastas)
Suiame1<- length(Suiame1_pos_fastas)
Suifus1<- length(Suifus1_pos_fastas)
#non_Suillus set
Rhivul1<- length(Rhivul1_pos_fastas)
Rhitru1<- length(Rhitru1_pos_fastas)
Amamu1<- length(Amamu1_pos_fastas)
Hebcy2<- length(Hebcy2_pos_fastas)
Lacbi2<- length(Lacbi2_pos_fastas)
Paxin1<- length(Paxin1_pos_fastas)
Pilcr1<- length(Pilcr1_pos_fastas)
Pismi1<- length(Pismi1_pos_fastas)
Sclci1<- length(Sclci1_pos_fastas)
Cananz1<- length(Cananz1_pos_fastas)
Gyrli1<- length(Gyrli1_pos_fastas)
Gaumor1<- length(Gaumor1_pos_fastas)
Hydru2<- length(Hydru2_pos_fastas)
Hyssto1<- length(Hyssto1_pos_fastas)
Lacam2<- length(Lacam2_pos_fastas)
Pisti1<- length(Pisti1_pos_fastas)
Rusbre1<- length(Rusbre1_pos_fastas)
Ruscom1<- length(Ruscom1_pos_fastas)
Rhisa1<- length(Rhisa1_pos_fastas)
Rhives1<- length(Rhives1_pos_fastas)
Rhivi1<- length(Rhivi1_pos_fastas)
Theter1<- length(Theter1_pos_fastas)
Thega1<- length(Thega1_pos_fastas)



totals.1<- data.frame(cbind(Suivar1, 
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
                            Suigr1,
                            Suidec1, 
                            Suicot1, 
                            Suicli1, 
                            Suibr2, 
                            Suibov1, 
                            Suiamp1, 
                            Suiame1,
                            Suifus1,
                            Rhivul1,
                            Rhitru1,
                            Amamu1,
                            Hebcy2,
                            Lacbi2,
                            Paxin1,
                            Pilcr1,
                            Pismi1,
                            Sclci1,
                            Cananz1,
                            Gyrli1,
                            Gaumor1,
                            Hydru2,
                            Hyssto1,
                            Lacam2,
                            Pisti1,
                            Rusbre1,
                            Ruscom1,
                            Rhisa1,
                            Rhives1,
                            Rhivi1,
                            Theter1,
                            Thega1),
                      row.names = "#GPCR's_w_GPCRHMM")


#get summary numbers for just the GPCRHMM and TMHMM parcing
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
Suigr1<- nrow(Suigr1_TMHMM_sm)
Suidec1<- nrow(Suidec1_TMHMM_sm)
Suicot1<- nrow(Suicot1_TMHMM_sm)
Suicli1<- nrow(Suicli1_TMHMM_sm)
Suibr2<- nrow(Suibr2_TMHMM_sm)
Suibov1<- nrow(Suibov1_TMHMM_sm)
Suiamp1<- nrow(Suiamp1_TMHMM_sm)
Suiame1<- nrow(Suiame1_TMHMM_sm)
Suifus1<- nrow(Suifus1_TMHMM_sm)
#non_Suillus set
Rhivul1<- nrow(Rhivul1_TMHMM_sm)
Rhitru1<- nrow(Rhitru1_TMHMM_sm)
Amamu1<- nrow(Amamu1_TMHMM_sm)
Hebcy2<- nrow(Hebcy2_TMHMM_sm)
Lacbi2<- nrow(Lacbi2_TMHMM_sm)
Paxin1<- nrow(Paxin1_TMHMM_sm)
Pilcr1<- nrow(Pilcr1_TMHMM_sm)
Pismi1<- nrow(Pismi1_TMHMM_sm)
Sclci1<- nrow(Sclci1_TMHMM_sm)
Cananz1<- nrow(Cananz1_TMHMM_sm)
Gyrli1<- nrow(Gyrli1_TMHMM_sm)
Gaumor1<- nrow(Gaumor1_TMHMM_sm)
Hydru2<- nrow(Hydru2_TMHMM_sm)
Hyssto1<- nrow(Hyssto1_TMHMM_sm)
Lacam2<- nrow(Lacam2_TMHMM_sm)
Pisti1<- nrow(Pisti1_TMHMM_sm)
Rusbre1<- nrow(Rusbre1_TMHMM_sm)
Ruscom1<- nrow(Ruscom1_TMHMM_sm)
Rhisa1<- nrow(Rhisa1_TMHMM_sm)
Rhives1<- nrow(Rhives1_TMHMM_sm)
Rhivi1<- nrow(Rhivi1_TMHMM_sm)
Theter1<- nrow(Theter1_TMHMM_sm)
Thega1<- nrow(Thega1_TMHMM_sm)



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
                            Suigr1,
                            Suidec1, 
                            Suicot1, 
                            Suicli1, 
                            Suibr2, 
                            Suibov1, 
                            Suiamp1, 
                            Suiame1,
                            Suifus1,
                            Rhivul1,
                            Rhitru1,
                            Amamu1,
                            Hebcy2,
                            Lacbi2,
                            Paxin1,
                            Pilcr1,
                            Pismi1,
                            Sclci1,
                            Cananz1,
                            Gyrli1,
                            Gaumor1,
                            Hydru2,
                            Hyssto1,
                            Lacam2,
                            Pisti1,
                            Rusbre1,
                            Ruscom1,
                            Rhisa1,
                            Rhives1,
                            Rhivi1,
                            Theter1,
                            Thega1),
                      row.names = "#GPCR's_w_GPCRHMM_and_TMHMM")

#get summary numbers for just the GPCRHMM + Phobius parcing 
Suivar1<- nrow(Suivar1_Phobius_sm)
Suitom1<- nrow(Suitom1_Phobius_sm)
Suisub1<- nrow(Suisub1_Phobius_sm)
Suisu1<- nrow(Suisu1_Phobius_sm)
Suipla1<- nrow(Suipla1_Phobius_sm)
Suipic1<- nrow(Suipic1_Phobius_sm)
Suipal1<- nrow(Suipal1_Phobius_sm)
Suiocc1<- nrow(Suiocc1_Phobius_sm)
Suilu4<- nrow(Suilu4_Phobius_sm)
Suilak1<- nrow(Suilak1_Phobius_sm)
Suihi1<- nrow(Suihi1_Phobius_sm)
Suigr1<- nrow(Suigr1_Phobius_sm)
Suidec1<- nrow(Suidec1_Phobius_sm)
Suicot1<- nrow(Suicot1_Phobius_sm)
Suicli1<- nrow(Suicli1_Phobius_sm)
Suibr2<- nrow(Suibr2_Phobius_sm)
Suibov1<- nrow(Suibov1_Phobius_sm)
Suiamp1<- nrow(Suiamp1_Phobius_sm)
Suiame1<- nrow(Suiame1_Phobius_sm)
Suifus1<- nrow(Suifus1_Phobius_sm)
#non_Suillus set
Rhivul1<- nrow(Rhivul1_Phobius_sm)
Rhitru1<- nrow(Rhitru1_Phobius_sm)
Amamu1<- nrow(Amamu1_Phobius_sm)
Hebcy2<- nrow(Hebcy2_Phobius_sm)
Lacbi2<- nrow(Lacbi2_Phobius_sm)
Paxin1<- nrow(Paxin1_Phobius_sm)
Pilcr1<- nrow(Pilcr1_Phobius_sm)
Pismi1<- nrow(Pismi1_Phobius_sm)
Sclci1<- nrow(Sclci1_Phobius_sm)
Cananz1<- nrow(Cananz1_Phobius_sm)
Gyrli1<- nrow(Gyrli1_Phobius_sm)
Gaumor1<- nrow(Gaumor1_Phobius_sm)
Hydru2<- nrow(Hydru2_Phobius_sm)
Hyssto1<- nrow(Hyssto1_Phobius_sm)
Lacam2<- nrow(Lacam2_Phobius_sm)
Pisti1<- nrow(Pisti1_Phobius_sm)
Rusbre1<- nrow(Rusbre1_Phobius_sm)
Ruscom1<- nrow(Ruscom1_Phobius_sm)
Rhisa1<- nrow(Rhisa1_Phobius_sm)
Rhives1<- nrow(Rhives1_Phobius_sm)
Rhivi1<- nrow(Rhivi1_Phobius_sm)
Theter1<- nrow(Theter1_Phobius_sm)
Thega1<- nrow(Thega1_Phobius_sm)


totals.3<- data.frame(cbind(Suivar1, 
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
                            Suigr1,
                            Suidec1, 
                            Suicot1, 
                            Suicli1, 
                            Suibr2, 
                            Suibov1, 
                            Suiamp1, 
                            Suiame1,
                            Suifus1,
                            Rhivul1,
                            Rhitru1,
                            Amamu1,
                            Hebcy2,
                            Lacbi2,
                            Paxin1,
                            Pilcr1,
                            Pismi1,
                            Sclci1,
                            Cananz1,
                            Gyrli1,
                            Gaumor1,
                            Hydru2,
                            Hyssto1,
                            Lacam2,
                            Pisti1,
                            Rusbre1,
                            Ruscom1,
                            Rhisa1,
                            Rhives1,
                            Rhivi1,
                            Theter1,
                            Thega1),
                      row.names = "#GPCR's_w_GPCRHMM_and_Phobius")

#write out totals to a file
results<- rbind(totals.1, totals.2, totals.3)

write.csv(results, quote = FALSE, file = "GPCR_totals_alloutgroups.csv")




#Get fastas of the positive phobius results to run in GPCR-CA and pull p-fams
Suivar1_pos_fastas_from_Phobius<- Suivar1_in[c(which(names(Suivar1_in) %in% Suivar1_Phobius_sm$SEQENCE_ID))]
Suitom1_pos_fastas_from_Phobius<- Suitom1_in[c(which(names(Suitom1_in) %in% Suitom1_Phobius_sm$SEQENCE_ID))]
Suisub1_pos_fastas_from_Phobius<- Suisub1_in[c(which(names(Suisub1_in) %in% Suisub1_Phobius_sm$SEQENCE_ID))]
Suisu1_pos_fastas_from_Phobius<- Suisu1_in[c(which(names(Suisu1_in) %in% Suisu1_Phobius_sm$SEQENCE_ID))]
Suipla1_pos_fastas_from_Phobius<- Suipla1_in[c(which(names(Suipla1_in) %in% Suipla1_Phobius_sm$SEQENCE_ID))]
Suipic1_pos_fastas_from_Phobius<- Suipic1_in[c(which(names(Suipic1_in) %in% Suipic1_Phobius_sm$SEQENCE_ID))]
Suipal1_pos_fastas_from_Phobius<- Suipal1_in[c(which(names(Suipal1_in) %in% Suipal1_Phobius_sm$SEQENCE_ID))]
Suiocc1_pos_fastas_from_Phobius<- Suiocc1_in[c(which(names(Suiocc1_in) %in% Suiocc1_Phobius_sm$SEQENCE_ID))]
Suilu4_pos_fastas_from_Phobius<- Suilu4_in[c(which(names(Suilu4_in) %in% Suilu4_Phobius_sm$SEQENCE_ID))]
Suilak1_pos_fastas_from_Phobius<- Suilak1_in[c(which(names(Suilak1_in) %in% Suilak1_Phobius_sm$SEQENCE_ID))]
Suihi1_pos_fastas_from_Phobius<- Suihi1_in[c(which(names(Suihi1_in) %in% Suihi1_Phobius_sm$SEQENCE_ID))]
Suigr1_pos_fastas_from_Phobius<- Suigr1_in[c(which(names(Suigr1_in) %in% Suigr1_Phobius_sm$SEQENCE_ID))]
Suidec1_pos_fastas_from_Phobius<- Suidec1_in[c(which(names(Suidec1_in) %in% Suidec1_Phobius_sm$SEQENCE_ID))]
Suicot1_pos_fastas_from_Phobius<- Suicot1_in[c(which(names(Suicot1_in) %in% Suicot1_Phobius_sm$SEQENCE_ID))]
Suicli1_pos_fastas_from_Phobius<- Suicli1_in[c(which(names(Suicli1_in) %in% Suicli1_Phobius_sm$SEQENCE_ID))]
Suibr2_pos_fastas_from_Phobius<- Suibr2_in[c(which(names(Suibr2_in) %in% Suibr2_Phobius_sm$SEQENCE_ID))]
Suibov1_pos_fastas_from_Phobius<- Suibov1_in[c(which(names(Suibov1_in) %in% Suibov1_Phobius_sm$SEQENCE_ID))]
Suiamp1_pos_fastas_from_Phobius<- Suiamp1_in[c(which(names(Suiamp1_in) %in% Suiamp1_Phobius_sm$SEQENCE_ID))]
Suiame1_pos_fastas_from_Phobius<- Suiame1_in[c(which(names(Suiame1_in) %in% Suiame1_Phobius_sm$SEQENCE_ID))]
Suifus1_pos_fastas_from_Phobius<- Suifus1_in[c(which(names(Suifus1_in) %in% Suifus1_Phobius_sm$SEQENCE_ID))]

#non-Suillus set
Rhivul1_pos_fastas_from_Phobius<- Rhivul1_in[c(which(names(Rhivul1_in) %in% Rhivul1_Phobius_sm$SEQENCE_ID))]
Rhitru1_pos_fastas_from_Phobius<- Rhitru1_in[c(which(names(Rhitru1_in) %in% Rhitru1_Phobius_sm$SEQENCE_ID))]
Amamu1_pos_fastas_from_Phobius<- Amamu1_in[c(which(names(Amamu1_in) %in% Amamu1_Phobius_sm$SEQENCE_ID))]
Hebcy2_pos_fastas_from_Phobius<- Hebcy2_in[c(which(names(Hebcy2_in) %in% Hebcy2_Phobius_sm$SEQENCE_ID))]
Lacbi2_pos_fastas_from_Phobius<- Lacbi2_in[c(which(names(Lacbi2_in) %in% Lacbi2_Phobius_sm$SEQENCE_ID))]
Paxin1_pos_fastas_from_Phobius<- Paxin1_in[c(which(names(Paxin1_in) %in% Paxin1_Phobius_sm$SEQENCE_ID))]
Pilcr1_pos_fastas_from_Phobius<- Pilcr1_in[c(which(names(Pilcr1_in) %in% Pilcr1_Phobius_sm$SEQENCE_ID))]
Pismi1_pos_fastas_from_Phobius<- Pismi1_in[c(which(names(Pismi1_in) %in% Pismi1_Phobius_sm$SEQENCE_ID))]
Sclci1_pos_fastas_from_Phobius<- Sclci1_in[c(which(names(Sclci1_in) %in% Sclci1_Phobius_sm$SEQENCE_ID))]
Cananz1_pos_fastas_from_Phobius<- Cananz1_in[c(which(names(Cananz1_in) %in% Cananz1_Phobius_sm$SEQENCE_ID))]
Gyrli1_pos_fastas_from_Phobius<- Gyrli1_in[c(which(names(Gyrli1_in) %in% Gyrli1_Phobius_sm$SEQENCE_ID))]
Gaumor1_pos_fastas_from_Phobius<- Gaumor1_in[c(which(names(Gaumor1_in) %in% Gaumor1_Phobius_sm$SEQENCE_ID))]
Hydru2_pos_fastas_from_Phobius<- Hydru2_in[c(which(names(Hydru2_in) %in% Hydru2_Phobius_sm$SEQENCE_ID))]
Hyssto1_pos_fastas_from_Phobius<- Hyssto1_in[c(which(names(Hyssto1_in) %in% Hyssto1_Phobius_sm$SEQENCE_ID))]
Lacam2_pos_fastas_from_Phobius<- Lacam2_in[c(which(names(Lacam2_in) %in% Lacam2_Phobius_sm$SEQENCE_ID))]
Pisti1_pos_fastas_from_Phobius<- Pisti1_in[c(which(names(Pisti1_in) %in% Pisti1_Phobius_sm$SEQENCE_ID))]
Rusbre1_pos_fastas_from_Phobius<- Rusbre1_in[c(which(names(Rusbre1_in) %in% Rusbre1_Phobius_sm$SEQENCE_ID))]
Ruscom1_pos_fastas_from_Phobius<- Ruscom1_in[c(which(names(Ruscom1_in) %in% Ruscom1_Phobius_sm$SEQENCE_ID))]
Rhisa1_pos_fastas_from_Phobius<- Rhisa1_in[c(which(names(Rhisa1_in) %in% Rhisa1_Phobius_sm$SEQENCE_ID))]
Rhives1_pos_fastas_from_Phobius<- Rhives1_in[c(which(names(Rhives1_in) %in% Rhives1_Phobius_sm$SEQENCE_ID))]
Rhivi1_pos_fastas_from_Phobius<- Rhivi1_in[c(which(names(Rhivi1_in) %in% Rhivi1_Phobius_sm$SEQENCE_ID))]
Theter1_pos_fastas_from_Phobius<- Theter1_in[c(which(names(Theter1_in) %in% Theter1_Phobius_sm$SEQENCE_ID))]
Thega1_pos_fastas_from_Phobius<- Thega1_in[c(which(names(Thega1_in) %in% Thega1_Phobius_sm$SEQENCE_ID))]


#print the fastas of the positive phobius results to run in GPCR-CA + pull pfams
write.fasta(Suivar1_pos_fastas_from_Phobius, names = names(Suivar1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suivar1_GPCRs_from_Phobius.fasta")
write.fasta(Suitom1_pos_fastas_from_Phobius, names = names(Suitom1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suitom1_GPCRs_from_Phobius.fasta")
write.fasta(Suisub1_pos_fastas_from_Phobius, names = names(Suisub1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisub1_GPCRs_from_Phobius.fasta")
write.fasta(Suisu1_pos_fastas_from_Phobius, names = names(Suisu1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisu1_GPCRs_from_Phobius.fasta")
write.fasta(Suipla1_pos_fastas_from_Phobius, names = names(Suipla1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipla1_GPCRs_from_Phobius.fasta")
write.fasta(Suipic1_pos_fastas_from_Phobius, names = names(Suipic1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipic1_GPCRs_from_Phobius.fasta")
write.fasta(Suipal1_pos_fastas_from_Phobius, names = names(Suipal1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipal1_GPCRs_from_Phobius.fasta")
write.fasta(Suiocc1_pos_fastas_from_Phobius, names = names(Suiocc1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiocc1_GPCRs_from_Phobius.fasta")
write.fasta(Suilu4_pos_fastas_from_Phobius, names = names(Suilu4_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilu4_GPCRs_from_Phobius.fasta")
write.fasta(Suilak1_pos_fastas_from_Phobius, names = names(Suilak1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilak1_GPCRs_from_Phobius.fasta")
write.fasta(Suihi1_pos_fastas_from_Phobius, names = names(Suihi1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suihi1_GPCRs_from_Phobius.fasta")
write.fasta(Suigr1_pos_fastas_from_Phobius, names = names(Suigr1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suigr1_GPCRs_from_Phobius.fasta")
write.fasta(Suidec1_pos_fastas_from_Phobius, names = names(Suidec1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suidec1_GPCRs_from_Phobius.fasta")
write.fasta(Suicot1_pos_fastas_from_Phobius, names = names(Suicot1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicot1_GPCRs_from_Phobius.fasta")
write.fasta(Suicli1_pos_fastas_from_Phobius, names = names(Suicli1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicli1_GPCRs_from_Phobius.fasta")
write.fasta(Suibr2_pos_fastas_from_Phobius, names = names(Suibr2_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibr2_GPCRs_from_Phobius.fasta")
write.fasta(Suibov1_pos_fastas_from_Phobius, names = names(Suibov1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibov1_GPCRs_from_Phobius.fasta")
write.fasta(Suiamp1_pos_fastas_from_Phobius, names = names(Suiamp1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiamp1_GPCRs_from_Phobius.fasta")
write.fasta(Suiame1_pos_fastas_from_Phobius, names = names(Suiame1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiame1_GPCRs_from_Phobius.fasta")
write.fasta(Suifus1_pos_fastas_from_Phobius, names = names(Suifus1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suifus1_GPCRs_from_Phobius.fasta")

#non_Suillus set
write.fasta(Rhivul1_pos_fastas_from_Phobius, names = names(Rhivul1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivul1_GPCRs_from_Phobius.fasta")
write.fasta(Rhitru1_pos_fastas_from_Phobius, names = names(Rhitru1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhitru1_GPCRs_from_Phobius.fasta")
write.fasta(Amamu1_pos_fastas_from_Phobius, names = names(Amamu1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Amamu1_GPCRs_from_Phobius.fasta")
write.fasta(Hebcy2_pos_fastas_from_Phobius, names = names(Hebcy2_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hebcy2_GPCRs_from_Phobius.fasta")
write.fasta(Lacbi2_pos_fastas_from_Phobius, names = names(Lacbi2_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacbi2_GPCRs_from_Phobius.fasta")
write.fasta(Paxin1_pos_fastas_from_Phobius, names = names(Paxin1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Paxin1_GPCRs_from_Phobius.fasta")
write.fasta(Pilcr1_pos_fastas_from_Phobius, names = names(Pilcr1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pilcr1_GPCRs_from_Phobius.fasta")
write.fasta(Pismi1_pos_fastas_from_Phobius, names = names(Pismi1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pismi1_GPCRs_from_Phobius.fasta")
write.fasta(Sclci1_pos_fastas_from_Phobius, names = names(Sclci1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Sclci1_GPCRs_from_Phobius.fasta")
write.fasta(Cananz1_pos_fastas_from_Phobius, names = names(Cananz1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Cananz1_GPCRs_from_Phobius.fasta")
write.fasta(Gyrli1_pos_fastas_from_Phobius, names = names(Gyrli1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gyrli1_GPCRs_from_Phobius.fasta")
write.fasta(Gaumor1_pos_fastas_from_Phobius, names = names(Gaumor1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Gaumor1_GPCRs_from_Phobius.fasta")
write.fasta(Hydru2_pos_fastas_from_Phobius, names = names(Hydru2_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hydru2_GPCRs_from_Phobius.fasta")
write.fasta(Hyssto1_pos_fastas_from_Phobius, names = names(Hyssto1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Hyssto1_GPCRs_from_Phobius.fasta")
write.fasta(Lacam2_pos_fastas_from_Phobius, names = names(Lacam2_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Lacam2_GPCRs_from_Phobius.fasta")
write.fasta(Pisti1_pos_fastas_from_Phobius, names = names(Pisti1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Pisti1_GPCRs_from_Phobius.fasta")
write.fasta(Rusbre1_pos_fastas_from_Phobius, names = names(Rusbre1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rusbre1_GPCRs_from_Phobius.fasta")
write.fasta(Ruscom1_pos_fastas_from_Phobius, names = names(Ruscom1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Ruscom1_GPCRs_from_Phobius.fasta")
write.fasta(Rhisa1_pos_fastas_from_Phobius, names = names(Rhisa1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhisa1_GPCRs_from_Phobius.fasta")
write.fasta(Rhives1_pos_fastas_from_Phobius, names = names(Rhives1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhives1_GPCRs_from_Phobius.fasta")
write.fasta(Rhivi1_pos_fastas_from_Phobius, names = names(Rhivi1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Rhivi1_GPCRs_from_Phobius.fasta")
write.fasta(Theter1_pos_fastas_from_Phobius, names = names(Theter1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Theter1_GPCRs_from_Phobius.fasta")
write.fasta(Thega1_pos_fastas_from_Phobius, names = names(Thega1_pos_fastas_from_Phobius), open = "w", nbchar = 60, as.string = FALSE, file.out = "Thega1_GPCRs_from_Phobius.fasta")

#run TOPCONS on GPCR lists 
#http://single.topcons.net///index.php
#import "all topologies" outputflle, save as a .csv file

#kill those factors
options(stringsAsFactors = FALSE)
#import files 
TOPCONS_Theter1<- read.table("TOPCONS_Theter1_GPCRS.txt", fill = TRUE, header = FALSE)

View(TOPCONS_Theter1)
#format TOPCONS files
Theter1.names<- TOPCONS_Theter1[grep(">", TOPCONS_Theter1$V1),]
Theter1.seq<- TOPCONS_Theter1[ - grep(">", TOPCONS_Theter1$V1),]
TOPCONS_Theter1<-cbind(Theter1.names, Theter1.seq)
TOPCONS_Theter1.2<- data.frame(gsub(">", "", TOPCONS_Theter1))
colnames(TOPCONS_Theter1.2)<- c("name", "seq")

#format the fasta files to match against
Theter1_fastas<-t(as.data.frame(Theter1_pos_fastas_from_Phobius))
Theter1_fastas2<- data.frame(cbind(row.names(Theter1_fastas), Theter1_fastas))
colnames(Theter1_fastas2)<- c("name", "seq")

#rename fastas 
seq.df<- data.frame(Theter1_fastas2[,2])
map.key<- data.frame(TOPCONS_Theter1.2[,2])
#get indices of "M" regions 
map.list.list<- apply(map.key, 2, function (x) str_locate_all(pattern ='M', x))
#get starts only
map.list.list<- map.list.list[[1]]
map.df<- mapply('[', map.list.list, TRUE, 1)

#itterate over all 
Theter1_cat_TM_list<- unlist(Map(function(x, y) paste(substring(x, unlist(y), 
                                                 unlist(y)), collapse=""), seq.df[[1]], map.df))

##re-attach sequence names and format back into fasta files. 
#get list for fasta headers
names<- TOPCONS_Theter1.2$name

#rename rows
Theter1_cat_TM_list<- data.frame(Theter1_cat_TM_list)
rownames(Theter1_cat_TM_list)<- TOPCONS_Theter1.2$name

#add the carrot back into the row names
rownames(Theter1_cat_TM_list) <- paste0(">",  rownames(Theter1_cat_TM_list))

#format and print this to a file 
write.table(Theter1_cat_TM_list, file = "Theter1_cat_TM_list.fasta", col.names=FALSE, quote = FALSE, sep = "\n")

#quick QC check - make sure the pro length is the same as the respective M lists. 
#str_count(map.key$TOPCONS_Theter1.2...2., "M")
#nchar(Theter1_cat_TM_list[1,1])
#nchar(Theter1_cat_TM_list[7,1])
#looks good. 

#align the outputfiles using 


#modify to parse each TM region seperately (so we can align them seprately)
#step 1, parse non-consecuative numbers into different bins in M_poistion

#first colapse M_position rows into a single row:
library("tidyr")
library("dplyr")
M_position<- data.frame(M_position)
#concatonate?
#M_position_colapse<-unite(M_position, new, 1:ncol(M_position), sep=" ", remove = TRUE)
#works, but is this really what you want? 

Breaks2<- cbind(apply(M_position, 1, function(x) apply(M_position, 2, function (y) c(0, which(diff(x) != 1), length(x)) )))
#almost . . . it's not giving different brake numbers for the different TMs though. . . 

#re-write to bin entries into diffrent lists
Breaks3<- cbind(apply(M_position, 1, function(x) apply(M_position, 2, function (x) c(0, which(diff(x) != 1), length(x)) )))


View(Breaks2)

sapply(seq(length(Breaks2) - 1), 
       function(i) Vec[(Breaks2[i] + 1):Breaks2[i+1]])



Vec<- c(1,2,3,4,7,8,9,10,12,13) 

class(Vec)
class(M_position_colapse[1,])

Breaks <- c(0, which(diff(Vec) != 1), length(Vec)) 

sapply(seq(length(Breaks) - 1), 
         function(i) Vec[(Breaks[i] + 1):Breaks[i+1]]) 





#######stats and graphs
#get two sets 
#transform
Results_t<-t(results)

#add a names col to search with 
names.list<- row.names(Results_t)
Results_t2<- data.frame(cbind(Results_t, names.list))

Results_Suillus<-data.frame(Results_t2[ grep("Sui", Results_t2[,4],), ])
Results_Suillus


Results_Other<-data.frame(Results_t2[ grep("Sui", Results_t2[,4],invert = TRUE), ])
Results_Other


#make repeating rows
Results_Suillus_sm<- data.frame(as.numeric(Results_Suillus[,3]))
Suillus<- cbind(Results_Suillus_sm, rep("a_Suillus", nrow(Results_Suillus_sm)))
colnames(Suillus)<- c("%", "group")


Results_Other_sm<- data.frame(as.numeric(Results_Other[,3]))
Other_ECM<- cbind(Results_Other_sm, rep("b_Other", nrow(Results_Other_sm)))
colnames(Other_ECM)<- c("%", "group")

GPCR_df_for_stats<- data.frame(rbind(Suillus, Other_ECM))


#######START HERE WITH GRAPHS######## 
med1= mean(as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "a_Suillus"])
med2= mean(as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "b_Other"])
par(mar = c(6.5, 8.5, 3, 3.5), mgp = c(6, 2.5, 0))
boxplot(GPCR_df_for_stats$X. ~ GPCR_df_for_stats$group)
stripchart(GPCR_df_for_stats$X. ~ GPCR_df_for_stats$group,
           vertical = TRUE,
           method = "jitter", jitter = 0.2, 
           pch = 16, 
           col = c("#F0502B",  "#F79552"),
           bg = rep(c("#F0502B", "#F79552")),
           cex.axis = 0.7,
           ylim=c(0,100), 
           ylab = "GPCR distribution", 
           axes = FALSE, 
           cex = 1.3, add = TRUE)
box()
axis(2)
mtext(text = c(round(med1,2), round(med2, 2)),side=1,at=c(1,2),line = 1, font = 3)
#segments(x0 = .7, y0 =  med1, x1 = 1.3, y1=med1, lwd = 2, col = "black" )
#segments(x0 = 1.7, y0 =  med2, x1 = 2.3, y1=med2, lwd = 2, col = "black" )

#stats on the above
#test normality 
shapiro.test(as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "a_Suillus"])
shapiro.test(as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "b_Other"])
#normal

#test for equal variance 
var.test((as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "a_Suillus"]), (as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "b_Other"]))
#variance is aprox. equal

#t-test 
t.test((as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "a_Suillus"]), (as.numeric(GPCR_df_for_stats$X.) [GPCR_df_for_stats$group == "b_Other"]))
#not different







#to class GPCR's

#read in the file 
Suivar_PCA_GPCR<- read.delim("Suivar_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suivar_PCA_GPCR_fam1<- data.frame(Suivar_PCA_GPCR[ grep("Family", Suivar_PCA_GPCR$X),])
colnames(Suivar_PCA_GPCR_fam1)<- "X"
Suivar_PCA_GPCR_fam2<- as.matrix(Suivar_PCA_GPCR_fam1[ grep("Sub", Suivar_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suivar_PCA_GPCR_fam2)<- "X"
Suivar_PCA_GPCR_fam3<- as.data.frame(Suivar_PCA_GPCR_fam2)
#get sub-family
Suivar_PCA_GPCR_subfam<- as.matrix(Suivar_PCA_GPCR_fam1[ grep("SubFamily", Suivar_PCA_GPCR_fam1$X,),])
Suivar_PCA_GPCR_subfam2<- as.data.frame(Suivar_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suivar_PCA_GPCR_subsubfam<- as.matrix(Suivar_PCA_GPCR_fam1[ grep("Sub-subFamily", Suivar_PCA_GPCR_fam1$X,),])
Suivar_PCA_GPCR_subsubfam2<-as.data.frame(Suivar_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suivar_PCA_GPCR_subtype<- as.matrix(Suivar_PCA_GPCR[ grep("Subtype", Suivar_PCA_GPCR$X,),])
Suivar_PCA_GPCR_subtype<- as.data.frame(Suivar_PCA_GPCR_subtype)
#bin
Suivar_GPCR_table_fam<- table(Suivar_PCA_GPCR_fam3)
Suivar_GPCR_table_subfam<- table(Suivar_PCA_GPCR_subfam2)
Suivar_GPCR_table_subsubfam<- table(Suivar_PCA_GPCR_subsubfam2)
Suivar_GPCR_table_subtype<- table(Suivar_PCA_GPCR_subtype)



#read in the file 
Suitom_PCA_GPCR<- read.delim("Suitom_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suitom_PCA_GPCR_fam1<- data.frame(Suitom_PCA_GPCR[ grep("Family", Suitom_PCA_GPCR$X),])
colnames(Suitom_PCA_GPCR_fam1)<- "X"
Suitom_PCA_GPCR_fam2<- as.matrix(Suitom_PCA_GPCR_fam1[ grep("Sub", Suitom_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suitom_PCA_GPCR_fam2)<- "X"
Suitom_PCA_GPCR_fam3<- as.data.frame(Suitom_PCA_GPCR_fam2)
#get sub-family
Suitom_PCA_GPCR_subfam<- as.matrix(Suitom_PCA_GPCR_fam1[ grep("SubFamily", Suitom_PCA_GPCR_fam1$X,),])
Suitom_PCA_GPCR_subfam2<- as.data.frame(Suitom_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suitom_PCA_GPCR_subsubfam<- as.matrix(Suitom_PCA_GPCR_fam1[ grep("Sub-subFamily", Suitom_PCA_GPCR_fam1$X,),])
Suitom_PCA_GPCR_subsubfam2<-as.data.frame(Suitom_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suitom_PCA_GPCR_subtype<- as.matrix(Suitom_PCA_GPCR[ grep("Subtype", Suitom_PCA_GPCR$X,),])
Suitom_PCA_GPCR_subtype<- as.data.frame(Suitom_PCA_GPCR_subtype)
#bin
Suitom_GPCR_table_fam<- table(Suitom_PCA_GPCR_fam3)
Suitom_GPCR_table_subfam<- table(Suitom_PCA_GPCR_subfam2)
Suitom_GPCR_table_subsubfam<- table(Suitom_PCA_GPCR_subsubfam2)
Suitom_GPCR_table_subtype<- table(Suitom_PCA_GPCR_subtype)


#read in the file 
Suisub_PCA_GPCR<- read.delim("Suisub_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suisub_PCA_GPCR_fam1<- data.frame(Suisub_PCA_GPCR[ grep("Family", Suisub_PCA_GPCR$X),])
colnames(Suisub_PCA_GPCR_fam1)<- "X"
Suisub_PCA_GPCR_fam2<- as.matrix(Suisub_PCA_GPCR_fam1[ grep("Sub", Suisub_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suisub_PCA_GPCR_fam2)<- "X"
Suisub_PCA_GPCR_fam3<- as.data.frame(Suisub_PCA_GPCR_fam2)
#get sub-family
Suisub_PCA_GPCR_subfam<- as.matrix(Suisub_PCA_GPCR_fam1[ grep("SubFamily", Suisub_PCA_GPCR_fam1$X,),])
Suisub_PCA_GPCR_subfam2<- as.data.frame(Suisub_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suisub_PCA_GPCR_subsubfam<- as.matrix(Suisub_PCA_GPCR_fam1[ grep("Sub-subFamily", Suisub_PCA_GPCR_fam1$X,),])
Suisub_PCA_GPCR_subsubfam2<-as.data.frame(Suisub_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suisub_PCA_GPCR_subtype<- as.matrix(Suisub_PCA_GPCR[ grep("Subtype", Suisub_PCA_GPCR$X,),])
Suisub_PCA_GPCR_subtype<- as.data.frame(Suisub_PCA_GPCR_subtype)
#bin
Suisub_GPCR_table_fam<- table(Suisub_PCA_GPCR_fam3)
Suisub_GPCR_table_subfam<- table(Suisub_PCA_GPCR_subfam2)
Suisub_GPCR_table_subsubfam<- table(Suisub_PCA_GPCR_subsubfam2)
Suisub_GPCR_table_subtype<- table(Suisub_PCA_GPCR_subtype)

#read in the file 
Suisu_PCA_GPCR<- read.delim("Suisu_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suisu_PCA_GPCR_fam1<- data.frame(Suisu_PCA_GPCR[ grep("Family", Suisu_PCA_GPCR$X),])
colnames(Suisu_PCA_GPCR_fam1)<- "X"
Suisu_PCA_GPCR_fam2<- as.matrix(Suisu_PCA_GPCR_fam1[ grep("Sub", Suisu_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suisu_PCA_GPCR_fam2)<- "X"
Suisu_PCA_GPCR_fam3<- as.data.frame(Suisu_PCA_GPCR_fam2)
#get sub-family
Suisu_PCA_GPCR_subfam<- as.matrix(Suisu_PCA_GPCR_fam1[ grep("SubFamily", Suisu_PCA_GPCR_fam1$X,),])
Suisu_PCA_GPCR_subfam2<- as.data.frame(Suisu_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suisu_PCA_GPCR_subsubfam<- as.matrix(Suisu_PCA_GPCR_fam1[ grep("Sub-subFamily", Suisu_PCA_GPCR_fam1$X,),])
Suisu_PCA_GPCR_subsubfam2<-as.data.frame(Suisu_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suisu_PCA_GPCR_subtype<- as.matrix(Suisu_PCA_GPCR[ grep("Subtype", Suisu_PCA_GPCR$X,),])
Suisu_PCA_GPCR_subtype<- as.data.frame(Suisu_PCA_GPCR_subtype)
#bin
Suisu_GPCR_table_fam<- table(Suisu_PCA_GPCR_fam3)
Suisu_GPCR_table_subfam<- table(Suisu_PCA_GPCR_subfam2)
Suisu_GPCR_table_subsubfam<- table(Suisu_PCA_GPCR_subsubfam2)
Suisu_GPCR_table_subtype<- table(Suisu_PCA_GPCR_subtype)

#read in the file 
Suipla_PCA_GPCR<- read.delim("Suipla_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suipla_PCA_GPCR_fam1<- data.frame(Suipla_PCA_GPCR[ grep("Family", Suipla_PCA_GPCR$X),])
colnames(Suipla_PCA_GPCR_fam1)<- "X"
Suipla_PCA_GPCR_fam2<- as.matrix(Suipla_PCA_GPCR_fam1[ grep("Sub", Suipla_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suipla_PCA_GPCR_fam2)<- "X"
Suipla_PCA_GPCR_fam3<- as.data.frame(Suipla_PCA_GPCR_fam2)
#get sub-family
Suipla_PCA_GPCR_subfam<- as.matrix(Suipla_PCA_GPCR_fam1[ grep("SubFamily", Suipla_PCA_GPCR_fam1$X,),])
Suipla_PCA_GPCR_subfam2<- as.data.frame(Suipla_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suipla_PCA_GPCR_subsubfam<- as.matrix(Suipla_PCA_GPCR_fam1[ grep("Sub-subFamily", Suipla_PCA_GPCR_fam1$X,),])
Suipla_PCA_GPCR_subsubfam2<-as.data.frame(Suipla_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suipla_PCA_GPCR_subtype<- as.matrix(Suipla_PCA_GPCR[ grep("Subtype", Suipla_PCA_GPCR$X,),])
Suipla_PCA_GPCR_subtype<- as.data.frame(Suipla_PCA_GPCR_subtype)
#bin
Suipla_GPCR_table_fam<- table(Suipla_PCA_GPCR_fam3)
Suipla_GPCR_table_subfam<- table(Suipla_PCA_GPCR_subfam2)
Suipla_GPCR_table_subsubfam<- table(Suipla_PCA_GPCR_subsubfam2)
Suipla_GPCR_table_subtype<- table(Suipla_PCA_GPCR_subtype)



#read in the file 
Suipic_PCA_GPCR<- read.delim("Suipic_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suipic_PCA_GPCR_fam1<- data.frame(Suipic_PCA_GPCR[ grep("Family", Suipic_PCA_GPCR$X),])
colnames(Suipic_PCA_GPCR_fam1)<- "X"
Suipic_PCA_GPCR_fam2<- as.matrix(Suipic_PCA_GPCR_fam1[ grep("Sub", Suipic_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suipic_PCA_GPCR_fam2)<- "X"
Suipic_PCA_GPCR_fam3<- as.data.frame(Suipic_PCA_GPCR_fam2)
#get sub-family
Suipic_PCA_GPCR_subfam<- as.matrix(Suipic_PCA_GPCR_fam1[ grep("SubFamily", Suipic_PCA_GPCR_fam1$X,),])
Suipic_PCA_GPCR_subfam2<- as.data.frame(Suipic_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suipic_PCA_GPCR_subsubfam<- as.matrix(Suipic_PCA_GPCR_fam1[ grep("Sub-subFamily", Suipic_PCA_GPCR_fam1$X,),])
Suipic_PCA_GPCR_subsubfam2<-as.data.frame(Suipic_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suipic_PCA_GPCR_subtype<- as.matrix(Suipic_PCA_GPCR[ grep("Subtype", Suipic_PCA_GPCR$X,),])
Suipic_PCA_GPCR_subtype<- as.data.frame(Suipic_PCA_GPCR_subtype)
#bin
Suipic_GPCR_table_fam<- table(Suipic_PCA_GPCR_fam3)
Suipic_GPCR_table_subfam<- table(Suipic_PCA_GPCR_subfam2)
Suipic_GPCR_table_subsubfam<- table(Suipic_PCA_GPCR_subsubfam2)
Suipic_GPCR_table_subtype<- table(Suipic_PCA_GPCR_subtype)



#read in the file 
Suipal_PCA_GPCR<- read.delim("Suipal_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suipal_PCA_GPCR_fam1<- data.frame(Suipal_PCA_GPCR[ grep("Family", Suipal_PCA_GPCR$X),])
colnames(Suipal_PCA_GPCR_fam1)<- "X"
Suipal_PCA_GPCR_fam2<- as.matrix(Suipal_PCA_GPCR_fam1[ grep("Sub", Suipal_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suipal_PCA_GPCR_fam2)<- "X"
Suipal_PCA_GPCR_fam3<- as.data.frame(Suipal_PCA_GPCR_fam2)
#get sub-family
Suipal_PCA_GPCR_subfam<- as.matrix(Suipal_PCA_GPCR_fam1[ grep("SubFamily", Suipal_PCA_GPCR_fam1$X,),])
Suipal_PCA_GPCR_subfam2<- as.data.frame(Suipal_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suipal_PCA_GPCR_subsubfam<- as.matrix(Suipal_PCA_GPCR_fam1[ grep("Sub-subFamily", Suipal_PCA_GPCR_fam1$X,),])
Suipal_PCA_GPCR_subsubfam2<-as.data.frame(Suipal_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suipal_PCA_GPCR_subtype<- as.matrix(Suipal_PCA_GPCR[ grep("Subtype", Suipal_PCA_GPCR$X,),])
Suipal_PCA_GPCR_subtype<- as.data.frame(Suipal_PCA_GPCR_subtype)
#bin
Suipal_GPCR_table_fam<- table(Suipal_PCA_GPCR_fam3)
Suipal_GPCR_table_subfam<- table(Suipal_PCA_GPCR_subfam2)
Suipal_GPCR_table_subsubfam<- table(Suipal_PCA_GPCR_subsubfam2)
Suipal_GPCR_table_subtype<- table(Suipal_PCA_GPCR_subtype)

#read in the file 
Suiocc_PCA_GPCR<- read.delim("Suiocc_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suiocc_PCA_GPCR_fam1<- data.frame(Suiocc_PCA_GPCR[ grep("Family", Suiocc_PCA_GPCR$X),])
colnames(Suiocc_PCA_GPCR_fam1)<- "X"
Suiocc_PCA_GPCR_fam2<- as.matrix(Suiocc_PCA_GPCR_fam1[ grep("Sub", Suiocc_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suiocc_PCA_GPCR_fam2)<- "X"
Suiocc_PCA_GPCR_fam3<- as.data.frame(Suiocc_PCA_GPCR_fam2)
#get sub-family
Suiocc_PCA_GPCR_subfam<- as.matrix(Suiocc_PCA_GPCR_fam1[ grep("SubFamily", Suiocc_PCA_GPCR_fam1$X,),])
Suiocc_PCA_GPCR_subfam2<- as.data.frame(Suiocc_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suiocc_PCA_GPCR_subsubfam<- as.matrix(Suiocc_PCA_GPCR_fam1[ grep("Sub-subFamily", Suiocc_PCA_GPCR_fam1$X,),])
Suiocc_PCA_GPCR_subsubfam2<-as.data.frame(Suiocc_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suiocc_PCA_GPCR_subtype<- as.matrix(Suiocc_PCA_GPCR[ grep("Subtype", Suiocc_PCA_GPCR$X,),])
Suiocc_PCA_GPCR_subtype<- as.data.frame(Suiocc_PCA_GPCR_subtype)
#bin
Suiocc_GPCR_table_fam<- table(Suiocc_PCA_GPCR_fam3)
Suiocc_GPCR_table_subfam<- table(Suiocc_PCA_GPCR_subfam2)
Suiocc_GPCR_table_subsubfam<- table(Suiocc_PCA_GPCR_subsubfam2)
Suiocc_GPCR_table_subtype<- table(Suiocc_PCA_GPCR_subtype)


#read in the file 
Suilu_PCA_GPCR<- read.delim("Suilu_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suilu_PCA_GPCR_fam1<- data.frame(Suilu_PCA_GPCR[ grep("Family", Suilu_PCA_GPCR$X),])
colnames(Suilu_PCA_GPCR_fam1)<- "X"
Suilu_PCA_GPCR_fam2<- as.matrix(Suilu_PCA_GPCR_fam1[ grep("Sub", Suilu_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suilu_PCA_GPCR_fam2)<- "X"
Suilu_PCA_GPCR_fam3<- as.data.frame(Suilu_PCA_GPCR_fam2)
#get sub-family
Suilu_PCA_GPCR_subfam<- as.matrix(Suilu_PCA_GPCR_fam1[ grep("SubFamily", Suilu_PCA_GPCR_fam1$X,),])
Suilu_PCA_GPCR_subfam2<- as.data.frame(Suilu_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suilu_PCA_GPCR_subsubfam<- as.matrix(Suilu_PCA_GPCR_fam1[ grep("Sub-subFamily", Suilu_PCA_GPCR_fam1$X,),])
Suilu_PCA_GPCR_subsubfam2<-as.data.frame(Suilu_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suilu_PCA_GPCR_subtype<- as.matrix(Suilu_PCA_GPCR[ grep("Subtype", Suilu_PCA_GPCR$X,),])
Suilu_PCA_GPCR_subtype<- as.data.frame(Suilu_PCA_GPCR_subtype)
#bin
Suilu_GPCR_table_fam<- table(Suilu_PCA_GPCR_fam3)
Suilu_GPCR_table_subfam<- table(Suilu_PCA_GPCR_subfam2)
Suilu_GPCR_table_subsubfam<- table(Suilu_PCA_GPCR_subsubfam2)
Suilu_GPCR_table_subtype<- table(Suilu_PCA_GPCR_subtype)

#read in the file 
Suilak_PCA_GPCR<- read.delim("Suilak_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suilak_PCA_GPCR_fam1<- data.frame(Suilak_PCA_GPCR[ grep("Family", Suilak_PCA_GPCR$X),])
colnames(Suilak_PCA_GPCR_fam1)<- "X"
Suilak_PCA_GPCR_fam2<- as.matrix(Suilak_PCA_GPCR_fam1[ grep("Sub", Suilak_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suilak_PCA_GPCR_fam2)<- "X"
Suilak_PCA_GPCR_fam3<- as.data.frame(Suilak_PCA_GPCR_fam2)
#get sub-family
Suilak_PCA_GPCR_subfam<- as.matrix(Suilak_PCA_GPCR_fam1[ grep("SubFamily", Suilak_PCA_GPCR_fam1$X,),])
Suilak_PCA_GPCR_subfam2<- as.data.frame(Suilak_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suilak_PCA_GPCR_subsubfam<- as.matrix(Suilak_PCA_GPCR_fam1[ grep("Sub-subFamily", Suilak_PCA_GPCR_fam1$X,),])
Suilak_PCA_GPCR_subsubfam2<-as.data.frame(Suilak_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suilak_PCA_GPCR_subtype<- as.matrix(Suilak_PCA_GPCR[ grep("Subtype", Suilak_PCA_GPCR$X,),])
Suilak_PCA_GPCR_subtype<- as.data.frame(Suilak_PCA_GPCR_subtype)
#bin
Suilak_GPCR_table_fam<- table(Suilak_PCA_GPCR_fam3)
Suilak_GPCR_table_subfam<- table(Suilak_PCA_GPCR_subfam2)
Suilak_GPCR_table_subsubfam<- table(Suilak_PCA_GPCR_subsubfam2)
Suilak_GPCR_table_subtype<- table(Suilak_PCA_GPCR_subtype)

#read in the file 
Suihi_PCA_GPCR<- read.delim("Suihi_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suihi_PCA_GPCR_fam1<- data.frame(Suihi_PCA_GPCR[ grep("Family", Suihi_PCA_GPCR$X),])
colnames(Suihi_PCA_GPCR_fam1)<- "X"
Suihi_PCA_GPCR_fam2<- as.matrix(Suihi_PCA_GPCR_fam1[ grep("Sub", Suihi_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suihi_PCA_GPCR_fam2)<- "X"
Suihi_PCA_GPCR_fam3<- as.data.frame(Suihi_PCA_GPCR_fam2)
#get sub-family
Suihi_PCA_GPCR_subfam<- as.matrix(Suihi_PCA_GPCR_fam1[ grep("SubFamily", Suihi_PCA_GPCR_fam1$X,),])
Suihi_PCA_GPCR_subfam2<- as.data.frame(Suihi_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suihi_PCA_GPCR_subsubfam<- as.matrix(Suihi_PCA_GPCR_fam1[ grep("Sub-subFamily", Suihi_PCA_GPCR_fam1$X,),])
Suihi_PCA_GPCR_subsubfam2<-as.data.frame(Suihi_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suihi_PCA_GPCR_subtype<- as.matrix(Suihi_PCA_GPCR[ grep("Subtype", Suihi_PCA_GPCR$X,),])
Suihi_PCA_GPCR_subtype<- as.data.frame(Suihi_PCA_GPCR_subtype)
#bin
Suihi_GPCR_table_fam<- table(Suihi_PCA_GPCR_fam3)
Suihi_GPCR_table_subfam<- table(Suihi_PCA_GPCR_subfam2)
Suihi_GPCR_table_subsubfam<- table(Suihi_PCA_GPCR_subsubfam2)
Suihi_GPCR_table_subtype<- table(Suihi_PCA_GPCR_subtype)

#read in the file 
Suigr_PCA_GPCR<- read.delim("Suigr_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suigr_PCA_GPCR_fam1<- data.frame(Suigr_PCA_GPCR[ grep("Family", Suigr_PCA_GPCR$X),])
colnames(Suigr_PCA_GPCR_fam1)<- "X"
Suigr_PCA_GPCR_fam2<- as.matrix(Suigr_PCA_GPCR_fam1[ grep("Sub", Suigr_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suigr_PCA_GPCR_fam2)<- "X"
Suigr_PCA_GPCR_fam3<- as.data.frame(Suigr_PCA_GPCR_fam2)
#get sub-family
Suigr_PCA_GPCR_subfam<- as.matrix(Suigr_PCA_GPCR_fam1[ grep("SubFamily", Suigr_PCA_GPCR_fam1$X,),])
Suigr_PCA_GPCR_subfam2<- as.data.frame(Suigr_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suigr_PCA_GPCR_subsubfam<- as.matrix(Suigr_PCA_GPCR_fam1[ grep("Sub-subFamily", Suigr_PCA_GPCR_fam1$X,),])
Suigr_PCA_GPCR_subsubfam2<-as.data.frame(Suigr_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suigr_PCA_GPCR_subtype<- as.matrix(Suigr_PCA_GPCR[ grep("Subtype", Suigr_PCA_GPCR$X,),])
Suigr_PCA_GPCR_subtype<- as.data.frame(Suigr_PCA_GPCR_subtype)
#bin
Suigr_GPCR_table_fam<- table(Suigr_PCA_GPCR_fam3)
Suigr_GPCR_table_subfam<- table(Suigr_PCA_GPCR_subfam2)
Suigr_GPCR_table_subsubfam<- table(Suigr_PCA_GPCR_subsubfam2)
Suigr_GPCR_table_subtype<- table(Suigr_PCA_GPCR_subtype)

#read in the file 
Suidec_PCA_GPCR<- read.delim("Suidec_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suidec_PCA_GPCR_fam1<- data.frame(Suidec_PCA_GPCR[ grep("Family", Suidec_PCA_GPCR$X),])
colnames(Suidec_PCA_GPCR_fam1)<- "X"
Suidec_PCA_GPCR_fam2<- as.matrix(Suidec_PCA_GPCR_fam1[ grep("Sub", Suidec_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suidec_PCA_GPCR_fam2)<- "X"
Suidec_PCA_GPCR_fam3<- as.data.frame(Suidec_PCA_GPCR_fam2)
#get sub-family
Suidec_PCA_GPCR_subfam<- as.matrix(Suidec_PCA_GPCR_fam1[ grep("SubFamily", Suidec_PCA_GPCR_fam1$X,),])
Suidec_PCA_GPCR_subfam2<- as.data.frame(Suidec_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suidec_PCA_GPCR_subsubfam<- as.matrix(Suidec_PCA_GPCR_fam1[ grep("Sub-subFamily", Suidec_PCA_GPCR_fam1$X,),])
Suidec_PCA_GPCR_subsubfam2<-as.data.frame(Suidec_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suidec_PCA_GPCR_subtype<- as.matrix(Suidec_PCA_GPCR[ grep("Subtype", Suidec_PCA_GPCR$X,),])
Suidec_PCA_GPCR_subtype<- as.data.frame(Suidec_PCA_GPCR_subtype)
#bin
Suidec_GPCR_table_fam<- table(Suidec_PCA_GPCR_fam3)
Suidec_GPCR_table_subfam<- table(Suidec_PCA_GPCR_subfam2)
Suidec_GPCR_table_subsubfam<- table(Suidec_PCA_GPCR_subsubfam2)
Suidec_GPCR_table_subtype<- table(Suidec_PCA_GPCR_subtype)

#read in the file 
Suicot_PCA_GPCR<- read.delim("Suicot_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suicot_PCA_GPCR_fam1<- data.frame(Suicot_PCA_GPCR[ grep("Family", Suicot_PCA_GPCR$X),])
colnames(Suicot_PCA_GPCR_fam1)<- "X"
Suicot_PCA_GPCR_fam2<- as.matrix(Suicot_PCA_GPCR_fam1[ grep("Sub", Suicot_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suicot_PCA_GPCR_fam2)<- "X"
Suicot_PCA_GPCR_fam3<- as.data.frame(Suicot_PCA_GPCR_fam2)
#get sub-family
Suicot_PCA_GPCR_subfam<- as.matrix(Suicot_PCA_GPCR_fam1[ grep("SubFamily", Suicot_PCA_GPCR_fam1$X,),])
Suicot_PCA_GPCR_subfam2<- as.data.frame(Suicot_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suicot_PCA_GPCR_subsubfam<- as.matrix(Suicot_PCA_GPCR_fam1[ grep("Sub-subFamily", Suicot_PCA_GPCR_fam1$X,),])
Suicot_PCA_GPCR_subsubfam2<-as.data.frame(Suicot_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suicot_PCA_GPCR_subtype<- as.matrix(Suicot_PCA_GPCR[ grep("Subtype", Suicot_PCA_GPCR$X,),])
Suicot_PCA_GPCR_subtype<- as.data.frame(Suicot_PCA_GPCR_subtype)
#bin
Suicot_GPCR_table_fam<- table(Suicot_PCA_GPCR_fam3)
Suicot_GPCR_table_subfam<- table(Suicot_PCA_GPCR_subfam2)
Suicot_GPCR_table_subsubfam<- table(Suicot_PCA_GPCR_subsubfam2)
Suicot_GPCR_table_subtype<- table(Suicot_PCA_GPCR_subtype)


#read in the file 
Suicli_PCA_GPCR<- read.delim("Suicli_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suicli_PCA_GPCR_fam1<- data.frame(Suicli_PCA_GPCR[ grep("Family", Suicli_PCA_GPCR$X),])
colnames(Suicli_PCA_GPCR_fam1)<- "X"
Suicli_PCA_GPCR_fam2<- as.matrix(Suicli_PCA_GPCR_fam1[ grep("Sub", Suicli_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suicli_PCA_GPCR_fam2)<- "X"
Suicli_PCA_GPCR_fam3<- as.data.frame(Suicli_PCA_GPCR_fam2)
#get sub-family
Suicli_PCA_GPCR_subfam<- as.matrix(Suicli_PCA_GPCR_fam1[ grep("SubFamily", Suicli_PCA_GPCR_fam1$X,),])
Suicli_PCA_GPCR_subfam2<- as.data.frame(Suicli_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suicli_PCA_GPCR_subsubfam<- as.matrix(Suicli_PCA_GPCR_fam1[ grep("Sub-subFamily", Suicli_PCA_GPCR_fam1$X,),])
Suicli_PCA_GPCR_subsubfam2<-as.data.frame(Suicli_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suicli_PCA_GPCR_subtype<- as.matrix(Suicli_PCA_GPCR[ grep("Subtype", Suicli_PCA_GPCR$X,),])
Suicli_PCA_GPCR_subtype<- as.data.frame(Suicli_PCA_GPCR_subtype)
#bin
Suicli_GPCR_table_fam<- table(Suicli_PCA_GPCR_fam3)
Suicli_GPCR_table_subfam<- table(Suicli_PCA_GPCR_subfam2)
Suicli_GPCR_table_subsubfam<- table(Suicli_PCA_GPCR_subsubfam2)
Suicli_GPCR_table_subtype<- table(Suicli_PCA_GPCR_subtype)


#read in the file 
Suibr_PCA_GPCR<- read.delim("Suibr_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suibr_PCA_GPCR_fam1<- data.frame(Suibr_PCA_GPCR[ grep("Family", Suibr_PCA_GPCR$X),])
colnames(Suibr_PCA_GPCR_fam1)<- "X"
Suibr_PCA_GPCR_fam2<- as.matrix(Suibr_PCA_GPCR_fam1[ grep("Sub", Suibr_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suibr_PCA_GPCR_fam2)<- "X"
Suibr_PCA_GPCR_fam3<- as.data.frame(Suibr_PCA_GPCR_fam2)
#get sub-family
Suibr_PCA_GPCR_subfam<- as.matrix(Suibr_PCA_GPCR_fam1[ grep("SubFamily", Suibr_PCA_GPCR_fam1$X,),])
Suibr_PCA_GPCR_subfam2<- as.data.frame(Suibr_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suibr_PCA_GPCR_subsubfam<- as.matrix(Suibr_PCA_GPCR_fam1[ grep("Sub-subFamily", Suibr_PCA_GPCR_fam1$X,),])
Suibr_PCA_GPCR_subsubfam2<-as.data.frame(Suibr_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suibr_PCA_GPCR_subtype<- as.matrix(Suibr_PCA_GPCR[ grep("Subtype", Suibr_PCA_GPCR$X,),])
Suibr_PCA_GPCR_subtype<- as.data.frame(Suibr_PCA_GPCR_subtype)
#bin
Suibr_GPCR_table_fam<- table(Suibr_PCA_GPCR_fam3)
Suibr_GPCR_table_subfam<- table(Suibr_PCA_GPCR_subfam2)
Suibr_GPCR_table_subsubfam<- table(Suibr_PCA_GPCR_subsubfam2)
Suibr_GPCR_table_subtype<- table(Suibr_PCA_GPCR_subtype)


#read in the file 
Suibov_PCA_GPCR<- read.delim("Suibov_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suibov_PCA_GPCR_fam1<- data.frame(Suibov_PCA_GPCR[ grep("Family", Suibov_PCA_GPCR$X),])
colnames(Suibov_PCA_GPCR_fam1)<- "X"
Suibov_PCA_GPCR_fam2<- as.matrix(Suibov_PCA_GPCR_fam1[ grep("Sub", Suibov_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suibov_PCA_GPCR_fam2)<- "X"
Suibov_PCA_GPCR_fam3<- as.data.frame(Suibov_PCA_GPCR_fam2)
#get sub-family
Suibov_PCA_GPCR_subfam<- as.matrix(Suibov_PCA_GPCR_fam1[ grep("SubFamily", Suibov_PCA_GPCR_fam1$X,),])
Suibov_PCA_GPCR_subfam2<- as.data.frame(Suibov_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suibov_PCA_GPCR_subsubfam<- as.matrix(Suibov_PCA_GPCR_fam1[ grep("Sub-subFamily", Suibov_PCA_GPCR_fam1$X,),])
Suibov_PCA_GPCR_subsubfam2<-as.data.frame(Suibov_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suibov_PCA_GPCR_subtype<- as.matrix(Suibov_PCA_GPCR[ grep("Subtype", Suibov_PCA_GPCR$X,),])
Suibov_PCA_GPCR_subtype<- as.data.frame(Suibov_PCA_GPCR_subtype)
#bin
Suibov_GPCR_table_fam<- table(Suibov_PCA_GPCR_fam3)
Suibov_GPCR_table_subfam<- table(Suibov_PCA_GPCR_subfam2)
Suibov_GPCR_table_subsubfam<- table(Suibov_PCA_GPCR_subsubfam2)
Suibov_GPCR_table_subtype<- table(Suibov_PCA_GPCR_subtype)

#read in the file 
Suiamp_PCA_GPCR<- read.delim("Suiamp_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suiamp_PCA_GPCR_fam1<- data.frame(Suiamp_PCA_GPCR[ grep("Family", Suiamp_PCA_GPCR$X),])
colnames(Suiamp_PCA_GPCR_fam1)<- "X"
Suiamp_PCA_GPCR_fam2<- as.matrix(Suiamp_PCA_GPCR_fam1[ grep("Sub", Suiamp_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suiamp_PCA_GPCR_fam2)<- "X"
Suiamp_PCA_GPCR_fam3<- as.data.frame(Suiamp_PCA_GPCR_fam2)
#get sub-family
Suiamp_PCA_GPCR_subfam<- as.matrix(Suiamp_PCA_GPCR_fam1[ grep("SubFamily", Suiamp_PCA_GPCR_fam1$X,),])
Suiamp_PCA_GPCR_subfam2<- as.data.frame(Suiamp_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suiamp_PCA_GPCR_subsubfam<- as.matrix(Suiamp_PCA_GPCR_fam1[ grep("Sub-subFamily", Suiamp_PCA_GPCR_fam1$X,),])
Suiamp_PCA_GPCR_subsubfam2<-as.data.frame(Suiamp_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suiamp_PCA_GPCR_subtype<- as.matrix(Suiamp_PCA_GPCR[ grep("Subtype", Suiamp_PCA_GPCR$X,),])
Suiamp_PCA_GPCR_subtype<- as.data.frame(Suiamp_PCA_GPCR_subtype)
#bin
Suiamp_GPCR_table_fam<- table(Suiamp_PCA_GPCR_fam3)
Suiamp_GPCR_table_subfam<- table(Suiamp_PCA_GPCR_subfam2)
Suiamp_GPCR_table_subsubfam<- table(Suiamp_PCA_GPCR_subsubfam2)
Suiamp_GPCR_table_subtype<- table(Suiamp_PCA_GPCR_subtype)

#read in the file 
Suiame_PCA_GPCR<- read.delim("Suiame_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Suiame_PCA_GPCR_fam1<- data.frame(Suiame_PCA_GPCR[ grep("Family", Suiame_PCA_GPCR$X),])
colnames(Suiame_PCA_GPCR_fam1)<- "X"
Suiame_PCA_GPCR_fam2<- as.matrix(Suiame_PCA_GPCR_fam1[ grep("Sub", Suiame_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Suiame_PCA_GPCR_fam2)<- "X"
Suiame_PCA_GPCR_fam3<- as.data.frame(Suiame_PCA_GPCR_fam2)
#get sub-family
Suiame_PCA_GPCR_subfam<- as.matrix(Suiame_PCA_GPCR_fam1[ grep("SubFamily", Suiame_PCA_GPCR_fam1$X,),])
Suiame_PCA_GPCR_subfam2<- as.data.frame(Suiame_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Suiame_PCA_GPCR_subsubfam<- as.matrix(Suiame_PCA_GPCR_fam1[ grep("Sub-subFamily", Suiame_PCA_GPCR_fam1$X,),])
Suiame_PCA_GPCR_subsubfam2<-as.data.frame(Suiame_PCA_GPCR_subsubfam)
#get subtype (where avail)
Suiame_PCA_GPCR_subtype<- as.matrix(Suiame_PCA_GPCR[ grep("Subtype", Suiame_PCA_GPCR$X,),])
Suiame_PCA_GPCR_subtype<- as.data.frame(Suiame_PCA_GPCR_subtype)
#bin
Suiame_GPCR_table_fam<- table(Suiame_PCA_GPCR_fam3)
Suiame_GPCR_table_subfam<- table(Suiame_PCA_GPCR_subfam2)
Suiame_GPCR_table_subsubfam<- table(Suiame_PCA_GPCR_subsubfam2)
Suiame_GPCR_table_subtype<- table(Suiame_PCA_GPCR_subtype)


###non-Suillus set
#read in the file 
Rhivul_PCA_GPCR<- read.delim("Rhivul_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Rhivul_PCA_GPCR_fam1<- data.frame(Rhivul_PCA_GPCR[ grep("Family", Rhivul_PCA_GPCR$X),])
colnames(Rhivul_PCA_GPCR_fam1)<- "X"
Rhivul_PCA_GPCR_fam2<- as.matrix(Rhivul_PCA_GPCR_fam1[ grep("Sub", Rhivul_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Rhivul_PCA_GPCR_fam2)<- "X"
Rhivul_PCA_GPCR_fam3<- as.data.frame(Rhivul_PCA_GPCR_fam2)
#get sub-family
Rhivul_PCA_GPCR_subfam<- as.matrix(Rhivul_PCA_GPCR_fam1[ grep("SubFamily", Rhivul_PCA_GPCR_fam1$X,),])
Rhivul_PCA_GPCR_subfam2<- as.data.frame(Rhivul_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Rhivul_PCA_GPCR_subsubfam<- as.matrix(Rhivul_PCA_GPCR_fam1[ grep("Sub-subFamily", Rhivul_PCA_GPCR_fam1$X,),])
Rhivul_PCA_GPCR_subsubfam2<-as.data.frame(Rhivul_PCA_GPCR_subsubfam)
#get subtype (where avail)
Rhivul_PCA_GPCR_subtype<- as.matrix(Rhivul_PCA_GPCR[ grep("Subtype", Rhivul_PCA_GPCR$X,),])
Rhivul_PCA_GPCR_subtype<- as.data.frame(Rhivul_PCA_GPCR_subtype)
#bin
Rhivul_GPCR_table_fam<- table(Rhivul_PCA_GPCR_fam3)
Rhivul_GPCR_table_subfam<- table(Rhivul_PCA_GPCR_subfam2)
Rhivul_GPCR_table_subsubfam<- table(Rhivul_PCA_GPCR_subsubfam2)
Rhivul_GPCR_table_subtype<- table(Rhivul_PCA_GPCR_subtype)

#read in the file 
Rhitru_PCA_GPCR<- read.delim("Rhitru_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Rhitru_PCA_GPCR_fam1<- data.frame(Rhitru_PCA_GPCR[ grep("Family", Rhitru_PCA_GPCR$X),])
colnames(Rhitru_PCA_GPCR_fam1)<- "X"
Rhitru_PCA_GPCR_fam2<- as.matrix(Rhitru_PCA_GPCR_fam1[ grep("Sub", Rhitru_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Rhitru_PCA_GPCR_fam2)<- "X"
Rhitru_PCA_GPCR_fam3<- as.data.frame(Rhitru_PCA_GPCR_fam2)
#get sub-family
Rhitru_PCA_GPCR_subfam<- as.matrix(Rhitru_PCA_GPCR_fam1[ grep("SubFamily", Rhitru_PCA_GPCR_fam1$X,),])
Rhitru_PCA_GPCR_subfam2<- as.data.frame(Rhitru_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Rhitru_PCA_GPCR_subsubfam<- as.matrix(Rhitru_PCA_GPCR_fam1[ grep("Sub-subFamily", Rhitru_PCA_GPCR_fam1$X,),])
Rhitru_PCA_GPCR_subsubfam2<-as.data.frame(Rhitru_PCA_GPCR_subsubfam)
#get subtype (where avail)
Rhitru_PCA_GPCR_subtype<- as.matrix(Rhitru_PCA_GPCR[ grep("Subtype", Rhitru_PCA_GPCR$X,),])
Rhitru_PCA_GPCR_subtype<- as.data.frame(Rhitru_PCA_GPCR_subtype)
#bin
Rhitru_GPCR_table_fam<- table(Rhitru_PCA_GPCR_fam3)
Rhitru_GPCR_table_subfam<- table(Rhitru_PCA_GPCR_subfam2)
Rhitru_GPCR_table_subsubfam<- table(Rhitru_PCA_GPCR_subsubfam2)
Rhitru_GPCR_table_subtype<- table(Rhitru_PCA_GPCR_subtype)

#read in the file 
Amamu_PCA_GPCR<- read.delim("Amamu_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Amamu_PCA_GPCR_fam1<- data.frame(Amamu_PCA_GPCR[ grep("Family", Amamu_PCA_GPCR$X),])
colnames(Amamu_PCA_GPCR_fam1)<- "X"
Amamu_PCA_GPCR_fam2<- as.matrix(Amamu_PCA_GPCR_fam1[ grep("Sub", Amamu_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Amamu_PCA_GPCR_fam2)<- "X"
Amamu_PCA_GPCR_fam3<- as.data.frame(Amamu_PCA_GPCR_fam2)
#get sub-family
Amamu_PCA_GPCR_subfam<- as.matrix(Amamu_PCA_GPCR_fam1[ grep("SubFamily", Amamu_PCA_GPCR_fam1$X,),])
Amamu_PCA_GPCR_subfam2<- as.data.frame(Amamu_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Amamu_PCA_GPCR_subsubfam<- as.matrix(Amamu_PCA_GPCR_fam1[ grep("Sub-subFamily", Amamu_PCA_GPCR_fam1$X,),])
Amamu_PCA_GPCR_subsubfam2<-as.data.frame(Amamu_PCA_GPCR_subsubfam)
#get subtype (where avail)
Amamu_PCA_GPCR_subtype<- as.matrix(Amamu_PCA_GPCR[ grep("Subtype", Amamu_PCA_GPCR$X,),])
Amamu_PCA_GPCR_subtype<- as.data.frame(Amamu_PCA_GPCR_subtype)
#bin
Amamu_GPCR_table_fam<- table(Amamu_PCA_GPCR_fam3)
Amamu_GPCR_table_subfam<- table(Amamu_PCA_GPCR_subfam2)
Amamu_GPCR_table_subsubfam<- table(Amamu_PCA_GPCR_subsubfam2)
Amamu_GPCR_table_subtype<- table(Amamu_PCA_GPCR_subtype)


#read in the file 
Hebcy_PCA_GPCR<- read.delim("Hebcy_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Hebcy_PCA_GPCR_fam1<- data.frame(Hebcy_PCA_GPCR[ grep("Family", Hebcy_PCA_GPCR$X),])
colnames(Hebcy_PCA_GPCR_fam1)<- "X"
Hebcy_PCA_GPCR_fam2<- as.matrix(Hebcy_PCA_GPCR_fam1[ grep("Sub", Hebcy_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Hebcy_PCA_GPCR_fam2)<- "X"
Hebcy_PCA_GPCR_fam3<- as.data.frame(Hebcy_PCA_GPCR_fam2)
#get sub-family
Hebcy_PCA_GPCR_subfam<- as.matrix(Hebcy_PCA_GPCR_fam1[ grep("SubFamily", Hebcy_PCA_GPCR_fam1$X,),])
Hebcy_PCA_GPCR_subfam2<- as.data.frame(Hebcy_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Hebcy_PCA_GPCR_subsubfam<- as.matrix(Hebcy_PCA_GPCR_fam1[ grep("Sub-subFamily", Hebcy_PCA_GPCR_fam1$X,),])
Hebcy_PCA_GPCR_subsubfam2<-as.data.frame(Hebcy_PCA_GPCR_subsubfam)
#get subtype (where avail)
Hebcy_PCA_GPCR_subtype<- as.matrix(Hebcy_PCA_GPCR[ grep("Subtype", Hebcy_PCA_GPCR$X,),])
Hebcy_PCA_GPCR_subtype<- as.data.frame(Hebcy_PCA_GPCR_subtype)
#bin
Hebcy_GPCR_table_fam<- table(Hebcy_PCA_GPCR_fam3)
Hebcy_GPCR_table_subfam<- table(Hebcy_PCA_GPCR_subfam2)
Hebcy_GPCR_table_subsubfam<- table(Hebcy_PCA_GPCR_subsubfam2)
Hebcy_GPCR_table_subtype<- table(Hebcy_PCA_GPCR_subtype)

#read in the file 
Lacbi_PCA_GPCR<- read.delim("Lacbi_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Lacbi_PCA_GPCR_fam1<- data.frame(Lacbi_PCA_GPCR[ grep("Family", Lacbi_PCA_GPCR$X),])
colnames(Lacbi_PCA_GPCR_fam1)<- "X"
Lacbi_PCA_GPCR_fam2<- as.matrix(Lacbi_PCA_GPCR_fam1[ grep("Sub", Lacbi_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Lacbi_PCA_GPCR_fam2)<- "X"
Lacbi_PCA_GPCR_fam3<- as.data.frame(Lacbi_PCA_GPCR_fam2)
#get sub-family
Lacbi_PCA_GPCR_subfam<- as.matrix(Lacbi_PCA_GPCR_fam1[ grep("SubFamily", Lacbi_PCA_GPCR_fam1$X,),])
Lacbi_PCA_GPCR_subfam2<- as.data.frame(Lacbi_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Lacbi_PCA_GPCR_subsubfam<- as.matrix(Lacbi_PCA_GPCR_fam1[ grep("Sub-subFamily", Lacbi_PCA_GPCR_fam1$X,),])
Lacbi_PCA_GPCR_subsubfam2<-as.data.frame(Lacbi_PCA_GPCR_subsubfam)
#get subtype (where avail)
Lacbi_PCA_GPCR_subtype<- as.matrix(Lacbi_PCA_GPCR[ grep("Subtype", Lacbi_PCA_GPCR$X,),])
Lacbi_PCA_GPCR_subtype<- as.data.frame(Lacbi_PCA_GPCR_subtype)
#bin
Lacbi_GPCR_table_fam<- table(Lacbi_PCA_GPCR_fam3)
Lacbi_GPCR_table_subfam<- table(Lacbi_PCA_GPCR_subfam2)
Lacbi_GPCR_table_subsubfam<- table(Lacbi_PCA_GPCR_subsubfam2)
Lacbi_GPCR_table_subtype<- table(Lacbi_PCA_GPCR_subtype)

#read in the file 
Paxin_PCA_GPCR<- read.delim("Paxin_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Paxin_PCA_GPCR_fam1<- data.frame(Paxin_PCA_GPCR[ grep("Family", Paxin_PCA_GPCR$X),])
colnames(Paxin_PCA_GPCR_fam1)<- "X"
Paxin_PCA_GPCR_fam2<- as.matrix(Paxin_PCA_GPCR_fam1[ grep("Sub", Paxin_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Paxin_PCA_GPCR_fam2)<- "X"
Paxin_PCA_GPCR_fam3<- as.data.frame(Paxin_PCA_GPCR_fam2)
#get sub-family
Paxin_PCA_GPCR_subfam<- as.matrix(Paxin_PCA_GPCR_fam1[ grep("SubFamily", Paxin_PCA_GPCR_fam1$X,),])
Paxin_PCA_GPCR_subfam2<- as.data.frame(Paxin_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Paxin_PCA_GPCR_subsubfam<- as.matrix(Paxin_PCA_GPCR_fam1[ grep("Sub-subFamily", Paxin_PCA_GPCR_fam1$X,),])
Paxin_PCA_GPCR_subsubfam2<-as.data.frame(Paxin_PCA_GPCR_subsubfam)
#get subtype (where avail)
Paxin_PCA_GPCR_subtype<- as.matrix(Paxin_PCA_GPCR[ grep("Subtype", Paxin_PCA_GPCR$X,),])
Paxin_PCA_GPCR_subtype<- as.data.frame(Paxin_PCA_GPCR_subtype)
#bin
Paxin_GPCR_table_fam<- table(Paxin_PCA_GPCR_fam3)
Paxin_GPCR_table_subfam<- table(Paxin_PCA_GPCR_subfam2)
Paxin_GPCR_table_subsubfam<- table(Paxin_PCA_GPCR_subsubfam2)
Paxin_GPCR_table_subtype<- table(Paxin_PCA_GPCR_subtype)

#read in the file 
Pilcr_PCA_GPCR<- read.delim("Pilcr_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Pilcr_PCA_GPCR_fam1<- data.frame(Pilcr_PCA_GPCR[ grep("Family", Pilcr_PCA_GPCR$X),])
colnames(Pilcr_PCA_GPCR_fam1)<- "X"
Pilcr_PCA_GPCR_fam2<- as.matrix(Pilcr_PCA_GPCR_fam1[ grep("Sub", Pilcr_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Pilcr_PCA_GPCR_fam2)<- "X"
Pilcr_PCA_GPCR_fam3<- as.data.frame(Pilcr_PCA_GPCR_fam2)
#get sub-family
Pilcr_PCA_GPCR_subfam<- as.matrix(Pilcr_PCA_GPCR_fam1[ grep("SubFamily", Pilcr_PCA_GPCR_fam1$X,),])
Pilcr_PCA_GPCR_subfam2<- as.data.frame(Pilcr_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Pilcr_PCA_GPCR_subsubfam<- as.matrix(Pilcr_PCA_GPCR_fam1[ grep("Sub-subFamily", Pilcr_PCA_GPCR_fam1$X,),])
Pilcr_PCA_GPCR_subsubfam2<-as.data.frame(Pilcr_PCA_GPCR_subsubfam)
#get subtype (where avail)
Pilcr_PCA_GPCR_subtype<- as.matrix(Pilcr_PCA_GPCR[ grep("Subtype", Pilcr_PCA_GPCR$X,),])
Pilcr_PCA_GPCR_subtype<- as.data.frame(Pilcr_PCA_GPCR_subtype)
#bin
Pilcr_GPCR_table_fam<- table(Pilcr_PCA_GPCR_fam3)
Pilcr_GPCR_table_subfam<- table(Pilcr_PCA_GPCR_subfam2)
Pilcr_GPCR_table_subsubfam<- table(Pilcr_PCA_GPCR_subsubfam2)
Pilcr_GPCR_table_subtype<- table(Pilcr_PCA_GPCR_subtype)


#read in the file 
Pismi_PCA_GPCR<- read.delim("Pismi_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Pismi_PCA_GPCR_fam1<- data.frame(Pismi_PCA_GPCR[ grep("Family", Pismi_PCA_GPCR$X),])
colnames(Pismi_PCA_GPCR_fam1)<- "X"
Pismi_PCA_GPCR_fam2<- as.matrix(Pismi_PCA_GPCR_fam1[ grep("Sub", Pismi_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Pismi_PCA_GPCR_fam2)<- "X"
Pismi_PCA_GPCR_fam3<- as.data.frame(Pismi_PCA_GPCR_fam2)
#get sub-family
Pismi_PCA_GPCR_subfam<- as.matrix(Pismi_PCA_GPCR_fam1[ grep("SubFamily", Pismi_PCA_GPCR_fam1$X,),])
Pismi_PCA_GPCR_subfam2<- as.data.frame(Pismi_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Pismi_PCA_GPCR_subsubfam<- as.matrix(Pismi_PCA_GPCR_fam1[ grep("Sub-subFamily", Pismi_PCA_GPCR_fam1$X,),])
Pismi_PCA_GPCR_subsubfam2<-as.data.frame(Pismi_PCA_GPCR_subsubfam)
#get subtype (where avail)
Pismi_PCA_GPCR_subtype<- as.matrix(Pismi_PCA_GPCR[ grep("Subtype", Pismi_PCA_GPCR$X,),])
Pismi_PCA_GPCR_subtype<- as.data.frame(Pismi_PCA_GPCR_subtype)
#bin
Pismi_GPCR_table_fam<- table(Pismi_PCA_GPCR_fam3)
Pismi_GPCR_table_subfam<- table(Pismi_PCA_GPCR_subfam2)
Pismi_GPCR_table_subsubfam<- table(Pismi_PCA_GPCR_subsubfam2)
Pismi_GPCR_table_subtype<- table(Pismi_PCA_GPCR_subtype)


#read in the file 
Sclci_PCA_GPCR<- read.delim("Sclci_PCA_GPCR.txt", col.names = "X", stringsAsFactors = FALSE)
#get family 
Sclci_PCA_GPCR_fam1<- data.frame(Sclci_PCA_GPCR[ grep("Family", Sclci_PCA_GPCR$X),])
colnames(Sclci_PCA_GPCR_fam1)<- "X"
Sclci_PCA_GPCR_fam2<- as.matrix(Sclci_PCA_GPCR_fam1[ grep("Sub", Sclci_PCA_GPCR_fam1$X, invert = TRUE),])
colnames(Sclci_PCA_GPCR_fam2)<- "X"
Sclci_PCA_GPCR_fam3<- as.data.frame(Sclci_PCA_GPCR_fam2)
#get sub-family
Sclci_PCA_GPCR_subfam<- as.matrix(Sclci_PCA_GPCR_fam1[ grep("SubFamily", Sclci_PCA_GPCR_fam1$X,),])
Sclci_PCA_GPCR_subfam2<- as.data.frame(Sclci_PCA_GPCR_subfam)
#get sub- sub family (where avail)
Sclci_PCA_GPCR_subsubfam<- as.matrix(Sclci_PCA_GPCR_fam1[ grep("Sub-subFamily", Sclci_PCA_GPCR_fam1$X,),])
Sclci_PCA_GPCR_subsubfam2<-as.data.frame(Sclci_PCA_GPCR_subsubfam)
#get subtype (where avail)
Sclci_PCA_GPCR_subtype<- as.matrix(Sclci_PCA_GPCR[ grep("Subtype", Sclci_PCA_GPCR$X,),])
Sclci_PCA_GPCR_subtype<- as.data.frame(Sclci_PCA_GPCR_subtype)
#bin
Sclci_GPCR_table_fam<- table(Sclci_PCA_GPCR_fam3)
Sclci_GPCR_table_subfam<- table(Sclci_PCA_GPCR_subfam2)
Sclci_GPCR_table_subsubfam<- table(Sclci_PCA_GPCR_subsubfam2)
Sclci_GPCR_table_subtype<- table(Sclci_PCA_GPCR_subtype)


####reformat data
#function to reformat GPCRs
reformat_GPCRs<- function(x) {
  t<- data.frame(x)
  t2<- t(t)
  colnames(t2)<-t2[1,]
  t3<- data.frame(t2)
  t4<- t3[2,]
}

options(stringsAsFactors = FALSE)

#run on each
Suivar1_GPCRs_rf<- reformat_GPCRs(Suivar_GPCR_table_fam)
Suitom1_GPCRs_rf<- reformat_GPCRs(Suitom_GPCR_table_fam)
Suisub1_GPCRs_rf<- reformat_GPCRs(Suisub_GPCR_table_fam)
Suisu1_GPCRs_rf<- reformat_GPCRs(Suisu_GPCR_table_fam)
Suipla1_GPCRs_rf<- reformat_GPCRs(Suipla_GPCR_table_fam)
Suipic1_GPCRs_rf<- reformat_GPCRs(Suipic_GPCR_table_fam)
Suipal1_GPCRs_rf<- reformat_GPCRs(Suipal_GPCR_table_fam)
Suiocc1_GPCRs_rf<- reformat_GPCRs(Suiocc_GPCR_table_fam)
Suilu4_GPCRs_rf<- reformat_GPCRs(Suilu_GPCR_table_fam)
Suilak1_GPCRs_rf<- reformat_GPCRs(Suilak_GPCR_table_fam)
Suihi1_GPCRs_rf<- reformat_GPCRs(Suihi_GPCR_table_fam)
Suigr1_GPCRs_rf<- reformat_GPCRs(Suigr_GPCR_table_fam)
Suidec1_GPCRs_rf<- reformat_GPCRs(Suidec_GPCR_table_fam)
Suicot1_GPCRs_rf<- reformat_GPCRs(Suicot_GPCR_table_fam)
Suicli1_GPCRs_rf<- reformat_GPCRs(Suicli_GPCR_table_fam)
Suibr2_GPCRs_rf<- reformat_GPCRs(Suibr_GPCR_table_fam)
Suibov1_GPCRs_rf<- reformat_GPCRs(Suibov_GPCR_table_fam)
Suiamp1_GPCRs_rf<- reformat_GPCRs(Suiamp_GPCR_table_fam)
Suiame1_GPCRs_rf<- reformat_GPCRs(Suiame_GPCR_table_fam)
#non-Suillus set
Rhivul1_GPCRs_rf<- reformat_GPCRs(Rhivul_GPCR_table_fam)
Rhitru1_GPCRs_rf<- reformat_GPCRs(Rhitru_GPCR_table_fam)
Amamu1_GPCRs_rf<- reformat_GPCRs(Amamu_GPCR_table_fam)
Hebcy2_GPCRs_rf<- reformat_GPCRs(Hebcy_GPCR_table_fam)
Lacbi2_GPCRs_rf<- reformat_GPCRs(Lacbi_GPCR_table_fam)
Paxin1_GPCRs_rf<- reformat_GPCRs(Paxin_GPCR_table_fam)
Pilcr1_GPCRs_rf<- reformat_GPCRs(Pilcr_GPCR_table_fam)
Pismi1_GPCRs_rf<- reformat_GPCRs(Pismi_GPCR_table_fam)
Sclci1_GPCRs_rf<- reformat_GPCRs(Sclci_GPCR_table_fam)


GPCR_class_counts_Suillus<- plyr::rbind.fill(Suivar1_GPCRs_rf,
                                             Suitom1_GPCRs_rf,
                                             Suisub1_GPCRs_rf,
                                             Suisu1_GPCRs_rf,
                                             Suipla1_GPCRs_rf,
                                             Suipic1_GPCRs_rf,
                                             Suipal1_GPCRs_rf,
                                             Suiocc1_GPCRs_rf,
                                             Suilu4_GPCRs_rf,
                                             Suilak1_GPCRs_rf,
                                             Suihi1_GPCRs_rf,
                                             Suigr1_GPCRs_rf,
                                             Suidec1_GPCRs_rf,
                                             Suicot1_GPCRs_rf,
                                             Suicli1_GPCRs_rf,
                                             Suibr2_GPCRs_rf,
                                             Suibov1_GPCRs_rf,
                                             Suiamp1_GPCRs_rf,
                                             Suiame1_GPCRs_rf)

row.names(GPCR_class_counts_Suillus)<- c("Suivar1", 
                                         "Suitom1",
                                         "Suisub1",
                                         "Suisu1",
                                         "Suipla1",
                                         "Suipic1",
                                         "Suipal1",
                                         "Suiocc1",
                                         "Suilu4",
                                         "Suilak1",
                                         "Suihi1",
                                         "Suigr1",
                                         "Suidec1",
                                         "Suicot1",
                                         "Suicli1",
                                         "Suibr2",
                                         "Suibov1",
                                         "Suiamp1",
                                         "Suiame1")




GPCR_class_counts_Other<- plyr::rbind.fill(Rhivul1_GPCRs_rf,
                                           Rhitru1_GPCRs_rf,
                                           Amamu1_GPCRs_rf,
                                           Hebcy2_GPCRs_rf,
                                           Lacbi2_GPCRs_rf,
                                           Paxin1_GPCRs_rf,
                                           Pilcr1_GPCRs_rf,
                                           Pismi1_GPCRs_rf,
                                           Sclci1_GPCRs_rf)


row.names(GPCR_class_counts_Other)<- c("Rhivul1",
                                       "Rhitru1",
                                       "Amamu1",
                                       "Hebcy2",
                                       "Lacbi2",
                                       "Paxin1",
                                       "Pilcr1",
                                       "Pismi1",
                                       "Sclci1")


#replace NAs with zeros

GPCR_class_counts_Suillus[is.na(GPCR_class_counts_Suillus)] <- 0
GPCR_class_counts_Other[is.na(GPCR_class_counts_Other)] <- 0



#get averages 
Suillus_mean_Class_A<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Class.A.Rhodopsin.like.)), 0)
Suillus_mean_Class_B<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Class.B.Secretin.like.)), 0)
Suillus_mean_Class_D<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Vomeronasal.receptors..V1R...V3R..)), 0)

Other_mean_Class_A<- round(mean(as.numeric(GPCR_class_counts_Other$Family.Class.A.Rhodopsin.like.)), 0)
Other_mean_Class_B<- round(mean(as.numeric(GPCR_class_counts_Other$Family.Class.B.Secretin.like.)), 0)
Other_mean_Class_D<- round(mean(as.numeric(GPCR_class_counts_Other$Family.Vomeronasal.receptors..V1R...V3R..)), 0)


#####make waffle plot 
##Suillus waffle
GPCR_dist_Suillus <- c(`Class A`=Suillus_mean_Class_A, `Class B`=Suillus_mean_Class_B, 
                       `Class D`=Suillus_mean_Class_D)
waffle(GPCR_dist_Suillus/1, rows=5, size=0.5, 
       colors=c("#405952", "#5677A1", "#9B9A79"), 
       title="GPCR distrbution Suillus")

##Other waffle
GPCR_dist_Other <- c(`Class A`=Other_mean_Class_A, `Class B`=Other_mean_Class_B, 
                     `Class D`=Other_mean_Class_D)
waffle(GPCR_dist_Other/1, rows=5, size=0.5, 
       colors=c("#405952", "#5677A1", "#9B9A79"), 
       title="GPCR distrbution Other", 
       pad = 4)



#####make plots to show dist. of data
#change to numeric
GPCR_class_counts_Suillus[] <- lapply(GPCR_class_counts_Suillus, function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
sapply(GPCR_class_counts_Suillus, class)

GPCR_class_counts_Other[] <- lapply(GPCR_class_counts_Other, function(x) {
  if(is.character(x)) as.numeric(as.character(x)) else x
})
sapply(GPCR_class_counts_Other, class)


#change to longform 
#first add repeat col (ours is only 1)
GPCR_class_counts_Suillus$id<- seq(1:nrow(GPCR_class_counts_Suillus))
GPCR_class_counts_Other$id<- seq(1:nrow(GPCR_class_counts_Other))
#make longform
GPCR_class_counts_Suillus_long <- GPCR_class_counts_Suillus %>% group_by(id) %>%
  gather(data = GPCR_class_counts_Suillus, id, Family.Class.A.Rhodopsin.like., Family.Class.B.Secretin.like., Family.Vomeronasal.receptors..V1R...V3R..)
GPCR_class_counts_Other_long <- GPCR_class_counts_Other %>% group_by(id) %>%
  gather(data = GPCR_class_counts_Other, id, Family.Class.A.Rhodopsin.like., Family.Class.B.Secretin.like., Family.Vomeronasal.receptors..V1R...V3R..)


#set col names
colnames(GPCR_class_counts_Suillus_long)<- c("treatment", "value")
colnames(GPCR_class_counts_Other_long)<- c("treatment", "value")

# make violin plot
#Suillus
with(GPCR_class_counts_Suillus_long, 
     vioplot( value[treatment=="Family.Class.A.Rhodopsin.like."], 
              value[treatment=="Family.Class.B.Secretin.like."], 
              value[treatment=="Family.Vomeronasal.receptors..V1R...V3R.."],  
              col=c("#405952", "#5677A1", "#9B9A79"), 
              names=c("Class A","Class B","Class D"),
              main = "Suillus"))


#Other
with(GPCR_class_counts_Other_long, 
     vioplot( value[treatment=="Family.Class.A.Rhodopsin.like."], 
              value[treatment=="Family.Class.B.Secretin.like."], 
              value[treatment=="Family.Vomeronasal.receptors..V1R...V3R.."],  
              col=c("#405952", "#5677A1", "#9B9A79"), 
              names=c("Class A","Class B","Class D"),
              main = "Other",
              ylim = c(0,70)))




#####stats and figures for Suillus by Host association
#add host association
#write.csv(x = GPCR_class_counts_Suillus, file = "Suillus_by_host_association_GPCRs.csv")
Suillus_by_host_association_GPCRs<-read.csv("Suillus_by_host_association_GPCRs.csv")

#reformat to long
#make longform
Suillus_by_host_association_long <- Suillus_by_host_association_GPCRs %>% group_by(id) %>%
  gather(data = Suillus_by_host_association_GPCRs, id, classA, classB, classD)



#set col names
colnames(Suillus_by_host_association_long)<- c("species", "host", "class", "counts")

#parse by host
Red_hosts_gpcrs<- Suillus_by_host_association_long[Suillus_by_host_association_long$host == "R",]
White_hosts_gpcrs<- Suillus_by_host_association_long[Suillus_by_host_association_long$host == "W",]
Larch_hosts_gpcrs<- Suillus_by_host_association_long[Suillus_by_host_association_long$host == "L",]

Suillus_redmean_Class_A1<-Red_hosts_gpcrs[Red_hosts_gpcrs$class == "classA",]
Suillus_redmean_Class_A<- round(mean(as.numeric(Suillus_redmean_Class_A1$counts)), 0)
Suillus_redmean_Class_B1<-Red_hosts_gpcrs[Red_hosts_gpcrs$class == "classB",]
Suillus_redmean_Class_B<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Class.B.Secretin.like.)), 0)
Suillus_redmean_Class_D1<-Red_hosts_gpcrs[Red_hosts_gpcrs$class == "classD",]
Suillus_redmean_Class_D<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Vomeronasal.receptors..V1R...V3R..)), 0)


Suillus_host_dist_red <- c(`Class A`=Suillus_redmean_Class_A, `Class B`=Suillus_redmean_Class_B, 
                           `Class D`=Suillus_redmean_Class_D)

waffle(Suillus_host_dist_red/1, rows=5, size=0.5, 
       colors=c("#405952", "#5677A1", "#9B9A79"), 
       title="GPCR distrbution Red", 
       pad = 4)


Suillus_Whitemean_Class_A1<-White_hosts_gpcrs[White_hosts_gpcrs$class == "classA",]
Suillus_Whitemean_Class_A<- round(mean(as.numeric(Suillus_Whitemean_Class_A1$counts)), 0)
Suillus_Whitemean_Class_B1<-White_hosts_gpcrs[White_hosts_gpcrs$class == "classB",]
Suillus_Whitemean_Class_B<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Class.B.Secretin.like.)), 0)
Suillus_Whitemean_Class_D1<-White_hosts_gpcrs[White_hosts_gpcrs$class == "classD",]
Suillus_Whitemean_Class_D<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Vomeronasal.receptors..V1R...V3R..)), 0)


Suillus_host_dist_White <- c(`Class A`=Suillus_Whitemean_Class_A, `Class B`=Suillus_Whitemean_Class_B, 
                             `Class D`=Suillus_Whitemean_Class_D)

waffle(Suillus_host_dist_White/1, rows=5, size=0.5, 
       colors=c("#405952", "#5677A1", "#9B9A79"), 
       title="GPCR distrbution White", 
       pad = 4)

Suillus_Larchmean_Class_A1<-Larch_hosts_gpcrs[Larch_hosts_gpcrs$class == "classA",]
Suillus_Larchmean_Class_A<- round(mean(as.numeric(Suillus_Larchmean_Class_A1$counts)), 0)
Suillus_Larchmean_Class_B1<-Larch_hosts_gpcrs[Larch_hosts_gpcrs$class == "classB",]
Suillus_Larchmean_Class_B<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Class.B.Secretin.like.)), 0)
Suillus_Larchmean_Class_D1<-Larch_hosts_gpcrs[Larch_hosts_gpcrs$class == "classD",]
Suillus_Larchmean_Class_D<- round(mean(as.numeric(GPCR_class_counts_Suillus$Family.Vomeronasal.receptors..V1R...V3R..)), 0)



Suillus_host_dist_Larch <- c(`Class A`=Suillus_Larchmean_Class_A, `Class B`=Suillus_Larchmean_Class_B, 
                             `Class D`=Suillus_Larchmean_Class_D)

waffle(Suillus_host_dist_Larch/1, rows=5, size=0.5, 
       colors=c("#405952", "#5677A1", "#9B9A79"), 
       title="GPCR distrbution Larch", 
       pad = 4)

#strip plot by GPCR type by host
#make a list with each host group in it
x <- list("Red"=Red_hosts_gpcrs$counts, "White"=White_hosts_gpcrs$counts, "Larch" = Larch_hosts_gpcrs$counts)


#by class
xclassA <- list("Red"= Suillus_redmean_Class_A1$counts, "White"= Suillus_Whitemean_Class_A1$counts, "Larch" =  Suillus_Larchmean_Class_A1$counts)
boxplot(xclassA, 
        border="#405952",
        outline = FALSE, boxlty = 0,
        whisklty = 0, staplelty = 0,
        ylim = c(0,80), 
        main="Total GPCRs by host") 
stripchart(xclassA,
           ylab="n GPCRs",
           method="jitter",
           col="#405952",
           vertical=TRUE,
           pch=19,
           add = TRUE
)

xclassB <- list("Red"= Suillus_redmean_Class_B1$counts, "White"= Suillus_Whitemean_Class_B1$counts, "Larch" =  Suillus_Larchmean_Class_B1$counts)
boxplot(xclassB, 
        border="#5677A1",
        outline = FALSE, boxlty = 0,
        whisklty = 0, staplelty = 0,
        add = TRUE) 

stripchart(xclassB,
           main="Total GPCRs by host",
           ylab="n GPCRs",
           method="jitter",
           col="#5677A1",
           vertical=TRUE,
           pch=19, 
           add = TRUE
)

xclassD <- list("Red"= Suillus_redmean_Class_D1$counts, "White"= Suillus_Whitemean_Class_D1$counts, "Larch" =  Suillus_Larchmean_Class_D1$counts)
boxplot(xclassD, 
        border="#9B9A79",
        outline = FALSE, boxlty = 0,
        whisklty = 0, staplelty = 0,
        add = TRUE) 
stripchart(xclassD,
           main="Total GPCRs by host",
           ylab="n GPCRs",
           method="jitter",
           col="#9B9A79",
           vertical=TRUE,
           pch=19, 
           add = TRUE
)

###make bar chart to represent totals in ea cat.
#add col to get GPCR sums
Suillus_by_host_association_GPCRs$GPCR_total<- rowSums(Suillus_by_host_association_GPCRs[, 2:4])

#means
red_set<- Suillus_by_host_association_GPCRs[Suillus_by_host_association_GPCRs$host == "R",]
white_set<- Suillus_by_host_association_GPCRs[Suillus_by_host_association_GPCRs$host == "W",]
larch_set<- Suillus_by_host_association_GPCRs[Suillus_by_host_association_GPCRs$host == "L",]

red_mean<-round(mean(red_set$GPCR_total),0)
white_mean<-round(mean(white_set$GPCR_total),0)
larch_mean<-round(mean(larch_set$GPCR_total),0)

GPCR_totals<- data.frame(rbind(red_mean, white_mean, larch_mean))
colnames(GPCR_totals)<- "counts"
GPCR_totals$host_sp<- row.names(GPCR_totals)

#make it into a table
GPCR_totals_long <- GPCR_totals[rep(seq(nrow(GPCR_totals)), GPCR_totals$counts),]
GPCR_table<- table(GPCR_totals_long)


#lets make this a stacked barplot by subtype 
red_mean<- data.frame(round(colMeans(red_set[2:4]), 0))
colnames(red_mean)<- "Red"
white_mean<- data.frame(round(colMeans(white_set[2:4]), 0))
colnames(white_mean)<- "White"
larch_mean<- data.frame(round(colMeans(larch_set[2:4]), 0))
colnames(larch_mean)<- "Larch"

stacked_set<-cbind(red_mean, white_mean, larch_mean)
stacked_set<-as.matrix(stacked_set)

barplot(stacked_set, main = "mean GPCRs by host", 
        col=c("#405952", "#5677A1", "#9B9A79"), 
        border="white", 
        space=0.04, 
        font.axis=2, 
        ylab = "n GPCRs", 
        ylim = c(0,80))








###### stats ######
#####ANOVA 

#type 2 ANOVA for Suillus vs. Other ECM 

#format input files
GPCR_class_counts_Other_long$Group<- rep("O", nrow(GPCR_class_counts_Other_long))
GPCR_class_counts_Suillus_long$Group<- rep("S", nrow(GPCR_class_counts_Suillus_long))

#bind them
class_counts_add_df<- rbind(GPCR_class_counts_Other_long, GPCR_class_counts_Suillus_long)

#make model
model.for.t2<-lm(class_counts_add_df$value ~ class_counts_add_df$treatment * class_counts_add_df$Group)

model.for.t3<- Anova(model.for.t2, type=2)
model.for.t3
summary(model.for.t3)
#significant 

#define interaction term
class_counts_add_df$Group_GPCR_interaction <- as.factor(interaction(class_counts_add_df$treatment, class_counts_add_df$Group))

#re-run model with interaction term
model.for.t4<- lm(class_counts_add_df$value ~ class_counts_add_df$Group_GPCR_interaction)

#take a look - lines in common mean not significantly different 
sidelines(pairwise(model.for.t4, class_counts_add_df$Group_GPCR_interaction,confidence = 0.95, type = "hsd"))

#full result with hsd -- note you can't get p-vals with hsd though with cfcdae. 
pairwise(model.for.t4, class_counts_add_df$Group_GPCR_interaction,confidence = 0.95, type = "hsd")

#run model for multiple t tests so that you can get p-vals
pairs_out<- pairwise(model.for.t4, class_counts_add_df$Group_GPCR_interaction,confidence = 0.95, type = "regwr")

sidelines(pairs_out)
#terpenes and Other are dignificantly different 

#get p-value using mulitple t tests and holm adjustment for multiple comparisons. 
p_vals<-with(class_counts_add_df, pairwise.t.test(class_counts_add_df$value,class_counts_add_df$Group_GPCR_interaction,
                                                               p.adjust.method="holm"))

p_vals





#type 2 ANOVA for Suillus by host association 
#subset df for anova
Suillus_by_host_association_long2<- Suillus_by_host_association_long[c(Suillus_by_host_association_long$host == "W" | Suillus_by_host_association_long$host == "R" | Suillus_by_host_association_long$host == "L"),]

#get means for each (for graph)
red_mean
sum(red_mean)

white_mean
sum(white_mean)

larch_mean
sum(larch_mean)

library("car")
#make model
model.for.t2<-lm(Suillus_by_host_association_long2$counts ~ Suillus_by_host_association_long2$class * Suillus_by_host_association_long2$host)

model.for.t3<- Anova(model.for.t2, type=2)

summary(model.for.t3)
#significant 

#run post hoc tests (multiple t test)
library("cfcdae")

#define interaction term
Suillus_by_host_association_long2$Group_GPCR_interaction <- as.factor(interaction(Suillus_by_host_association_long2$class, Suillus_by_host_association_long2$host))

#re-run model with interaction term
model.for.t4<- lm(Suillus_by_host_association_long2$counts ~ Suillus_by_host_association_long2$Group_GPCR_interaction)

#take a look - lines in common mean not significantly different 
sidelines(pairwise(model.for.t4, Suillus_by_host_association_long2$Group_GPCR_interaction,confidence = 0.95, type = "hsd"))

#full result with hsd -- note you can't get p-vals with hsd though with cfcdae. 
pairwise(model.for.t4, Suillus_by_host_association_long2$Group_GPCR_interaction,confidence = 0.95, type = "hsd")

#run model for multiple t tests so that you can get p-vals
pairs_out<- pairwise(model.for.t4, Suillus_by_host_association_long2$Group_GPCR_interaction,confidence = 0.95, type = "regwr")

sidelines(pairs_out)
#terpenes and Other are dignificantly different 

#get p-value using mulitple t tests and holm adjustment for multiple comparisons. 
p_vals<-with(Suillus_by_host_association_long2, pairwise.t.test(Suillus_by_host_association_long2$counts,Suillus_by_host_association_long2$Group_GPCR_interaction,
                                                         p.adjust.method="holm"))

