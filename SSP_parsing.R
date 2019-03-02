#parse SSP's
#load libraries
library("seqinr")
library("data.table")
#data was generated using SignalP on the command line, followed by tmhmm on the "mature sequence" output form Signal P, also on the command line

setwd("~/Desktop/Project_Suillus_comp_genomics/R")
#read in the input files
Suivar1_TMHMM<- data.frame(read.csv("Suivar1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suitom1_TMHMM<- data.frame(read.csv("Suitom1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suisub1_TMHMM<- data.frame(read.csv("Suisub1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suisu1_TMHMM<- data.frame(read.csv("Suisu1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipla1_TMHMM<- data.frame(read.csv("Suipla1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipic1_TMHMM<- data.frame(read.csv("Suipic1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suipal1_TMHMM<- data.frame(read.csv("Suipal1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiocc1_TMHMM<- data.frame(read.csv("Suiocc1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suilu4_TMHMM<- data.frame(read.csv("Suilu4_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suilak1_TMHMM<- data.frame(read.csv("Suilak1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suihi1_TMHMM<- data.frame(read.csv("Suihi1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suigr1_TMHMM<- data.frame(read.csv("Suigr1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suidec1_TMHMM<- data.frame(read.csv("Suidec1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suicot1_TMHMM<- data.frame(read.csv("Suicot1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suicli1_TMHMM<- data.frame(read.csv("Suicli1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suibr2_TMHMM<- data.frame(read.csv("Suibr2_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suibov1_TMHMM<- data.frame(read.csv("Suibov1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiamp1_TMHMM<- data.frame(read.csv("Suiamp1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))
Suiame1_TMHMM<- data.frame(read.csv("Suiame1_TMHMM_SSP_output.txt", header = FALSE, col.names = "header"))

View(Suivar1_TMHMM)

#isolate relevent lines 
#get only lines starting with a hash
Suivar1_TMHMM_sm<- data.frame(Suivar1_TMHMM[ grep("# jgi", Suivar1_TMHMM$header),])
Suitom1_TMHMM_sm<- data.frame(Suitom1_TMHMM[ grep("# jgi", Suitom1_TMHMM$header),])
Suisub1_TMHMM_sm<- data.frame(Suisub1_TMHMM[ grep("# jgi", Suisub1_TMHMM$header),])
Suisu1_TMHMM_sm<- data.frame(Suisu1_TMHMM[ grep("# jgi", Suisu1_TMHMM$header),])
Suipla1_TMHMM_sm<- data.frame(Suipla1_TMHMM[ grep("# jgi", Suipla1_TMHMM$header),])
Suipic1_TMHMM_sm<- data.frame(Suipic1_TMHMM[ grep("# jgi", Suipic1_TMHMM$header),])
Suipal1_TMHMM_sm<- data.frame(Suipal1_TMHMM[ grep("# jgi", Suipal1_TMHMM$header),])
Suiocc1_TMHMM_sm<- data.frame(Suiocc1_TMHMM[ grep("# jgi", Suiocc1_TMHMM$header),])
Suilu4_TMHMM_sm<- data.frame(Suilu4_TMHMM[ grep("# jgi", Suilu4_TMHMM$header),])
Suilak1_TMHMM_sm<- data.frame(Suilak1_TMHMM[ grep("# jgi", Suilak1_TMHMM$header),])
Suihi1_TMHMM_sm<- data.frame(Suihi1_TMHMM[ grep("# jgi", Suihi1_TMHMM$header),])
Suigr1_TMHMM_sm<- data.frame(Suigr1_TMHMM[ grep("# jgi", Suigr1_TMHMM$header),])
Suidec1_TMHMM_sm<- data.frame(Suidec1_TMHMM[ grep("# jgi", Suidec1_TMHMM$header),])
Suicot1_TMHMM_sm<- data.frame(Suicot1_TMHMM[ grep("# jgi", Suicot1_TMHMM$header),])
Suicli1_TMHMM_sm<- data.frame(Suicli1_TMHMM[ grep("# jgi", Suicli1_TMHMM$header),])
Suibr2_TMHMM_sm<- data.frame(Suibr2_TMHMM[ grep("# jgi", Suibr2_TMHMM$header),])
Suibov1_TMHMM_sm<- data.frame(Suibov1_TMHMM[ grep("# jgi", Suibov1_TMHMM$header),])
Suiamp1_TMHMM_sm<- data.frame(Suiamp1_TMHMM[ grep("# jgi", Suiamp1_TMHMM$header),])
Suiame1_TMHMM_sm<- data.frame(Suiame1_TMHMM[ grep("# jgi", Suiame1_TMHMM$header),])



#combine the size and predicted TMD# by columns 
Suivar1_TMHMM_sm.by.col<- data.frame(cbind(Suivar1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suivar1_TMHMM_sm[c(FALSE, TRUE), ]))
Suitom1_TMHMM_sm.by.col<- data.frame(cbind(Suitom1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suitom1_TMHMM_sm[c(FALSE, TRUE), ]))
Suisub1_TMHMM_sm.by.col<- data.frame(cbind(Suisub1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suisub1_TMHMM_sm[c(FALSE, TRUE), ]))
Suisu1_TMHMM_sm.by.col<- data.frame(cbind(Suisu1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suisu1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipla1_TMHMM_sm.by.col<- data.frame(cbind(Suipla1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suipla1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipic1_TMHMM_sm.by.col<- data.frame(cbind(Suipic1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suipic1_TMHMM_sm[c(FALSE, TRUE), ]))
Suipal1_TMHMM_sm.by.col<- data.frame(cbind(Suipal1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suipal1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiocc1_TMHMM_sm.by.col<- data.frame(cbind(Suiocc1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suiocc1_TMHMM_sm[c(FALSE, TRUE), ]))
Suilu4_TMHMM_sm.by.col<- data.frame(cbind(Suilu4_TMHMM_sm[c(TRUE, FALSE), ],
                                Suilu4_TMHMM_sm[c(FALSE, TRUE), ]))
Suilak1_TMHMM_sm.by.col<- data.frame(cbind(Suilak1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suilak1_TMHMM_sm[c(FALSE, TRUE), ]))
Suihi1_TMHMM_sm.by.col<- data.frame(cbind(Suihi1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suihi1_TMHMM_sm[c(FALSE, TRUE), ]))
Suigr1_TMHMM_sm.by.col<- data.frame(cbind(Suigr1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suigr1_TMHMM_sm[c(FALSE, TRUE), ]))
Suidec1_TMHMM_sm.by.col<- data.frame(cbind(Suidec1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suidec1_TMHMM_sm[c(FALSE, TRUE), ]))
Suicot1_TMHMM_sm.by.col<- data.frame(cbind(Suicot1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suicot1_TMHMM_sm[c(FALSE, TRUE), ]))
Suicli1_TMHMM_sm.by.col<- data.frame(cbind(Suicli1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suicli1_TMHMM_sm[c(FALSE, TRUE), ]))
Suibr2_TMHMM_sm.by.col<- data.frame(cbind(Suibr2_TMHMM_sm[c(TRUE, FALSE), ],
                                Suibr2_TMHMM_sm[c(FALSE, TRUE), ]))
Suibov1_TMHMM_sm.by.col<- data.frame(cbind(Suibov1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suibov1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiamp1_TMHMM_sm.by.col<- data.frame(cbind(Suiamp1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suiamp1_TMHMM_sm[c(FALSE, TRUE), ]))
Suiame1_TMHMM_sm.by.col<- data.frame(cbind(Suiame1_TMHMM_sm[c(TRUE, FALSE), ],
                                Suiame1_TMHMM_sm[c(FALSE, TRUE), ]))


#only proteins with no TMD
Suivar1_TMHMM_no_TMD<- data.frame(Suivar1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suivar1_TMHMM_sm.by.col$X2), ])
Suitom1_TMHMM_no_TMD<- data.frame(Suitom1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suitom1_TMHMM_sm.by.col$X2),])
Suisub1_TMHMM_no_TMD<- data.frame(Suisub1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suisub1_TMHMM_sm.by.col$X2),])
Suisu1_TMHMM_no_TMD<- data.frame(Suisu1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suisu1_TMHMM_sm.by.col$X2),])
Suipla1_TMHMM_no_TMD<- data.frame(Suipla1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipla1_TMHMM_sm.by.col$X2),])
Suipic1_TMHMM_no_TMD<- data.frame(Suipic1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipic1_TMHMM_sm.by.col$X2),])
Suipal1_TMHMM_no_TMD<- data.frame(Suipal1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suipal1_TMHMM_sm.by.col$X2),])
Suiocc1_TMHMM_no_TMD<- data.frame(Suiocc1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiocc1_TMHMM_sm.by.col$X2),])
Suilu4_TMHMM_no_TMD<- data.frame(Suilu4_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suilu4_TMHMM_sm.by.col$X2),])
Suilak1_TMHMM_no_TMD<- data.frame(Suilak1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suilak1_TMHMM_sm.by.col$X2),])
Suihi1_TMHMM_no_TMD<- data.frame(Suihi1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suihi1_TMHMM_sm.by.col$X2),])
Suigr1_TMHMM_no_TMD<- data.frame(Suigr1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suigr1_TMHMM_sm.by.col$X2),])
Suidec1_TMHMM_no_TMD<- data.frame(Suidec1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suidec1_TMHMM_sm.by.col$X2),])
Suicot1_TMHMM_no_TMD<- data.frame(Suicot1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suicot1_TMHMM_sm.by.col$X2),])
Suicli1_TMHMM_no_TMD<- data.frame(Suicli1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suicli1_TMHMM_sm.by.col$X2),])
Suibr2_TMHMM_no_TMD<- data.frame(Suibr2_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suibr2_TMHMM_sm.by.col$X2),])
Suibov1_TMHMM_no_TMD<- data.frame(Suibov1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suibov1_TMHMM_sm.by.col$X2),])
Suiamp1_TMHMM_no_TMD<- data.frame(Suiamp1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiamp1_TMHMM_sm.by.col$X2),])
Suiame1_TMHMM_no_TMD<- data.frame(Suiame1_TMHMM_sm.by.col[ grep("Number of predicted TMHs:  0", Suiame1_TMHMM_sm.by.col$X2),])

#function to isolate the aa length
isolate_aa <- function(in_df) { 
a<- data.frame(str_split(in_df$X1, "Length:"))
b<- a[2,]
c<- data.table::transpose(b)
out_df<- cbind(in_df, c)
return= out_df
}

#run the function
Suivar1_df<- isolate_aa(Suivar1_TMHMM_no_TMD)
Suitom1_df<- isolate_aa(Suitom1_TMHMM_no_TMD)
Suisub1_df<- isolate_aa(Suisub1_TMHMM_no_TMD)
Suisu1_df<- isolate_aa(Suisu1_TMHMM_no_TMD)
Suipla1_df<- isolate_aa(Suipla1_TMHMM_no_TMD)
Suipic1_df<- isolate_aa(Suipic1_TMHMM_no_TMD)
Suipal1_df<- isolate_aa(Suipal1_TMHMM_no_TMD)
Suiocc1_df<- isolate_aa(Suiocc1_TMHMM_no_TMD)
Suilu4_df<- isolate_aa(Suilu4_TMHMM_no_TMD)
Suilak1_df<- isolate_aa(Suilak1_TMHMM_no_TMD)
Suihi1_df<- isolate_aa(Suihi1_TMHMM_no_TMD)
Suigr1_df<- isolate_aa(Suigr1_TMHMM_no_TMD)
Suidec1_df<- isolate_aa(Suidec1_TMHMM_no_TMD)
Suicot1_df<- isolate_aa(Suicot1_TMHMM_no_TMD)
Suicli1_df<- isolate_aa(Suicli1_TMHMM_no_TMD)
Suibr2_df<- isolate_aa(Suibr2_TMHMM_no_TMD)
Suibov1_df<- isolate_aa(Suibov1_TMHMM_no_TMD)
Suiamp1_df<- isolate_aa(Suiamp1_TMHMM_no_TMD)
Suiame1_df<- isolate_aa(Suiame1_TMHMM_no_TMD)

#gotta convert the aa number into a a numeric type or R think's it's a char. string.
Suivar1_df.num<-transform(Suivar1_df, V1 = as.numeric(V1))
Suitom1_df.num<-transform(Suitom1_df, V1 = as.numeric(V1))
Suisub1_df.num<-transform(Suisub1_df, V1 = as.numeric(V1))
Suisu1_df.num<-transform(Suisu1_df, V1 = as.numeric(V1))
Suipla1_df.num<-transform(Suipla1_df, V1 = as.numeric(V1))
Suipic1_df.num<-transform(Suipic1_df, V1 = as.numeric(V1))
Suipal1_df.num<-transform(Suipal1_df, V1 = as.numeric(V1))
Suiocc1_df.num<-transform(Suiocc1_df, V1 = as.numeric(V1))
Suilu4_df.num<-transform(Suilu4_df, V1 = as.numeric(V1))
Suilak1_df.num<-transform(Suilak1_df, V1 = as.numeric(V1))
Suihi1_df.num<-transform(Suihi1_df, V1 = as.numeric(V1))
Suigr1_df.num<-transform(Suigr1_df, V1 = as.numeric(V1))
Suidec1_df.num<-transform(Suidec1_df, V1 = as.numeric(V1))
Suicot1_df.num<-transform(Suicot1_df, V1 = as.numeric(V1))
Suicli1_df.num<-transform(Suicli1_df, V1 = as.numeric(V1))
Suibr2_df.num<-transform(Suibr2_df, V1 = as.numeric(V1))
Suibov1_df.num<-transform(Suibov1_df, V1 = as.numeric(V1))
Suiamp1_df.num<-transform(Suiamp1_df, V1 = as.numeric(V1))
Suiame1_df.num<-transform(Suiame1_df, V1 = as.numeric(V1))

#only proteins < 300 aa in length
Suivar1_TMHMM_no_TMD_300<- Suivar1_df.num[Suivar1_df.num$V1 < 300, ]
Suitom1_TMHMM_no_TMD_300<- Suitom1_df.num[Suitom1_df.num$V1 < 300, ]
Suisub1_TMHMM_no_TMD_300<- Suisub1_df.num[Suisub1_df.num$V1 < 300, ]
Suisu1_TMHMM_no_TMD_300<- Suisu1_df.num[Suisu1_df.num$V1 < 300, ]
Suipla1_TMHMM_no_TMD_300<- Suipla1_df.num[Suipla1_df.num$V1 < 300, ]
Suipic1_TMHMM_no_TMD_300<- Suipic1_df.num[Suipic1_df.num$V1 < 300, ]
Suipal1_TMHMM_no_TMD_300<- Suipal1_df.num[Suipal1_df.num$V1 < 300, ]
Suiocc1_TMHMM_no_TMD_300<- Suiocc1_df.num[Suiocc1_df.num$V1 < 300, ]
Suilu4_TMHMM_no_TMD_300<- Suilu4_df.num[Suilu4_df.num$V1 < 300, ]
Suilak1_TMHMM_no_TMD_300<- Suilak1_df.num[Suilak1_df.num$V1 < 300, ]
Suihi1_TMHMM_no_TMD_300<- Suihi1_df.num[Suihi1_df.num$V1 < 300, ]
Suigr1_TMHMM_no_TMD_300<- Suigr1_df.num[Suigr1_df.num$V1 < 300, ]
Suidec1_TMHMM_no_TMD_300<- Suidec1_df.num[Suidec1_df.num$V1 < 300, ]
Suicot1_TMHMM_no_TMD_300<- Suicot1_df.num[Suicot1_df.num$V1 < 300, ]
Suicli1_TMHMM_no_TMD_300<- Suicli1_df.num[Suicli1_df.num$V1 < 300, ]
Suibr2_TMHMM_no_TMD_300<- Suibr2_df.num[Suibr2_df.num$V1 < 300, ]
Suibov1_TMHMM_no_TMD_300<- Suibov1_df.num[Suibov1_df.num$V1 < 300, ]
Suiamp1_TMHMM_no_TMD_300<- Suiamp1_df.num[Suiamp1_df.num$V1 < 300, ]
Suiame1_TMHMM_no_TMD_300<- Suiame1_df.num[Suiame1_df.num$V1 < 300, ]


#get #SSP's per genome
#get summary numbers for each of these and attach them to the previous output file. 
Suivar1<- nrow(Suivar1_TMHMM_no_TMD_300)
Suitom1<- nrow(Suitom1_TMHMM_no_TMD_300)
Suisub1<- nrow(Suisub1_TMHMM_no_TMD_300)
Suisu1<- nrow(Suisu1_TMHMM_no_TMD_300)
Suipla1<- nrow(Suipla1_TMHMM_no_TMD_300)
Suipic1<- nrow(Suipic1_TMHMM_no_TMD_300)
Suipal1<- nrow(Suipal1_TMHMM_no_TMD_300)
Suiocc1<- nrow(Suiocc1_TMHMM_no_TMD_300)
Suilu4<- nrow(Suilu4_TMHMM_no_TMD_300)
Suilak1<- nrow(Suilak1_TMHMM_no_TMD_300)
Suihi1<- nrow(Suihi1_TMHMM_no_TMD_300)
Suigr1<- nrow(Suigr1_TMHMM_no_TMD_300)
Suidec1<- nrow(Suidec1_TMHMM_no_TMD_300)
Suicot1<- nrow(Suicot1_TMHMM_no_TMD_300)
Suicli1<- nrow(Suicli1_TMHMM_no_TMD_300)
Suibr2<- nrow(Suibr2_TMHMM_no_TMD_300)
Suibov1<- nrow(Suibov1_TMHMM_no_TMD_300)
Suiamp1<- nrow(Suiamp1_TMHMM_no_TMD_300)
Suiame1<- nrow(Suiame1_TMHMM_no_TMD_300)


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
                      row.names = "#SSPs_signalP,TMHMM,lt_300aa")
totals.2[1,]

#write out totals to a file
#results<- rbind(totals, totals.2[1,])
write.csv(totals.2, quote = FALSE, file = "SSP_totals.csv")


#now parce out the SSP's from the original aa sequences. 

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


#function to isolate the isolate_ID
isolate_ID <- function(in_df) { 
  a<- data.frame(strsplit(as.character(in_df[,1]), split=" "))
  b<- data.frame(a[1,])
  c<- data.table::transpose(b)
  return= c
}


  
#run function to isolate IDs  
Suivar1_subset <- data.frame(substr(Suivar1_TMHMM_no_TMD_300$X1, 3,nchar(Suivar1_TMHMM_no_TMD_300$X1)))
Suivar1.2_subset<- isolate_ID(Suivar1_subset) 
Suivar1.3_subset<- data.frame(lapply(Suivar1.2_subset, gsub, pattern=' ', replacement=''))

Suitom1_subset <- data.frame(substr(Suitom1_TMHMM_no_TMD_300$X1, 3,nchar(Suitom1_TMHMM_no_TMD_300$X1)))
Suitom1.2_subset<- isolate_ID(Suitom1_subset) 
Suitom1.3_subset<- data.frame(lapply(Suitom1.2_subset, gsub, pattern=' ', replacement=''))

Suisub1_subset <- data.frame(substr(Suisub1_TMHMM_no_TMD_300$X1, 3,nchar(Suisub1_TMHMM_no_TMD_300$X1)))
Suisub1.2_subset<- isolate_ID(Suisub1_subset) 
Suisub1.3_subset<- data.frame(lapply(Suisub1.2_subset, gsub, pattern=' ', replacement=''))

Suisu1_subset <- data.frame(substr(Suisu1_TMHMM_no_TMD_300$X1, 3,nchar(Suisu1_TMHMM_no_TMD_300$X1)))
Suisu1.2_subset<- isolate_ID(Suisu1_subset) 
Suisu1.3_subset<- data.frame(lapply(Suisu1.2_subset, gsub, pattern=' ', replacement=''))

Suipla1_subset <- data.frame(substr(Suipla1_TMHMM_no_TMD_300$X1, 3,nchar(Suipla1_TMHMM_no_TMD_300$X1)))
Suipla1.2_subset<- isolate_ID(Suipla1_subset) 
Suipla1.3_subset<- data.frame(lapply(Suipla1.2_subset, gsub, pattern=' ', replacement=''))

Suipic1_subset <- data.frame(substr(Suipic1_TMHMM_no_TMD_300$X1, 3,nchar(Suipic1_TMHMM_no_TMD_300$X1)))
Suipic1.2_subset<- isolate_ID(Suipic1_subset) 
Suipic1.3_subset<- data.frame(lapply(Suipic1.2_subset, gsub, pattern=' ', replacement=''))

Suipal1_subset <- data.frame(substr(Suipal1_TMHMM_no_TMD_300$X1, 3,nchar(Suipal1_TMHMM_no_TMD_300$X1)))
Suipal1.2_subset<- isolate_ID(Suipal1_subset) 
Suipal1.3_subset<- data.frame(lapply(Suipal1.2_subset, gsub, pattern=' ', replacement=''))

Suiocc1_subset <- data.frame(substr(Suiocc1_TMHMM_no_TMD_300$X1, 3,nchar(Suiocc1_TMHMM_no_TMD_300$X1)))
Suiocc1.2_subset<- isolate_ID(Suiocc1_subset) 
Suiocc1.3_subset<- data.frame(lapply(Suiocc1.2_subset, gsub, pattern=' ', replacement=''))

Suilu4_subset <- data.frame(substr(Suilu4_TMHMM_no_TMD_300$X1, 3,nchar(Suilu4_TMHMM_no_TMD_300$X1)))
Suilu4.2_subset<- isolate_ID(Suilu4_subset) 
Suilu4.3_subset<- data.frame(lapply(Suilu4.2_subset, gsub, pattern=' ', replacement=''))

Suilak1_subset <- data.frame(substr(Suilak1_TMHMM_no_TMD_300$X1, 3,nchar(Suilak1_TMHMM_no_TMD_300$X1)))
Suilak1.2_subset<- isolate_ID(Suilak1_subset) 
Suilak1.3_subset<- data.frame(lapply(Suilak1.2_subset, gsub, pattern=' ', replacement=''))

Suihi1_subset <- data.frame(substr(Suihi1_TMHMM_no_TMD_300$X1, 3,nchar(Suihi1_TMHMM_no_TMD_300$X1)))
Suihi1.2_subset<- isolate_ID(Suihi1_subset) 
Suihi1.3_subset<- data.frame(lapply(Suihi1.2_subset, gsub, pattern=' ', replacement=''))

Suigr1_subset <- data.frame(substr(Suigr1_TMHMM_no_TMD_300$X1, 3,nchar(Suigr1_TMHMM_no_TMD_300$X1)))
Suigr1.2_subset<- isolate_ID(Suigr1_subset) 
Suigr1.3_subset<- data.frame(lapply(Suigr1.2_subset, gsub, pattern=' ', replacement=''))

Suidec1_subset <- data.frame(substr(Suidec1_TMHMM_no_TMD_300$X1, 3,nchar(Suidec1_TMHMM_no_TMD_300$X1)))
Suidec1.2_subset<- isolate_ID(Suidec1_subset) 
Suidec1.3_subset<- data.frame(lapply(Suidec1.2_subset, gsub, pattern=' ', replacement=''))

Suicot1_subset <- data.frame(substr(Suicot1_TMHMM_no_TMD_300$X1, 3,nchar(Suicot1_TMHMM_no_TMD_300$X1)))
Suicot1.2_subset<- isolate_ID(Suicot1_subset) 
Suicot1.3_subset<- data.frame(lapply(Suicot1.2_subset, gsub, pattern=' ', replacement=''))

Suicli1_subset <- data.frame(substr(Suicli1_TMHMM_no_TMD_300$X1, 3,nchar(Suicli1_TMHMM_no_TMD_300$X1)))
Suicli1.2_subset<- isolate_ID(Suicli1_subset) 
Suicli1.3_subset<- data.frame(lapply(Suicli1.2_subset, gsub, pattern=' ', replacement=''))

Suibr2_subset <- data.frame(substr(Suibr2_TMHMM_no_TMD_300$X1, 3,nchar(Suibr2_TMHMM_no_TMD_300$X1)))
Suibr2.2_subset<- isolate_ID(Suibr2_subset) 
Suibr2.3_subset<- data.frame(lapply(Suibr2.2_subset, gsub, pattern=' ', replacement=''))

Suibov1_subset <- data.frame(substr(Suibov1_TMHMM_no_TMD_300$X1, 3,nchar(Suibov1_TMHMM_no_TMD_300$X1)))
Suibov1.2_subset<- isolate_ID(Suibov1_subset) 
Suibov1.3_subset<- data.frame(lapply(Suibov1.2_subset, gsub, pattern=' ', replacement=''))

Suiamp1_subset <- data.frame(substr(Suiamp1_TMHMM_no_TMD_300$X1, 3,nchar(Suiamp1_TMHMM_no_TMD_300$X1)))
Suiamp1.2_subset<- isolate_ID(Suiamp1_subset) 
Suiamp1.3_subset<- data.frame(lapply(Suiamp1.2_subset, gsub, pattern=' ', replacement=''))

Suiame1_subset <- data.frame(substr(Suiame1_TMHMM_no_TMD_300$X1, 3,nchar(Suiame1_TMHMM_no_TMD_300$X1)))
Suiame1.2_subset<- isolate_ID(Suiame1_subset) 
Suiame1.3_subset<- data.frame(lapply(Suiame1.2_subset, gsub, pattern=' ', replacement=''))


#get fastas of only positive hits for EffectorP analysis 
Suivar1_SSP_fastas<- Suivar1_in[c(which(names(Suivar1_in) %in% Suivar1.3_subset$V1))]
Suitom1_SSP_fastas<- Suitom1_in[c(which(names(Suitom1_in) %in% Suitom1.3_subset$V1))]
Suisub1_SSP_fastas<- Suisub1_in[c(which(names(Suisub1_in) %in% Suisub1.3_subset$V1))]
Suisu1_SSP_fastas<- Suisu1_in[c(which(names(Suisu1_in) %in% Suisu1.3_subset$V1))]
Suipla1_SSP_fastas<- Suipla1_in[c(which(names(Suipla1_in) %in% Suipla1.3_subset$V1))]
Suipic1_SSP_fastas<- Suipic1_in[c(which(names(Suipic1_in) %in% Suipic1.3_subset$V1))]
Suipal1_SSP_fastas<- Suipal1_in[c(which(names(Suipal1_in) %in% Suipal1.3_subset$V1))]
Suiocc1_SSP_fastas<- Suiocc1_in[c(which(names(Suiocc1_in) %in% Suiocc1.3_subset$V1))]
Suilu4_SSP_fastas<- Suilu4_in[c(which(names(Suilu4_in) %in% Suilu4.3_subset$V1))]
Suilak1_SSP_fastas<- Suilak1_in[c(which(names(Suilak1_in) %in% Suilak1.3_subset$V1))]
Suihi1_SSP_fastas<- Suihi1_in[c(which(names(Suihi1_in) %in% Suihi1.3_subset$V1))]
Suigr1_SSP_fastas<- Suigr1_in[c(which(names(Suigr1_in) %in% Suigr1.3_subset$V1))]
Suidec1_SSP_fastas<- Suidec1_in[c(which(names(Suidec1_in) %in% Suidec1.3_subset$V1))]
Suicot1_SSP_fastas<- Suicot1_in[c(which(names(Suicot1_in) %in% Suicot1.3_subset$V1))]
Suicli1_SSP_fastas<- Suicli1_in[c(which(names(Suicli1_in) %in% Suicli1.3_subset$V1))]
Suibr2_SSP_fastas<- Suibr2_in[c(which(names(Suibr2_in) %in% Suibr2.3_subset$V1))]
Suibov1_SSP_fastas<- Suibov1_in[c(which(names(Suibov1_in) %in% Suibov1.3_subset$V1))]
Suiamp1_SSP_fastas<- Suiamp1_in[c(which(names(Suiamp1_in) %in% Suiamp1.3_subset$V1))]
Suiame1_SSP_fastas<- Suiame1_in[c(which(names(Suiame1_in) %in% Suiame1.3_subset$V1))]

#make sure nothing screwy went on and the fasta are the correct lenght
#n SSP's 
totals.2[1,]
# length of a few files
length(Suivar1_SSP_fastas)
length(Suitom1_SSP_fastas)
length(Suiamp1_SSP_fastas)
length(Suisu1_SSP_fastas)
# numbers match. 

#print fastas to files to run EffectorP
write.fasta(Suivar1_SSP_fastas, names = names(Suivar1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suivar1_SSPs_for_EffectorP.fasta")
write.fasta(Suitom1_SSP_fastas, names = names(Suitom1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suitom1_SSPs_for_EffectorP.fasta")
write.fasta(Suisub1_SSP_fastas, names = names(Suisub1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisub1_SSPs_for_EffectorP.fasta")
write.fasta(Suisu1_SSP_fastas, names = names(Suisu1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suisu1_SSPs_for_EffectorP.fasta")
write.fasta(Suipla1_SSP_fastas, names = names(Suipla1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipla1_SSPs_for_EffectorP.fasta")
write.fasta(Suipic1_SSP_fastas, names = names(Suipic1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipic1_SSPs_for_EffectorP.fasta")
write.fasta(Suipal1_SSP_fastas, names = names(Suipal1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suipal1_SSPs_for_EffectorP.fasta")
write.fasta(Suiocc1_SSP_fastas, names = names(Suiocc1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiocc1_SSPs_for_EffectorP.fasta")
write.fasta(Suilu4_SSP_fastas, names = names(Suilu4_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilu4_SSPs_for_EffectorP.fasta")
write.fasta(Suilak1_SSP_fastas, names = names(Suilak1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suilak1_SSPs_for_EffectorP.fasta")
write.fasta(Suihi1_SSP_fastas, names = names(Suihi1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suihi1_SSPs_for_EffectorP.fasta")
write.fasta(Suigr1_SSP_fastas, names = names(Suigr1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suigr1_SSPs_for_EffectorP.fasta")
write.fasta(Suidec1_SSP_fastas, names = names(Suidec1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suidec1_SSPs_for_EffectorP.fasta")
write.fasta(Suicot1_SSP_fastas, names = names(Suicot1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicot1_SSPs_for_EffectorP.fasta")
write.fasta(Suicli1_SSP_fastas, names = names(Suicli1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suicli1_SSPs_for_EffectorP.fasta")
write.fasta(Suibr2_SSP_fastas, names = names(Suibr2_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibr2_SSPs_for_EffectorP.fasta")
write.fasta(Suibov1_SSP_fastas, names = names(Suibov1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suibov1_SSPs_for_EffectorP.fasta")
write.fasta(Suiamp1_SSP_fastas, names = names(Suiamp1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiamp1_SSPs_for_EffectorP.fasta")
write.fasta(Suiame1_SSP_fastas, names = names(Suiame1_SSP_fastas), open = "w", nbchar = 60, as.string = FALSE, file.out = "Suiame1_SSPs_for_EffectorP.fasta")


########
#now run EffectorP version2 on the web interface

#######


#read in EffectorP output files
Suivar1_effectors<- seqinr::read.fasta(file = "Suivar1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suitom1_effectors<- seqinr::read.fasta(file = "Suitom1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisub1_effectors<- seqinr::read.fasta(file = "Suisub1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suisu1_effectors<- seqinr::read.fasta(file = "Suisu1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipla1_effectors<- seqinr::read.fasta(file = "Suipla1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipic1_effectors<- seqinr::read.fasta(file = "Suipic1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suipal1_effectors<- seqinr::read.fasta(file = "Suipal1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiocc1_effectors<- seqinr::read.fasta(file = "Suiocc1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilu4_effectors<- seqinr::read.fasta(file = "Suilu4_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suilak1_effectors<- seqinr::read.fasta(file = "Suilak1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suihi1_effectors<- seqinr::read.fasta(file = "Suihi1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suigr1_effectors<- seqinr::read.fasta(file = "Suigr1_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suidec1_effectors<- seqinr::read.fasta(file = "Suidec1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicot1_effectors<- seqinr::read.fasta(file = "Suicot1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suicli1_effectors<- seqinr::read.fasta(file = "Suicli1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibr2_effectors<- seqinr::read.fasta(file = "Suibr2_EffectorCandidates.fasta", 
                                      seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suibov1_effectors<- seqinr::read.fasta(file = "Suibov1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiamp1_effectors<- seqinr::read.fasta(file = "Suiamp1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)
Suiame1_effectors<- seqinr::read.fasta(file = "Suiame1_EffectorCandidates.fasta", 
                                       seqtype = "AA",as.string = TRUE, set.attributes = FALSE)



#get totals
#total number of putative effectors
Suivar1<- length(Suivar1_effectors)
Suitom1<- length(Suitom1_effectors)
Suisub1<- length(Suisub1_effectors)
Suisu1<- length(Suisu1_effectors)
Suipla1<- length(Suipla1_effectors)
Suipic1<- length(Suipic1_effectors)
Suipal1<- length(Suipal1_effectors)
Suiocc1<- length(Suiocc1_effectors)
Suilu4<- length(Suilu4_effectors)
Suilak1<- length(Suilak1_effectors)
Suihi1<- length(Suihi1_effectors)
Suigr1<- length(Suigr1_effectors)
Suidec1<- length(Suidec1_effectors)
Suicot1<- length(Suicot1_effectors)
Suicli1<- length(Suicli1_effectors)
Suibr2<- length(Suibr2_effectors)
Suibov1<- length(Suibov1_effectors)
Suiamp1<- length(Suiamp1_effectors)
Suiame1<- length(Suiame1_effectors)

total_effectors<- data.frame(cbind(Suivar1, 
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
                                   Suiame1),
                             row.names = "n_putative_effectors_from_EffectorP")



##SSP's per genome
Suivar1<- nrow(Suivar1_TMHMM_no_TMD_300)
Suitom1<- nrow(Suitom1_TMHMM_no_TMD_300)
Suisub1<- nrow(Suisub1_TMHMM_no_TMD_300)
Suisu1<- nrow(Suisu1_TMHMM_no_TMD_300)
Suipla1<- nrow(Suipla1_TMHMM_no_TMD_300)
Suipic1<- nrow(Suipic1_TMHMM_no_TMD_300)
Suipal1<- nrow(Suipal1_TMHMM_no_TMD_300)
Suiocc1<- nrow(Suiocc1_TMHMM_no_TMD_300)
Suilu4<- nrow(Suilu4_TMHMM_no_TMD_300)
Suilak1<- nrow(Suilak1_TMHMM_no_TMD_300)
Suihi1<- nrow(Suihi1_TMHMM_no_TMD_300)
Suigr1<- nrow(Suigr1_TMHMM_no_TMD_300)
Suidec1<- nrow(Suidec1_TMHMM_no_TMD_300)
Suicot1<- nrow(Suicot1_TMHMM_no_TMD_300)
Suicli1<- nrow(Suicli1_TMHMM_no_TMD_300)
Suibr2<- nrow(Suibr2_TMHMM_no_TMD_300)
Suibov1<- nrow(Suibov1_TMHMM_no_TMD_300)
Suiamp1<- nrow(Suiamp1_TMHMM_no_TMD_300)
Suiame1<- nrow(Suiame1_TMHMM_no_TMD_300)


total_SSPs<- data.frame(cbind(Suivar1, 
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
                            Suiame1),
                      row.names = "#SSPs_signalP,TMHMM,lt_300aa")


#Proteins per genome
Suivar1<- length(Suivar1_in)
Suitom1<- length(Suitom1_in)
Suisub1<- length(Suisub1_in)
Suisu1<- length(Suisu1_in)
Suipla1<- length(Suipla1_in)
Suipic1<- length(Suipic1_in)
Suipal1<- length(Suipal1_in)
Suiocc1<- length(Suiocc1_in)
Suilu4<- length(Suilu4_in)
Suilak1<- length(Suilak1_in)
Suihi1<- length(Suihi1_in)
Suigr1<- length(Suigr1_in)
Suidec1<- length(Suidec1_in)
Suicot1<- length(Suicot1_in)
Suicli1<- length(Suicli1_in)
Suibr2<- length(Suibr2_in)
Suibov1<- length(Suibov1_in)
Suiamp1<- length(Suiamp1_in)
Suiame1<- length(Suiame1_in)


total_proteins<-data.frame(cbind(Suivar1, 
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
                                 Suiame1),
                           row.names = "n_putative_proteins_from_gene_cat")  

#bind results and transform
results<- rbind(total_proteins[1,], total_SSPs[1,], total_effectors[1,])
results.1<- t(results)

#get % SSP's out of total genes
percent_SSPs_out_of_total_genes<- (results.1[,2] / results.1[,1])*100
percent_SSPs_out_of_total_genes<- signif(percent_SSPs_out_of_total_genes, digits = 3)
#get % of Effectors out of total SSP's 
percent_effectors_out_of_total_genes<- (results.1[,3] / results.1[,1])*100
percent_effectors_out_of_total_genes<- signif(percent_effectors_out_of_total_genes, digits = 2)
#get % of Effectors out of total genes
percent_effectors_out_of_SSPs<- (results.1[,3] / results.1[,2])*100
percent_effectors_out_of_SSPs<-signif(percent_effectors_out_of_SSPs, digits = 4)

#slap it together
results_with_percentages<- cbind(results.1, percent_SSPs_out_of_total_genes, percent_effectors_out_of_total_genes, percent_effectors_out_of_SSPs)
  
#write it up
write.csv(results, quote = FALSE, file = "SSP_and_Effector_totals.csv")

