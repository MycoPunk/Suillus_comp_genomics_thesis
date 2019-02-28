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

cbind(df[c(TRUE, FALSE), ],
      df[c(FALSE, TRUE), ])




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
Suisub_df<- isolate_aa(Suisub1_TMHMM_no_TMD)
Suisu1_df<- isolate_aa(Suisu1_TMHMM_no_TMD)
Suipla1_df<- isolate_aa(Suipla1_TMHMM_no_TMD)
Suipic1_df<- isolate_aa(Suipic1_TMHMM_no_TMD)
Suipal_df<- isolate_aa(Suipal1_TMHMM_no_TMD)
Suiocc_df<- isolate_aa(Suiocc1_TMHMM_no_TMD)
Suilu4_df<- isolate_aa(Suilu4_TMHMM_no_TMD)
Suilak1_df<- isolate_aa(Suilak1_TMHMM_no_TMD)
Suihi1_df<- isolate_aa(Suihi1_TMHMM_no_TMD)
Suigr1_df<- isolate_aa(Suigr1_TMHMM_no_TMD)
Suidec1_df<- isolate_aa(Suidec1_TMHMM_no_TMD)
Suicot1_df<- isolate_aa(Suicot1_TMHMM_no_TMD)
Suicli_df<- isolate_aa(Suicli1_TMHMM_no_TMD)
Suiabr2_df<- isolate_aa(Suibr2_TMHMM_no_TMD)
Suibov1_df<- isolate_aa(Suibov1_TMHMM_no_TMD)
Suiamp_df<- isolate_aa(Suiamp1_TMHMM_no_TMD)
Suiame1_df<- isolate_aa(Suiame1_TMHMM_no_TMD)

#gotta convert the aa number into a a numeric type or R think's it's a char. string.
Suivar1_df.num<-transform(Suivar1_df, V1 = as.numeric(V1))






#only proteins < 300 aa in length
Suivar1_TMHMM_no_TMD_300<- Suivar1_df.num[Suivar1_df.num$V1 < 300, ]


#get #SSP's per genome
nrow(Suivar1_TMHMM_no_TMD_300)



