library(ggplot2)
library(tidyr)
library(dplyr)
library(ggpubr)
library(stringr)
library(reshape2)
library(gridExtra)

######################################################
###########Read SARS-CoV-2 genomic data from gff
#Read annotation and get data
GFFDf <- read.csv("./NC_045512.2.gff3", skip = 2, header = F, sep="\t")
geneDf<-subset(GFFDf,GFFDf$V3 == "gene")
geneDfA<-separate(data = geneDf, col = V9, into = c("A1","A2","gene","A4","A5","A6","A7"), sep = ";")
geneDfA$gene <- gsub('Name=', '', geneDfA$gene)


########################################################
#######Read insertions data
#Read insertions info
GISAIDins<-read.csv("./GISAID_msa_0617_ab_insertions.csv", sep = "\t",
                    header = F)

###!!!!This is the part that got updated after the revision!!!!
#Exclude sequences with multiple insertions
FilterDfCounts<-GISAIDins%>% group_by(V1) %>% count()
FilterDf<-subset(FilterDfCounts, FilterDfCounts$n <3)
GISAIDins_noMulti<-subset(GISAIDins,GISAIDins$V1 %in% FilterDf$V1)

#filter all non ATGC
INS_ATGC<-subset(GISAIDins_noMulti, grepl("^(a|t|g|c)+$",GISAIDins_noMulti$V5))
length(unique(INS_ATGC$V1))

###New in this version
#First let's subset all that multiple by three
INS_ATGC_three<-subset(INS_ATGC, INS_ATGC$V4 %% 3 == 0)
INS_ATGC_notthree<-subset(INS_ATGC, ((INS_ATGC$V4 %% 3 > 0)& (INS_ATGC$V4 > 1)))
#check basic stat
INS_ATGC_notthree %>% group_by(V4) %>% count()
INS_ATGC_notthree %>% group_by(V2) %>% count()

INS_ATGC_Save<-rbind(INS_ATGC_three, INS_ATGC_notthree)

####
#Save output
write.table(INS_ATGC_Save,"./GISAID_msa0617_insertions_filtered.csv", row.names = F, sep = "\t",
            quote = F)

##Extract insertion events
GISAIDins_noN_counts<-INS_ATGC_Save %>% group_by(V2,V3,V4,V5) %>% count()

#more summary
uniqueEventsNotThree<-subset(GISAIDins_noN_counts, GISAIDins_noN_counts$V4 %% 3 >0)
uniqueEventsNotThree %>% group_by(V2) %>% count()

