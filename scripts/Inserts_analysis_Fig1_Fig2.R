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
GFFDf <- read.csv("../data/NC_045512.2.gff3", skip = 2, header = F, sep="\t")
geneDf<-subset(GFFDf,GFFDf$V3 == "gene")
geneDfA<-separate(data = geneDf, col = V9, into = c("A1","A2","gene","A4","A5","A6","A7"), sep = ";")
geneDfA$gene <- gsub('Name=', '', geneDfA$gene)

ymin<-c(0,-1)
ymax<-c(1,0)
geneDfA$ymin<-rep(ymin, 6,length.out = length(geneDfA$V1))
geneDfA$ymax<-rep(ymax, 6,length.out = length(geneDfA$V1))


##################################################
##Draw SARS-CoV-2 genome

#Plot for figure 1
SARSCoVGenome <- ggplot(data = geneDfA)+
  geom_rect(aes(xmin = V4, xmax = V5, ymin = ymin, ymax = ymax,
                fill = gene),
            color = "black",
            size =0.2)+
  geom_rect(xmin=0,xmax=29903,ymin=-0.05,ymax=0.05)+
  geom_text(aes(x=V4+(V5-V4)/2, 
                y = rep(c(1.3,-1.3),6,length.out=(length(geneDfA$V1))),
                label = gene), 
            size =3.5,
            angle = 90, 
            hjust= rep(c(0,1),6,length.out=(length(geneDfA$V1))))+
  ylim(c(-2.6,+2.7))+
  scale_fill_manual(values =c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
                              "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                              "#d9d9d9", "#bc80bd", "#ccebc5"))+
  scale_x_continuous(expand = c(0.01, 0), limits = c(0,30000))+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(color = "white", size =10),
        axis.ticks.y = element_line(color = "white"),
        plot.margin=unit(c(0.2,1,0,1), "cm"))
SARSCoVGenome
#########
#Flat plot
#add flat x and y
geneDfA$ymin2<-rep(-0.3, length.out = length(geneDfA$V1))
geneDfA$ymax2<-rep(0, length.out = length(geneDfA$V1))
#
SARSCoVGenome_plain<- ggplot(data = geneDfA)+
  geom_rect(xmin=0,xmax=29903,ymin=-0.3,ymax=0,
            color = "black",
            size =0.2)+
  geom_rect(aes(xmin = V4, xmax = V5, ymin = ymin2, ymax = ymax2,
                fill = gene),
            color = "black",
            size =0.2)+
  # geom_text(aes(x=V4+(V5-V4)/2, 
  #               y = -0.5,
  #               label = gene), 
  #           size =3.5,
  #           angle = 90, 
  #           hjust= rep(c(0,1),6,length.out=(length(geneDfA$V1))))+
  scale_y_continuous(expand = c(0.01, 0))+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  scale_fill_manual(values =c("#8dd3c7", "#ffffb3", "#bebada", "#fb8072",
                              "#80b1d3", "#fdb462", "#b3de69", "#fccde5",
                              "#d9d9d9", "#bc80bd", "#ccebc5"))+
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_blank(),
        axis.line = element_blank(),
        axis.text.y = element_text(color = "white", size =10),
        axis.text.x= element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(color = "white"),
        plot.margin=unit(c(-0.1,1,0,1), "cm"))
SARSCoVGenome_plain

########################################################
########################################################
#######Read insertions data
#Read insertions info
GISAIDins<-read.csv("../data/GISAID_full_msa.insertions.out.csv", sep = "\t",
                    header = F)
#Exclude sequences with multiple insertions
FilterDfCounts<-GISAIDins%>% group_by(V1) %>% count()
FilterDf<-subset(FilterDfCounts, FilterDfCounts$n <3)
GISAIDins_noMulti<-subset(GISAIDins,GISAIDins$V1 %in% FilterDf$V1)

#filter all non ATGC
INS_ATGC<-subset(GISAIDins_noMulti, grepl("^(a|t|g|c)+$",GISAIDins_noMulti$V5))
length(unique(INS_ATGC$V1))

# INS_ATGC_Save<-subset(INS_ATGC, INS_ATGC$V4 %% 3 == 0)
# ####
# #Save output
# write.table(INS_ATGC_Save,"../data/GISAID_full_msa_insertions_filtered.csv", row.names = F, sep = "\t")

##Extract insertion events
GISAIDins_noN_counts<-INS_ATGC %>% group_by(V2,V3,V4,V5) %>% count()

#And prepare data for plotting
Insertions_noN_for_plotting<-GISAIDins_noN_counts[,c("V3","V4","V5","n")]
Insertions_noN_for_plotting$xmin<-Insertions_noN_for_plotting$V3 -40
Insertions_noN_for_plotting$xmax<-Insertions_noN_for_plotting$V3 +40
Insertions_noN_for_plotting$fill<-ifelse(Insertions_noN_for_plotting$n == 1, 
                                         "Singleton", "Multi")
#Calculate the multiplication
Insertions_noN_for_plotting$multi<-Insertions_noN_for_plotting$V4%%3

#Subset only those that multiply by 3
Threemulpl_insertions<-subset(Insertions_noN_for_plotting, Insertions_noN_for_plotting$multi ==0)
names(Threemulpl_insertions)<-c("Position","Length","Insertion", "Times", "xmin","xmax","fill","multi")
#correct wrong position in one case
Threemulpl_insertions$Position[Threemulpl_insertions$Position == 28028]<-28031

#Read data about monophyly and mechanisms
PhylogenyInsertions<-read.csv("../data/Insertions_phylogeny_all.csv", sep ="\t", header = T,
                              stringsAsFactors = F)

#Merge dataframes and add statuses
Insertions_plotting_merged<-merge(Threemulpl_insertions,PhylogenyInsertions, all=T)

Insertions_plotting_merged[is.na(Insertions_plotting_merged)]<-""

###Remove insertions that were not confirmed by reads
Insertions_plotting_merged_f<-subset(Insertions_plotting_merged, Insertions_plotting_merged$Confirmed_by_reads !="N")

Insertions_plotting_merged_f$Fill<-ifelse(Insertions_plotting_merged_f$Confirmed_by_reads == "Y",
                                        "Confirmed by reads", 
                                        ifelse(Insertions_plotting_merged_f$USHER == "Yes",
                                               "Monophyletic",
                                               ifelse(Insertions_plotting_merged_f$fill == "Multi", "Multiple genomes", "Singleton")))
Insertions_plotting_merged_f$Length2<-ifelse(Insertions_plotting_merged_f$Length < 30, 
                                           Insertions_plotting_merged_f$Length, 35)
# ####
# #export final table
# write.table(Insertions_plotting_merged_f, "./Inserts_141_all_features_export.csv", row.names = F, sep = "\t")

##############Plot insertions

#Plot length distribution (Figure 1a)
Insertions_length_distribution<- ggplot(data=Insertions_plotting_merged_f)+
  geom_bar(aes(x = Length,y = ..count..),
           fill = "#08306b")+
  theme_classic()+
  scale_x_continuous(limits = c(0,200), expand = c(0,0),
                     breaks= c(3,9,15,21,27,seq(50, 200, 50)))+
  scale_y_continuous(expand = c(0.01,0),
                     breaks = c(1,seq(5,80,5)))+
  theme(plot.margin=unit(c(0.2,1,0.1,1), "cm"))
Insertions_length_distribution

#Plot inserts distribution along genome (Figure 1b)
ThreemulInsertionsPlot<-ggplot(data=Insertions_plotting_merged_f)+
  geom_vline(aes(xintercept = Position), linetype = "dashed", color = "#d9d9d9",
             alpha=0.8)+
  geom_jitter(aes(x = Position,
                  y = Length2,
                  fill = Fill),
              color = "black",
              shape = 25, 
              size =5, 
              height = 0.5,
              alpha =0.7)+
  scale_fill_manual(values = c("#00441b","#41ab5d", "#c7e9c0", "#bababa"), name = "")+
  scale_y_continuous(breaks=c(0,3,6,9,12,15,18,21,24,27))+
  ylab("insertion length")+
  scale_x_continuous(expand= c(0.01,0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  geom_hline(yintercept = 7, color = "#a50f15")+
  theme_classic()+
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.y = element_text(size =10),
        plot.margin=unit(c(0.3,1,0.3,1), "cm"))
ThreemulInsertionsPlot

#############################################
#######Calculate nucleotide composition

#Subset inserts of the length 3
Three_insertions<-subset(Insertions_plotting_merged_f, Insertions_plotting_merged_f$Length < 4)
#Subset inserts of the length 6
Six_insertions<-subset(Insertions_plotting_merged_f, Insertions_plotting_merged_f$Length == 6)
##################################
#Subset long insertions
Long_insertions<-subset(Insertions_plotting_merged_f, Insertions_plotting_merged_f$Length > 8)
#############################
#Calculate nucleotide composition
nucl_composition<- function(df) {
  tempdf<-df[,c("Position","Insertion")]
  tempdf$a <- str_count(tempdf$Insertion, "a")
  tempdf$t <- str_count(tempdf$Insertion, "t")
  tempdf$g <- str_count(tempdf$Insertion, "g")
  tempdf$c <- str_count(tempdf$Insertion, "c")
  out<-colSums(tempdf[,c("a","t","g","c")])
  out_final<-unlist(c(df[1,2],out))
  return(out_final)
}
##############
Three_ins_nuc_counts<-nucl_composition(Three_insertions)
Six_ins_nuc_counts<-nucl_composition(Six_insertions)
Long_ins_nuc_counts<-nucl_composition(Long_insertions)
#Information for the whole genome is precalculated from the reference strain
Whole<-c(0,8954,9594,5863,5492)
##Merge results into dataframe
ATGC_countsDF<-as.data.frame(t(data.frame(Three_ins_nuc_counts,
                                          Six_ins_nuc_counts,Long_ins_nuc_counts, Whole)))

data<-ATGC_countsDF[,2:5]
ATGC_portionsDF<-as.data.frame(t(apply(data,1, function(x) x*100/sum(x))))
ATGC_portionsDF$length<-factor(c("3 nucl","6 nucl","> 8 nucl","genome"),
                               levels = c("Short c","3 nucl","6 nucl","> 8 nucl","genome"))
ATGC_portion_melted<-reshape2::melt(ATGC_portionsDF, id = "length")

##Plot results
NucleotidesCounts<-ggplot(data=ATGC_portion_melted, aes(x =length,y = value, fill=variable))+
  geom_bar(stat="identity",position=position_dodge(), color = "black", size =0.1)+
  scale_fill_manual(name = "Nucleotide", labels = c("A","T","G","C"),
                    values = c("#4daf4a","#e41a1c","#fecc5c","#377eb8"))+
  ylab("%")+
  xlab("")+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,1,0.3,1), "cm"))
NucleotidesCounts

####Supplementary figure
#Same as Fig1b but only for experimentally verified
ExperimentallyVerified<-subset(Insertions_plotting_merged_f, 
                               Insertions_plotting_merged_f$Fill == "Confirmed by reads")
Three_insertions_e<-subset(ExperimentallyVerified, ExperimentallyVerified$Length < 4)
Six_insertions_e<-subset(ExperimentallyVerified, ExperimentallyVerified$Length == 6)
Long_insertions_e<-subset(ExperimentallyVerified, ExperimentallyVerified$Length > 8)

Three_ins_nuc_counts_e<-nucl_composition(Three_insertions_e)
Six_ins_nuc_counts_e<-nucl_composition(Six_insertions_e)
Long_ins_nuc_counts_e<-nucl_composition(Long_insertions_e)

ATGC_countsDF_e<-as.data.frame(t(data.frame(Three_ins_nuc_counts_e,
                                          Six_ins_nuc_counts_e,Long_ins_nuc_counts_e, Whole)))

data_e<-ATGC_countsDF_e[,2:5]
ATGC_portionsDF_e<-as.data.frame(t(apply(data_e,1, function(x) x*100/sum(x))))
ATGC_portionsDF_e$length<-factor(c("3 nucl","6 nucl","> 8 nucl","genome"),
                               levels = c("Short c","3 nucl","6 nucl","> 8 nucl","genome"))
ATGC_portion_e_melted<-reshape2::melt(ATGC_portionsDF_e, id = "length")

#
NucleotidesCounts_e<-ggplot(data=ATGC_portion_e_melted, aes(x =length,y = value, fill=variable))+
  geom_bar(stat="identity",position=position_dodge(), color = "black", size =0.1)+
  scale_fill_manual(name = "Nucleotide", labels = c("A","T","G","C"),
                    values = c("#4daf4a","#e41a1c","#fecc5c","#377eb8"))+
  ylab("%")+
  xlab("")+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,1,0.3,1), "cm"))
NucleotidesCounts_e

#GC content in all inserts + PRRA
#Add PRRA to dataset
PRRA<-data_frame(Length = 12, Insertion = "cctcggcgggca", Fill = "PRRA")

SupplementaryNucCount<-Insertions_plotting_merged_f[,c("Length", "Insertion", "Fill")]
SupplementaryNucCountPRRA<-rbind(SupplementaryNucCount, PRRA)
SupplementaryNucCountPRRA$Length2<-ifelse((SupplementaryNucCountPRRA$Length > 8), "long",SupplementaryNucCountPRRA$Length)

SupplementaryNucCountPRRA$GC <- str_count(SupplementaryNucCountPRRA$Insertion,"g") +
  str_count(SupplementaryNucCountPRRA$Insertion,"c")
SupplementaryNucCountPRRA$AT <- str_count(SupplementaryNucCountPRRA$Insertion,"a") +
  str_count(SupplementaryNucCountPRRA$Insertion,"t")
SupplementaryNucCountPRRA$GCp <- SupplementaryNucCountPRRA$GC*100/(SupplementaryNucCountPRRA$GC +
                                                                     SupplementaryNucCountPRRA$AT)

Sup_GC_content<-ggplot(data = SupplementaryNucCountPRRA)+
  geom_histogram(aes(y=GCp,
                     fill = Fill))+
  facet_wrap(~Length2, scales = 'free_x', ncol =3)+
  scale_fill_manual(values = c("#00441b","#41ab5d", "#c7e9c0", "#3182bd","#bababa"))+
  ylab("% GC")+
  theme_bw()+
  theme(legend.position = "top")
Sup_GC_content

SupFigure1<-ggarrange(NucleotidesCounts_e,
                      Sup_GC_content,
                      ncol =1,
                      labels = c("a","b"))
SupFigure1
ggsave("../figures/Supplementary_Figure1.svg", plot= SupFigure1,
       width = 18, height = 18, dpi =300, units ="cm")

#########################################################
#########################################################
#Read RNA structure and transcriptome data

#Kim data from Kim 2020
kimdata<-read.csv("../data/external/Kim_et_al_2020_mmc2.csv", header = T, sep="\t")[,c(1:3)]
colnames(kimdata) <- c("five", "three", "count")

#calculate number of reads for each bin along genome
kimdata_melted<-melt(kimdata, id.vars = "count")


require(Hmisc)
intervals<-seq(0,29930, 100)

require(data.table)
groups_cut1 <- data.table(bin = levels(cut(kimdata_melted[, "value"],breaks=intervals, 
                                           labels = intervals[2:length(intervals)])))
groups_cut2<-groups_cut1
groups_cut1$variable<-rep("three", length(groups_cut1$bin))
groups_cut2$variable<-rep("five", length(groups_cut2$bin))
groups_cut<-rbind(groups_cut1, groups_cut2)

kimdata_melted<-transform(kimdata_melted, bin = cut(value, breaks = intervals,
                                                    labels = intervals[2:length(intervals)]))

##Sum up read counts for each bin
kimdata_bins_pre<-kimdata_melted %>% group_by(bin,variable) %>% summarise(totVal=sum(count))
kimdata_bins<-merge(kimdata_bins_pre, groups_cut, all = TRUE)
kimdata_bins$logcount<-log(kimdata_bins$totVal)
kimdata_bins[is.na(kimdata_bins)]<- 0
kimdata_bins$bin_num<-as.numeric(kimdata_bins$bin)*100

###Precalculated data from Huston et al., 2021

Structure_Huston<-read.csv("../data/external/SARS-CoV-2_Full_Length_Secondary_Structure_Map.ct",header =F, sep ="")[-1,]
Structure_Huston$stem<-ifelse(Structure_Huston$V5 != 0, "stem", "loop")

intcount<-0
intID<-c(0)
for (i in 2:nrow(Structure_Huston)) {
  if(Structure_Huston[i-1,7] == "loop") {
    if (Structure_Huston[i,7] == "loop") {
      intID<-c(intID,0)
    }
    else {
      intcount<-intcount +1
      intID<-c(intID,intcount)
    }
  }
  else {
    if (Structure_Huston[i,7] == "stem") {
      intID<-c(intID,intcount)
    }
    else {
      intID<-c(intID,0)
    }
  }
}
Structure_Huston$intID<-intID
####################

foldtop<-Structure_Huston[,c("V6","stem","intID")] %>% filter(intID >0) %>% group_by(intID) %>% top_n(1,V6)
foldbottom<-Structure_Huston[,c("V6","stem","intID")] %>% filter(intID >0) %>% group_by(intID) %>% top_n(-1,V6)
#
foldplotting<-merge(foldbottom, foldtop, by= c('intID','stem'))

##Plot both datasets together (Figure 1c)
RNAstrAndsgRNAs_together<-ggplot(data=foldplotting)+
  # geom_rect(aes(xmin=V6.x,
  #               xmax=V6.y,
  #               ymin=0,
  #               ymax= 16, 
  #               fill = ""),
  #           alpha =0.8,
  #           size =0.1)+
  geom_line(data=kimdata_bins,
            aes(x = bin_num,
                y = logcount,
                colour = variable),
            size =1)+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0,30000),breaks=seq(0,30000,2000), name ="")+
  scale_y_continuous(expand = c(0.01, 0), name = "log(read count)")+
  scale_colour_brewer(palette="Dark2", name = "Kim et al., 2020", labels = c("5' end", "3' end"))+
  #scale_fill_manual(values = c("#ffeda0"), name = "Huston et al., 2021" , labels = c("paired bases"))+
  theme_classic()
RNAstrAndsgRNAs_together

#####################
#Calculate correlation

#Bin all insertions
All_ins<-Insertions_plotting_merged_f[1:4]
All_ins<-transform(All_ins, bin = cut(Position, breaks = intervals,
                                      labels = intervals[2:length(intervals)]))
Ins_bins_pre<- All_ins %>% group_by(bin) %>% summarise(count = n())
Ins_bins<-merge(Ins_bins_pre, groups_cut1[,c(1)], all = TRUE)
Ins_bins[is.na(Ins_bins)]<- 0

#Kimdata is already in bins: kimdata_bins
kimdata_bins_three<-subset(kimdata_bins, kimdata_bins$variable == "three")

#Bin well-structured regions

#correlation
cor.test(Ins_bins$count,kimdata_bins_three$totVal)

##################
####Distance to the closest junction
AllJunctions<-unique(sort(kimdata_melted[,"value"]))
Insertions_All<-Insertions_plotting_merged_f[,1]

#Real data
ClosestJunction<-sapply(Insertions_All, function(Insertions_All, AllJunctions) {
  AllJunctions[which.min(abs(Insertions_All-AllJunctions))]}, AllJunctions)
DistanceRealData<-abs(Insertions_All-ClosestJunction)

#Simulated data
SimulationDistanceToJuction<-vector()
for (j in seq(1,1000)){
  InsPos<-floor(runif(length(Insertions_All),min=1, max = 29903))
  Closest<-sapply(InsPos, function(InsPos, AllJunctions) {
    AllJunctions[which.min(abs(InsPos-AllJunctions))]}, AllJunctions)
  Dist<-abs(InsPos-Closest)
  SimulationDistanceToJuction<-c(SimulationDistanceToJuction,Dist)
}

#pvalue
wilcox.test(DistanceRealData,SimulationDistanceToJuction)
#plot distance to closest junction
DistanceToJunction<-ggplot()+
  geom_histogram(aes(DistanceRealData,
                    after_stat(count*100/sum(count))),
                fill = "#41ae76",
                binwidth =5,
                size=1)+
  geom_histogram(aes(SimulationDistanceToJuction,
                    after_stat(count*100/sum(count))),
                fill = "#525252",
                binwidth = 5, alpha =0.7)+
  annotate("text", x=250, y=10.0,label = "p-value == 1.1^-{12}", parse =T)+
  scale_x_continuous(expand = c(0.01, 0),name ="distance to closest junction (nt)",
                     limits=c(-5,1000))+
  scale_y_continuous(expand = c(0.01, 0), name = "%")+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,0.5,0.3,0.2), "cm"))
DistanceToJunction

####Whether insertion is in stem or in loop

countStems<-function(InsList) {
  stemcount<-0
  for(elem in InsList) {
    checkdf<-subset(foldplotting,
                    foldplotting$V6.x <= elem &
                      foldplotting$V6.y >= elem)
    stemcount<-stemcount +nrow(checkdf)
  }
  return(stemcount)
}
#Real data in loops
StemsRealData<-countStems(Insertions_All)

#Simulation
SimulationStems<-vector()
for (j in seq(1,1000)){
  InsPos<-floor(runif(length(Insertions_All),min=1, max = 29903))
  SimulationStems<-c(SimulationStems, countStems(InsPos))
}

#p-value
pvalueStems<-length(SimulationStems[SimulationStems < StemsRealData])/1000
#plot
InsertionsInStemsPlot<-ggplot()+
  geom_histogram(aes(SimulationStems,
                 after_stat(count*100/sum(count))),
                 fill = "#525252",
                 binwidth = 1, alpha =0.5)+
  geom_vline(xintercept =StemsRealData, color = "#2b8cbe",size=2)+
  scale_x_continuous(expand = c(0.01, 0),name ="Number of insertions in stems")+
  scale_y_continuous(expand = c(0.01, 0), name = "%")+
  annotate("text", x=70, y=5.7,label = paste("p-value = ",pvalueStems))+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,0.2,0.3,0.2), "cm"))
InsertionsInStemsPlot

###Supplementary Figure 4
#Same as 1e and 1f but for 56 high-confident inserts
HighlyConfidentInserts<-subset(Insertions_plotting_merged_f, (Insertions_plotting_merged_f$Length <7 &
                                                                Insertions_plotting_merged_f$Confirmed_by_reads == "Y") |
                                 (Insertions_plotting_merged_f$Length > 8 &
                                    Insertions_plotting_merged_f$Fill != "Singleton"))
Insertions_All_HC<-HighlyConfidentInserts[,1]

ClosestJunction_HC<-sapply(Insertions_All_HC, function(Insertions_All_HC, AllJunctions) {
  AllJunctions[which.min(abs(Insertions_All_HC-AllJunctions))]}, AllJunctions)
DistanceRealData_HC<-abs(Insertions_All_HC-ClosestJunction)

#Simulation
SimulationDistanceToJuction_HC<-vector()
for (j in seq(1,1000)){
  InsPos<-floor(runif(length(Insertions_All_HC),min=1, max = 29903))
  Closest<-sapply(InsPos, function(InsPos, AllJunctions) {
    AllJunctions[which.min(abs(InsPos-AllJunctions))]}, AllJunctions)
  Dist<-abs(InsPos-Closest)
  SimulationDistanceToJuction_HC<-c(SimulationDistanceToJuction_HC,Dist)
}
wilcox.test(DistanceRealData_HC,SimulationDistanceToJuction_HC)
#plot distance to closest junction
DistanceToJunction_HC<-ggplot()+
  geom_histogram(aes(DistanceRealData_HC,
                     after_stat(count*100/sum(count))),
                 fill = "#41ae76",
                 binwidth =5,
                 size=1)+
  geom_histogram(aes(SimulationDistanceToJuction_HC,
                     after_stat(count*100/sum(count))),
                 fill = "#525252",
                 binwidth = 5, alpha =0.7)+
  annotate("text", x=250, y=10.0,label = "p-value == 8.7^-{5}", parse =T)+
  scale_x_continuous(expand = c(0.01, 0),name ="distance to closest junction (nt)",
                     limits=c(-5,1000))+
  scale_y_continuous(expand = c(0.01, 0), name = "%")+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,0.5,0.3,0.2), "cm"))
DistanceToJunction_HC

##
StemsRealData_HC<-countStems(Insertions_All_HC)

#Simulation
SimulationStems_HC<-vector()
for (j in seq(1,1000)){
  InsPos<-floor(runif(length(Insertions_All_HC),min=1, max = 29903))
  SimulationStems_HC<-c(SimulationStems_HC, countStems(InsPos))
}

#p-value
pvalueStems_HC<-length(SimulationStems_HC[SimulationStems_HC < StemsRealData_HC])/1000
#plot
InsertionsInStemsPlot_HC<-ggplot()+
  geom_histogram(aes(SimulationStems_HC,
                     after_stat(count*100/sum(count))),
                 fill = "#525252",
                 binwidth = 1, alpha =0.5)+
  geom_vline(xintercept =StemsRealData_HC, color = "#2b8cbe",size=2)+
  scale_x_continuous(expand = c(0.01, 0),name ="Number of insertions in stems")+
  scale_y_continuous(expand = c(0.01, 0), name = "%")+
  annotate("text", x=20, y=5.7,label = paste("p-value = ",pvalueStems_HC))+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        plot.margin=unit(c(0.3,0.2,0.3,0.2), "cm"))
InsertionsInStemsPlot_HC

SupplFig4<-ggarrange(DistanceToJunction_HC,
                     InsertionsInStemsPlot_HC,
                     ncol =2,
                     labels = c("a","b"))
SupplFig4

ggsave("../figures/Supplementary_Figure4.svg", plot= SupplFig4,
       width = 28, height = 15, dpi =300, units ="cm")
#################################
###Build whole Figure 1
#middle aligned panel
Figure1_bottom<-ggarrange(SARSCoVGenome,
                        ThreemulInsertionsPlot,
                        RNAstrAndsgRNAs_together,
                        ncol =1,
                        align = "v",
                        labels = c("c","","d"),
                        heights = c(0.55,1,0.4),
                        common.legend = F,
                        legend = "right")
#Build top and bottom parts
Figure1_top<-ggarrange(NULL,Insertions_length_distribution, NULL,
                       labels=c("a","",""),
                       ncol = 3,
                       widths = c(0.1,1,0.25))
Figure1_middle<-ggarrange(NULL,NucleotidesCounts,NULL,
                          ncol =3,
                          widths = c(0.3,1,0.3),
                          labels = c("","b",""))
Figure1_Simulations<-ggarrange(DistanceToJunction,
                               InsertionsInStemsPlot,
                               ncol =2,
                               widths = c(1, 0.6),
                               labels = c("e","f"))
#Combine all in one figure
Figure1_full<-ggarrange(Figure1_top,
                        Figure1_middle,
                        Figure1_bottom,
                        Figure1_Simulations,
                        ncol =1,
                        heights = c(0.3,0.25, 1, 0.3))
Figure1_full
#Save figure
ggsave("../figures/Figure1.svg", plot= Figure1_full,
       width = 26, height = 34.5, dpi =600, units ="cm")


###############################
#Investigation of insertion mechanisms
######

#For 25 long insertions map their origin where known

Long_insertions_verified<-subset(Long_insertions, Long_insertions$Fill != "Singleton")
Long_insertions_verified<-separate(data = Long_insertions_verified, col = Origin_position, into = c("Origin_position", "OP2"), sep = ";")

Mechanism_plotDf<-Long_insertions_verified %>% 
  group_by(Position) %>% dplyr::mutate(num = 1:n())

##Plot origin
Mechanism_plot<-ggplot(data=Mechanism_plotDf)+
  geom_vline(aes(xintercept = Position), linetype = "dashed", color = "#d9d9d9",
             alpha=0.5)+
  geom_point(aes(x = Position,
                 size =Length2,
                 y = num,
                 fill = Fill),
             shape = 25,
             color = "black",
             alpha =0.7)+
  scale_size(breaks = c(9,12,15,21),labels = c(9,12,15,21), range =c(5,10), name = "Length")+
  scale_fill_manual(values = c("#00441b","#41ab5d", "#c7e9c0"), name = "")+
  scale_colour_manual(values = c("#969696","#8c510a"), name = "Strand",
                      labels = c("direct", "compliment"))+
  scale_y_continuous(expand = c(0, 0), breaks=c(0,1,2,3),
                     limits = c (0.5,4.5))+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  geom_curve(data = subset(Mechanism_plotDf, Mechanism_plotDf$Origin_position !="") %>% 
               filter(as.numeric(Origin_position) < Position),
             aes(xend= Position, y = num-0.5, x = as.numeric(Origin_position), 
                 yend = num-0.5, color = Mechanism),
             curvature = -0.4,
             lineend = "square",
             arrow = arrow(length = unit(0.1,"cm")))+
  geom_curve(data = subset(Mechanism_plotDf, Mechanism_plotDf$Origin_position !="") %>% 
               filter(as.numeric(Origin_position) > Position),
             aes(xend= Position, y = num-0.5, x = as.numeric(Origin_position), 
                 yend = num-0.5, color = Mechanism),
             curvature = 0.4,
             lineend = "square",
             arrow = arrow(length = unit(0.1,"cm")))+
  theme_classic()+
  theme(#legend.position = "none",
    axis.title = element_blank(),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    plot.margin=unit(c(1,1,0,1), "cm"))

Mechanism_plot

##################################
##Combine with sgRNAs junction data from Kim et al., 2020

template_switch<-subset(Long_insertions_verified, Long_insertions_verified$Origin_position !="")[,c("Position", "Origin_position")]
temp1<-transform(template_switch, Position = as.numeric(Position))
template_switch_num<-transform(temp1, Origin_position = as.numeric(Origin_position))
template_switch_num$five<-with(template_switch_num, pmin(Position,Origin_position))
template_switch_num$three<-with(template_switch_num, pmax(Position,Origin_position))

nh<-100 #Length of neighbourhood to look at

template_switch_num$five_min<-template_switch_num$five - nh
template_switch_num$five_max<-template_switch_num$five + nh
template_switch_num$three_min<-template_switch_num$three - nh
template_switch_num$three_max<-template_switch_num$three + nh

Conditions_list<-as.list(as.data.frame(t(template_switch_num[,c(5:8)]),nrows =8))
#filter the dataset
supporting_reads<-subset(kimdata,
                         kimdata$five >= Conditions_list[[1]][1] &
                           kimdata$five <= Conditions_list[[1]][2] &
                           kimdata$three >= Conditions_list[[1]][3] &
                           kimdata$three <= Conditions_list[[1]][4])
originalDataConfirmed<-0
for (i in 2:length(Conditions_list)) {
  tmp<-subset(kimdata,
              kimdata$five >= Conditions_list[[i]][1] &
                kimdata$five <= Conditions_list[[i]][2] &
                kimdata$three >= Conditions_list[[i]][3] &
                kimdata$three <= Conditions_list[[i]][4])
  if(dim(tmp)[1] != 0)
  {
    originalDataConfirmed<-originalDataConfirmed+1
  }
  supporting_reads<-rbind(supporting_reads,tmp)
}
originalDataConfirmed

###Plot junctions information

supporting_reads_plot<-ggplot()+
  geom_curve(data=kimdata,
             aes(x= five, y = 0.1, 
                 xend = three,yend = 0.1, size = count),
             color ="grey",
             alpha =0.5,
             curvature = 0.3)+
  scale_size(range = c(0.1,2))+
  geom_curve(data=supporting_reads,
             aes(x= five, y = 0.1, 
                 xend = three,yend = 0.1),
             color = "#54278f",
             size =0.5,
             curvature = 0.3)+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0,30000),breaks=seq(0,30000,2000))+
  scale_y_continuous(expand = c(0.01, 0), limits = c(-0.5,0.1))+
  theme_classic()+
  theme(#legend.position = "none",
    axis.title = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(color = "white", size =10),
    axis.ticks.y = element_line(color = "white"),
    plot.margin=unit(c(0,1,1,1), "cm"))
supporting_reads_plot

###########
#Permutation test
###########
Original_data<-template_switch_num[,c(1,2)]
Original_data$experiment<-rep(0, length(Original_data$Position))

##
#Scenario 1: only origin position is selected randomly

results<-vector()
for (j in seq(1,1000)){
  Origin_position<-floor(runif(8,min=1, max = 29903))
  tmp<-data.frame(Position = Original_data$Position,
                  Origin_position,
                  experiment = rep(j, length(Origin_position)))

  tmp$five<-with(tmp, pmin(Position,Origin_position))
  tmp$three<-with(tmp, pmax(Position,Origin_position))
  tmp$five_min<-tmp$five - nh
  tmp$five_max<-tmp$five + nh
  tmp$three_min<-tmp$three - nh
  tmp$three_max<-tmp$three + nh

  tmp_list<-as.list(as.data.frame(t(tmp[,c(6:9)]),nrows =8))
  check <-0
  for (i in 1:length(Conditions_list)) {
    checkdf<-subset(kimdata,
                kimdata$five >= tmp_list[[i]][1] &
                  kimdata$five <= tmp_list[[i]][2] &
                  kimdata$three >= tmp_list[[i]][3] &
                  kimdata$three <= tmp_list[[i]][4])
    if (dim(checkdf)[1] != 0)
    {
      check<-check+1
    }
  }
  results<-append(results,check)

}
results

##
#Scenario 2: both origin position and insertion position are selected randomly

results_both<-vector()
for (j in seq(1,1000)){
  Origin_position<-floor(runif(8,min=1, max = 29903))
  Position<-floor(runif(8,min=1, max = 29903))
  tmp<-data.frame(Position,
                  Origin_position)

  tmp$five<-with(tmp, pmin(Position,Origin_position))
  tmp$three<-with(tmp, pmax(Position,Origin_position))
  tmp$five_min<-tmp$five - nh
  tmp$five_max<-tmp$five + nh
  tmp$three_min<-tmp$three - nh
  tmp$three_max<-tmp$three + nh

  tmp_list<-as.list(as.data.frame(t(tmp[,c(5:8)]),nrows =8))
  check <-0
  for (i in 1:length(Conditions_list)) {
    checkdf<-subset(kimdata,
                    kimdata$five >= tmp_list[[i]][1] &
                      kimdata$five <= tmp_list[[i]][2] &
                      kimdata$three >= tmp_list[[i]][3] &
                      kimdata$three <= tmp_list[[i]][4])
    if (dim(checkdf)[1] != 0)
    {
      check<-check+1
    }
  }
  results_both<-append(results_both,check)

}
results_both

##########
#Plot permutations
pvalue_one_end <- length(results[results == 2])/1000
pvalue_both_ends <- length(results[results_both == 2])/1000

#Plot scenario 1
Permutations_one<- ggplot() + 
  aes(results)+ 
  geom_histogram(binwidth=1, colour="black",
                 fill = "#08306b")+
  geom_vline(xintercept =2, color = "#980043",
             size =2)+
  annotate("text", x = 1.2 , y =750, label = paste("p-value = ",pvalue_one_end))+
  theme_classic()+
  ggtitle("1000 permutations of origin site")+
  theme(axis.text = element_text(size =10),
        title = element_text(size =10),
        axis.title.x= element_blank())

#Plot scenario 2
Permutations_both<-ggplot() + 
  aes(results_both)+ 
  geom_histogram(binwidth=1, colour="black",
                 fill = "#08306b")+
  geom_vline(xintercept =2, color = "#980043",
             size =2)+
  ggtitle("1000 permutations of both sites")+
  annotate("text", x = 1.2 , y =750, label = paste("p-value = ",pvalue_both_ends))+
  theme_classic()+
  theme(axis.text = element_text(size =10),
        title = element_text(size =10),
        axis.title.x= element_blank())

#############
#Build Figure 2

Figure2_top<-ggarrange(Mechanism_plot,
                       SARSCoVGenome_plain,
                       supporting_reads_plot,
                       ncol =1,
                       align = "v",
                       heights = c(0.5,0.1,0.55),
                       common.legend = T,
                       legend = "top")
Figure2_bottom<-ggarrange(NULL, Permutations_one, NULL, Permutations_both, NULL,
                                widths = c(0.3,1, 0.25, 1,0.3),
                                ncol = 5,
                                labels = c("","b","","c",""))
Figure2_full<-ggarrange(Figure2_top,
          Figure2_bottom,
          labels=c("a",""),
          heights = c(1,0.3),
          ncol =1)
Figure2_full

#Save Figure
ggsave("../figures/Figure2.svg", plot= Figure2_full,
       width = 30, height = 20, dpi =300, units ="cm")



