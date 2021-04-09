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
PhylogenyInsertions<-read.csv("../data/Insertions_phylogeny_long.csv", sep ="\t", header = T,
                              stringsAsFactors = F)

#Merge dataframes and add statuses
Insertions_plotting_merged<-merge(Threemulpl_insertions,PhylogenyInsertions, all=T)

Insertions_plotting_merged[is.na(Insertions_plotting_merged)]<-""
Insertions_plotting_merged$Fill<-ifelse(Insertions_plotting_merged$Confirmed_by_reads_or_multiple == "Y",
                                        "Confirmed by reads", 
                                        ifelse(Insertions_plotting_merged$USHER == "Yes",
                                               "Monophyletic",
                                               ifelse(Insertions_plotting_merged$fill == "Multi", "Multiple genomes", "Singleton")))
Insertions_plotting_merged$Length2<-ifelse(Insertions_plotting_merged$Length < 30, 
                                           Insertions_plotting_merged$Length, 35)

##############Plot insertions

#Plot length distribution (Figure 1a)
Insertions_length_distribution<- ggplot(data=Insertions_plotting_merged)+
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
ThreemulInsertionsPlot<-ggplot(data=Insertions_plotting_merged)+
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
Three_insertions<-subset(Insertions_plotting_merged, Insertions_plotting_merged$Length < 4)
#Subset inserts of the length 6
Six_insertions<-subset(Insertions_plotting_merged, Insertions_plotting_merged$Length == 6)
##################################
#Subset long insertions
Long_insertions<-subset(Insertions_plotting_merged, Insertions_plotting_merged$Length > 8)
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
ATGC_countsDF<-as.data.frame(t(data.frame(Three_ins_nuc_counts,Six_ins_nuc_counts,Long_ins_nuc_counts, Whole)))

data<-ATGC_countsDF[,2:5]
ATGC_portionsDF<-as.data.frame(t(apply(data,1, function(x) x*100/sum(x))))
ATGC_portionsDF$length<-factor(c("3 nucl","6 nucl","> 8 nucl","genome"),
                               levels = c("3 nucl","6 nucl","> 8 nucl","genome"))
ATGC_portion_melted<-melt(ATGC_portionsDF, id = "length")

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

WellFolded<-read.csv("../data/external/well_folded_regions_Huston_v2.csv",header =T, sep ='\t')

foldtop<-WellFolded[,c("X","intervals","idx")] %>% group_by(intervals) %>% top_n(1,X)
foldbottom<-WellFolded[,c("X","intervals","idx")] %>% group_by(intervals) %>% top_n(-1,X)
#
foldplotting<-merge(foldbottom, foldtop, by= c('intervals','idx'))

##Plot both datasets together (Figure 1c)
RNAstrAndsgRNAs_together<-ggplot(data=foldplotting)+
  geom_rect(aes(xmin=X.x,
                xmax=X.y,
                ymin=0,
                ymax= 16, 
                fill = ""),
            alpha =0.8,
            size =0.1)+
  geom_line(data=kimdata_bins,
            aes(x = bin_num,
                y = logcount,
                colour = variable),
            size =1)+
  scale_x_continuous(expand = c(0.01, 0),limits = c(0,30000),breaks=seq(0,30000,2000), name ="")+
  scale_y_continuous(expand = c(0.01, 0), name = "log(read count)")+
  scale_colour_brewer(palette="Dark2", name = "Kim et al., 2020", labels = c("5' end", "3' end"))+
  scale_fill_manual(values = c("#ffeda0"), name = "Huston et al., 2021" , labels = c("well-folded region"))+
  theme_classic()
RNAstrAndsgRNAs_together

#####################
#Calculate correlation

#Bin all insertions
All_ins<-Insertions_plotting_merged[1:4]
All_ins<-transform(All_ins, bin = cut(Position, breaks = intervals,
                                      labels = intervals[2:length(intervals)]))
Ins_bins_pre<- All_ins %>% group_by(bin) %>% summarise(count = n())
Ins_bins<-merge(Ins_bins_pre, groups_cut1[,c(1)], all = TRUE)
Ins_bins[is.na(Ins_bins)]<- 0

#Kimdata is already in bins: kimdata_bins
kimdata_bins_three<-subset(kimdata_bins, kimdata_bins$variable == "three")

#correlation
cor.test(Ins_bins$count,kimdata_bins_three$totVal)

#################################
###Build whole Figure 1
#middle aligned panel
Figure1_middle<-ggarrange(SARSCoVGenome,
                        ThreemulInsertionsPlot,
                        RNAstrAndsgRNAs_together,
                        ncol =1,
                        align = "v",
                        labels = c("b","","c"),
                        heights = c(0.5,1,0.4),
                        common.legend = F,
                        legend = "right")
#Build top and bottom parts
Figure1_top<-ggarrange(NULL,Insertions_length_distribution, NULL,
                       labels=c("a","",""),
                       ncol = 3,
                       widths = c(0.1,1,0.25))
Figure1_bottom<-ggarrange(NULL,NucleotidesCounts, NULL,
                          ncol =3,
                          widths = c(0.5,1,0.5),
                          labels = c("","d",""))
#Combine all in one figure
Figure1_full<-ggarrange(Figure1_top,
                        Figure1_middle,
                        Figure1_bottom,
                        ncol =1,
                        heights = c(0.3,1,0.25))
Figure1_full
#Save figure
ggsave("../figures/Figure1.svg", plot= Figure1_full,
       width = 26, height = 31, dpi =300, units ="cm")


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
                 yend = num-0.5, color = Mechanism.Igor.),
             curvature = -0.4,
             lineend = "square",
             arrow = arrow(length = unit(0.1,"cm")))+
  geom_curve(data = subset(Mechanism_plotDf, Mechanism_plotDf$Origin_position !="") %>% 
               filter(as.numeric(Origin_position) > Position),
             aes(xend= Position, y = num-0.5, x = as.numeric(Origin_position), 
                 yend = num-0.5, color = Mechanism.Igor.),
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
#Build Figure 3

Figure3_top<-ggarrange(Mechanism_plot,
                       SARSCoVGenome_plain,
                       supporting_reads_plot,
                       ncol =1,
                       align = "v",
                       heights = c(0.5,0.1,0.55),
                       common.legend = T,
                       legend = "top")
Figure3_bottom<-ggarrange(NULL, Permutations_one, NULL, Permutations_both, NULL,
                                widths = c(0.3,1, 0.25, 1,0.3),
                                ncol = 5,
                                labels = c("","b","","c",""))
Figure3_full<-ggarrange(Figure3_top,
          Figure3_bottom,
          labels=c("a",""),
          heights = c(1,0.3),
          ncol =1)
Figure3_full

#Save Figure
ggsave("../figures/Figure3.svg", plot= Figure3_full,
       width = 30, height = 20, dpi =300, units ="cm")



