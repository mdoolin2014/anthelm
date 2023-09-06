#Script by M. Doolin for 16S microbiome data from demultiplexed, paired-end reads.
#The purpose of this script is to run 16S data from the Anthelmintic project testing
# the impacts of antiparasitic drugs on the microbiome. 

# Updated from March 2022 GTDB run, now running in May 2023. 
# Adapted from previous versions from myself, Dylan Klure (2019) and from the tutorials: 
# https://benjjneb.github.io/dada2/tutorial.html
# https://joey711.github.io/phyloseq/
# Comments are my own.

#SessionInfo 5May2023 (and still true 2Aug23):
##R version 4.3.0 (2023-04-21)
##Platform: aarch64-apple-darwin20 (64-bit)
##Running under: macOS Monterey 12.6.1


#Imports created using qiime2-2022.11 on the Utah CHPC.
# Used deblur denoising, classified with V4 region of silva.



##### Important packages #####
#required packages for running the whole 16S DADA2 pipeline: 
library(Rcpp)
library(ggplot2)
library(ape)
library(gridExtra)
library(SummarizedExperiment)
library(dada2) 
library(S4Vectors)
library(stats4)
library(IRanges)
library(XVector)
library(RSQLite)
library(ShortRead)
library(Biostrings)
library(phyloseq)
library(microbiome)
library(RColorBrewer)
library(ggbreak)
library(ggbeeswarm)
library(reshape2)


##### Import qiime artifacts into R to make a physeq object #####
library(qiime2R) #package version 0.99.6


# reads, as well as singletons and doubletons, in qiime. Can go ahead and make ps1.
ps1 <- qza_to_phyloseq(features="qiime_outputs/final-table.qza", 
                       taxonomy="qiime_outputs/taxonomy.qza", 
                       tree ="qiime_outputs/tree_root.qza",
                       metadata="~/Documents/Anthelm_16S/metadata/AnthelmMetadatBySample_Novaseq.tsv")

ps1 #This has all mouse samples, the mock community, and the blanks.
##otu_table()   OTU Table:         [ 7289 taxa and 83 samples ]
##sample_data() Sample Data:       [ 83 samples by 16 sample variables ]
##tax_table()   Taxonomy Table:    [ 7289 taxa by 7 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 7289 tips and 7103 internal nodes ]

View(sample_data(ps1))


#Adding a few things to the ps1 sample data. 
ps1@sam_data[["sum"]] <- sample_sums(ps1)
sam <- data.frame(rownames(sample_data(ps1)))
ps1@sam_data[["sample"]] <- rownames(sample_data(ps1))
ASVs <- estimate_richness(ps1, measures="Observed")
sample_data(ps1) <- cbind(sample_data(ps1), ASVs)
ps1@sam_data[["treatment"]] <- factor(ps1@sam_data[["treatment"]], levels=c("control", "ALB", "MTZ"))
ps1@sam_data[["TD"]] <- factor(ps1@sam_data[["TD"]], levels=c("control_D0", "control_D7", "control_D21", 
                                                   "ALB_D0", "ALB_D7", "ALB_D21", 
                                                   "MTZ_D0", "MTZ_D7", "MTZ_D21"))
ps1@sam_data[["day"]] <- factor(ps1@sam_data[["day"]], levels=c("D0", "D7", "D21"))
View(sample_data(ps1))


#   ##### Rarefied physeq objects #####
#Rarefy.
psrare <- rarefy_even_depth(psexp, rngseed=10, sample.size=20000)
##phyloseq-class experiment-level object
##otu_table()   OTU Table:         [ 798 taxa and 72 samples ]
##sample_data() Sample Data:       [ 72 samples by 19 sample variables ]
##tax_table()   Taxonomy Table:    [ 798 taxa by 7 taxonomic ranks ]
##phy_tree()    Phylogenetic Tree: [ 798 tips and 794 internal nodes ]


##### Modeling BODY WEIGHT over time #####
#   ##### All samples, with lmer #####

library(lme4)
library(lmerTest)

dfrare$day_num <- as.numeric(dfrare$day_num)

#Create the models for all 3 metrics you want. 
wtmod <- lmer(body_wt ~ treatment + day_num + sex + treatment*day_num + (1|mouse), data=dfrare, REML=TRUE) 
wtmod
anova(wtmod)
ranova(wtmod)


#   ##### Comparing D0 body wt #####

mod <- aov(body_wt ~ treatment, data=dfrareD0)
summary(mod)


#.  ##### Weight plots #####

#.      ##### Plotting D0 weights. #####
D0wt <- ggplot(aes(x=treatment, y=body_wt, color=treatment), data=dfrareD0) +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(aes(colour=treatment), alpha=0.9, size=10, cex=4) +
  theme_bw() +
  scale_color_manual(values=c("#367F6A", "#C15E37", "#506796")) +
  theme(text=element_text(size=20), legend.position="none", 
        axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("D0 body mass, Anthelm")) + 
  ylab("Body mass (g)") + 
  xlab("Treatment Group") +
  ylim(17.5,30)

#
#.      ##### Plot change in weights between pretreat and recovery. #####
wt_spread$treatment <- factor(wt_spread$treatment, levels=c("control", "ALB", "MTZ"))
p1 <- ggplot(aes(x=treatment, y=prerecov, colour=treatment), data=wt_spread) +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(aes(colour=treatment), alpha=0.9, size=10, cex=4) +
  #geom_text(mapping=aes(label=mouse.y), hjust=-0.2, vjust=0) +
  scale_color_manual(values=c("#367F6A", "#C15E37", "#506796")) +
  theme_bw() +
  theme(text=element_text(size=20), legend.position="none", axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("Change in body mass pretreatment to recovery")) + 
  ylab("Change in Body Mass (g)") + 
  xlab("Treatment")

#

##### Alpha diversity#####

#      ##### ggplot diversity plotting #####

#         ##### plotting D0vsD7 #####

dfrareD07$TD <- factor(dfrareD07$TD, levels=c("control_D0", "control_D7", "ALB_D0", "ALB_D7",  "MTZ_D0", "MTZ_D7"))
min(dfrareD07$Observed) 
max(dfrareD07$Observed)
#Plot with ggplot. Observed range is 275-409
asvplot <- ggplot(dfrareD07, aes(x=TD, y=Observed, colour=TD)) +
  theme_bw() +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(aes(colour=TD), alpha=0.9, size=8, cex=4) +
  #geom_text(mapping = aes(label = mouse), size = 4, vjust = 2) +  
  scale_color_manual(values=c("#6bc0a7", "#367F6A", "#f98e63", "#C15E37", "#8ba0ca", "#506796")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5), legend.position="none",
        text=element_text(size=15)) +
  ggtitle(paste0("Anthelm 16S ASV richness D0 D7")) + 
  ylab("Observed ASVs") + 
  xlab("Treatment") +
  ylim(270,440) +
  facet_wrap(~treatment, scales="free_x")
#

#         ##### And to just look at D21 by itself.  #####
min(dfrareD21$Observed) 
max(dfrareD21$Observed) #range 303-403
D21asv <- ggplot(dfrareD21, aes(x=treatment, y=Observed, colour=treatment)) +
  theme_bw() +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(aes(colour=treatment), alpha=0.9, size=10, cex=4) +
  scale_color_manual(values=c("#1c4539", "#70351e", "#2c3a57")) +
  theme(legend.position="none", 
        text=element_text(size=24)) +
  ggtitle(paste0("Anthelm 16S ASV richness D21")) + 
  ylab("Observed ASVs") + 
  xlab("Treatment") +
  ylim(300, 420)



#   ##### Mixed Effects Models for repeated measures to include treatment, day, and individual. #####

#      ##### Mixed effects just including individual over time within each group #####
#https://biostats.w.uib.no/post-hoc-tests-multiple-comparisons-in-linear-mixed-effect-models/
library(lme4)
library(lmerTest) #this has the random effect quantification. 

#Create the models for all 3 metrics you want. 
obsmod <- lmer(Observed ~ day + (1|mouse), data=dfrarecont) 
shanmod <- lmer(Shannon ~ day + (1|mouse), data=dfrarecont) 
anova(shanmod)
anova(obsmod)

#   ##### Plotting change in alpha diversity between time points #####

#      ##### Looking at alpha diversity pre-to-posttreatment in individuals #####

obs <- ggplot(data=dfrare, aes(x=day, y=Observed, group=mouse, color=treatment)) +
  geom_line() + 
  scale_color_brewer(palette="Set2") +
  geom_point(aes(size=8)) + ##facet_wrap(facets="cage") +
  ggtitle("ASV Richness by individual") +
  facet_grid(~treatment) +
  theme_bw() +
  theme(legend.position="none", text=element_text(size=15))




##### Beta diversity#####

#   ##### Betadisper #####
library(vegan) #version 2.5-7

#Before deciding which beta div analysis, run a dispersion test to see if appropriate for parametrics
x <- distance(psraremtz, "bray")
y <- as(sample_data(psraremtz), "data.frame")
groups <- y[["day"]]
mod <- betadisper(x, groups)
mod
anova(mod)

plot(mod)  
plot(TukeyHSD(mod))


#   ##### Distance metrics and PCoA plots #####
pstrans1 <- transform_sample_counts(psrareD07, function(otu) otu/sum(otu))
#Run distance metrics.
pstrans1 <- transform_sample_counts(psrareD0, function(otu) otu/sum(otu)) 
pstrans2 <- transform_sample_counts(psrareD7, function(otu) otu/sum(otu)) 
pstrans3 <- transform_sample_counts(psrareD21, function(otu) otu/sum(otu)) 
#Bray ordinations
ps.ordi1 <- ordinate(pstrans1, "PCoA", "bray")
ps.ordi2 <- ordinate(pstrans2, "PCoA", "bray")
ps.ordi3 <- ordinate(pstrans3, "PCoA", "bray")


#      ##### PCOA of whole dataset #####

brayplot <- plot_ordination(psrare, ps.ordi1) +
  geom_point(aes(color=TD, shape=treatment), size =10) + 
  theme_bw() +
  theme(text=element_text(size=20), legend.position="none") + 
  scale_shape_manual(values=c(16,17,15)) +
  scale_color_manual(values=c("#6bc0a7", "#367F6A", "#f98e63","#C15E37", "#8ba0ca", "#506796")) +
  #scale_color_manual(values=c("#6bc0a7", "#367F6A","#1c4539", "#f98e63","#C15E37","#70351e", "#8ba0ca", "#506796", "#2c3a57")) +
  #scale_fill_manual(palette="Greys") +
  #geom_text(mapping = aes(label = mouse), size = 4, vjust = 2) +
  ggtitle("Bray PCoA -- psrareD07") +  
  facet_wrap(~day) +
  stat_ellipse(aes(color=TD), type="t", level=0.9)



#      ##### Plot PCoAs within a treatment #####

pstrans1 <- transform_sample_counts(psrarectrlD07, function(otu) otu/sum(otu)) 
pstrans2 <- transform_sample_counts(psrarealbD07, function(otu) otu/sum(otu)) 
pstrans3 <- transform_sample_counts(psraremtzD07, function(otu) otu/sum(otu)) 
#Bray ordinations
ps.ordi1 <- ordinate(pstrans1, "PCoA", "bray")
ps.ordi2 <- ordinate(pstrans2, "PCoA", "bray")
ps.ordi3 <- ordinate(pstrans3, "PCoA", "bray")

#         ##### D0vsD7 PCOA within a treatment #####

ctrlplot <- plot_ordination(psrarectrlD07, ps.ordi1) +
  geom_point(aes(color=day), size =10) + 
  theme_bw() +
  theme(text=element_text(size=20)) + 
  #scale_shape_manual(values=c(21,22,24)) +
  scale_color_manual(values=c("#6bc0a7", "#367F6A")) +
  #scale_fill_manual(palette="Greys") +
  #geom_text(mapping = aes(label = mouse), size = 4, vjust = 2) +
  ggtitle("jaccard PCoA -- psrarectrl") +
  stat_ellipse(aes(color=day), type="t", level=0.9)


#   ##### Statistical tests on beta diversity #####
library(pairwiseAdonis)

#      ##### PERMANOVAs on psrare for day and treatment #####
metadata <- as(sample_data(psrare), "data.frame") #make metadata into a df.
perms <- with(metadata, how(nperm = 5000, blocks = cage)) #Note increasing number of permutations makes results more consistent.
d <- phyloseq::distance(psrare, method="wunifrac")
mod <- adonis2(d ~ treatment + day + treatment * day, data = metadata, permutations = perms)
mod

#Gives pairwise outputs results
metadata <- as(sample_data(psrare), "data.frame") 
pair.mod<-pairwise.adonis2(d ~ treatment + day + treatment * day, data=metadata, strata="cage", nperm=5000)
pair.mod

##mod <- adonis2(phyloseq::distance(psrare, method="wunifrac") ~ treatment + day + treatment * day, data = metadata, permutations = perms)
##               Df  SumOfSqs        R2        F   Pr(>F)    
##treatment      2  0.09619 0.03121  2.9415 0.0002 ***
##day            2  1.85295 0.60114 56.6623 0.0002 ***
##treatment:day  4  0.10313 0.03346  1.5768 0.1548    
##Residual      63  1.03010 0.33419                   
##Total         71  3.08237 1.00000 

##mod <- adonis2(phyloseq::distance(psrare, method="unifrac") ~ treatment + day + treatment * day, data = metadata, permutations = perms)
##             Df SumOfSqs      R2      F Pr(>F)  
##treatment      2   0.1412 0.04459 1.7053 0.0002 ***
##day            2   0.2261 0.07143 2.7319 0.0002 ***
##treatment:day  4   0.1911 0.06035 1.1540 0.0122 *  
##Residual      63   2.6075 0.82363                  
##Total         71   3.1659 1.00000   

##mod <- adonis2(phyloseq::distance(psrare, method="bray") ~ treatment + day + treatment * day, data = metadata, permutations = perms)
##              Df SumOfSqs      R2      F Pr(>F)  
##treatment      2   0.4917 0.05512  2.8224 0.00020 ***
##day            2   2.4928 0.27941 14.3083 0.00020 ***
##treatment:day  4   0.4492 0.05035  1.2892 0.08458 .  
##Residual      63   5.4879 0.61512                    
##Total         71   8.9216 1.00000      

##mod <- adonis2(phyloseq::distance(psrare, method="jaccard") ~ treatment + day + treatment * day, data = metadata, permutations = perms)
##              Df SumOfSqs      R2      F Pr(>F)    
##treatment      2   0.7624 0.05018 2.3236 0.0002 ***
##day            2   3.1884 0.20987 9.7174 0.0002 ***
##treatment:day  4   0.9058 0.05962 1.3803 0.0218 *  
##Residual      63  10.3357 0.68032                  
##Total         71  15.1923 1.00000      


#      ##### PERMANOVAs on other ps within each treatment #####

metadata <- as(sample_data(psrarealb), "data.frame") #make metadata into a df.
perms <- with(metadata, how(nperm = 5000, blocks = cage)) 
mod <- adonis2(phyloseq::distance(psrarealb.core, method="wunifrac") ~ day, data = metadata,
               permutations = perms)
mod

#Gives pairwise outputs 
metadata <- as(sample_data(psraremtz), "data.frame") 
d <- phyloseq::distance(psraremtz, method="jaccard")
pair.mod<-pairwise.adonis2(d ~ day, data=metadata, strata="cage", nperm=5000)
pair.mod

#Testing for significance across days, within each of the treatments. 
#psrarectrl by day, blocked by cage. 
#Unifrac (R2=0.11188 F=1.3228  p=0.002**) with significant diffs D0-D7 (p=0.026), D7-D21 (p=0.002), D0-D21 (p<0.002)
#WUnifrac (R2=0.67862 F=22.171  p=2e-04***) with significant diffs D7-D21 (p=0.003), D0-D21 (p=0.010), not D0-D7.
#Bray (R2=0.31146, F=4.7496, p=2e-04***) with significant diffs D0-D7 (p=0.035), D7-D21 (p<0.002), D0-D21 (p<0.002).
#Jaccard (R2=0.25044 F=3.5083  p=2e-04***) also significant diffs D0-D7 (p=0.033), D7-D21 (p<0.002), D0-D21 (p<0.002).

#psrarealb by day, blocked by cage.
#Unifrac (R2=0.10986 F=1.2958 p=0.008598**) significant diffs D7-D21, D0-D21, not D0-D7.
#WUnifrac (R2=0.61809, F=16.993, p=2e-04***) significant diffs D7-D21, D0-D21, not D0-D7.
#Bray (R2=0.32051 F=4.9527 p=2e-04***) significant diffs D7-D21, D0-D21, not D0-D7.
#Jaccard (R2=0.25817 F=3.6542 p=2e-04***) significant diffs D7-D21, D0-D21, not D0-D7.

#psraremtz
#Unifrac (R2=0.18516 F=2.386  p=2e-04***) with signif diffs D0-D7 (p=0.004), D0-D21 (p=0.0004), D7-D21 (p=0.001). 
#WUnifrac (R2=0.65454 F=19.894 p=2e-04***) with signif diffs D0-D21 (p=0.004), D7-D21 (p=0.002), but not D0-D7. 
#Bray (R2=0.41673 F=7.5018 p=2e-04 ***) with signif diffs D0-D7 (p=0.004), D0-D21 (p=0.003), D7-D21 (p=0.001). 
#Jaccard (R2=0.34372 F=5.4992 p=2e-04 ***) with signif diffs D0-D7 (p=0.004), D0-D21 (p=0.004), D7-D21 (p=0.001). 


#Diversity metrics for core microbial communities:
#psrarectrl.core is significant for day for all metrics except for unifrac (p=0.122)
#psrarealb.core is significant for day for all metrics except for unifrac (p=0.4141)
#psraremtz.core is significant for day for all metrics except unifrac (p=0.771) 


#      ##### PERMANOVAs on other ps within each day without cage #####

#Without including cage, just testing the effect of treatment group at each dat. 

metadata <- as(sample_data(psrareD21), "data.frame") #make metadata into a df.
perms <- with(metadata, how(nperm = 1000)) 
mod <- adonis2(phyloseq::distance(psrareD21, method="wunifrac") ~ treatment, data = metadata,
               permutations = perms)
mod

#psrareD0 by treatment. 
#Unifrac (R2=0.06778 F=0.7634 p=0.9451)  
#WUnifrac (R2=0.06903 F=0.7785 0.6064)
#Bray (R2=0.08101 F=0.9256 p=0.5035)

#psrareD7 by treatment. 
#Unifrac (R2=0.15568 F=1.936 p=0.002997) 
#WUnifrac (R2=0.2138 F=2.8554 p=0.007992)
#Bray (R2=0.21731 F=2.9153 p=0.000999) 

#psrareD21 by treatment. 
#Unifrac (R2=0.11619 F=1.3803 p=0.03996) 
#WUnifrac (R2=0.2096 F=2.7844 p=0.02597)
#Bray (R2=0.15586 F=1.9387 p=0.03397) 
#Df 2,23

#Which are different?
library(pairwiseAdonis)
#Gives pairwise outputs results
metadata <- as(sample_data(psrareD7), "data.frame") 
d <- phyloseq::distance(psrareD7, method="wunifrac")
pair.mod<-pairwise.adonis2(d ~ treatment, data=metadata, nperm=1000)
pair.mod

#D7 MTZ different from the other two based on Bray, Unifrac and WUnifrac. 
#D21, MTZ differnt from the other two based on Unifrac. ALB diff from control in wunifrac and diff
# from both other groups based on Bray. 

#
#   ##### Pulling out distances #####

#Pulling out beta diversity distances. 
#Using this tutorial https://www.biostars.org/p/452069/
library(reshape2)

#Run the code to make the distance matrix 
U<-phyloseq::distance(psrare, "unifrac")

#pull the distances into usable form. 
U.m <- melt(as.matrix(U))
U.m = U.m %>% #Takes out the same exact sample compared to itself, which all give 0 distance. 
  filter(as.character(Var1) != as.character(Var2)) %>%
  mutate_if(is.factor, as.character)
sd = dfrare %>%
  select(sample, day_num, mouse, treatment) %>% #If adding mouse back in, need to add here too.
  mutate_if(is.factor,as.character)
colnames(sd) = c("Var1", "day_num", "mouse", "treatment") #and also need to add here.
U.sd = left_join(U.m, sd, by = "Var1")
colnames(sd) = c("Var2", "day_num", "mouse", "treatment")
U.sd = left_join(U.sd, sd, by = "Var2")
U.sd <- subset(U.sd, U.sd$mouse.x==U.sd$mouse.y)

U.sd$treatment.x <- factor(U.sd$treatment.x, levels=c("control","ALB","MTZ"))
U.sd$day_num.y <- as.numeric(U.sd$day_num.y)

U.sd.D0D7 <- subset(U.sd, U.sd$day_num.x=="0" & U.sd$day_num.y=="7")
U.sd.D7D21 <- subset(U.sd, U.sd$day_num.x=="7" & U.sd$day_num.y=="21")
U.sd.D0D21 <- subset(U.sd, U.sd$day_num.x=="0" & U.sd$day_num.y=="21")

U.bind <- rbind(U.sd.D0D21, U.sd.D0D7, U.sd.D7D21)
U.bind$diff <- paste(U.bind$day_num.x, U.bind$day_num.y, sep=" vs ")
nrow(U.bind)

summary(subU.sd)



#   ##### Comparing beta diversity change across treatment groups #####
#      ##### Plotting distance change grouped by treatments excluding D0-D21 #####
punif <-  ggplot(U.sd.D0721, aes(x = diff, y = value, colour=treatment.x)) +
  theme_bw() + 
  theme(text=element_text(size=15), legend.position="none") +
  geom_boxplot(outlier.shape=NA) +
  geom_beeswarm(aes(colour=treatment.x), alpha=0.9, size=6, cex=3) +
  #geom_text(mapping=aes(label=mouse.y), hjust=-0.2, vjust=0) +
  scale_color_manual(values=c("#367F6A", "#C15E37", "#506796")) +
  theme(axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle(paste0("Distance Metric = Unifrac, Anthelm psrare")) + 
  ylab("Unifrac distance") + 
  xlab("treatment") +
  facet_grid(~treatment.x)


#   ##### Kruskal-Wallace on D0 vs. D7 to look at whether there are differences by group #####
head(J.sd.D0D7)

TC.U <- kruskal.test(value~treatment.y, data=U.sd.D0D7)

TC.U
#


##### DeSeq2 for 16S Differential Abundance #####

library(phyloseq)
library(DESeq2) 
library(microbiome)
library(ashr)
library(ape)


#   ##### Run DESeq #####

#Choose your physeq object:
obj <- psctrlD021

dsps = phyloseq_to_deseq2(obj, ~ day)
dsps = DESeq(dsps, test="Wald", fitType="parametric", sfType="poscounts")
res = results(dsps, cooksCutoff = FALSE)
resultsNames(dsps)   
#Add correct coefficient based on resultsNames
resLFC <- lfcShrink(dsps, coef="day_D21_vs_D0", type="ashr") 
summary(resLFC)

#Create a summary table with whatever p-adj you want (using alpha). Can change to see more results. 
resOrdered <- resLFC[order(resLFC$log2FoldChange),]

#All outcomes, for volcano plots. 
deseq.all=cbind(as(resOrdered, "data.frame"), as(tax_table(obj)[rownames(resOrdered), ], "matrix"))


#      ##### Grab signif results and write csv. #####

#Now to only grab the significant results. 
alpha = 0.05
sigtab = resOrdered[which(resOrdered$padj < alpha), ]
sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(obj)[rownames(sigtab), ], "matrix"))
sigtab$name <- paste(sigtab$Family, sep=" ", rownames(sigtab))
dim(sigtab)
#View(sigtab)


write.csv(sigtab, "DESeq/csv_files/psalb_D0vsD21_DESeq2.csv")

#

#      ##### Horizontal barplot like Lefse #####
#By day
mtzD07tab
mtsD721tab
albD07tab
albD721tab
ctrlD07tab
ctrlD721tab

#colors
c0 <- "#6bc0a7"
c7 <- "#367F6A"
c21 <- "#1c4539"
a0 <- "#f98e63"
a7 <- "#C15E37"
a21 <- "#70351e"
m0 <- "#8ba0ca"
m7 <- "#506796"
m21 <- "#2c3a57"
NS <- "#878787"

#Which table?
tab <- psctrlD021tab
#min(psctrlD021tab$log2FoldChange)
#max(psctrlD021tab$log2FoldChange)
#-24.5, 5.3


#Grab significant 
alpha = 0.05
sigtab = tab[which(tab$padj < alpha), ]
sigtab$name <- paste(sigtab$Family, sep=" ", rownames(sigtab))
dim(sigtab)
#


#         ##### For D0 vs. D7 values #####
ggp <- ggplot(sigtab, aes(x=reorder(name, -log2FoldChange), log2FoldChange,fill=log2FoldChange >0)) +  
  geom_bar(stat = "identity", color="#000000") +
  scale_fill_manual(name = 'log2FoldChange > 0', values = setNames(c(m7, m0),c(T, F))) + #first color is the positive color
  theme_bw() +
  ggtitle("mtzD0vsD7 DESeq2 with ashr shrinkage") +
  ylim(-25, 25) +
  theme(legend.position = "none")

ggctrlD0D7
ggalbD0D7
ggmtzD0D7 <- ggp +  coord_flip()

grid.arrange(ggctrlD0D7, ggalbD0D7, ncol=2)



#         ##### Comparing treatments within the same time point #####

#By treatment. Some did not have any sig. diff, so don't have them listed here. 
ca_D0tab #Did not create a plot bc only 1 ASV and have volcano. 
mc_D7tab
ma_D7tab
ca_D7tab
ps_mc_D21tab

ggp <- ggplot(sigtab, aes(x=reorder(name, -log2FoldChange), log2FoldChange,fill=log2FoldChange >0)) +  
  geom_bar(stat = "identity", color="#000000") +
  scale_fill_manual(name = 'log2FoldChange > 0', values = setNames(c(m7, c7),c(T, F))) + #first color is the positive color
  theme_bw() +
  ggtitle("Control vs. MTZ D7 DESeq2 with ashr shrinkage") +
  #ylim(-25, 25) +
  theme(legend.position = "none")

ggp.ps_mc_D21tab <- ggp + coord_flip()


ggp.mc_D7tab 
ggp.ma_D7tab
ggp.ca_D7tab
ggp.ps_mc_D21tab

grid.arrange(ggp.ca_D7tab, ggp.ma_D7tab, ggp.mc_D7tab,ncol=3)
#

