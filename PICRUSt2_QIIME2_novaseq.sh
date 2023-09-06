#Qiime2 picrust script
#Running based on tutorial on q2-picrust2 github.
https://github.com/picrust/picrust2/wiki/q2-picrust2-Tutorial
#If, for some reason, the plugin isn't working, you can use picrust outside of
# qiime as well.
https://github.com/picrust/picrust2/wiki/PICRUSt2-Tutorial-(v2.4.2)
# This is more manipulatable but just not as stream-lined.
# Right now, my goal is to just get some results with the typical pipeline, not
# trying to mess with any of the settings.

#Activate qiime environment.
conda activate qiime2-2022.11

conda install q2-picrust2=2023.2 -c conda-forge -c bioconda -c gavinmdouglas
#Grabbing picrust2 plug-in directly from github while I was in my qiime environment.
#They only have the most up to date picrust2 plug-in available here. Even though
# I don't have the most recent qiime version, it still seems like it would work...

qiime picrust2 full-pipeline --help
#See all of the possible inputs and calls.

#My qiime outputs are all in one folder in my Documents. Update directory.
cd ~/Documents/Anthelm_16S/GTDB_run2/

#You can look at your biom file in terminal if you want
biom head -i psrare.biom

#Probably have created final table as a phyloseq object. To convert to qiime
# file, first create a biom file in R and have your qiime2 metadata and repseq ready.
# Convert the biom file into a .qza qiime frequency table file in qiime:
qiime tools import \
  --input-path psrare.biom \
  --type 'FeatureTable[Frequency]' \
  --input-format BIOMV100Format \
  --output-path psrare-table.qza 

#The newer versions of qiime2, including 2021.4 (what I'm using) wants a
# BIOMV210Format as default, so must specify it's a BIOMV100Format

qiime picrust2 full-pipeline \
  --i-table psrare-table.qza \
  --i-seq ~/Documents/Mar2023Run/qiime_outputs/repseq_nochim.qza \
  --output-dir picrust2_psrare \
  --p-placement-tool epa-ng \
  --p-threads 4 \
  --p-hsp-method pic \
  --p-max-nsti 2 \
  --verbose

#Placement tool and hidden state prediction (hsp) methods are up for debate.
# Might want to try different combinations to see how much results change.

#With verbose call, it'll tell you any issues at each step. This includes:
# Seqs that don't align well to reference
# Whether there were any ASVs that don't make the cut for max-NSTI.

cd /Users/mdoolin/Documents/Anthelm_16S/Novaseq_run/picrust2_psrare
qiime feature-table summarize \
   --i-table pathway_abundance.qza \
   --o-visualization pathway_abundance.qzv

qiime feature-table summarize \
   --i-table ec_metagenome.qza \
   --o-visualization ec_metagenome.qzv

qiime feature-table summarize \
   --i-table ko_metagenome.qza \
   --o-visualization ko_metagenome.qzv

#Convert any of these into a visualizations.
qiime tools view pathway_abundance.qzv
qiime tools view ko_metagenome.qzv
qiime tools view ec_metagenome.qzv

#These visualizations aren't super informative. They're just summarizing the
# features in each sample. Do something else.


#Create visualizations with the metadata included.
#Make sure metadata file is a .tsv, in the right format for sample-id, and in the
# working directory.


#Get the sampling depth from the minimum frequency summary page on the visualization.
qiime diversity core-metrics \
   --i-table pathway_abundance.qza \
   --p-sampling-depth 968867 \
   --m-metadata-file ../psrare_metadat.tsv \
   --output-dir pathway_coremetrics \
   --p-n-jobs 4

qiime diversity core-metrics \
   --i-table ko_metagenome.qza \
   --p-sampling-depth 8160147 \
   --m-metadata-file ../psrare_metadat.tsv \
   --output-dir ko_coremetrics \
   --p-n-jobs 4

qiime diversity core-metrics \
    --i-table ec_metagenome.qza \
    --p-sampling-depth 4599893 \
    --m-metadata-file ../psrare_metadat.tsv \
    --output-dir ec_coremetrics \
    --p-n-jobs 4

qiime tools view pathway_coremetrics/bray_curtis_emperor.qzv
qiime tools view ko_coremetrics/bray_curtis_emperor.qzv
qiime tools view ec_coremetrics/bray_curtis_emperor.qzv


#Exporting things to use in R. First KO's. 
#Export the .qza as a biom file.
qiime tools export \
  --input-path ko_metagenome.qza \
  --output-path ko-biom-output

#Change the format of the feature-table.biom file
biom convert \
   -i ko-biom-output/feature-table.biom \
   -o kometagenome_biom.tsv \
   --to-tsv

qiime tools export \
  --input-path KO_coremetrics/bray_curtis_distance_matrix.qza \
  --output-path KO_coremetrics/bray_curtis_distance_matrix

qiime tools export \
  --input-path KO_coremetrics/jaccard_distance_matrix.qza \
  --output-path KO_coremetrics/jaccard_distance_matrix

#Now to do the pathway abundance too. 
#Export the .qza as a biom file.
qiime tools export \
  --input-path pathway_abundance.qza \
  --output-path pathabund-biom-output

#Change the format of the feature-table.biom file
biom convert \
   -i pathabund-biom-output/feature-table.biom \
   -o pathabund_biom.tsv \
   --to-tsv

qiime tools export \
  --input-path pathway_coremetrics/bray_curtis_distance_matrix.qza \
  --output-path pathway_coremetrics/bray_curtis_distance_matrix

qiime tools export \
  --input-path pathway_coremetrics/jaccard_distance_matrix.qza \
  --output-path pathway_coremetrics/jaccard_distance_matrix


#Trying to figure out how to find the file with all of the metadata that tells
# what each ortholog is. Maybe this is tied within the biom file? But it doesn't
# export with the Biom file.

cd ~/Documents/Anthelm_16S/GTDB_run2/picrust2_psraremtz/ko_coremetrics

qiime metadata tabulate \
   --m-input-file observed_features_vector.qza \
   --o-visualization ko_observed_features_vector.qzv

#Get back to where you want to be, more general folder.
cd ~/Documents/Anthelm_16S/GTDB_run2



#But what if we want to use subsetted data? Option 1: Input subsetted phyloseq object
# from R to QIIME.


############## Grab a phyloseq object to bring into Qiime2 (also in an R script)############


#In R, run these things (same script in in 1U4URerun_16S_FromQiime2.R script)
library(biomformat) #version 1.22.0


#Create the taxa table for use as FeatureTable[Sequence] (I think)
tax <- as(tax_table(psexp4000), "matrix")
tax_cols <- colnames(tax)
tax <- as.data.frame(tax)
tax$taxonomy <- do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co] <- NULL
write.table(tax, "psexp4000_tax.txt", quote=FALSE, col.names=FALSE, sep="\t")

#Create your OTU table for the biom file, as FeatureTable[Frequency]
otu <- t(as(otu_table(psexp4000),"matrix"))
otu_biom <- make_biom(data=otu)
write_biom(otu_biom, "psexp4000_biom.biom")
#Could also export this as a text file, but biom file from phyloseq object
# seems easier.

#Now write out your metadata for this phyloseq object.
write.table(sample_data(psexp4000), "psexp4000_metadat.txt", sep="\t",
            row.names=FALSE, col.names=TRUE, quote=FALSE)

# Ok, now go to qiime script saved in
#  ~/Documents/1U4U_16S/Rerun_scripts/PICRUSt2_QIIME2.sh
# for the rest of the magic, importing into qiime.



################## Back into qiime script now.  ######
#Some nice code for manipulation of which samples are in your qiime objects.


#If starting from R outputs and wanting to run this, you can convert the pieces
# from phyloseq into a qiime objext.

cd ~/Documents/Anthelm_16S/qiime2_GTDB/

#Import feature table as biom file.
qiime tools import \
--input-path psrare_biom.biom \
--type 'FeatureTable[Frequency]' \
--input-format BIOMV100Format \
--output-path psrare-feature-table.qza

#Import taxonomy table.
qiime tools import \
--type 'FeatureData[Taxonomy]' \
--input-format HeaderlessTSVTaxonomyFormat \
--input-path psrare_tax.txt \
--output-path psrare-taxonomy.qza

###Nevermind the 20 lines above this. Can come back if interested.
#This didn't work for picrust full pipeline because you need to have a representative
# sequence table. So back to working from the qiime artifact final-table.qza.

#We can rarefy and subset to look at samples of interest. Easier to do this all in
# qiime if the intial object was made in qiime. Otherwise, need to figure out how
# to pull representative sequences out of a phyloseq object, which would require
# breaking down the phylseq object to pull out samples, which is tough based on
# existing forum answers.

#Make sure inputs are all in there.
cd ~/Documents/1U4U_16S/qiime2_rerun_results/

#Making a rarefied table.
qiime feature-table rarefy \
  --i-table pstmp-exp-table.qza \
  --p-sampling-depth 9000 \
  --o-rarefied-table psrare-table.qza

#Pull out the metadata sheet from the rarefied table.

qiime metadata tabulate \
   --m-input-file rare4000-alltimepoints-table.qza \
   --o-visualization rare4000-alltimepoints-table.qzv



#Open the visualization and download the tsv from the qiime2 viewer.
#Copy the first column and make sure it's got the same name as your master metadata
# for your sample id's. For this, mine needs to be "sample-id".
#In R, merge based on the sample names in that sample-id column you downloaded from
# the rarefied qiime artifact.



#######Jumping into R:
#setwd("~/Documents/1U4U_16S/qiime2_rerun_results")
#make sure your master metadata is in the right folder.


#library(readr)
#setwd("~/Documents/1U4U_16S/qiime2_rerun_results")
#master <- read_tsv("MapBySampleID_1U4U_final_wRerun.tsv")
#names <- read_tsv("rare4000-alltimepoints-names.tsv")

#merged <- merge(master, names, by="sample-id")
#nrow(merged)

#Cool, looks like that worked since the final merged df has the right number of rows
# for the rarefied data instead of the unrarefied data.
#Now subset the merged dataset to make it only have experimental samples, and
# to only have the important samples. This way, I will easily be able to subset qiime
# artifacts for picrust to mimic the dfexp4000 and dfimp4000 samples.

#q2_exp4000_metadat <- subset(merged, !day %in% c("none","pre"))

#q2_imp4000_metadat <- subset(merged, day %in% c("D-5", "D14", "D41"))
#nrow(q2_exp4000_metadat)
#unique(q2_exp4000_metadat$day)

#Ok great, now write out those metadata files.

#write_tsv(q2_exp4000_metadat, "q2_exp4000_metadat.tsv")




############# Back To Qiime ###############
#Now subset the rarefied table to look only at timepoints I'm interested in.
# Do this in qiime, with the subsetted and exported metadata from the R code.

cd ~/Documents/1U4U_16S/qiime2_rerun_results

#Make an artifact that is the same as psexp4000. (Excludes pre-samples of LO mice.)
qiime feature-table filter-samples \
  --i-table final-table.qza \
  --m-metadata-file psrare_metadat.tsv \
  --o-filtered-table pstmp-exp-table.qza

#Make an artifact that is the same as psimp4000. (Includes only D-5, D14, D41 samples.)
#And also will make the hi and lo versions. Have made updated metadata files for
# these.
qiime feature-table filter-samples \
  --i-table pstmp-table.qza \
  --m-metadata-file psrare_metadat.tsv \
  --o-filtered-table imp4000hi-table.qza

#If you want to export as a biom file, do this:
qiime tools export \
  --input-path imp4000hi-table.qza \
  --output-path Biom_files/imp4000hi-table.biom
#Warning that it stays as named "feature-table.biom" so will haev to update names
# or update output folders if saving several.
