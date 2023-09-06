#!/bin/bash
#SBATCH --account=penguin-np
#SBATCH --partition=penguin-shared-np
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=250000
#SBATCH --time=56:00:00
#SBATCH -J Mar2023WholeRun
#SBATCH -o /uufs/chpc.utah.edu/common/home/dearing-group1/Mar2023Run/Preprocess16S.outerror

###### Submitted:  Job ID: 7763209, Submitted at 10:00 on 3May2023. ######
###### Ran without incident, but initially only gave 36h for 192 samples. Needed more time, and
###### Resubmitted starting at Deblur step with excess time: Job ID: 7769096, Submitted at 10:00 on 3May2023. ######
###### Ran successfully, with deblur step taking almost 

#This script is for 16S sequences generated on the Illumina NovaSeq with some comparative seqs from 
# the previous MiSeq run to compare. Using Deblur but could use DADA2 if Deblur errors. 

# Load qiime2. Source it in your miniconda3 env
module use ~/MyModules
module load miniconda3/latest
source activate qiime2-2022.11

SCRATCH=/scratch/general/vast/u1212572
WRKDIR=/uufs/chpc.utah.edu/common/home/dearing-group1/Mar2023Run
MANIFEST=${WRKDIR}/Mar2023_SeqRun_Manifest.csv
CLASSIFIER=${HOME}/16S/metadata/silva138-1_515f-806r-average-classifier_fromZenodo.qza
METADATbyGNomexID=${WRKDIR}/MapByGNomexID_wcomparisons.tsv
METADATbySampleID=${WRKDIR}/MapBySampleID_wcomparisons.tsv

mkdir ${WRKDIR}/qiime_outputs



#Make sure you're working from SCRATCH and then exporting to WRKDIR.
cd ${SCRATCH}


##### 1. IMPORT DATA
#Beware of which file type you're using. May give you an error for either
# tab-delimited or with extra spaces.
qiime tools import \
 --type 'SampleData[PairedEndSequencesWithQuality]' \
 --input-path ${MANIFEST} \
 --output-path seqs_import.qza \
 --input-format PairedEndFastqManifestPhred33

# Visualize the reads
qiime demux summarize \
 --i-data seqs_import.qza \
 --o-visualization seqs_import.qzv

cp seqs_import.qz[av] ${WRKDIR}/qiime_outputs


##### 2. REMOVE PRIMERS
#Run cutadapt on demuxed reads
#Updated to have correct earth microbiome primer sequences for this run. For 1U4U and other
# runs with Zac, have used Takahashi that are V3-V4 instead of just V4. 
qiime cutadapt trim-paired \
   --i-demultiplexed-sequences seqs_import.qza \
   --p-cores 8 \
   --p-front-f GTGYCAGCMGCCGCGGTAA \
   --p-front-r GGACTACNVGGGTWTCTAAT \
   --p-discard-untrimmed TRUE \
   --p-error-rate 0 \
   --o-trimmed-sequences seqs_trim.qza \
   --verbose


# Visualize the reads
qiime demux summarize \
 --i-data seqs_trim.qza \
 --o-visualization seqs_trim.qzv

cp cutadapt\seqs_trim.qzv ${WRKDIR}/qiime_outputs

##### 3. JOIN PAIRED END SEQUENCES
#We will use deblur to denoise reads instead of DADA2. For deblur, you join the
# seqs first, and then denoise. In DADA2, it's the opposite.
#Joining pairs. we're allowing maxdiffs to be pretty high bc of long seqs.
qiime vsearch merge-pairs \
  --i-demultiplexed-seqs seqs_trim.qza \
  --o-merged-sequences seqs_trim_join.qza \
  --p-minmergelen 248 \
  --p-minovlen 5 \
  --p-maxdiffs 2 \
  --p-allowmergestagger \
  --p-threads 8 \
  --verbose

qiime demux summarize \
  --i-data seqs_trim_join.qza \
  --o-visualization seqs_trim_join.qzv


cp seqs_trim_join.qz[av] ${WRKDIR}/qiime_outputs



##### 4. QUALITY FILTER AFTER JOINING.
#Can filter quality score before or after joining, but I'm choosing to do it afterward.
#If before, just remove the "joined" part from the qiime quality-filter q-score-joined
# and make sure you have the right input object.
qiime quality-filter q-score \
 --i-demux seqs_trim_join.qza \
 --p-min-quality 4 \
 --o-filtered-sequences seqs_trim_join_filt.qza \
 --o-filter-stats seqs_trim_join_filt_stats.qza

#view your q-score output
qiime metadata tabulate \
  --m-input-file seqs_trim_join_filt_stats.qza \
  --o-visualization seqs_trim_join_filt_stats.qzv

#Look at sequences after denoising and filtering.
qiime demux summarize \
 --i-data seqs_trim_join_filt.qza \
 --o-visualization seqs_trim_join_filt.qzv

cp seqs_trim_join_filt_stats.qzv ${WRKDIR}/qiime_outputs
cp seqs_trim_join_filt.qzv ${WRKDIR}/qiime_outputs



##### 5. DENOISE WITH DEBLUR.
#If you're doing just V4, think about trimming to 250 bp, but for V3-V4, trim to
# 400-ish. Zac seems to trim to 392 for his things. Do this.
#If you're on CHPC, don't be shy about jobs-to-start up to the number of cores you're using. 
# Will run in parallel and be MUCH faster!
qiime deblur denoise-16S \
  --i-demultiplexed-seqs seqs_trim_join_filt.qza \
  --p-trim-length 250 \
  --p-jobs-to-start 16 \
  --o-table deblur_table.qza \
  --o-representative-sequences deblur_repseq.qza \
  --o-stats deblur_stats.qza

qiime deblur visualize-stats \
  --i-deblur-stats deblur_stats.qza \
  --o-visualization deblur_stats.qzv
  
  
qiime feature-table summarize \
  --i-table deblur_table.qza \
  --o-visualization deblur_table.qzv
  
  
cp deblur_stats.qzv ${WRKDIR}/qiime_outputs
cp deblur_table.qzv ${WRKDIR}/qiime_outputs


##### 6. FILTER OUT CHIMERAS
#Filter out chimeras from your first table. Using methods to keep borderline chimeras.
qiime vsearch uchime-denovo \
  --i-table deblur_table.qza \
  --i-sequences deblur_repseq.qza \
  --output-dir uchime_out

qiime feature-table filter-features \
  --i-table deblur_table.qza \
  --m-metadata-file uchime_out/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-table table_nochim.qza

qiime feature-table filter-seqs \
  --i-data deblur_repseq.qza \
  --m-metadata-file uchime_out/chimeras.qza \
  --p-exclude-ids \
  --o-filtered-data repseq_nochim.qza

#Visualize the table and sequences.
qiime feature-table summarize \
  --i-table table_nochim.qza \
  --o-visualization table_nochim.qzv
  
qiime feature-table tabulate-seqs \
  --i-data repseq_nochim.qza \
  --o-visualization repseq_nochim.qzv


cp repseq_nochim.qz[av] ${WRKDIR}/qiime_outputs
cp table_nochim.qz[av] ${WRKDIR}/qiime_outputs



##### 7. BUILD A PHYLOGENY
#Now, build phylogeny. This is important for doing any UniFRAC calculations
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences repseq_nochim.qza \
  --o-alignment aligned_repseq.qza \
  --o-masked-alignment masked_aligned_repseq.qza \
  --o-tree tree_unroot.qza \
  --o-rooted-tree tree_root.qza
  #Here we made both a rooted and unrooted tree for these taxa.

cp tree_unroot.qza ${WRKDIR}/qiime_outputs
cp tree_root.qza ${WRKDIR}/qiime_outputs


##### 8. CLASSIFY YOUR READS BY A CLASSIFIER YOU'VE TRAINED.
#Ok now we're actually going to classify by the CLASSIFIER
qiime feature-classifier classify-sklearn \
  --i-classifier ${CLASSIFIER} \
  --i-reads repseq_nochim.qza \
  --o-classification taxonomy.qza \
  --p-n-jobs 8

#Make a visualizer after classification
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

cp taxonomy.qz[av] ${WRKDIR}/qiime_outputs


##### 9. FURTHER FILTERING OF YOUR TABLE
#Filter out any features (aka ASVs) that are singletons or doubletons.
qiime feature-table filter-features \
    --i-table table_nochim.qza \
    --p-min-samples 3 \
    --p-min-frequency 10 \
    --o-filtered-table doubleton-lowfreq-filt-table.qza
#Added 2 commands here, so they much have at least 3 total reads and be seen at
# least a total of 10 times across all samples.

#Filter out mitochondria and chloroplast
qiime taxa filter-table \
    --i-table doubleton-lowfreq-filt-table.qza \
    --i-taxonomy taxonomy.qza \
    --p-exclude mitochondria,chloroplast \
    --o-filtered-table final-table-byGNomexID.qza



##### 10. RENAME SAMPLES BASED ON METADATA COLUMN.
  #Rename samples based on sample id instead of run ID for future labeling/figures.
  #Will then have to make sure that calling samples in the future is done by new
  # sample name.
qiime feature-table rename-ids \
  --i-table final-table-byGNomexID.qza \
  --m-metadata-file ${METADATbyGNomexID} \
  --m-metadata-column corrected_name \
  --o-renamed-table final-table.qza
    #Note: If adding metadata in qiime2, it doesn't like using "yes" or "no", so make sure
    # that you don't use only those words as a binary category.

qiime metadata tabulate \
  --m-input-file final-table.qza \
  --o-visualization final-table.qzv

cp final-table.qz[av] ${WRKDIR}/qiime_outputs


#Make a stacked bar to visualize the outputs. I'll do everything else in R.
qiime taxa barplot \
  --i-table final-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file ${METADATbySampleID} \
  --o-visualization table_taxbarplot.qzv

cp table_taxbarplot.qzv ${WRKDIR}/qiime_outputs




#When you know your script will run correctly, clean everything up from your
# scratch space. Don't do this if you're trying to track down where something
# is erroring.
#rm *.qz[av]




#Locally just going to separate qiime files into mouse and woodrat samples to then look at the 
# taxa barplots in qiime, with the interactive features. 
cd /Users/mdoolin/Documents/Mar2023Run
#OR Subset your feature table based on a category in your whole.
qiime feature-table filter-samples \
   --i-table qiime_outputs/final-table.qza \
   --m-metadata-file metadata/MapBySampleID_wo_comparisons.tsv \
   --p-where "[correct_type]!='mouse'" \
   --o-filtered-table qiime_outputs/coprophagy-table.qza

cd qiime_outputs
qiime taxa barplot \
  --i-table AnthelmandBlanks-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /Users/mdoolin/Documents/Anthelm_16S/metadata/AnthelmMetadatBySample_Novaseq.tsv \
  --o-visualization Anthelm_taxbarplot.qzv
  
  qiime taxa barplot \
  --i-table coprophagy-table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file /Users/mdoolin/Documents/Mar2023Run/metadata/MapBySampleID_WoodratsOnly.txt \
  --o-visualization coprophagy_taxbarplot.qzv


