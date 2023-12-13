######### Qiime2 workflow from demultiplexed seqs ########## 

# Since we have so many files, we will use a manifest file to imoprt them
# We only use forward reads, as read quality is too poor to overlap forward and reverse reads

qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path sample_manifest2.tsv \
  --output-path single-end-demux.qza \
  --input-format SingleEndFastqManifestPhred33V2

# [this takes a while]

qiime demux summarize \
--i-data single-end-demux.qza \
--o-visualization single-end-demux.qzv

# trim primers usuing cutadapt 27F 534R

# trim primer sequence
qiime cutadapt trim-single \
--i-demultiplexed-sequences single-end-demux.qza \
--p-adapter CTGTCTCTTAT \
--p-front AGAGTTTGATYMTGGCTCAG \
--p-discard-untrimmed \
--o-trimmed-sequences single-end-trimmed.qza \
--verbose > cutadapt-log.txt


qiime demux summarize \
--i-data single-end-trimmed.qza \
--o-visualization single-end-trimmed.qzv

qiime dada2 denoise-single \
--i-demultiplexed-seqs single-end-trimmed.qza \
--p-trunc-len 128 \
--p-chimera-method consensus \
--output-dir DADA2_denoising_output \
--verbose


# visualise dada2 output

# summarise features remaining after dada2 denoising
qiime feature-table summarize \
--i-table DADA2_denoising_output/table.qza \
--o-visualization DADA2_denoising_output/table.qzv

# rep sequences 
qiime feature-table tabulate-seqs \
--i-data DADA2_denoising_output/representative_sequences.qza \
--o-visualization DADA2_denoising_output/representative_sequences.qzv

qiime metadata tabulate \
--m-input-file DADA2_denoising_output/denoising_stats.qza \
--o-visualization DADA2_denoising_output/denoising_stats.qzv


# Run using Alex Poret's classifier

qiime feature-classifier classify-sklearn \
--i-reads DADA2_denoising_output/representative_sequences.qza \
--i-classifier alex_class/classifier.qza \
--output-dir alex_classifier_output \
--verbose

# Visualise
qiime metadata tabulate \
--m-input-file alex_classifier_output/classification.qza \
--o-visualization alex_classifier_output/alex_classification.qzv


#output data aggregated to genus level

qiime taxa collapse \
  --i-table DADA2_denoising_output/table.qza \
  --i-taxonomy alex_classifier_output/classification.qza \
  --p-level 6 --o-collapsed-table alex_classifier_output/collapsed_taxtable.qza


# make an accesisble object we can download 
qiime metadata tabulate \
--m-input-file alex_classifier_output/collapsed_taxtable.qza \
--o-visualization alex_classifier_output/collapsed_taxtable.qzv
