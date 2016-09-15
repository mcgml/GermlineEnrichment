#!/bin/bash -euxfo pipefail
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_LANE
version="dev"

# Directory structure required for pipeline
#
# /data
# └── results
#     └── seqId
#         ├── panel1
#         │   ├── sample1
#         │   ├── sample2
#         │   └── sample3
#         └── panel2
#             ├── sample1
#             ├── sample2
#             └── sample3
#
# Script 2 runs in panel folder

#load sample & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#Analyse patterns of covariation in the sequence dataset
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I RealignedBams.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_recal_data.table \
-ip "$padding" \
-nct 12 \
-dt NONE

#Do a second pass to analyze covariation remaining after recalibration
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR "$seqId"_recal_data.table \
-I RealignedBams.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_post_recal_data.table \
-nct 12 \
-ip "$padding" \
-dt NONE

#Generate BQSR plots
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$seqId"_recal_data.table \
-after "$seqId"_post_recal_data.table \
-plots "$seqId"_recalibration_plots.pdf \
-csv "$seqId"_recalibration.csv \
-dt NONE

### Clean up ###

#run script 3
find -mindepth 1 -maxdepth 1 -type d -exec bash -c "cd {} && qsub 3_GermlineEnrichment.sh" \;