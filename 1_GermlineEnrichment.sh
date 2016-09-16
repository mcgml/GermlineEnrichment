#!/bin/bash
set -euo pipefail
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
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
# Script 1 runs in sample folder, requires uBAM

#load sample & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Preprocessing ###
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
INTERLEAVE=true \
NON_PF=true | \
/share/apps/cutadapt-distros/cutadapt-1.10/bin/cutadapt \
-a "$read1Adapter" \
-A "$read2Adapter" \
--interleaved \
--minimum-length 30 \
--quiet \ 
- | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem -M -t 7 -p /state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta /dev/stdin | \  
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MergeBamAlignment \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \ 
OUTPUT="$seqId"_"$sampleId"_clean.bam \
R=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
CREATE_INDEX=true \
ADD_MATE_CIGAR=true \
CLIP_OVERLAPPING_READS=true \
INCLUDE_SECONDARY_ALIGNMENTS=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
ATTRIBUTES_TO_RETAIN=XS

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_clean.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-ip "$padding" \
-nt 12 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realign.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam \
-compress 0 \
-dt NONE

### Clean up ###
rm "$seqId"_"$sampleId"_rmdup.bam "$seqId"_"$sampleId"_rmdup.bai

#create BAM list for script 2
find $PWD -name "$seqId"_"$sampleId"_realigned.bam >> ../RealignedBams.list

#check if all BAMs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../RealignedBams.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineEnrichment.sh .. && cd .. && qsub 2_GermlineEnrichment.sh
fi