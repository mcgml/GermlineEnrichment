#!/bin/bash -e
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
# Script 1 runs in sample folder

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#map one or more lanes of data for a single sample
for fastqPair in $(ls "$sampleId"_*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do
        
    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2*fastq.gz)
    
    #Trim adapters and remove short reads
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -o $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') \
    -p $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') \
    --minimum-length 30 \
    "$read1Fastq" \
    "$read2Fastq"

    #Align reads to reference genome, sort by coordinate and convert to BAM
    /share/apps/bwa-distros/bwa-0.7.15/bwa mem \
    -M \
    -R '@RG\tID:'"$seqId"_"$laneId"_"$sampleId"'\tSM:'"$sampleId"'\tPL:ILLUMINA\tLB:'"$worklistId"_"$panel"_"$sampleId"'\tPU:'"$seqId"_"$laneId" \
    -t 8 \
    /state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
    $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools sort -m8G -@4 -l0 -o "$seqId"_"$sampleId"_"$laneId"_sorted.bam

done

#merge mulitple lanes
/share/apps/samtools-distros/samtools-1.3.1/samtools merge \
-@ 12 \
-u "$seqId"_"$sampleId"_all_sorted.bam \
"$seqId"_"$sampleId"_*_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_all_sorted.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
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
rm -r tmp
rm *fastq *sorted.bam *rmdup.bam *rmdup.bai

#create BAM list for script 2
find $PWD -name "$seqId"_"$sampleId"_realigned.bam >> ../RealignedBams.list

#check if all BAMs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../RealignedBams.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineEnrichment.sh .. && cd .. && qsub 2_GermlineEnrichment.sh
fi