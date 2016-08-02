#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
#PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

#Description: Germline Illumina TruSight
#Author: Matthew Lyon, All Wales Medical Genetics Lab
#Status: Development
#Mode: BY_SAMPLE
version="1.0.0"

#TODO file staging
#TODO file piping
#TODO MD5 input files
#TODO QC
#TODO multithread

#load sample variables
. variables

#trim adapters and remove short reads
/share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
-a CTGTCTCTTATACACATCT \
-A CTGTCTCTTATACACATCT \
-o "$seqId"_"$sampleId"_R1_trimmed.fastq \
-p "$seqId"_"$sampleId"_R2_trimmed.fastq \
"$read1Fastq" \
"$read2Fastq"

#Align reads to reference genome
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-R '@RG\tID:'"$seqId"_"$laneId"_"$sampleId"'\tSM:'"$sampleId"'\tPL:ILLUMINA\tLB:'"$worksheetId_$sampleId" \
/data/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
"$seqId"_"$sampleId"_R1_trimmed.fastq "$seqId"_"$sampleId"_R2_trimmed.fastq \
| /share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$seqId"_"$sampleId"_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_sorted.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_dupMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.intervals \
-L "$version"/"$bedFileName" \
-ip 100 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realigned.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam \
-compress 0 \
-dt NONE

#Analyze patterns of covariation in the sequence dataset
/share/apps/jre-distros/jre1.8.0_71/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_realigned.bam \
-L "$version"/"$bedFileName" \
-o "$seqId"_"$sampleId"_recal_data.table \
-ip 200 \
-dt NONE

#Do a second pass to analyze covariation remaining after recalibration
/share/apps/jre-distros/jre1.8.0_71/bin/java -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-I "$seqId"_"$sampleId"_realigned.bam \
-L "$version"/"$bedFileName" \
-o "$seqId"_"$sampleId"_post_recal_data.table \
-ip 200 \
-dt NONE

#Generate before/after plots
/share/apps/jre-distros/jre1.8.0_71/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$seqId"_"$sampleId"_recal_data.table \
-after "$seqId"_"$sampleId"_post_recal_data.table \
-plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
-L "$version"/"$bedFileName" \
-ip 200 \
-dt NONE

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_71/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR "$seqId"_"$sampleId"_post_recal_data.table \
-o "$seqId"_"$sampleId"_recal.bam \
-compress 0 \
-dt NONE

#fix bam
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -n -l0 "$seqId"_"$sampleId"_recal.bam | 
/share/apps/samtools-distros/samtools-1.3.1/samtools fixmates |
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -o "$seqId"_"$sampleId".bam

#variant calling with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_71/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-I "$seqId"_"$sampleId".bam \
-L "$version"/"$bedFileName" \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
--emitRefConfidence GVCF \
-dt NONE

#clean up
rm -r tmp
rm "$seqId"_"$sampleId"_R?_trimmed.fastq
rm "$seqId"_"$sampleId"_rmdup.ba?
rm "$seqId"_"$sampleId"_sorted.ba?
rm "$seqId"_"$sampleId"_recal.ba?
rm "$seqId"_"$sampleId"_realigned.ba?