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
--minimum-length 36 \
-o "$runId"_"$sampleId"_R1_trimmed.fastq \
-p "$runId"_"$sampleId"_R2_trimmed.fastq \
"$read1Fastq" \
"$read2Fastq"

#Align reads to reference genome
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-R '@RG\tID:'"$runId"_"$lane"_"$sampleId"'\tSM:'"$sampleId"'\tPL:ILLUMINA\tLB:'"$experimentName" \
/data/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
"$runId"_"$sampleId"_R1_trimmed.fastq "$runId"_"$sampleId"_R2_trimmed.fastq \
| /share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$runId"_"$sampleId"_sorted.bam

#Mark duplicate reads
/usr/java/jdk1.8.0/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$runId"_"$sampleId"_sorted.bam \
OUTPUT="$runId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$runId"_"$sampleId"_dupMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Analyze patterns of covariation in the sequence dataset
/usr/java/jdk1.8.0/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$runId"_"$sampleId"_rmdup.bam \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/"$version"/"$bedFileName" \
-o "$runId"_"$sampleId"_recal_data.table \
-ip 200 \
-dt NONE

#Do a second pass to analyze covariation remaining after recalibration
/usr/java/jdk1.8.0/bin/java -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-BQSR "$runId"_"$sampleId"_recal_data.table \
-I "$runId"_"$sampleId"_rmdup.bam \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/"$version"/"$bedFileName" \
-o "$runId"_"$sampleId"_post_recal_data.table \
-ip 200 \
-dt NONE

#Generate before/after plots
/usr/java/jdk1.8.0/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$runId"_"$sampleId"_recal_data.table \
-after "$runId"_"$sampleId"_post_recal_data.table \
-plots "$runId"_"$sampleId"_recalibration_plots.pdf \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/"$version"/"$bedFileName" \
-ip 200 \
-dt NONE

#Apply the recalibration to your sequence data
/usr/java/jdk1.8.0/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$runId"_"$sampleId"_rmdup.bam \
-BQSR "$runId"_"$sampleId"_post_recal_data.table \
-o "$runId"_"$sampleId".bam \
-compress 0 \
-dt NONE

#variant calling with Haplotypecaller
/usr/java/jdk1.8.0/bin/java -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-I "$runId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/"$version"/"$bedFileName" \
-o "$runId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
--emitRefConfidence GVCF \
-dt NONE