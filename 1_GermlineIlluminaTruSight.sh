#!/bin/bash
#PBS -l walltime=04:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Illumina TruSight Pipeline. Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="dev"

#TODO file staging
#TODO file piping
#TODO MD5 input files
#TODO QC
#TODO multithread
#TODO handle multiple lanes
#TODO fastqc summary
#TODO insert size

#load sample variables
. variables

#split BED files by contig
grep -E ^"[1-22]\t" "$version"/"$bedFileName" > autosomal.bed
grep -E ^"X\t" "$version"/"$bedFileName" > x.bed
grep -E ^"Y\t" "$version"/"$bedFileName" > y.bed

#trim adapters and remove short reads
/share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
-a CTGTCTCTTATACACATCT \
-A CTGTCTCTTATACACATCT \
-o "$seqId"_"$sampleId"_R1_trimmed.fastq \
-p "$seqId"_"$sampleId"_R2_trimmed.fastq \
--minimum-length 30 \
"$read1Fastq" \
"$read2Fastq"

#fastqc
/share/apps/fastqc-distros/fastqc_v0.11.5/fastqc --extract "$seqId"_"$sampleId"_R1_trimmed.fastq
/share/apps/fastqc-distros/fastqc_v0.11.5/fastqc --extract "$seqId"_"$sampleId"_R2_trimmed.fastq

#Align reads to reference genome, sort by coordinate and convert to BAM
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-R '@RG\tID:'"$seqId"_"$laneNo"_"$sampleId"'\tSM:'"$sampleId"'\tPL:ILLUMINA\tLB:'"$worksheetId_$sampleId" \
/data/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
"$seqId"_"$sampleId"_R1_trimmed.fastq "$seqId"_"$sampleId"_R2_trimmed.fastq | \
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$seqId"_"$sampleId"_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_sorted.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-L "$version"/"$bedFileName" \
-ip 200 \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realign.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam \
-compress 0 \
-dt NONE

#Analyse patterns of covariation in the sequence dataset
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$seqId"_"$sampleId"_recal_data.table \
-after "$seqId"_"$sampleId"_post_recal_data.table \
-plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
-L "$version"/"$bedFileName" \
-ip 200 \
-dt NONE

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR "$seqId"_"$sampleId"_post_recal_data.table \
-o "$seqId"_"$sampleId"_recal.bam \
-compress 0 \
-dt NONE

#fix mate information
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -n -l0 "$seqId"_"$sampleId"_recal.bam |
/share/apps/samtools-distros/samtools-1.3.1/samtools fixmate - - |
/share/apps/samtools-distros/samtools-1.3.1/samtools sort -o "$seqId"_"$sampleId".bam
/share/apps/samtools-distros/samtools-1.3.1/samtools index "$seqId"_"$sampleId".bam
mv "$seqId"_"$sampleId".bam.bai "$seqId"_"$sampleId".bai

#variant calling with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-I "$seqId"_"$sampleId".bam \
-L "$version"/"$bedFileName" \
-o "$seqId"_"$sampleId".g.vcf \
-bamout "$seqId"_"$sampleId"_hc.bam \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
--emitRefConfidence GVCF \
-dt NONE

#structural variant calling with pindel
echo -e "$seqId"_"$sampleId".bam"\t"300"\t""$sampleId" > pindel.txt
/share/apps/pindel-distros/pindel-0.2.5b8/pindel \
-f /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-i pindel.txt \
-c ALL \
-T 12 \
--max_range_index 6 \
-o "$seqId"_"$sampleId"_pindel

#convert pindel output to VCF format
/share/apps/pindel-distros/pindel-0.2.5b8/pindel2vcf \
-P "$seqId"_"$sampleId"_pindel \
-r /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-R human_g1k_v37 \
-d none \
-e 3 \
--min_size 50 \
--min_coverage 10 \
-v "$seqId"_"$sampleId"_pindel.vcf

#extract 1kg on-target autosomal snps
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T SelectVariants \
--variant /data/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-o 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-env \
-ef \
-L autosomal.bed \
-dt NONE

#calculate dna contamination
/share/apps/verifyBamID-distros/verifyBamID-1.1.1/bin/verifyBamID \
--bam "$seqId"_"$sampleId".bam \
--vcf 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
--out "$seqId"_"$sampleId"_contamination \
--maxDepth 1000 \
--precise \
--ignoreRG \
--verbose

#calculate mean coverage for autosomes
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T DepthOfCoverage \
-o autosome \
-L autosome.bed \
-I "$seqId"_"$sampleId".bam \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
--omitLocusTable \
-dt NONE

#calculate mean coverage for X chrom
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T DepthOfCoverage \
-o X \
-L x.bed \
-I "$seqId"_"$sampleId".bam \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
--omitLocusTable \
-dt NONE

#calculate mean coverage for Y chrom
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T DepthOfCoverage \
-o y \
-L y.bed \
-I "$seqId"_"$sampleId".bam \
--omitDepthOutputAtEachBase \
--omitIntervalStatistics \
--omitLocusTable \
-dt NONE

#extract mean autosome, x and y coverage & calculate gender
autoCount=$(head -n2 autosome.sample_summary | tail -n1 | cut -s -f3)
xCount=$(head -n2 x.sample_summary | tail -n1 | cut -s -f3)
yCount=$(head -n2 y.sample_summary | tail -n1 | cut -s -f3)

#generate per-base coverage
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L "$version"/"$bedFileName" \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
-ct 30 \
-dt NONE

#calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator2.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
"$version"/"$geneListFileName" \
/data/db/human/refseq/refseq.gtf > "$seqId"_"$sampleId"_PercentageCoverage.txt

#TODO collect QC metrics for sample


#create VCFs list for filtering
find `pwd` -name "$seqId"_"$sampleId".g.vcf >> ../VCFsforFiltering.list

#check if all VCFs are written
if [ ${#AnalysisDirs[@]} -eq $(wc -l ../VCFsforFiltering.list | cut -f1 -d' ') ]
  then
      cd ..
      qsub 2_GermlineIlluminaTruSight.sh
fi

#clean up
rm -r tmp
rm "$seqId"_"$sampleId"_R?_trimmed.fastq
rm "$seqId"_"$sampleId"_rmdup.ba?
rm "$seqId"_"$sampleId"_sorted.ba?
rm "$seqId"_"$sampleId"_recal.ba?
rm "$seqId"_"$sampleId"_realigned.ba?
rm autosome.sample_statistics x.sample_statistics y.sample_statistics
rm autosomal.bed x.bed y.bed
rm pindel.txt
rm "$seqId"_"$sampleId"_R?_trimmed_fastqc.html
rm "$seqId"_"$sampleId"_R?_trimmed_fastqc.zip