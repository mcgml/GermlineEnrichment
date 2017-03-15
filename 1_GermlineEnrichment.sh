#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.6.0"

# Script 1 runs in sample folder, requires fastq files split by lane
#TODO exlude NTC

countQCFlagFails() {
    #count how many core FASTQC tests failed
    grep -E "Basic Statistics|Per base sequence quality|Per tile sequence quality|Per sequence quality scores|Per base N content" "$1" | \
    grep -v ^PASS | \
    grep -v ^WARN | \
    wc -l | \
    sed 's/^[[:space:]]*//g'
}

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#record FASTQC pass/fail
rawSequenceQuality=PASS

#convert FASTQ to uBAM & add RGIDs
for fastqPair in $(ls "$sampleId"_S*.fastq.gz | cut -d_ -f1-3 | sort | uniq); do

    #parse fastq filenames
    laneId=$(echo "$fastqPair" | cut -d_ -f3)
    read1Fastq=$(ls "$fastqPair"_R1_*fastq.gz)
    read2Fastq=$(ls "$fastqPair"_R2_*fastq.gz)

    #trim adapters
    /share/apps/cutadapt-distros/cutadapt-1.9.1/bin/cutadapt \
    -a "$read1Adapter" \
    -A "$read2Adapter" \
    -m 35 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #TODO FASTQ screen

    #convert fastq to ubam
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar FastqToSam \
    F1="$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    F2="$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    O="$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    QUALITY_FORMAT=Standard \
    READ_GROUP_NAME="$seqId"_"$laneId"_"$sampleId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" \
    PLATFORM_UNIT="$seqId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir

    #fastqc
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R1.fastq
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc -d /state/partition1/tmpdir --threads 12 --extract "$seqId"_"$sampleId"_"$laneId"_R2.fastq

    #check FASTQC output
    if [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R1_fastqc/summary.txt) -gt 0 ] || [ $(countQCFlagFails "$seqId"_"$sampleId"_"$laneId"_R2_fastqc/summary.txt) -gt 0 ]; then
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq *_fastqc.zip *_fastqc.html

done

#merge lane bams
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
INTERLEAVE=true \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
COMPRESSION_LEVEL=0 \
TMP_DIR=/state/partition1/tmpdir | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-t 12 \
-p \
/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
/dev/stdin | \
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MergeBamAlignment \
EXPECTED_ORIENTATIONS=FR \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
REFERENCE_SEQUENCE=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
PAIRED_RUN=true \
SORT_ORDER="coordinate" \
IS_BISULFITE_SEQUENCE=false \
ALIGNED_READS_ONLY=false \
CLIP_ADAPTERS=false \
MAX_RECORDS_IN_RAM=2000000 \
ADD_MATE_CIGAR=true \
MAX_INSERTIONS_OR_DELETIONS=-1 \
PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
UNMAP_CONTAMINANT_READS=false \
CLIP_OVERLAPPING_READS=true \
ALIGNER_PROPER_PAIR_FLAGS=false \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \
TMP_DIR=/state/partition1/tmpdir

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_aligned.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/state/partition1/tmpdir

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx24g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip "$padding" \
-dt NONE

#Realign around indels
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T IndelRealigner \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-targetIntervals "$seqId"_"$sampleId"_realign.intervals \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realigned.bam \
-dt NONE

if [ "$includeBQSR" = true ] ; then

    echo "Performing BQSR ..."

    #Analyse patterns of covariation in the sequence dataset
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx6g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -I "$seqId"_"$sampleId"_realigned.bam \
    -L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    -o "$seqId"_"$sampleId"_recal_data.table \
    -ip "$padding" \
    -nct 12 \
    -dt NONE

    #Do a second pass to analyze covariation remaining after recalibration
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx6g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T BaseRecalibrator \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
    -knownSites /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
    -BQSR "$seqId"_"$sampleId"_recal_data.table \
    -I "$seqId"_"$sampleId"_realigned.bam \
    -L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    -o "$seqId"_"$sampleId"_post_recal_data.table \
    -nct 12 \
    -ip "$padding" \
    -dt NONE

    #Generate BQSR plots
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T AnalyzeCovariates \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -before "$seqId"_"$sampleId"_recal_data.table \
    -after "$seqId"_"$sampleId"_post_recal_data.table \
    -plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
    -csv "$seqId"_"$sampleId"_recalibration.csv \
    -dt NONE

    #Apply the recalibration to your sequence data
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx12g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T PrintReads \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -I "$seqId"_"$sampleId"_realigned.bam \
    -BQSR "$seqId"_"$sampleId"_recal_data.table \
    -o "$seqId"_"$sampleId".bam \
    -dt NONE

else
    
    echo "Skipping BQSR ..."

    cp "$seqId"_"$sampleId"_realigned.bam "$seqId"_"$sampleId".bam
    cp "$seqId"_"$sampleId"_realigned.bai "$seqId"_"$sampleId".bai

fi

### Variant calling ###

#SNPs and Indels GVCF with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 1 \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
--emitRefConfidence GVCF \
-dt NONE

### QC ###

#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar BedToIntervalList \
I=/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
O="$panel"_ROI.interval_list \
SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict 

#TODO mapping rate


#Calculate insert size: fragmentation performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_InsertMetrics.txt \
H="$seqId"_"$sampleId"_InsertMetrics.pdf

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_HsMetrics.txt \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$panel"_ROI.interval_list \
TARGET_INTERVALS="$panel"_ROI.interval_list

#Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx12g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 10 \
--omitIntervalStatistics \
--omitLocusTable \
-ct "$minimumCoverage" \
-nt 12 \
-dt NONE

#tabix index the per-base coverage file
awk -F'[\t|:]' '{if(NR>1) print $1"\t"$2"\t"$3}' "$seqId"_"$sampleId"_DepthOfCoverage | \
/share/apps/htslib-distros/htslib-1.3.1/bgzip > "$seqId"_"$sampleId"_DepthOfCoverage.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -b2 "$seqId"_"$sampleId"_DepthOfCoverage.gz

#Calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.2/CoverageCalculator-2.0.2.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_genes.txt \
/state/partition1/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p5 \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#sort BED and add file prefix
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort -i "$sampleId"_gaps.bed > "$seqId"_"$sampleId"_Coverage_Gaps.bed

#Extract 1kg autosomal snps for contamination analysis
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T SelectVariants \
--variant /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-o 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-env \
-ef \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-XL X -XL Y -XL MT \
-dt NONE

#Calculate dna contamination: sample-to-sample contamination
/share/apps/verifyBamID-distros/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID \
--vcf 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
--bam "$seqId"_"$sampleId".bam \
--out "$seqId"_"$sampleId"_Contamination \
--verbose \
--ignoreRG \
--chip-none \
--minMapQ 20 \
--maxDepth 1000 \
--precise

#Gather QC metrics TODO at/gc dropout
meanInsertSize=$(head -n8 "$seqId"_"$sampleId"_InsertMetrics.txt | tail -n1 | cut -s -f5) #mean insert size
sdInsertSize=$(head -n8 "$seqId"_"$sampleId"_InsertMetrics.txt | tail -n1 | cut -s -f6) #insert size standard deviation
duplicationRate=$(head -n8 "$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt | tail -n1 | cut -s -f9) #The percentage of mapped sequence that is marked as duplicate.
totalReads=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_HsMetrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases. NB BQSR requires >= 100M, ideally >= 1B
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection
freemix=$(tail -n1 "$seqId"_"$sampleId"_Contamination.selfSM | cut -s -f7) #percentage DNA contamination. Should be <= 0.02

#gender analysis using Y chrom coverage
#todo gender analysis without y coverage
awk '{if ($1 == "Y") print $0}' /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed > Y.bed

if [ $(wc -l Y.bed |cut -d' ' -f1) -gt 0 ] && [ $(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" 'BEGIN{printf "%3.0f", meanOnTargetCoverage}') -gt 10 ]; then

    #calc Y coverage
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx12g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T DepthOfCoverage \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -o "$seqId"_"$sampleId"_Y \
    --omitDepthOutputAtEachBase \
    --omitIntervalStatistics \
    --omitLocusTable \
    -L Y.bed \
    -I "$seqId"_"$sampleId".bam \
    --countType COUNT_FRAGMENTS \
    --minMappingQuality 20 \
    -dt NONE

    #extract Y mean coverage
    meanYCov=$(head -n2 "$seqId"_"$sampleId"_Y.sample_summary | tail -n1 | cut -s -f3)
    calcGender=$(awk -v meanOnTargetCoverage="$meanOnTargetCoverage" -v meanYCov="$meanYCov" 'BEGIN {if (meanYCov > 10 && (meanYCov / meanOnTargetCoverage) > 0.1){print "MALE"} else if (meanYCov < 10 && (meanYCov / meanOnTargetCoverage) < 0.1) {print "FEMALE" } else {print "UNKNOWN"} }')

    #clean up
    rm "$seqId"_"$sampleId"_Y.*

else
    calcGender="UNKNOWN"
fi

#Print QC metrics
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tPctTargetBasesCt\tMeanOnTargetCoverage\tGender\tEstimatedContamination\tMeanInsertSize\tSDInsertSize" > "$seqId"_"$sampleId"_QC.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$calcGender\t$freemix\t$meanInsertSize\t$sdInsertSize" >> "$seqId"_"$sampleId"_QC.txt

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Germline,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=GermlineEnrichment,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",PctTargetBasesCt\="$pctTargetBasesCt",EstimatedContamination\="$freemix",GenotypicGender\="$calcGender",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteVcfFilePath\=$(dirname $PWD)/"$seqId"_filtered_meta_annotated.vcf,RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> > "$seqId"_"$sampleId"_meta.txt

#Create PED file
#TSV Format: Family_ID, Individual_ID, Paternal_ID, Maternal_ID, Sex (1=male; 2=female; 0=unknown), Phenotype (Description or 1=unaffected, 2=affected, 0=missing). Missing data is 0
if [ ! -z ${familyId-} ]; then echo -ne "$familyId\t" > "$sampleId"_pedigree.ped; else echo -ne "0\t" > "$sampleId"_pedigree.ped; fi
echo -ne "$sampleId\t" >> "$sampleId"_pedigree.ped
if [ ! -z ${paternalId-} ]; then echo -ne "$paternalId\t" >> "$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$sampleId"_pedigree.ped; fi
if [ ! -z ${maternalId-} ]; then echo -ne "$maternalId\t" >> "$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$sampleId"_pedigree.ped; fi
if [ ! -z ${gender-} ]; then echo -ne "$gender\t" >> "$sampleId"_pedigree.ped; else echo -ne "0\t" >> "$sampleId"_pedigree.ped; fi
if [ ! -z ${phenotype-} ]; then echo -e "$phenotype" >> "$sampleId"_pedigree.ped; else echo -e "2" >> "$sampleId"_pedigree.ped; fi

cat "$sampleId"_pedigree.ped >> ../"$seqId"_pedigree.ped

### Clean up ###

#create final file lists
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../GVCFs.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../BAMs.list

#filter low coverage samples 
if [ $(echo "$meanOnTargetCoverage" | awk '{if ($1 > 50) print "true"; else print "false"}') = true ]; then
    find $PWD -name "$seqId"_"$sampleId".bam >> ../HighCoverageBams.list
fi

#delete unused files
rm "$seqId"_"$sampleId"*unaligned.bam "$seqId"_"$sampleId"_rmdup.bam "$seqId"_"$sampleId"_rmdup.bai "$seqId"_"$sampleId"_realigned.bam 
rm "$seqId"_"$sampleId"_realigned.bai 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf Y.bed "$panel"_ROI.interval_list
rm 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf.idx "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_aligned.bai
rm "$seqId"_"$sampleId"_Contamination.log "$sampleId"_gaps.bed "$seqId"_"$sampleId"_DepthOfCoverage "$seqId"_"$sampleId"_DepthOfCoverage.sample_statistics

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineEnrichment.sh .. && cd .. && qsub 2_GermlineEnrichment.sh
fi