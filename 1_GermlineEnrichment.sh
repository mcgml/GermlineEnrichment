#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.0.6"

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
# Script 1 runs in sample folder, requires fastq files split by lane

#TODO bin gender score

phoneTrello() {
    #Call trello API
    /share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2" #seqId & message
}

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

phoneTrello "$seqId" "Starting GermlineEnrichment analysis for $sampleId ..."

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
    -m 30 \
    -o "$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    -p "$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    "$read1Fastq" \
    "$read2Fastq"

    #convert fastq to ubam
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar FastqToSam \
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
        phoneTrello "$seqId" "$sampleId has failed raw sequence quality check for $laneId"
        rawSequenceQuality=FAIL
    fi

    #clean up
    rm "$seqId"_"$sampleId"_"$laneId"_R1.fastq "$seqId"_"$sampleId"_"$laneId"_R2.fastq

done

#merge lane bams
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar MergeSamFiles \
$(ls "$seqId"_"$sampleId"_*_unaligned.bam | sed 's/^/I=/' | tr '\n' ' ') \
SORT_ORDER=queryname \
ASSUME_SORTED=true \
VALIDATION_STRINGENCY=SILENT \
USE_THREADING=true \
MAX_RECORDS_IN_RAM=2000000 \
TMP_DIR=/state/partition1/tmpdir \
O="$seqId"_"$sampleId"_unaligned.bam

#uBam2fq, map & MergeBamAlignment
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar SamToFastq \
I="$seqId"_"$sampleId"_unaligned.bam \
FASTQ=/dev/stdout \
INTERLEAVE=true \
NON_PF=true \
MAX_RECORDS_IN_RAM=2000000 \
VALIDATION_STRINGENCY=SILENT \
TMP_DIR=/state/partition1/tmpdir | \
/share/apps/bwa-distros/bwa-0.7.15/bwa mem \
-M \
-t 10 \
-p \
/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
/dev/stdin | \
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar MergeBamAlignment \
EXPECTED_ORIENTATIONS=FR \
ATTRIBUTES_TO_RETAIN=X0 \
ALIGNED_BAM=/dev/stdin \
UNMAPPED_BAM="$seqId"_"$sampleId"_unaligned.bam \
OUTPUT="$seqId"_"$sampleId"_aligned.bam \
R=/state/partition1/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
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
ATTRIBUTES_TO_RETAIN=XS \
INCLUDE_SECONDARY_ALIGNMENTS=true \
CREATE_INDEX=true \
TMP_DIR=/state/partition1/tmpdir

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_aligned.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
MAX_RECORDS_IN_RAM=2000000 \
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
-nt 12 \
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
    -nct 12 \
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
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
-bamout "$seqId"_"$sampleId"_HC.bam \
--emitRefConfidence GVCF \
-dt NONE

### QC ###

#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar BedToIntervalList \
I=/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
O="$panel"_ROI.interval_list \
SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict 

#Calculate insert size: fragmentation performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_insert_metrics.txt \
H="$seqId"_"$sampleId"_insert_metrics.pdf

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.0/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
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
-ct "$minimumCoverage" \
-nt 12 \
-dt NONE

#Calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.2/CoverageCalculator-2.0.2.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_genes.txt \
/state/partition1/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p5 \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gender analysis using off-targed reads
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools slop \
-i /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-g /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai \
-b 300 > padded.bed

grep -P "^Y\t" /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | awk '{print $1"\t"0"\t"$2}' > Y.bed
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
-a Y.bed \
-b padded.bed \
-b /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/PAR.bed > Y.off.bed

grep -P "^X\t" /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | awk '{print $1"\t"0"\t"$2}' > X.bed
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
-a X.bed \
-b padded.bed \
-b /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/PAR.bed > X.off.bed

chromYCount=$(/share/apps/samtools-distros/samtools-1.3.1/samtools view \
-c \
-f0x0002 \
-F0x0100 \
-F0x0400 \
-q60 \
-L Y.off.bed \
"$seqId"_"$sampleId".bam)

chromXCount=$(/share/apps/samtools-distros/samtools-1.3.1/samtools view \
-c \
-f0x0002 \
-F0x0100 \
-F0x0400 \
-q60 \
-L X.off.bed \
"$seqId"_"$sampleId".bam)

gender=$(echo "print ($chromYCount / $(awk '{n+= $3-$2} END {print n}' Y.off.bed)) / ($chromXCount / $(awk '{n+= $3-$2} END {print n}' X.off.bed))" | perl)

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
-nt 12 \
-dt NONE

#Calculate dna contamination: sample-to-sample contamination
/share/apps/verifyBamID-distros/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID \
--vcf 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
--bam "$seqId"_"$sampleId".bam \
--out "$seqId"_"$sampleId"_contamination \
--verbose \
--ignoreRG \
--chip-none \
--minMapQ 20 \
--maxDepth 1000 \
--precise

#Gather QC metrics
meanInsertSize=$(head -n8 "$seqId"_"$sampleId"_insert_metrics.txt | tail -n1 | cut -s -f5) #mean insert size
sdInsertSize=$(head -n8 "$seqId"_"$sampleId"_insert_metrics.txt | tail -n1 | cut -s -f6) #insert size standard deviation
duplicationRate=$(head -n8 "$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt | tail -n1 | cut -s -f9) #The percentage of mapped sequence that is marked as duplicate.
totalReads=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f6) #The total number of reads in the SAM or BAM file examine.
pctSelectedBases=$(head -n8 "$seqId"_"$sampleId"_hs_metrics.txt | tail -n1 | cut -s -f19) #On+Near Bait Bases / PF Bases Aligned.
totalTargetedUsableBases=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f2) #total number of usable bases. NB BQSR requires >= 100M, ideally >= 1B
meanOnTargetCoverage=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f3) #avg usable coverage
pctTargetBasesCt=$(head -n2 $seqId"_"$sampleId"_DepthOfCoverage".sample_summary | tail -n1 | cut -s -f7) #percentage panel covered with good enough data for variant detection
freemix=$(tail -n1 "$seqId"_"$sampleId"_contamination.selfSM | cut -s -f7) #percentage DNA contamination. Should be <= 0.02

#Print QC metricss
echo -e "TotalReads\tRawSequenceQuality\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tpctTargetBasesCt\tMeanOnTargetCoverage\tGender\tFreemix\tMeanInsertSize\tSDInsertSize" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$rawSequenceQuality\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$gender\t$freemix\t$meanInsertSize\t$sdInsertSize" >> "$seqId"_"$sampleId"_qc.txt

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",WorklistId\="$worklistId",SeqId\="$seqId",Panel\="$panel",PipelineName\=GermlineEnrichment,PipelineVersion\="$version",RawSequenceQuality\="$rawSequenceQuality",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",pctTargetBasesCt\="$pctTargetBasesCt",Freemix\="$freemix",Gender\="$gender",TotalTargetedUsableBases\="$totalTargetedUsableBases",RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> > "$seqId"_"$sampleId"_meta.txt

### Clean up ###

#create final file lists
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../GVCFs.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../FinalBams.list

#delete unused files
rm "$seqId"_"$sampleId"*unaligned.bam "$seqId"_"$sampleId"_rmdup.bam "$seqId"_"$sampleId"_rmdup.bai "$seqId"_"$sampleId"_realigned.bam 
rm "$seqId"_"$sampleId"_realigned.bai X.off.bed X.bed Y.off.bed Y.bed 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf
rm 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf.idx padded.bed "$seqId"_"$sampleId"_aligned.bam "$seqId"_"$sampleId"_aligned.bai

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineEnrichment.sh .. && cd .. && qsub 2_GermlineEnrichment.sh
fi