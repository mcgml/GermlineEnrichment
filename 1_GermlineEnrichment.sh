#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="dev"

#TODO optimise pindel filter + ?LCR filter
#TODO ?BQSR on all samples for a single lane

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
    /data/db/human/mappers/b37/bwa/human_g1k_v37.fasta \
    $(echo "$read1Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') $(echo "$read2Fastq" | sed 's/\.fastq\.gz/_trimmed\.fastq/g') | \
    /share/apps/samtools-distros/samtools-1.3.1/samtools sort -l0 -o "$seqId"_"$sampleId"_"$laneId"_sorted.bam

done

#merge mulitple lanes
/share/apps/samtools-distros/samtools-1.3.1/samtools merge \
-@ 8 \
-u "$seqId"_"$sampleId"_all_sorted.bam \
"$seqId"_"$sampleId"_*_sorted.bam

#Mark duplicate reads
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar MarkDuplicates \
INPUT="$seqId"_"$sampleId"_all_sorted.bam \
OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
METRICS_FILE="$seqId"_"$sampleId"_MarkDuplicatesMetrics.txt \
CREATE_INDEX=true \
COMPRESSION_LEVEL=0

#Identify regions requiring realignment
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T RealignerTargetCreator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-known /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-known /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_rmdup.bam \
-o "$seqId"_"$sampleId"_realign.intervals \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-ip "$padding" \
-nt 8 \
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
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T BaseRecalibrator \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-knownSites /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/1000G_phase1.indels.b37.vcf \
-knownSites /data/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-I "$seqId"_"$sampleId"_realigned.bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_"$sampleId"_recal_data.table \
-ip "$padding" \
-nct 8 \
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
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_"$sampleId"_post_recal_data.table \
-nct 8 \
-ip "$padding" \
-dt NONE

#Generate BQSR plots
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T AnalyzeCovariates \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-before "$seqId"_"$sampleId"_recal_data.table \
-after "$seqId"_"$sampleId"_post_recal_data.table \
-plots "$seqId"_"$sampleId"_recalibration_plots.pdf \
-csv "$seqId"_"$sampleId"_recalibration.csv \
-dt NONE

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR "$seqId"_"$sampleId"_recal_data.table \
-o "$seqId"_"$sampleId".bam \
-compress 0 \
-nct 8 \
-dt NONE

### Variant calling ###

#SNPs and Indels with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
--emitRefConfidence GVCF \
--bamOutput "$seqId"_"$sampleId"_HC.bam \
-dt NONE

#Structural variants with pindel using padded BED file
echo -e "$seqId"_"$sampleId".bam"\t$expectedInsertSize\t""$sampleId" > pindel.txt

/share/apps/bedtools-distros/bedtools-2.24.0/bin/bedtools slop \
-i /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-g /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai \
-b 400 | /share/apps/bedtools-distros/bedtools-2.24.0/bin/bedtools merge > padded.bed

/share/apps/pindel-distros/pindel-0.2.5b8/pindel \
-f /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-j padded.bed \
-i pindel.txt \
-T 8 \
--max_range_index 6 \
--anchor_quality 20 \
--report_interchromosomal_events \
-L "$seqId"_"$sampleId"_pindel.log \
-o "$seqId"_"$sampleId"_pindel

#Convert pindel output to VCF format and filter calls
/share/apps/pindel-distros/pindel-0.2.5b8/pindel2vcf \
-P "$seqId"_"$sampleId"_pindel \
-r /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-R human_g1k_v37 \
-d none \
-e 3 \
--min_size 50 \
--min_coverage 10 \
-v "$seqId"_"$sampleId"_pindel.vcf

### QC ###

#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar BedToIntervalList \
I=/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
O="$bedFileName".interval_list \
SD=/data/db/human/gatk/2.8/b37/human_g1k_v37.dict

#Fastqc: raw sequence quality
echo -e "File\tBasicStatistics\tPerBaseSequenceQuality\tPerTileSequenceQuality\tPerSequenceQualityScores\tPerBaseNContent\tOverrepresentedSequences\tAdapterContent" > "$seqId"_"$sampleId"_fastqc.txt
for i in $(ls *_trimmed.fastq); do
    /share/apps/fastqc-distros/fastqc_v0.11.5/fastqc --threads 8 --extract "$i"
    fastqcFolder=$(echo "$i" | sed 's/\.fastq/_fastqc/g')
    
    basicStats=$(head -n1 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    perBaseSeqQuality=$(head -n2 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    perTileSeqQuality=$(head -n3 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    perSeqQualityScore=$(head -n4 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    perBaseNContent=$(head -n7 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    overRepresentedSeq=$(head -n10 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    adapterContent=$(head -n11 "$fastqcFolder"/summary.txt | tail -n1 |cut -s -f1)
    
    echo -e "$i\t$basicStats\t$perBaseSeqQuality\t$perTileSeqQuality\t$perSeqQualityScore\t$perBaseNContent\t$overRepresentedSeq\t$adapterContent" >> "$seqId"_"$sampleId"_fastqc.txt 
done

#Calculate insert size: fragmentation performance
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_insert_metrics.txt \
H="$seqId"_"$sampleId"_insert_metrics.pdf

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$bedFileName".interval_list \
TARGET_INTERVALS="$bedFileName".interval_list

#Generate per-base/per-target coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 10 \
--omitIntervalStatistics \
-ct "$minimumCoverage" \
-nt 8 \
-dt NONE

#Calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.0/CoverageCalculator-2.0.0.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_genes.txt \
/data/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p"$spliceSitePadding" \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gender analysis using off-targed reads
grep -P "^Y\t" /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | awk '{print $1"\t"0"\t"$2}' > Y.bed
/share/apps/bedtools-distros/bedtools2/bin/bedtools subtract \
-a Y.bed \
-b padded.bed \
-b /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/PAR.bed > Y.off.bed

grep -P "^X\t" /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | awk '{print $1"\t"0"\t"$2}' > X.bed
/share/apps/bedtools-distros/bedtools2/bin/bedtools subtract \
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
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T SelectVariants \
--variant /data/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-o 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-env \
-ef \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-XL X -XL Y \
-nt 8 \
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

#Print QC metrics
echo -e "TotalReads\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tpctTargetBasesCt\tMeanOnTargetCoverage\tGender\tFreemix\tMeanInsertSize\tSDInsertSize" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$gender\t$freemix\t$meanInsertSize\t$sdInsertSize" >> "$seqId"_"$sampleId"_qc.txt

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",WorklistId\="$worklistId",SeqId\="$seqId",Panel\="$panel",PipelineName\=GermlineEnrichment,PipelineVersion\="$version",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",pctTargetBasesCt\="$pctTargetBasesCt",Freemix\="$freemix",Gender\="$gender",RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> > "$seqId"_"$sampleId"_meta.txt

### Clean up ###
rm -r tmp

#create BAM and VCF list for script 2
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../VCFsforFiltering.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../FinalBAMs.list

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../VCFsforFiltering.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    echo -e "seqId=$seqId\npanel=$panel" > ../variables
    cp 2_GermlineEnrichment.sh .. && cd .. && qsub 2_GermlineEnrichment.sh
fi