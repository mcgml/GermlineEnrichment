#!/bin/bash -e
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR
mkdir /state/partition1/tmpdir/"$JOB_ID"

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
# Script 3 runs in sample folder

#TODO - BUG - fastqc log file not updated

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Preprocessing ###

#Apply the recalibration to your sequence data
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T PrintReads \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-I "$seqId"_"$sampleId"_realigned.bam \
-BQSR ../"$seqId"_recal_data.table \
-o "$seqId"_"$sampleId".bam \
-nct 12 \
-dt NONE

### Variant calling ###

#SNPs and Indels GVCF with Haplotypecaller
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T HaplotypeCaller \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_"$sampleId".g.vcf \
--genotyping_mode DISCOVERY \
-stand_emit_conf 10 \
-stand_call_conf 30 \
--emitRefConfidence GVCF \
-dt NONE

### QC ###

#Convert BED to interval_list for later
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar BedToIntervalList \
I=/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
O="$bedFileName".interval_list \
SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectInsertSizeMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_insert_metrics.txt \
H="$seqId"_"$sampleId"_insert_metrics.pdf

#HsMetrics: capture & pooling performance
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.5.0/picard.jar CollectHsMetrics \
I="$seqId"_"$sampleId".bam \
O="$seqId"_"$sampleId"_hs_metrics.txt \
R=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
BAIT_INTERVALS="$bedFileName".interval_list \
TARGET_INTERVALS="$bedFileName".interval_list

#Generate per-base coverage: variant detection sensitivity
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx12g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T DepthOfCoverage \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_"$sampleId"_DepthOfCoverage \
-I "$seqId"_"$sampleId".bam \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
--countType COUNT_FRAGMENTS \
--minMappingQuality 20 \
--minBaseQuality 10 \
--omitIntervalStatistics \
-ct "$minimumCoverage" \
-nt 12 \
-dt NONE

#Calculate gene percentage coverage
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx8g -jar /data/diagnostics/apps/CoverageCalculator-2.0.0/CoverageCalculator-2.0.0.jar \
"$seqId"_"$sampleId"_DepthOfCoverage \
/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_genes.txt \
/state/partition1/db/human/refseq/ref_GRCh37.p13_top_level.gff3 \
-p"$spliceSitePadding" \
-d"$minimumCoverage" \
> "$seqId"_"$sampleId"_PercentageCoverage.txt

#Gender analysis using off-targed reads
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools slop \
-i /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
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
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir/"$JOB_ID"/tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T SelectVariants \
--variant /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-o 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
-selectType SNP \
-restrictAllelesTo BIALLELIC \
-env \
-ef \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-XL X -XL Y \
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

#Print QC metrics
echo -e "TotalReads\tTotalTargetUsableBases\tDuplicationRate\tPctSelectedBases\tpctTargetBasesCt\tMeanOnTargetCoverage\tGender\tFreemix\tMeanInsertSize\tSDInsertSize" > "$seqId"_"$sampleId"_qc.txt
echo -e "$totalReads\t$totalTargetedUsableBases\t$duplicationRate\t$pctSelectedBases\t$pctTargetBasesCt\t$meanOnTargetCoverage\t$gender\t$freemix\t$meanInsertSize\t$sdInsertSize" >> "$seqId"_"$sampleId"_qc.txt

#print metaline for final VCF
echo \#\#SAMPLE\=\<ID\="$sampleId",WorklistId\="$worklistId",SeqId\="$seqId",Panel\="$panel",PipelineName\=GermlineEnrichment,PipelineVersion\="$version",MeanInsertSize\="$meanInsertSize",SDInsertSize\="$sdInsertSize",DuplicationRate\="$duplicationRate",TotalReads\="$totalReads",PctSelectedBases\="$pctSelectedBases",MeanOnTargetCoverage\="$meanOnTargetCoverage",pctTargetBasesCt\="$pctTargetBasesCt",Freemix\="$freemix",Gender\="$gender",RemoteBamFilePath\=$(find $PWD -type f -name "$seqId"_"$sampleId".bam)\> > "$seqId"_"$sampleId"_meta.txt

### Clean up ###
rm -r /state/partition1/tmpdir/"$JOB_ID"
rm *realigned.bam *realigned.bai X.off.bed X.bed Y.off.bed Y.bed 1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf \
1kg_highconfidence_autosomal_ontarget_monoallelic_snps.vcf.idx padded.bed 

#create file lists for script 4
find $PWD -name "$seqId"_"$sampleId".g.vcf >> ../GVCFs.list
find $PWD -name "$seqId"_"$sampleId".bam >> ../FinalBams.list

#check if all VCFs are written
if [ $(find .. -maxdepth 1 -mindepth 1 -type d | wc -l | sed 's/^[[:space:]]*//g') -eq $(sort ../GVCFs.list | uniq | wc -l | sed 's/^[[:space:]]*//g') ]; then
    cp 4_GermlineEnrichment.sh .. && cd .. && qsub 4_GermlineEnrichment.sh
fi