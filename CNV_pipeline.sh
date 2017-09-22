#!/bin/bash
#PBS -l walltime=40:00:00
#PBS -l ncpus=12
set -euo pipefail
#PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_SAMPLE
version="1.8.9"

# Script runs in panel folder and requires final Bams

#load sample & pipeline variables
. *.variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

convertVcfToBed(){
    /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
    -T VariantsToTable \
    -R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -V "$1" \
    -F CHROM -F POS -F END -F FILTER \
    -o $(echo "$1" | sed 's/\.vcf/\.txt/g') \
    -raw \
    -dt NONE

    awk -F"\t" 'NR>1 {print $1"\t"$2-1"\t"$3"\t"$4}' $(echo "$1" | sed 's/\.vcf/\.txt/g') > $(echo "$1" | sed 's/\.vcf/\.bed/g')
    rm $(echo "$1" | sed 's/\.vcf/\.txt/g')
}

### Prepare BED for CNV dosage ###

#counts reads per region
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T DiagnoseTargets \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
$(sed 's/^/-I /' HighCoverageBams.list | tr '\n' ' ') \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 230 \
--bad_mate_status_threshold 0.05 \
--coverage_status_threshold 0.33 \
--minimum_base_quality 10 \
--minimum_coverage 10 \
--maximum_coverage 500 \
-o "$panel"_CNV_b37.vcf \
-dt NONE

#calculate average mappability of regions
convertVcfToBed "$panel"_CNV_b37.vcf
awk -F"\t" '{print "chr"$1"\t"$2"\t"$3"\tr"NR}' "$panel"_CNV_b37.bed > "$panel"_CNV_hg19.bed
/share/apps/bigWigAverageOverBed-distros/bigWigAverageOverBed \
/data/db/human/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig \
-bedOut="$panel"_CNV_hg19_wgEncodeCrgMapabilityAlign100mer.bed "$panel"_CNV_hg19.bed "$panel"_CNV_hg19_wgEncodeCrgMapabilityAlign100mer.txt
sed 's/^chr//g' "$panel"_CNV_hg19_wgEncodeCrgMapabilityAlign100mer.bed | \
/share/apps/htslib-distros/htslib-1.4.1/bgzip -c > "$panel"_CNV_b37_wgEncodeCrgMapabilityAlign100mer.bed.gz
/share/apps/htslib-distros/htslib-1.4.1/tabix -p bed "$panel"_CNV_b37_wgEncodeCrgMapabilityAlign100mer.bed.gz

#annotate with mappability score
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools annotate \
-a "$panel"_CNV_b37_wgEncodeCrgMapabilityAlign100mer.bed.gz \
-h <(echo "##INFO=<ID=wgEncodeCrgMapabilityAlign100mer,Number=1,Type=Float,Description=\"Average wgEncodeCrgMapabilityAlign100mer score\">") \
-c CHROM,FROM,TO,-,wgEncodeCrgMapabilityAlign100mer \
"$panel"_CNV_b37.vcf > "$panel"_CNV_b37_wgEncodeCrgMapabilityAlign100mer.vcf

#calculate % overlap with repeat masker
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools coverage \
-a "$panel"_CNV_b37.bed \
-b /data/db/human/rmsk/rmsk_b37_sorted.bed | awk -F"\t" '{print $1"\t"$2"\t"$3"\tr"NR"\t"$8}' | \
/share/apps/htslib-distros/htslib-1.4.1/bgzip -c > "$panel"_CNV_b37_repeatmasker.bed.gz
/share/apps/htslib-distros/htslib-1.4.1/tabix -p bed "$panel"_CNV_b37_repeatmasker.bed.gz

#annotate with % overlap with repeat masker
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools annotate \
-a "$panel"_CNV_b37_repeatmasker.bed.gz \
-h <(echo "##INFO=<ID=RepeatMaskerOverlap,Number=1,Type=Float,Description=\"Pct overlap with RepeatMasker\">") \
-c CHROM,FROM,TO,-,RepeatMaskerOverlap \
"$panel"_CNV_b37_wgEncodeCrgMapabilityAlign100mer.vcf > "$panel"_CNV_b37_RepeatMaskerOverlap.vcf

#calculate GC of unpadded ROI
rm -f "$panel"_CNV_b37_GC.bed
cut -s -f1-3 "$panel"_CNV_b37.bed | while read line; do
    
    #extract targets from padded ROI
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools intersect \
    -a /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    -b <(echo "$line") | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools nuc \
    -bed - \
    -fi /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta | \
    awk -vl="$line" -F"\t" 'NR > 1 {j++;n+=$5} END {print l"\t\t"(n/j)}' >> "$panel"_CNV_b37_GC.bed

done

#compress and index GC bed
/share/apps/htslib-distros/htslib-1.4.1/bgzip -f "$panel"_CNV_b37_GC.bed
/share/apps/htslib-distros/htslib-1.4.1/tabix -p bed "$panel"_CNV_b37_GC.bed.gz

#annotate GC of unpadded ROI
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools annotate \
-a "$panel"_CNV_b37_GC.bed.gz \
-h <(echo "##INFO=<ID=UnpaddedGCContent,Number=1,Type=Float,Description=\"GC content of the unpadded interval\">") \
-c CHROM,FROM,TO,-,UnpaddedGCContent \
"$panel"_CNV_b37_RepeatMaskerOverlap.vcf > "$panel"_CNV_b37_GC.vcf

#annotate features with mapping quality
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$panel"_CNV_b37_GC.vcf \
$(sed 's/^/-I /' HighCoverageBams.list | tr '\n' ' ') \
-o "$panel"_CNV_b37_annotated.vcf \
-A MappingQualityZero -A RMSMappingQuality -A LowMQ \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 230 \
-dt NONE

#filter regions
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "IGC < 0.15 || UnpaddedGCContent < 0.15" \
--filterName "LowGC" \
--filterExpression "IGC > 0.85 || UnpaddedGCContent > 0.85" \
--filterName "HighGC" \
--filterExpression "wgEncodeCrgMapabilityAlign100mer < 0.85" \
--filterName "wgEncodeCrgMapability" \
--filterExpression "LowMQ[1] > 0.1" \
--filterName "LowMQ" \
--filterExpression "CHROM == 'X' || CHROM == 'Y' || CHROM == 'MT'" \
--filterName "NonDiploidContig" \
--filterExpression "RepeatMaskerOverlap > 0.25" \
--filterName "RepeatMaskerOverlap" \
-V "$panel"_CNV_b37_annotated.vcf \
-o "$panel"_CNV_b37_annotated_unmasked.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 230 \
-dt NONE

#mask HiSeq high depth
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$panel"_CNV_b37_annotated_unmasked.vcf \
-o "$panel"_CNV_b37_annotated_hiseq_masked.vcf \
--mask /data/db/human/hiSeqDepth/hiSeqDepthTop5Pct_b37_sorted.bed.gz \
--maskName hiSeqDepthTop5Pct \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 230 \
-dt NONE

#mask segmental duplications
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$panel"_CNV_b37_annotated_hiseq_masked.vcf \
-o "$panel"_CNV_b37_annotated_filtered.vcf \
--mask /data/db/human/genomicSuperDups/genomicSuperDups_b37.bed.gz \
--maskName genomicSuperDups \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 230 \
-dt NONE

#convert final vcf to bed
convertVcfToBed "$panel"_CNV_b37_annotated_filtered.vcf
grep PASS "$panel"_CNV_b37_annotated_filtered.bed | \
awk -F"\t" '{print $1"\t"$2"\t"$3"\tr"NR}' > "$panel"_CNV_b37_annotated_filtered_pass.bed

### CNV calling ###

#call CNVs using read depth
/share/apps/R-distros/R-3.3.1/bin/Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/ExomeDepth.R \
-b HighCoverageBams.list \
-f /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-r "$panel"_CNV_b37_annotated_filtered_pass.bed \
2>&1 | tee ExomeDepth.log

#print ExomeDepth metrics
echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_ExomeDepth_Metrics.txt
paste HighCoverageBams.list \
<(grep "Number of counted fragments" ExomeDepth.log | cut -d' ' -f6) \
<(grep "Correlation between reference and tests count" ExomeDepth.log | cut -d' ' -f8) >> "$seqId"_ExomeDepth_Metrics.txt

#add VCF headers and index
for vcf in $(ls *cnv.vcf); do
    /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-2.10.10/picard.jar UpdateVcfSequenceDictionary \
    I="$vcf" \
    O=$(echo "$vcf" | sed 's/\.vcf/_header\.vcf/g') \
    SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict
done

#annotate vcf calls with allele balance & filter
#annotate cnvs with probe GC

### SV calling ###

#Structural variant calling with Manta
while read bam; do
    /share/apps/manta-distros/manta-1.1.0.centos5_x86_64/bin/configManta.py \
    --bam "$bam" \
    --referenceFasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    --exome \
    --runDir manta
    manta/runWorkflow.py \
    --quiet \
    -m local \
    -j 12

    gzip -dc manta/results/variants/diploidSV.vcf.gz > $(echo "$bam" | sed 's/\.bam/_sv\.vcf/g')
    rm -r manta
done < HighCoverageBams.list

### interpretation ###

#decipher (DDG2P)
#dgv
#ISCA (clingen hap/trip insufficency score)
#OMIM
#pubmed (gene/exon/chromosomal loc)

#clean up
#TODO