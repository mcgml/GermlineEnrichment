#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="1.2.3"

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
# Script 2 runs in panel folder, requires final Bams &gVCFs

addMetaDataToVCF(){
    output=$(echo "$1" | sed 's/\.vcf/_meta\.vcf/g')
    grep '^##' "$1" > "$output"
    for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$1"); do
        cat "$sample"/"$seqId"_"$sample"_meta.txt >> "$output"
    done
    grep -v '^##' "$1" >> "$output"
}

annotateVCF(){
    #annotate VCF
    perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    --verbose \
    --no_progress \
    --everything \
    --fork 12 \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file "$1" \
    --format vcf \
    --output_file "$2" \
    --force_overwrite \
    --no_stats \
    --cache \
    --dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --no_intergenic \
    --offline \
    --cache_version 86 \
    --allele_number \
    --no_escape \
    --shift_hgvs 1 \
    --vcf \
    --refseq
    
    #check VEP has produced annotated VCF
    if [ ! -e "$2" ]; then
        cp "$1" "$2"
    fi
}

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_variants.vcf \
-dt NONE

#Annotate with low complexity region length using mdust
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools annotate \
-a /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \
-c CHROM,FROM,TO,LCRLen \
-h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
-o "$seqId"_variants.lcr.vcf \
"$seqId"_variants.vcf

#Select SNPs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.lcr.vcf \
-selectType SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_snps.vcf \
-dt NONE

#Filter SNPs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_snps.vcf \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 60.0" \
--filterName "FS" \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSum" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSum" \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select INDELs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.lcr.vcf \
-selectType INDEL \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_indels.vcf \
-dt NONE

#Filter INDELs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_indels.vcf \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 200.0" \
--filterName "FS" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSum" \
--filterExpression "LCRLen > 8" \
--filterName "LowComplexity" \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_indels_filtered.vcf \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_indels_filtered.vcf \
-o "$seqId"_filtered.vcf \
-nt 12 \
-genotypeMergeOptions UNSORTED \
-dt NONE

#Add VCF meta data to final VCF
addMetaDataToVCF "$seqId"_filtered.vcf

### ROH, SV & CNV analysis ###

#Structural variant calling with Manta
/share/apps/manta-distros/manta-1.0.3.centos5_x86_64/bin/configManta.py \
$(sed 's/^/--bam /' Bams.list | tr '\n' ' ') \
--referenceFasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--exome \
--runDir manta
manta/runWorkflow.py \
--quiet \
-m local \
-j 12

#unzip VCF
gzip -dc manta/results/variants/diploidSV.vcf.gz > "$seqId"_sv_filtered.vcf

#Add VCF meta data to SV VCF
addMetaDataToVCF "$seqId"_sv_filtered.vcf

#make CNV target BED file - sort, merge, increase bins to min 160bp, remove extreme GC & poor mappability bins
awk '{ if ($1 > 0 && $1 < 23) print $1"\t"$2"\t"$3 }' /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort -faidx /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
awk '{ len=$3-($2+1); if (len < 160) { slop=(160-len)/2; adjStart=$2-slop; adjEnd=$3+slop; printf "%s\t%.0f\t%.0f\n", $1,adjStart,adjEnd; } else {print $1"\t"$2"\t"$3} }' | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools nuc -fi /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta -bed - | \
awk '{ if ($5 >= 0.1 && $5 <= 0.9) print "chr"$1,$2,$3,$1"-"$2"-"$3 }' | tr ' ' '\t' > "$panel"_ROI_b37_window_gc.bed
/share/apps/bigWigAverageOverBed-distros/bigWigAverageOverBed /data/db/human/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig "$panel"_ROI_b37_window_gc.bed "$panel"_ROI_b37_window_gc_mappability.txt
awk '{ if ($5 > 0.95 && $6 > 0.95) print $1"\ttarget_"NR }' "$panel"_ROI_b37_window_gc_mappability.txt | tr '-' '\t' > "$panel"_ROI_b37_CNV.bed

#call CNVs using read depth
/share/apps/R-distros/R-3.3.1/bin/Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/ExomeDepth.R \
-b HighCoverageBams.list \
-f /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-r "$panel"_ROI_b37_CNV.bed \
2>&1 | tee ExomeDepth.log

#convert ExomeDepth output to BED format & move to sample folder
for i in $(ls X*bam.txt); do
    filename=$(echo "$i" | cut -c2- | sed 's/\./-/g' | sed 's/-bam-txt//g')
    sample=$(echo "$filename" | cut -d_ -f5)
    
    #Chrom,Start,Stop,Call;DQ,BF
    grep -v start "$i" | sed '/^$/d' | awk '{print $7"\t"$5-1"\t"$6"\t"$3";"$12"\t"$9}' > "$sample"/"$filename"_cnv.bed
    grep -v start "$i" | sed '/^$/d' | awk '{print $7"\t"$5-1"\t"$6"\t"$3}' > "$filename"_cnv.bed

    #annotate bed
    perl /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/variant_effect_predictor.pl \
    --verbose \
    --no_progress \
    --numbers \
    --symbol \
    --no_intergenic \
    --species homo_sapiens \
    --assembly GRCh37 \
    --input_file "$filename"_cnv.bed \
    --output_file "$sample"/"$filename"_cnv_annotated.txt \
    --force_overwrite \
    --no_stats \
    --cache \
    --dir /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --fasta /share/apps/vep-distros/ensembl-tools-release-86/scripts/variant_effect_predictor/annotations \
    --offline \
    --cache_version 86 \
    --fields Location,Allele,SYMBOL,Feature,EXON,INTRON \
    --refseq

    rm "$i"
    rm "$filename"_cnv.bed
done

#print ExomeDepth metrics
echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_exomedepth.metrics.txt
paste HighCoverageBams.list \
<(grep "Number of counted fragments" ExomeDepth.log | cut -d' ' -f6) \
<(grep "Correlation between reference and tests count" ExomeDepth.log | cut -d' ' -f8) >> "$seqId"_exomedepth.metrics.txt

#identify runs of homozygosity
/share/apps/htslib-distros/htslib-1.3.1/bgzip -c "$seqId"_filtered_meta.vcf > "$seqId"_filtered_meta.vcf.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -p vcf "$seqId"_filtered_meta.vcf.gz
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_filtered_meta.vcf); do
    /share/apps/bcftools-distros/bcftools-dev/bcftools roh -O r -s "$sample" -R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed "$seqId"_filtered_meta.vcf.gz | \
    grep -v '^#' | awk '{print $3"\t"$4-1"\t"$5"\t\t"$8}' > "$sample"/"$seqId"_"$sample"_roh.bed
done

### Annotation & Reporting ###

#annotate with VEP
annotateVCF "$seqId"_filtered_meta.vcf "$seqId"_filtered_meta_annotated.vcf
annotateVCF "$seqId"_sv_filtered_meta.vcf "$seqId"_sv_filtered_meta_annotated.vcf

#index annotated VCFs
/share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$seqId"_filtered_meta_annotated.vcf
/share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$seqId"_sv_filtered_meta_annotated.vcf

#write SNV & Indel dataset to table
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -jar /data/diagnostics/apps/VCFParse/VCFParse-1.2.0/VCFParse.jar \
-V "$seqId"_filtered_meta_annotated.vcf \
-O snv_indel \
-K

#write SV dataset to table
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -jar /data/diagnostics/apps/VCFParse/VCFParse-1.2.0/VCFParse.jar \
-V "$seqId"_sv_filtered_meta_annotated.vcf \
-O sv \
-K

### QC ###

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variant_evaluation.txt \
--eval:"$seqId" "$seqId"_filtered_meta_annotated.vcf \
--comp:omni2.5 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
--comp:hapmap3.3 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-nt 12 \
-dt NONE

### Clean up ###

#delete unused files
rm -r manta
rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx "$seqId"_variants.lcr.vcf "$seqId"_variants.lcr.vcf.idx "$panel"_ROI_b37_window_gc_mappability.txt
rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx "$seqId"_indels.vcf igv.log
rm "$seqId"_indels.vcf.idx "$seqId"_indels_filtered.vcf "$seqId"_indels_filtered.vcf.idx "$seqId"_filtered.vcf "$seqId"_filtered.vcf.idx
rm "$seqId"_filtered_meta.vcf.gz "$seqId"_filtered_meta.vcf.gz.tbi ExomeDepth.log GVCFs.list HighCoverageBams.list "$seqId"_sv_filtered.vcf "$panel"_ROI_b37_window_gc.bed 
rm "$seqId"_filtered_meta.vcf "$seqId"_sv_filtered_meta.vcf BAMs.list