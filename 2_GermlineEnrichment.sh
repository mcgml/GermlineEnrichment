#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="dev"

#TODO add ExomeDepth output to VCF for DB import
#TODO SNPRelate
#TODO PCA for ancestry
#TODO replace HRun
#TODO filter on repetitive sequences


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
# Script 2 runs in panel folder

phoneTrello()
{
    /share/apps/node-distros/node-v0.12.7-linux-x64/bin/node \
    /data/diagnostics/scripts/TrelloAPI.js \
    "$1" "$2"
}

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx16g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V VCFsforFiltering.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_variants.vcf \
--annotation TandemRepeatAnnotator \
--annotation HomopolymerRun \
--reference_window_stop 20 \
-nt 8 \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx16g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
-selectType SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_snps.vcf \
-nt 8 \
-dt NONE

#Filter SNPs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_snps.vcf \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 60.0" \
--filterName "FS" \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "SOR > 3.0" \
--filterName "SOR" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSum" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSum" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select INDELs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx16g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
-selectType INDEL \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_indels.vcf \
-nt 8 \
-dt NONE

#Filter INDELs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_indels.vcf \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 200.0" \
--filterName "FS" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSum" \
--filterExpression "InbreedingCoeff < -0.6" \
--filterName "InbreedingCoeff" \
--filterExpression "SOR > 10.0" \
--filterName "SOR" \
--filterExpression "HRun > 8" \
--filterName "HRun" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_indels_filtered.vcf \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_indels_filtered.vcf \
-o "$seqId"_variants_filtered.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-nt 8 \
-genotypeMergeOptions UNSORTED \
-dt NONE

#Filter low GQ calls
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants_filtered.vcf \
--genotypeFilterExpression "GQ < 20.0" \
--genotypeFilterName "LowGQ" \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
-o "$seqId"_genotypes_filtered.vcf \
-dt NONE

#Add VCF meta data
grep '^##' "$seqId"_genotypes_filtered.vcf > "$seqId"_genotypes_filtered_meta.vcf
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_genotypes_filtered.vcf); do
    cat "$sample"/"$seqId"_"$sample"_meta.txt >> "$seqId"_genotypes_filtered_meta.vcf
done
grep -v '^##' "$seqId"_genotypes_filtered.vcf >> "$seqId"_genotypes_filtered_meta.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variant_evaluation.txt \
--eval "$seqId"_genotypes_filtered_meta.vcf \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-nt 8 \
-dt NONE

### ROH and CNV analysis ###

#identify runs of homozygosity
/share/apps/htslib-distros/htslib-1.3.1/bgzip -c "$seqId"_genotypes_filtered.vcf > "$seqId"_genotypes_filtered.vcf.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -p vcf "$seqId"_genotypes_filtered.vcf.gz
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_genotypes_filtered.vcf); do
    /share/apps/bcftools-distros/bcftools-1.3.1/bcftools roh -R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed -s "$sample" "$seqId"_genotypes_filtered.vcf.gz | \
    grep -v '^#' | perl /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/bcftools_roh_range.pl | grep -v '#' | awk '{print $1"\t"$2-1"\t"$3"\t"$5}' > "$sample"/"$sample"_roh.bed
done

#Identify CNVs using read-depth
#grep -P '^[1-22]' /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed > autosomal.bed
#/share/apps/R-distros/R-3.3.1/bin/Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/ExomeDepth.R \
#-b FinalBAMs.list \
#-f /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
#-r autosomal.bed \
#2>&1 | tee log.txt

#print ExomeDepth metrics
echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_exomedepth.metrics.txt
paste FinalBAMs.list \
<(grep "Number of counted fragments" log.txt | cut -d' ' -f6) \
<(grep "Correlation between reference and tests count" log.txt | cut -d' ' -f8) \
>> "$seqId"_exomedepth.metrics.txt

### Clean up ###
rm -r tmp

#log run complete
phoneTrello "$seqId" "Analysis complete"