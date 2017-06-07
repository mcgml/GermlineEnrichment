#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="1.8.9"

# Script 2 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Build the SNP recalibration model
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_variants.vcf \
-resource:hapmap,known=false,training=true,truth=true,prior=15.0 /state/partition1/db/human/gatk/2.8/b37/hapmap_3.3.b37.vcf \
-resource:omni,known=false,training=true,truth=true,prior=12.0 /state/partition1/db/human/gatk/2.8/b37/1000G_omni2.5.b37.vcf \
-resource:1000G,known=false,training=true,truth=false,prior=10.0 /state/partition1/db/human/gatk/2.8/b37/1000G_phase1.snps.high_confidence.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-an QD \
-an FS \
-an SOR \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode SNP \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-rscriptFile "$seqId"_SNP_plots.R \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-dt NONE

#Apply the desired level of recalibration to the SNPs in the call set
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_variants.vcf \
-mode SNP \
--ts_filter_level 99.0 \
-recalFile "$seqId"_SNP.recal \
-tranchesFile "$seqId"_SNP.tranches \
-o "$seqId"_recalibrated_snps_raw_indels.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-dt NONE

#Build the Indel recalibration model
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantRecalibrator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-resource:mills,known=false,training=true,truth=true,prior=12.0 /state/partition1/db/human/gatk/2.8/b37/Mills_and_1000G_gold_standard.indels.b37.vcf \
-resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-an QD \
-an FS \
-an SOR \
-an MQRankSum \
-an ReadPosRankSum \
-an InbreedingCoeff \
-mode INDEL \
-tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \
--maxGaussians 4 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-rscriptFile "$seqId"_INDEL_plots.R \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-dt NONE

#Apply the desired level of recalibration to the Indels in the call set
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T ApplyRecalibration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-input "$seqId"_recalibrated_snps_raw_indels.vcf \
-mode INDEL \
--ts_filter_level 99.0 \
-recalFile "$seqId"_INDEL.recal \
-tranchesFile "$seqId"_INDEL.tranches \
-o "$seqId"_recalibrated_variants.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-dt NONE

#Filter genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_recalibrated_variants.vcf \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
--setFilteredGtToNocall \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_combined_filtered_100pad.vcf \
-dt NONE

#restrict to ROI but retain overlapping indels
/share/apps/htslib-distros/htslib-1.4/bgzip "$seqId"_combined_filtered_100pad.vcf
/share/apps/htslib-distros/htslib-1.4/tabix -p vcf "$seqId"_combined_filtered_100pad.vcf.gz
/share/apps/bcftools-distros/bcftools-1.4/bcftools view \
-R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
"$seqId"_combined_filtered_100pad.vcf.gz > "$seqId"_combined_filtered.vcf