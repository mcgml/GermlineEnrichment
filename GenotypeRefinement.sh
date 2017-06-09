#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="2.0.0"

# Script 3 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

#Apply only family priors to a callset
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T CalculateGenotypePosteriors \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_combined_filtered_100pad.vcf \
--skipPopulationPriors \
-ped "$seqId"_pedigree.ped \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_combined_filtered_100pad_GCP.vcf \
-dt NONE

#filter genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_combined_filtered_100pad_GCP.vcf \
-ped "$seqId"_pedigree.ped \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_combined_filtered_100pad_GCP_filtered.vcf \
-dt NONE

#phase genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T PhaseByTransmission \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_combined_filtered_100pad_GCP_filtered.vcf \
-ped "$seqId"_pedigree.ped \
-o "$seqId"_combined_filtered_100pad_GCP_filtered_phased.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-dt NONE