#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Illumina TruSight Pipeline (Paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="dev"

#TODO ExomeDepth
#TODO add metadata to VCF
#TODO collect QC data
#TODO het vs hom ratio, TsTv, % dbSNP etc

#load run variables
. variables

#merge calls from HC VCFs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V VCFsforFiltering.list \
-L "$version"/"$bedFileName" \
-o "$seqId"_Variants.vcf \
-nt 4 \
-dt NONE

#filter calls
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx2g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-T VariantFiltration \
-o "$seqId"_Variants_Filtered.vcf \
--variant "$seqId"_Variants.vcf \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
-L "$version"/"$bedFileName" \
-dt NONE

#identify runs of homozygosity
/share/apps/bcftools-distros/bcftools-1.2/bcftools roh \
-R "$version"/"$bedFileName" \
"$seqId"_Variants_Filtered.vcf

#log run complete
/share/apps/node-distros/node-v0.12.7-linux-x64/bin/node \
/data/diagnostics/scripts/TrelloAPI.js \
"$seqId" "$worksheetId"