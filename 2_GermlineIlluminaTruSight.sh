#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

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

### Joint variant calling and filtering ###

#Merge calls from HC VCFs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V VCFsforFiltering.list \
-L "$version"/"$bedFileName" \
-o "$seqId"_variants.vcf \
-nt 4 \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \ 
-T SelectVariants \ 
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variants.vcf \
-selectType SNP \ 
-L "$version"/"$bedFileName" \
-o "$seqId"_snps.vcf \
-nt 4 \
-dt NONE

#Filter SNPs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_snps.vcf 
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
-L "$version"/"$bedFileName" \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select INDELs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \ 
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
-selectType INDEL \
-L "$version"/"$bedFileName" \
-o "$seqId"_indels.vcf \
-nt 4 \
-dt NONE

#Filter INDELs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_indels.vcf
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 200.0" \
--filterName "FS" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSum" \
-L "$version"/"$bedFileName" \
-o "$seqId"_indels_filtered.vcf \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_indels_filtered.vcf \
-o "$seqId"_variants_filtered.vcf \
-L "$version"/"$bedFileName" \
-dt NONE

#Add metadata to filtered VCF
#TODO

### ROH and CNV analysis ###

#ExomeDepth
#TODO

#identify runs of homozygosity
/share/apps/htslib-distros/htslib-1.3.1/bgzip -c "$seqId"_variants_filtered.vcf > "$seqId"_variants_filtered.vcf.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -p vcf "$seqId"_variants_filtered.vcf.gz
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_variants_filtered.vcf); do
    /share/apps/bcftools-distros/bcftools-1.3.1/bcftools roh -R "$version"/"$bedFileName" -s "$sample" "$seqId"_variants_filtered.vcf.gz | \
    grep -v '^#' | perl bcftools_roh_range.pl | awk '{print $1'\t'$2-1'\t'$3'\t'$5}' > "$sample"/"$sample"_roh.bed
done

### QC ###

##TODO variantevaluation

### Clean up ###
#rm "$seqId"_variants_filtered.vcf.gz
#rm "$seqId"_variants_filtered.vcf.gz.tbi
#rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx
#rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx
#rm "$seqId"_indels.vcf "$seqId"_indels.vcf.idx
#rm "$seqId"_indels_filtered.vcf "$seqId"_indels_filtered.vcf.idx
#rm "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx

#log run complete
#/share/apps/node-distros/node-v0.12.7-linux-x64/bin/node \
#/data/diagnostics/scripts/TrelloAPI.js \
#"$seqId" "$worksheetId"