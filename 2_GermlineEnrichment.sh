#!/bin/bash
#PBS -l walltime=08:00:00
#PBS -l ncpus=12
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
#cd $PBS_O_WORKDIR

#Description: Germline Illumina TruSight Pipeline (Paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="dev"

#TODO add ExomeDepth output to VCF for DB import
#TODO SNPRelate
#TODO PCA for ancestry

#load run variables
. variables

### Joint variant calling and filtering ###

#Merge calls from HC VCFs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx8g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V VCFsforFiltering.list \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
-o "$seqId"_variants.vcf \
-nt 4 \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \ 
-T SelectVariants \ 
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variants.vcf \
-selectType SNP \ 
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
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
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select INDELs
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \ 
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
-selectType INDEL \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
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
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
-o "$seqId"_indels_filtered.vcf \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_indels_filtered.vcf \
-o "$seqId"_variants_filtered.vcf \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
-dt NONE

#Add metadata to filtered VCF
grep '^##' "$seqId"_variants_filtered.vcf > "$seqId"_variants_filtered_meta.vcf
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_variants_filtered.vcf); do

        #load variables into scope
        . "$sample"/variables

        #add metadata
        echo \#\#SAMPLE\=\<ID\="$sampleId",Tissue\=Blood,WorklistId\="$worklistId",SeqId\="$seqId",Assay\="$panel",PipelineName\=GermlineIlluminaTruSight,PipelineVersion\="$version",RemoteBamFilePath\=$(find $PWD "$sample"/"$seqId"_"$sampleId"_haplotypecaller.bam),RemoteVcfFilePath=$(find $PWD "$seqId"_variants_filtered_meta.vcf)\> >> "$seqId"_variants_filtered_meta.vcf
done
grep -v '^##' "$seqId"_variants_filtered.vcf >> "$seqId"_variants_filtered_meta.vcf

#validate and index final VCF
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants_filtered_meta.vcf \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-dt NONE

### ROH and CNV analysis ###

#Identify CNVs using read-depth
grep -P '^[1-22]' /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed > autosomal.bed
find $PWD -type f -name "*.bam" -not -name "haplotypecaller.bam" > final_bams.list
/share/apps/R-distros/R-3.3.1/bin/Rscript ExomeDepth.R \
-b final_bams.list \
-f /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-r autosomal.bed \
2>&1 | tee log.txt

#print ExomeDepth metrics
correlations=$(grep "Correlation between reference and tests count" log.txt | cut -d' ' -f8)
fragments=$(grep "Number of counted fragments" log.txt | cut -d' ' -f6)
echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_exomedepth.metrics.txt
paste "$RunID"_gt_100x_bams.list $(echo "$fragments") $(echo "$correlations") >> "$seqId"_exomedepth.metrics.txt

#identify runs of homozygosity
/share/apps/htslib-distros/htslib-1.3.1/bgzip -c "$seqId"_variants_filtered.vcf > "$seqId"_variants_filtered.vcf.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -p vcf "$seqId"_variants_filtered.vcf.gz
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_variants_filtered.vcf); do
    /share/apps/bcftools-distros/bcftools-1.3.1/bcftools roh -R /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed -s "$sample" "$seqId"_variants_filtered.vcf.gz | \
    grep -v '^#' | perl bcftools_roh_range.pl | awk '{print $1'\t'$2-1'\t'$3'\t'$5}' > "$sample"/"$sample"_roh.bed
done

### QC ###

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_71/bin/java -Djava.io.tmpdir=tmp -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variant_evaluation.txt \
--eval:"$seqId" "$seqId"_variants_filtered.vcf \
--dbsnp /data/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-L /data/diagnostics/pipelines/GermlineIlluminaTruSight/GermlineIlluminaTruSight-"$version"/"$panel"/"$panel"_Design.bed \
-dt NONE

### Clean up ###
#rm "$seqId"_variants_filtered.vcf.gz
#rm "$seqId"_variants_filtered.vcf.gz.tbi
#rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx
#rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx
#rm "$seqId"_indels.vcf "$seqId"_indels.vcf.idx
#rm "$seqId"_indels_filtered.vcf "$seqId"_indels_filtered.vcf.idx
#rm "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx
#rm "$seqId"_variants_filtered.vcf "$seqId"_variants_filtered.vcf.idx

#log run complete
#/share/apps/node-distros/node-v0.12.7-linux-x64/bin/node \
#/data/diagnostics/scripts/TrelloAPI.js \
#"$seqId" "$worksheetId"