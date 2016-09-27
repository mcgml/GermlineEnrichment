#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="dev"

#TODO SNPRelate
#TODO PCA for ancestry
#TODO upgrade dbSNP for variant evaluation
#TODO process ED output
#TODO annotate output
#TODO optimise ED

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

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--dbsnp /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_variants.vcf \
-nt 12 \
-dt NONE

#Annotate with low complexity region length using mdust
/share/apps/bcftools-distros/bcftools-1.3.1/bcftools annotate \
-a /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.mdust.v34.lpad1.bed.gz \
-c CHROM,FROM,TO,LCRLen \
-h <(echo '##INFO=<ID=LCRLen,Number=1,Type=Integer,Description="Overlapping mdust low complexity region length (mask cutoff: 34)">') \
-o "$seqId"_variants.lcr.vcf \
"$seqId"_variants.vcf

#Select SNPs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.lcr.vcf \
-selectType SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_snps.vcf \
-nt 12 \
-dt NONE

#Filter SNPs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_snps.vcf \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 60.0" \
--filterName "FS" \
--filterExpression "SOR > 4.0" \
--filterName "SOR" \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "MQRankSum < -12.5" \
--filterName "MQRankSum" \
--filterExpression "ReadPosRankSum < -8.0" \
--filterName "ReadPosRankSum" \
--genotypeFilterExpression "DP < 20" \
--filterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--filterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select INDELs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.lcr.vcf \
-selectType INDEL \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_indels.vcf \
-nt 12 \
-dt NONE

#Filter INDELs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_indels.vcf \
--filterExpression "QD < 3.5" \
--filterName "QD" \
--filterExpression "FS > 170.0" \
--filterName "FS" \
--filterExpression "SOR > 8.0" \
--filterName "SOR" \
--filterExpression "MQ < 40.0" \
--filterName "MQ" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSum" \
--filterExpression "InbreedingCoeff < -0.7" \
--filterName "InbreedingCoeff" \
--filterExpression "LCRLen > 8" \
--filterName "LowComplexity" \
--genotypeFilterExpression "DP < 20" \
--filterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--filterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-o "$seqId"_indels_filtered.vcf \
-dt NONE

#Structural variant calling with Manta
/share/apps/manta-distros/manta-1.0.0.centos5_x86_64/bin/configManta.py \
$(sed 's/^/--bam /' FinalBams.list | tr '\n' ' ') \
--referenceFasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--exome \
--runDir manta
manta/runWorkflow.py \
-m local \
-j 12

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_indels_filtered.vcf \
--variant manta/results/variants/diploidSV.vcf.gz \
-o "$seqId"_variants_filtered.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-nt 12 \
-genotypeMergeOptions UNSORTED \
-dt NONE

#Add VCF meta data to final VCF
grep '^##' "$seqId"_variants_filtered.vcf > "$seqId"_filtered_meta.vcf
for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_variants_filtered.vcf); do
    cat "$sample"/"$seqId"_"$sample"_meta.txt >> "$seqId"_filtered_meta.vcf
done
grep -v '^##' "$seqId"_variants_filtered.vcf >> "$seqId"_filtered_meta.vcf

### QC ###

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.6.0/GenomeAnalysisTK.jar \
-T VariantEval \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-o "$seqId"_variant_evaluation.txt \
--eval "$seqId"_filtered_meta.vcf \
--dbsnp /state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.vcf \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
-nt 12 \
-dt NONE

### ROH & CNV analysis ###

#Identify autosomal CNVs using read-depth
awk '{if ($1 > 0 && $1 < 23) print $1"\t"$2"\t"$3"\tbin"NR}' \
/data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed \
> autosomal.bed

/share/apps/R-distros/R-3.3.1/bin/Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/ExomeDepth.R \
-b FinalBams.list \
-f /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-r autosomal.bed \
2>&1 | tee ExomeDepth.log

#print ExomeDepth metrics
echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_exomedepth.metrics.txt
paste FinalBams.list \
<(grep "Number of counted fragments" ExomeDepth.log | cut -d' ' -f6) \
<(grep "Correlation between reference and tests count" ExomeDepth.log | cut -d' ' -f8) \
>> "$seqId"_exomedepth.metrics.txt

#identify runs of homozygosity
/share/apps/htslib-distros/htslib-1.3.1/bgzip -c "$seqId"_filtered_meta.vcf > "$seqId"_filtered_meta.vcf.gz
/share/apps/htslib-distros/htslib-1.3.1/tabix -p vcf "$seqId"_filtered_meta.vcf.gz

for sample in $(/share/apps/bcftools-distros/bcftools-1.3.1/bcftools query -l "$seqId"_filtered_meta.vcf); do
    /share/apps/bcftools-distros/bcftools-1.3.1/bcftools roh -R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI.bed -s "$sample" "$seqId"_filtered_meta.vcf.gz | \
    grep -v '^#' | perl /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/bcftools_roh_range.pl | grep -v '#' | awk '{print $1"\t"$2-1"\t"$3"\t"$5}' > "$sample"/"$sample"_roh.bed
done

### Funcitonal annotation ###

#annotate VCF with VEP
perl /share/apps/vep-distros/ensembl-tools-release-85/scripts/variant_effect_predictor/variant_effect_predictor.pl \
--cache \
--fasta /share/apps/vep-distros/ensembl-tools-release-85/scripts/variant_effect_predictor/annotations/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \
--dir /share/apps/vep-distros/ensembl-tools-release-85/scripts/variant_effect_predictor/annotations \
--no_progress \
--everything \
--fork 12 \
--species homo_sapiens \
--assembly GRCh37 \
--format vcf \
--no_stats \
--offline \
--refseq \
--allele_number \
--no_escape \
--shift_hgvs 1 \
--vcf \
--no_intergenic \
-i "$seqId"_filtered_meta.vcf \
-o "$seqId"_filtered_meta_annotated.vcf

### Clean up ###

#delete unused files
rm -r manta
rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx "$seqId"_variants.lcr.vcf "$seqId"_variants.lcr.vcf.idx \
"$seqId"_snps.vcf "$seqId"_snps.vcf.idx "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx "$seqId"_indels.vcf \
"$seqId"_indels.vcf.idx "$seqId"_indels_filtered.vcf "$seqId"_indels_filtered.vcf.idx "$seqId"_variants_filtered.vcf \
"$seqId"_variants_filtered.vcf.idx "$seqId"_genotypes_filtered.vcf "$seqId"_genotypes_filtered.vcf.idx "$seqId"_filtered_meta.vcf.gz \
"$seqId"_filtered_meta.vcf.gz.tbi autosomal.bed ExomeDepth.log GVCFs.list FinalBams.list

#log with Trello
/share/apps/node-distros/node-v4.4.7-linux-x64/bin/node \
/data/diagnostics/scripts/TrelloAPI.js \
"$seqId" "Pipeline complete. Ready for bioinformatics check"