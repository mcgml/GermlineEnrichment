#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="2.4.1"

# Script 2 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

addMetaDataToVCF(){
    output=$(echo "$1" | sed 's/\.vcf/_meta\.vcf/g')
    grep '^##' "$1" > "$output"
    for sample in $(/share/apps/bcftools-distros/bcftools-1.4.1/bcftools query -l "$1"); do
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
    --minimal \
    --refseq
    
    #check VEP has produced annotated VCF
    if [ ! -e "$2" ]; then
        cp "$1" "$2"
    fi

    #index annotated VCF
    /share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$2"
}

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
-selectType SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_snps.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Filter SNPs
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
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
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_snps_filtered.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Select non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
--selectTypeToExclude SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_non_snps.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Filter non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_non_snps.vcf \
--filterExpression "QUAL < 30.0" \
--filterName "LowQual" \
--filterExpression "QD < 2.0" \
--filterName "QD" \
--filterExpression "FS > 200.0" \
--filterName "FS" \
--filterExpression "SOR > 10.0" \
--filterName "SOR" \
--filterExpression "ReadPosRankSum < -20.0" \
--filterName "ReadPosRankSum" \
--filterExpression "InbreedingCoeff != 'NaN' && InbreedingCoeff < -0.8" \
--filterName "InbreedingCoeff" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_non_snps_filtered.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_non_snps_filtered.vcf \
-o "$seqId"_variants_filtered.vcf \
-genotypeMergeOptions UNSORTED \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Apply only family priors to a callset
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T CalculateGenotypePosteriors \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants_filtered.vcf \
--skipPopulationPriors \
-ped "$seqId"_pedigree.ped \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants_filtered_genotypes_refined.vcf \
-dt NONE

#filter genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants_filtered_genotypes_refined.vcf \
-ped "$seqId"_pedigree.ped \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants_filtered_genotypes_filtered.vcf \
-dt NONE

#Add VCF meta data to final VCFs
addMetaDataToVCF "$seqId"_variants_filtered_genotypes_filtered.vcf

#annotate with VEP
annotateVCF "$seqId"_variants_filtered_genotypes_filtered_meta.vcf "$seqId"_variants_filtered_genotypes_filtered_meta_vep.vcf

#add gnomad allele frequencies
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants_filtered_genotypes_filtered_meta_vep.vcf \
-A PossibleDeNovo \
--resource:GNOMAD_2.0.1_Genome_chr1 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.1.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr2 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.2.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr3 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.3.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr4 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.4.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr5 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.5.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr6 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.6.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr7 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.7.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr8 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.8.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr9 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.9.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr10 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.10.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr11 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.11.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr12 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.12.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr13 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.13.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr14 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.14.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr15 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.15.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr16 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.16.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr17 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.17.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr18 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.18.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr19 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.19.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr20 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.20.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr21 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.21.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chr22 /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.22.vcf.gz \
--resource:GNOMAD_2.0.1_Genome_chrX /data/db/human/gnomad/gnomad.genomes.r2.0.1.sites.X.vcf.gz \
--resource:GNOMAD_2.0.1_Exome /data/db/human/gnomad/gnomad.exomes.r2.0.1.sites.vcf.gz \
--resource:spidex /data/db/human/spidex/spidex_public_noncommercial_v1_0.b37.vcf.gz \
--resource:mcap /data/db/human/mcap/mcap_v1_0.b37.vcf.gz \
-E GNOMAD_2.0.1_Genome_chr1.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr2.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr3.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr4.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr5.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr6.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr7.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr8.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr9.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr10.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr11.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr12.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr13.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr14.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr15.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr16.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr17.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr18.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr19.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr20.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr21.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chr22.AF_POPMAX \
-E GNOMAD_2.0.1_Genome_chrX.AF_POPMAX \
-E GNOMAD_2.0.1_Exome.AF_POPMAX \
-E spidex.dpsi_max_tissue \
-E spidex.dpsi_zscore \
-E mcap.mcap \
-ped "$seqId"_pedigree.ped \
--resourceAlleleConcordance \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_filtered_annotated_padded.vcf.gz \
-dt NONE

#restrict variants to ROI but retain overlapping indels
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools view \
-R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
"$seqId"_filtered_annotated_padded.vcf.gz > "$seqId"_filtered_annotated_roi.vcf

#validate final VCF
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.8.0/GenomeAnalysisTK.jar \
-T ValidateVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_filtered_annotated_roi.vcf \
--reference_window_stop 300 \
-dt NONE

#report variants to text
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx48g -jar /data/diagnostics/apps/VariantReporterSpark/VariantReporterSpark-1.3.0/VariantReporterSpark.jar \
-V "$seqId"_filtered_annotated_roi.vcf \
-P "$seqId"_pedigree.ped \
-T 12 \
-N

### CNV analysis ###

#check one or more samples have high coverage
if [[ -e "HighCoverageBams.list" ]] && [[ $(wc -l "HighCoverageBams.list" | awk '{print $1}') -gt 4 ]]; then

    #make CNV bed
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools slop \
    -i /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    -g /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai \
    -b 250 | \
    grep -v ^X| \
    grep -v ^Y | \
    grep -v ^MT | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort \
    -faidx /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
    /share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools subtract \
    -A -a - -b /data/db/human/genomicSuperDups/genomicSuperDups_b37.bed.gz | \
    awk -F"\t" '{print $1"\t"$2"\t"$3"\tr"NR}' > "$panel"_ROI_b37_CNV.bed

    #call CNVs using read depth
    /share/apps/R-distros/R-3.3.1/bin/Rscript /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/ExomeDepth.R \
    -b HighCoverageBams.list \
    -f /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -r "$panel"_ROI_b37_CNV.bed \
    2>&1 | tee ExomeDepth.log

    #print ExomeDepth metrics
    echo -e "BamPath\tFragments\tCorrelation" > "$seqId"_ExomeDepth_Metrics.txt
    paste HighCoverageBams.list \
    <(grep "Number of counted fragments" ExomeDepth.log | cut -d' ' -f6) \
    <(grep "Correlation between reference and tests count" ExomeDepth.log | cut -d' ' -f8) >> "$seqId"_ExomeDepth_Metrics.txt

    #add CNV vcf headers and move to sample folder
    for vcf in $(ls *_cnv.vcf); do

        prefix=$(echo "$vcf" | sed 's/\.vcf//g')
        sampleId=$(/share/apps/bcftools-distros/bcftools-1.4.1/bcftools query -l "$vcf")

        #add VCF headers
        /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-2.12.2/picard.jar UpdateVcfSequenceDictionary \
        I="$vcf" \
        O="$sampleId"/"$seqId"_"$sampleId"_cnv.vcf \
        SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

        #gzip and tabix
        /share/apps/htslib-distros/htslib-1.4.1/bgzip "$sampleId"/"$seqId"_"$sampleId"_cnv.vcf
        /share/apps/htslib-distros/htslib-1.4.1/tabix -p vcf "$sampleId"/"$seqId"_"$sampleId"_cnv.vcf.gz

    done

    #move cnv metrics to sample folder
    for i in $(ls *cnv.txt); do
        mv "$i" $(echo "$i" | cut -d_ -f1)/"$seqId"_"$i";
    done

fi

### QC ###

#relatedness test
/share/apps/vcftools-distros/vcftools-0.1.14/build/bin/vcftools \
--relatedness2 \
--out "$seqId"_relatedness \
--vcf "$seqId"_variants_filtered.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-2.12.2/picard.jar CollectVariantCallingMetrics \
INPUT="$seqId"_filtered_annotated_padded.vcf.gz \
OUTPUT="$seqId"_CollectVariantCallingMetrics.txt \
DBSNP=/state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
THREAD_COUNT=8

### Clean up ###

#delete unused files
rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx "$seqId"_non_snps.vcf
rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx 
rm "$seqId"_non_snps.vcf.idx "$seqId"_non_snps_filtered.vcf "$seqId"_non_snps_filtered.vcf.idx
rm GVCFs.list igv.log BAMs.list variables
rm "$seqId"_variants_filtered_genotypes_filtered_meta.vcf "$seqId"_variants_filtered_genotypes_filtered_meta_vep.vcf
rm "$seqId"_variants_filtered_genotypes_filtered_meta_vep.vcf.idx "$seqId"_variants_filtered_genotypes_filtered.vcf
rm "$seqId"_variants_filtered_genotypes_filtered.vcf.idx
rm "$seqId"_variants_filtered_genotypes_refined.vcf "$seqId"_variants_filtered_genotypes_refined.vcf.idx
rm "$seqId"_variants_filtered.vcf "$seqId"_variants_filtered.vcf.idx
rm -f ExomeDepth.log HighCoverageBams.list "$seqId"_*_cnv.vcf