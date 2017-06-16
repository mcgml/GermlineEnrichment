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
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
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
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_non_snps_filtered.vcf \
-o "$seqId"_combined_filtered_100pad.vcf \
-genotypeMergeOptions UNSORTED \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-ped "$seqId"_pedigree.ped \
-dt NONE

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

#phase genotypes for trios
if [ $(awk '$3 != 0 && $4 != 0' "$seqId"_pedigree.ped | wc -l | sed 's/^[[:space:]]*//g') -gt 0 ]; then
    /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
    -T PhaseByTransmission \
    -R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
    -V "$seqId"_combined_filtered_100pad_GCP.vcf \
    -ped "$seqId"_pedigree.ped \
    -o "$seqId"_combined_filtered_100pad_GCP_phased.vcf \
    --DeNovoPrior 0.000001 \
    -L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
    -ip 100 \
    -mvf "$seqId"_MendelianViolations.txt \
    -dt NONE
else 
    #skip phasing if no trios are present
    cp "$seqId"_combined_filtered_100pad_GCP.vcf "$seqId"_combined_filtered_100pad_GCP_phased.vcf
    cp "$seqId"_combined_filtered_100pad_GCP.vcf.idx "$seqId"_combined_filtered_100pad_GCP_phased.vcf.idx
fi

#filter genotypes
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantFiltration \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_combined_filtered_100pad_GCP_phased.vcf \
-ped "$seqId"_pedigree.ped \
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
--setFilteredGtToNocall \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf \
-dt NONE

### ROH, SV & CNV analysis ###

#identify runs of homozygosity

#bgzip vcf and index with tabix
/share/apps/htslib-distros/htslib-1.4.1/bgzip -c "$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf > "$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf.gz
/share/apps/htslib-distros/htslib-1.4.1/tabix -p vcf "$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf.gz

#Annotate VCF with the AFs and run ROH
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools annotate \
-c CHROM,POS,REF,ALT,AF1KG \
-h /state/partition1/db/human/roh/AFs.tab.gz.hdr \
-a /state/partition1/db/human/roh/AFs.tab.gz \
"$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf.gz | \
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools roh \
-R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
--AF-tag AF1KG \
-M 100 \
-m /state/partition1/db/human/roh/genetic-map/genetic_map_chr{CHROM}_combined_b37.txt \
-o "$seqId"_ROH.txt

#Structural variant calling with Manta
/share/apps/manta-distros/manta-1.1.0.centos5_x86_64/bin/configManta.py \
$(sed 's/^/--bam /' HighCoverageBams.list | tr '\n' ' ') \
--referenceFasta /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--exome \
--runDir manta
manta/runWorkflow.py \
--quiet \
-m local \
-j 12

gzip -dc manta/results/variants/diploidSV.vcf.gz > "$seqId"_sv_filtered.vcf

#make CNV target BED file - sort, merge, increase bins to min 160bp, remove extreme GC & poor mappability bins
awk '{ if ($1 > 0 && $1 < 23) print $1"\t"$2"\t"$3 }' /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools sort -faidx /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta.fai | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
awk '{ len=$3-($2+1); if (len < 160) { slop=(160-len)/2; adjStart=$2-slop; adjEnd=$3+slop; printf "%s\t%.0f\t%.0f\n", $1,adjStart,adjEnd; } else {print $1"\t"$2"\t"$3} }' | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools merge | \
/share/apps/bedtools-distros/bedtools-2.26.0/bin/bedtools nuc -fi /data/db/human/gatk/2.8/b37/human_g1k_v37.fasta -bed - | \
awk '{ if ($5 >= 0.1 && $5 <= 0.9) print "chr"$1,$2,$3,$1"-"$2"-"$3 }' | tr ' ' '\t' > "$panel"_ROI_b37_window_gc.bed
/share/apps/bigWigAverageOverBed-distros/bigWigAverageOverBed /data/db/human/wgEncodeMapability/wgEncodeCrgMapabilityAlign100mer.bigWig "$panel"_ROI_b37_window_gc.bed "$panel"_ROI_b37_window_gc_mappability.txt
awk '{ if ($5 > 0.9 && $6 > 0.9) print $1"\ttarget_"NR }' "$panel"_ROI_b37_window_gc_mappability.txt | tr '-' '\t' > "$panel"_ROI_b37_CNV.bed

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

### Annotation & Reporting ###

#restrict variants to ROI but retain overlapping indels
/share/apps/bcftools-distros/bcftools-1.4.1/bcftools view \
-R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
"$seqId"_combined_filtered_100pad_GCP_phased_gtfiltered.vcf.gz > "$seqId"_combined_filtered.vcf

#Add VCF meta data to final VCFs
addMetaDataToVCF "$seqId"_combined_filtered.vcf
addMetaDataToVCF "$seqId"_sv_filtered.vcf

#annotate with VEP
annotateVCF "$seqId"_combined_filtered_meta.vcf "$seqId"_filtered_meta_annotated.vcf
annotateVCF "$seqId"_sv_filtered_meta.vcf "$seqId"_sv_filtered_meta_annotated.vcf

#add gnomad allele frequencies
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T VariantAnnotator \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_filtered_meta_annotated.vcf \
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
--resourceAlleleConcordance \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_filtered_meta_annotated_gnomad.vcf \
-dt NONE

#add CNV vcf headers, metadata and annotate
for vcf in $(ls *_cnv.vcf); do

    prefix=$(echo "$vcf" | sed 's/\.vcf//g')
    sampleId=$(/share/apps/bcftools-distros/bcftools-1.4.1/bcftools query -l "$vcf")

    #add VCF headers
    /share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar UpdateVcfSequenceDictionary \
    I="$vcf" \
    O="$prefix"_header.vcf \
    SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

    #add metadata, annotate & index
    addMetaDataToVCF "$prefix"_header.vcf
    annotateVCF "$prefix"_header_meta.vcf "$prefix"_meta_annotated.vcf

    #move files to sampleId folder
    mv "$prefix"_meta_annotated.vcf* "$sampleId"

    #delete unused files
    rm "$vcf" "$prefix"_header.vcf "$prefix"_header_meta.vcf
done

### QC ###

#relatedness test
/share/apps/vcftools-distros/vcftools-0.1.14/build/bin/vcftools \
--relatedness2 \
--out "$seqId"_relatedness \
--vcf "$seqId"_combined_filtered_100pad.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_131/bin/java -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx2g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectVariantCallingMetrics \
INPUT="$seqId"_filtered_meta_annotated_gnomad.vcf \
OUTPUT="$seqId"_CollectVariantCallingMetrics.txt \
DBSNP=/state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
THREAD_COUNT=4

### Clean up ###

#delete unused files
rm -r manta
rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx "$panel"_ROI_b37_window_gc_mappability.txt  "$seqId"_combined_filtered_meta.vcf
rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx "$seqId"_non_snps.vcf igv.log
rm "$seqId"_non_snps.vcf.idx "$seqId"_non_snps_filtered.vcf "$seqId"_non_snps_filtered.vcf.idx "$seqId"_combined_filtered_meta.vcf.gz
rm ExomeDepth.log GVCFs.list HighCoverageBams.list "$seqId"_sv_filtered.vcf "$panel"_ROI_b37_window_gc.bed
rm "$seqId"_sv_filtered_meta.vcf BAMs.list variables "$seqId"_combined_filtered.vcf "$seqId"_combined_filtered_meta.vcf.gz.tbi 
rm "$seqId"_combined_filtered_100pad.vcf "$seqId"_combined_filtered_100pad.vcf.idx "$seqId"_combined_filtered_100pad_GCP.vcf
rm "$seqId"_combined_filtered_100pad_GCP.vcf.idx "$seqId"_combined_filtered_100pad_GCP_filtered.vcf
rm "$seqId"_combined_filtered_100pad_GCP_filtered.vcf.idx igv.log