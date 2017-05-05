#!/bin/bash
#PBS -l walltime=12:00:00
#PBS -l ncpus=12
set -euo pipefail
PBS_O_WORKDIR=(`echo $PBS_O_WORKDIR | sed "s/^\/state\/partition1//" `)
cd $PBS_O_WORKDIR

#Description: Germline Enrichment Pipeline (Illumina paired-end). Not for use with other library preps/ experimental conditions.
#Author: Matt Lyon, All Wales Medical Genetics Lab
#Mode: BY_COHORT
version="1.8.7"

# Script 2 runs in panel folder, requires final Bams, gVCFs and a PED file
# Variant filtering assumes non-related samples. If familiy structures are known they MUST be provided in the PED file

#load run & pipeline variables
. variables
. /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel".variables

addMetaDataToVCF(){
    output=$(echo "$1" | sed 's/\.vcf/_meta\.vcf/g')
    grep '^##' "$1" > "$output"
    for sample in $(/share/apps/bcftools-distros/bcftools-1.4/bcftools query -l "$1"); do
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

    #index annotated VCF
    /share/apps/igvtools-distros/igvtools_2.3.75/igvtools index "$2"
}

makeCNVBed(){
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
}

### Joint variant calling and filtering ###

#Joint genotyping
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T GenotypeGVCFs \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V GVCFs.list \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_variants.vcf \
-ped "$seqId"_pedigree.ped \
-dt NONE

#Select SNPs
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
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
--setFilteredGtToNocall \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-o "$seqId"_snps_filtered.vcf \
-dt NONE

#Select non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx16g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T SelectVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
-V "$seqId"_variants.vcf \
--selectTypeToExclude SNP \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_non_snps.vcf \
-dt NONE

#Filter non-snps (INDEL, MIXED, MNP, SYMBOLIC, NO_VARIATION)
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
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
--genotypeFilterExpression "DP < 10" \
--genotypeFilterName "LowDP" \
--genotypeFilterExpression "GQ < 20" \
--genotypeFilterName "LowGQ" \
--setFilteredGtToNocall \
-L /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed \
-ip 100 \
-o "$seqId"_non_snps_filtered.vcf \
-dt NONE

#Combine filtered VCF files
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx4g -jar /share/apps/GATK-distros/GATK_3.7.0/GenomeAnalysisTK.jar \
-T CombineVariants \
-R /state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.fasta \
--variant "$seqId"_snps_filtered.vcf \
--variant "$seqId"_non_snps_filtered.vcf \
-o "$seqId"_combined_filtered.vcf \
-genotypeMergeOptions UNSORTED \
-dt NONE

#Add VCF meta data to final VCF
addMetaDataToVCF "$seqId"_combined_filtered.vcf

#bgzip vcf and index with tabix
/share/apps/htslib-distros/htslib-1.4/bgzip -c "$seqId"_combined_filtered_meta.vcf > "$seqId"_combined_filtered_meta.vcf.gz
/share/apps/htslib-distros/htslib-1.4/tabix -p vcf "$seqId"_combined_filtered_meta.vcf.gz

### ROH, SV & CNV analysis ###

#identify runs of homozygosity
for sample in $(/share/apps/bcftools-distros/bcftools-1.4/bcftools query -l "$seqId"_combined_filtered_meta.vcf); do
    /share/apps/bcftools-distros/bcftools-1.4/bcftools roh -O r -s "$sample" -R /data/diagnostics/pipelines/GermlineEnrichment/GermlineEnrichment-"$version"/"$panel"/"$panel"_ROI_b37.bed "$seqId"_combined_filtered_meta.vcf.gz | \
    grep -v '^#' | awk '{print $3"\t"$4-1"\t"$5"\t\t"$8}' > "$sample"/"$seqId"_"$sample"_roh.bed
done

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

#unzip VCF
gzip -dc manta/results/variants/diploidSV.vcf.gz > "$seqId"_sv_filtered.vcf

#Add VCF meta data to SV VCF
addMetaDataToVCF "$seqId"_sv_filtered.vcf

#make CNV target BED file
makeCNVBed

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

#annotate CNV calls with number of het calls
for vcf in $(ls *_cnv.vcf); do

    prefix=$(echo "$vcf" | sed 's/\.vcf//g')
    sampleId=$(/share/apps/bcftools-distros/bcftools-1.4/bcftools query -l "$vcf")

    #add VCF headers
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar UpdateVcfSequenceDictionary \
    I="$vcf" \
    O="$prefix"_header.vcf \
    SD=/state/partition1/db/human/gatk/2.8/b37/human_g1k_v37.dict

    #add metadata, annotate & index
    addMetaDataToVCF "$prefix"_header.vcf
    annotateVCF "$prefix"_header_meta.vcf "$prefix"_meta_annotated.vcf

    #write SNV & Indel dataset to table
    /share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -jar /data/diagnostics/apps/VCFParse/VCFParse-1.2.5/VCFParse.jar \
    -V "$prefix"_meta_annotated.vcf \
    -O "$seqId"_cnv \
    -K

    #move files to sampleId folder
    mv "$prefix"_meta_annotated.vcf* "$sampleId"

    #delete unused files
    rm "$vcf" "$prefix"_header.vcf "$prefix"_header_meta.vcf    
done

#annotate with VEP
annotateVCF "$seqId"_combined_filtered_meta.vcf "$seqId"_filtered_meta_annotated.vcf
annotateVCF "$seqId"_sv_filtered_meta.vcf "$seqId"_sv_filtered_meta_annotated.vcf

#write SNV & Indel dataset to table
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -jar /data/diagnostics/apps/VCFParse/VCFParse-1.2.5/VCFParse.jar \
-V "$seqId"_filtered_meta_annotated.vcf \
-O "$seqId"_snv-indel \
-K

#write SV dataset to table
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -jar /data/diagnostics/apps/VCFParse/VCFParse-1.2.5/VCFParse.jar \
-V "$seqId"_sv_filtered_meta_annotated.vcf \
-O "$seqId"_sv \
-K

#move variant reports to sampleId folder
for i in $(ls *_VariantReport.txt);do
    mv "$i" $(echo "$i" | sed 's/_VariantReport.txt//g' | sed 's/.*_//g')
done

### QC ###

#relatedness test
/share/apps/vcftools-distros/vcftools-0.1.14/build/bin/vcftools \
--relatedness2 \
--out "$seqId"_relatedness \
--vcf "$seqId"_filtered_meta_annotated.vcf

#Variant Evaluation
/share/apps/jre-distros/jre1.8.0_101/bin/java -Djava.io.tmpdir=/state/partition1/tmpdir -Xmx8g -jar /share/apps/picard-tools-distros/picard-tools-2.8.3/picard.jar CollectVariantCallingMetrics \
INPUT="$seqId"_filtered_meta_annotated.vcf \
OUTPUT="$seqId"_CollectVariantCallingMetrics.txt \
DBSNP=/state/partition1/db/human/gatk/2.8/b37/dbsnp_138.b37.excluding_sites_after_129.vcf \
THREAD_COUNT=4

### Clean up ###

#delete unused files
rm -r manta
rm "$seqId"_variants.vcf "$seqId"_variants.vcf.idx "$panel"_ROI_b37_window_gc_mappability.txt  "$seqId"_combined_filtered_meta.vcf
rm "$seqId"_snps.vcf "$seqId"_snps.vcf.idx "$seqId"_snps_filtered.vcf "$seqId"_snps_filtered.vcf.idx "$seqId"_non_snps.vcf igv.log
rm "$seqId"_non_snps.vcf.idx "$seqId"_non_snps_filtered.vcf "$seqId"_non_snps_filtered.vcf.idx "$seqId"_combined_filtered_meta.vcf.gz
rm ExomeDepth.log GVCFs.list HighCoverageBams.list "$seqId"_sv_filtered.vcf "$seqId"_combined_filtered.vcf.idx "$panel"_ROI_b37_window_gc.bed
rm "$seqId"_sv_filtered_meta.vcf BAMs.list variables "$seqId"_combined_filtered.vcf "$seqId"_combined_filtered_meta.vcf.gz.tbi