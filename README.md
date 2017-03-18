# GermlineEnrichment

Diagnostic NGS pipeline for SNPs/Indels/CNVs/SVs/LOH from germline panel/exome data (Illumina paired-end)

Requires variables files.

Launch with qsub 1_GermlineEnrichment.sh in the sample directory.

Ouputs:

-BAM
-VCF
-QC.txt
-Coverage
-Releatedness
-Variant metrics

Caveats:

-BQSR requires at least 100M bases (after filtering) to create a model. Roughly, it shouldnt be used for designs less than 0.5Mb.
-Script 2 requires PED file. By default one is created assuming all unrelated samples.
-Relatedness metrics:  First-degree relatives are ~0.25, and 2nd-degree ~0.125, and 3rd degree 0.0625. Unrelated patients can reach ~0.04.
-Expected variant metrics:
Sequencing Type   Variants    TiTv
WGS               ~4.4M       2.0-2.1
WES               ~41k        3.0-3.3
If your TiTv Ratio is too low, your callset likely has more false positives.
Source: http://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset


<table>
    <th>
        <td>Indel frequency</td>
        <td>Insertion to deletion ratio</td>
    </th>
    <tr>
        <td>Common</td>
        <td>~1</td>
    </tr>
        <tr>
        <td>Rare</td>
        <td>0.2-0.5</td>
    </tr>
</table>