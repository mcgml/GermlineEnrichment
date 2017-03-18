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

<h1>Relatedness</h1>
<table>
    <tr>
        <td>Same sample</td>
        <td>1st degree</td>
        <td>2nd degree</td>
        <td>3rd degree</td>
        <td>Unrelated</td>
    </tr>
    <tr>
        <td>~0.5</td>
        <td>~0.25</td>
        <td>~0.125</td>
        <td>~0.0625</td>
        <td>0-0.04</td>
    </tr>
</table>

<h1>Expected variant metrics</h1>
<h2>Source: http://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset</h2>
<h3>Variants</h3>
<table>
    <tr>
        <td>Type</td>
        <td>Variants</td>
        <td>TiTv</td>
    </tr>
    <tr>
        <td>WGS</td>
        <td>~4.4M</td>
        <td>2.0-2.1</td>
    </tr>
    <tr>
        <td>WES</td>
        <td>~41k</td>
        <td>3.0-3.3</td>
    </tr>
</table>
<h4>If your TiTv Ratio is too low, your callset likely has more false positives.</h4>

<table>
    <tr>
        <td>Indel frequency</td>
        <td>Insertion to deletion ratio</td>
    </tr>
    <tr>
        <td>Common</td>
        <td>~1</td>
    </tr>
    <tr>
        <td>Rare</td>
        <td>0.2-0.5</td>
    </tr>
</table>