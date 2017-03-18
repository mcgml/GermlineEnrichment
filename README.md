<h1>GermlineEnrichment</h1>
<h2>Description<h2>
<p>Diagnostic NGS pipeline for SNPs/Indels/CNVs/SVs/LOH from germline panel/exome data (Illumina paired-end)</p>
<p>Requires variables files. See https://github.com/mcgml/MakeVariableFiles</p>
<p>Launch with qsub 1_GermlineEnrichment.sh in the sample directory. Assumes Torque/PBS is installed</p>
<h2>Caveats</h2>
<ul>
  <li>BQSR requires at least 100M bases post filtering to create an accurate model. Roughly, it shouldnt be used for designs less than 0.5Mb.</li>
  <li>Script 2 requires PED file. By default one is created assuming all unrelated samples. Downstream filtering assumes samples are unrelated unless specified in the PED</li>
</ul>
<h2>Ouputs<h2>
<ul>
  <li>BAM alignement</li>
  <li>VCF files</li>
  <li>QC metrics</li>
  <li>Tabix indexed coverage per base</li>
</ul>
<h2>Relatedness</h2>
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
<h2>Expected variant metrics</h2>
<h4>Variants</h4>
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
<p>If your TiTv Ratio is too low, your callset likely has more false positives.</p>
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
<p>A significant deviation from the expected ratios listed in the table above could indicate a bias resulting from artifactual variants.</p>
<h5>Source: http://gatkforums.broadinstitute.org/gatk/discussion/6308/evaluating-the-quality-of-a-variant-callset</h5>