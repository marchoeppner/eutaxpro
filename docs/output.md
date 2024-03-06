# Outputs 

## Reports

<details markdown=1>
<summary>reports</summary>

- `name_of_pipeline_run`.taxonomy_by_sample.tsv: A table with accumulated results - one row per sample using the following format:

```
sample  reads   hits
SampleA 12678   Sus scrofa:75.5,Bos taurus:24.5
```

where hits are a sorted list of identified taxa and their respective percentages of the total read count. 

</details>

## Quality control

<details markdown=1>
<summary>MultiQC</summary>

- MultiQC/`name_of_pipeline_run`_multiqc_report.html: A graphical and interactive report of various QC steps and results

</details>

## Raw outputs

<details markdown=1>
<summary>OTUs</summary>

- `name_of_pipeline_run`.usearch_global.tsv - the Number of reads mapping against each respective OTU, per sample
- `name_of_pipeline_run`.precluster.fasta - the final set of OTUs in FASTA format

</details>

<details markdown=1>
<summary>vsearch</summary>

This folder contains the various intermediate processing outputs and is mostly there for debugging purposes. 

</summary>