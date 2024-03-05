# Usage information

This is not a full release. Please note that some things may not work as intended yet. 

# Running the pipeline

A basic execution of the pipeline looks as follows:

a) Without a site-specific config file

```
nextflow run marchoeppner/eutaxpro -profile standard,singularity --input samples.csv --reference_base /path/to/references --run_name pipeline-test
```
where `path_to_references` corresponds to the location in which you have [installed](installation.md) the pipeline references. 

In this example, the pipeline will assume it runs on a single computer with the singularity container engine available. Other options to provision software are:

`-profile standard,docker` 

`-profile standard,podman` 

`-profile standard,conda` 

b) with a site-specific config file

```
nextflow run marchoeppner/gmo-check -profile lsh --input samples.csv --genome tomato --run_name pipeline-text
```

In this example, both `--reference_base` and the choice of software provisioning are already set in the local configuration and don't have to provided as command line argument. 

# Options

## `--input samplesheet.csv` [default = null]

This pipeline expects a CSV-formatted sample sheet to properly pull various meta data through the processes. The required format looks as follows:

```
sample_id,platform,R1,R2
S100,ILLUMINA,/home/marc/projects/gaba/data/S100_R1.fastq.gz,/home/marc/projects/gaba/data/S100_R2.fastq.gz
```

If the pipeline sees more than one set of reads for a given sample ID, it will concatenate them automatically at the appropriate time. 

## `--primer_set default` [default = "default"]

The name of the pre-configured primer set to use for read clipping. At the moment, only one set is available which corresponds to the §64 guide lines of the German BVL guide lines under L00.00-184. 

Available options:

- default

## `--run_name Fubar` [default = null]

A mandatory name for this run, to be included with the result files. 

## `--email me@google.com` [ default = null]

An email address to which the MultiQC report is send after pipeline completion. This requires for the executing system to have [sendmail](https://rimuhosting.com/support/settingupemail.jsp?mta=sendmail) configured. 


## `--reference_base` [default = null ]

The location of where the pipeline references are installed on your system. This will typically be pre-set in your site-specific config file and is only needed when you run without one. 

## `--outdir results` [default = results]

The location where the results are stored. Usually this will be `results`in the location from where you run the nextflow process. However, this option also accepts any other path in your file system(s). 

# Specialist options

Only change these if you have a good reason to do so. 

## `--vsearch_min_cov` [ default = 5 ]
The minimum amount of coverage required for an OTU to be created from the read data. 

## `--vsearch_cluster_id` [ default = 5 ]
The percentage similarity for ASUs to be collapsed into OTUs. If you set this to 100, ASUs will not be collapsed at all, which will generate a higher resolution call set at the cost of added noise. 
