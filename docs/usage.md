# Usage information

This is not a full release. Please note that some things may not work as intended yet. 

## Running the pipeline

Please see our [installation guide](installation.md) to learn how to set up this pipeline first. 

A basic execution of the pipeline looks as follows:

a) Without a site-specific config file

```
nextflow run marchoeppner/eutaxpro -profile standard,singularity --input samples.csv --reference_base /path/to/references --run_name pipeline-test
```
where `path_to_references` corresponds to the location in which you have [installed](installation.md) the pipeline references (this can be omitted to trigger an on-the-fly temporary installation, but is not recommended in production). 

In this example, the pipeline will assume it runs on a single computer with the singularity container engine available. Available options to provision software are:

`-profile standard,singularity`

`-profile standard,docker` 

`-profile standard,podman` 

`-profile standard,conda` 

b) with a site-specific config file

```
nextflow run marchoeppner/eutaxpro -profile lsh --input samples.csv --run_name pipeline-test 
```

In this example, both `--reference_base` and the choice of software provisioning are already set in the local configuration and don't have to provided as command line argument. 

## Options

### `--input samplesheet.csv` [default = null]

This pipeline expects a CSV-formatted sample sheet to properly pull various meta data through the processes. The required format looks as follows:

```
sample_id,platform,R1,R2
S100,ILLUMINA,/home/marc/projects/gaba/data/S100_R1.fastq.gz,/home/marc/projects/gaba/data/S100_R2.fastq.gz
```

If the pipeline sees more than one set of reads for a given sample ID, it will concatenate them automatically at the appropriate time. 

Allowed platforms are:

* ILLUMINA (expecting PE Illumina reads)
* NANOPORE (expecting ONT reads in fastq format)
* PACBIO (expecting Pacbio CCS reads in fastq format)
* TORRENT (expecting single-end IonTorrent reads in fastq format)

Note that only Illumina processing is currently enabled - the rest is "coming eventually". 

### `--primer_set par64_illumina` [default = "par64_illumina"]

The name of the pre-configured primer set to use for read clipping. At the moment, only one set is available which corresponds to the §64 German BVL guide lines L00.00-184. More sets will be added over time.

Available options:

- par64_illumina (German BVL L00.00-184)

Alternatively, you can specify your own primers as described in the following.

### `--primers primers.txt` [ default = null ]

If you wish to use a set of primers not already configured for this pipeline, you can provide it with this option. You will also have to specify which mitochondrial gene this primer set is targeting using the `--gene` option described elsewhere. 

This text file will be read by [Ptrimmer](https://pubmed.ncbi.nlm.nih.gov/31077131/) to remove PCR primers from the adapter-clipped reads. Please see the Ptrimmer [documentation](https://github.com/DMU-lilab/pTrimmer) on how to create such a config file or look at the [example](../assets/ptrimmer/par64_illumina.txt) included with this pipeline. 

### `--gene` [default = null]

If you do not use a pre-configured primer set, you will also need to tell the pipeline which mitochondrial gene you are targeting. Available options are:

- srrna
- lrrna
- co1
- cytb

### `--run_name Fubar` [default = null]

A mandatory name for this run, to be included with the result files. 

### `--email me@google.com` [ default = null]

An email address to which the MultiQC report is send after pipeline completion. This requires for the executing system to have [sendmail](https://rimuhosting.com/support/settingupemail.jsp?mta=sendmail) configured. 


### `--reference_base` [default = null ]

The location of where the pipeline references are installed on your system. This will typically be pre-set in your site-specific config file and is only needed when you run without one. 

This option can be ommitted to trigger an on-the-fly temporary installation in the work directory. This is however not recommended as it creates unecessary traffic for the hoster of the references. See our [installation guide](installation.md) to learn how to install the references permanently on your system.

### `--outdir results` [default = results]

The location where the results are stored. Usually this will be `results`in the location from where you run the nextflow process. However, this option also accepts any other path in your file system(s). 

## Specialist options

Only change these if you have a good reason to do so. 

### `--vsearch_min_cov` [ default = 5 ]
The minimum amount of coverage required for an OTU to be created from the read data. 

### `--vsearch_cluster_id` [ default = 5 ]
The percentage similarity for ASUs to be collapsed into OTUs. If you set this to 100, ASUs will not be collapsed at all, which will generate a higher resolution call set at the cost of added noise. 
