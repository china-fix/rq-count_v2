# rq-count_v2

***

# Introduction

***
q-count v2 is a Snakemake pipeline designed to quickly process raw FASTQ files generated by TraDis-like amplicon sequencing data. Briefly, for each sample, raw data from the MiSeq sequencer underwent initial preprocessing using FASP v0.24.1 to obtain high-quality clean reads—essential for downstream analysis. These clean reads were then aligned to the corresponding reference genome using the Burrows-Wheeler Aligner (BWA) v0.7.19. Specifically, reads from target backbone bacterium-derived mutants were aligned to the corresponding backbone bacterial reference genome. Next, samtools v1.21 was used to calculate and generate a table file documenting the read-mapping depth values at each base position in the reference genome. These computational steps are orchestrated in this Snakemake pipeline.

## Key Features

***

* **Snakemake Script:** The provided Snakemake script efficiently handles reads pre-processing and mapping, ensuring a smooth and organized workflow.
* **Python Scripts:** The "scripts" folder houses Python scripts that enable users to extract mapping information and calculate reads count.

# Dependencies

***

1. [python 3.11.3](https://www.python.org/) and following librarires are required;

* [Biopython 1.79](https://biopython.org/)
* [pandas 1.4.2](https://pandas.pydata.org/)
* [scikit-learn 1.1.1](https://scikit-learn.org/)

2. [snakemake 7.25.3](https://snakemake.github.io/);
3. [conda 23.5.0](https://docs.conda.io/);


# Quick Start (run it on Linux system)

***

## Installation


```
$ git clone https://github.com/china-fix/rq-count_v2.git
```

* [ ]  export the rq-count_v2 folder to PATH (optional);
* [ ]  install conda and create a new snakemake env;

```
$ conda install -n base -c conda-forge mamba
$ conda activate base
$ conda activate snakemake
$ snakemake --help
```

* [ ]  install python biopython pandas scikit-learn

```
$ conda install python biopython pandas scikit-learn -n snakemake
```

## Usage

### Raw reads pre-processing and mapping

* [ ] In your working folder, establish a structured directory with the following arrangement:

1. The primary directory is named "in."
2. Inside the "in" directory, create a subfolder named "REF."In this "REF" subfolder, you should place the reference genome in FASTA format. Ensure that the reference file is named with a ".REF" extension. For example, you can name it "reference\_genome.fasta.REF."
3. Alongside the "REF" subfolder, place your raw sequencing reads data files. Each sample should consist of one file, typically denoted as "\_.fq.gz" .
    For instance, if you have two samples, "sample1" and "sample2," you would place their corresponding read files as "your\_raw\_reads\_sample1\.fq.gz," and "your\_raw\_reads\_sample2\.fq.gz."
    This organized folder structure will facilitate efficient data management and analysis for your sequencing project.

```.
└── in
    ├── REF
    │   └── your_reference_fasta_file.REF
    ├── your_raw_reads_sample1.fq.gz
    └── your_raw_reads_sample2.fq.gz
```

* [ ] run the command in your working direcotry

```
snakemake -s Snakefile -c 8 --use-conda
```

* [ ] New folder named `report` will be created in your working directory, the `.tab` file recorded the mapping depth 


# Citation

....

# Acknowledgments

