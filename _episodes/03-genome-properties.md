---
title: "Exploring input data and genome characteristics"
teaching: 0
exercises: 0
questions:
- "How do the characteristics of short-read and long-read datasets differ?"
- "What underlying characteristics of a genome may complicate genome assembly?"
objectives:
- "Build an understanding of how input data characteristics and underlying genome properties can impact genome assembly"
keypoints:
- "While long-read data are beneficial for producing high-quality genomes, short-read data can help understand the underlying properties of the genome of the focal species."
- "It is important to consider characteristics such as genome size, ploidy, and heterozygosity before assembling a genome, as this may influence program choice, parameter-setting, and overall assembly quality."
---

{% include links.md %}

## 3.1 Setup

First, let's move into the `genome_assembly` directory for today's workshop. Check out the directory structure, what data has been provided, and make a new directory for results outputs.

```
cd ~/obss_2022/genome_assembly/
ls 
ls data/
mkdir results/
```

What can you tell about the data? {% comment %}  looking for: .fastq, size, Illumina and ONT. You'll note that only one fastq file - which contains both forward and reverse reads. {% endcomment %} 

## 3.2 Input data quality assessment

Our raw data has been pre-processed, so let's first get an overview of the quality and metrics of these processed data. For our Nanopore long-read data, there are specific tools that can assess quality and other metrics. One example is the program `NanoStat`. 

As we saw earlier, we can use SLURM to submit scripts to the queue, which allows us to do multiple processes simultaneously. Let's navigate to the `scripts` directory and make a SLURM script to assess our Nanopore data.

```
cd ~/obss_2022/genome_assembly/scripts/
nano nanostat.sl
```

Copy the script below into the new file. Here you want to assess only the dataset that you are using in your assembly, so you will need to replace the `*` in the input and output filenames to one of A-E. 

```
#!/bin/bash -e

#SBATCH --account=nesi02659
#SBATCH --job-name=nanostat
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=12:00
#SBATCH --mem=1G

module purge
module load NanoStat 

cd ~/obss_2022/genome_assembly/data/

# replace * with one of A-E
NanoStat --fastq all_trimmed_ont_*.fastq -t1 > ~/obss_2022/genome_assembly/results/Nanostat_all_trimmed_ont_*
```

One benefit of using NeSI is that many commonly-used programs are pre-installed via modules. We use `module purge` to clear the module space before loading the modules we need for the job using the `module load` command. 

Save the script and run it using `sbatch nanostat.sl`.

While NanoStat runs, we can use FastQC to assess the Illumina short read data in the same way that we used it yesterday. 

```
module load FastQC

cd results/

fastqc ~/obss_2022/genome_assembly/data/All_trimmed_illumina.fastq
```

Once `FastQC` has finished, check your queue to see whether NanoStat is still running, using `squeue -u <nesi.id>`. If it's finished, use `less` to check your `.err` and `.out` logs for any issues.

Now take a look at the results for `FastQC` and `NanoStat`, and discuss the overall metrics and quality with your neighbour. How do our short-read and long-read data sets differ from one another? 

## 3.3 K-mer counting

Before assembling a genome, it is helpful to assess the characteristics of the input data. These data can tell us useful information about expected genome size (if we don't already have an estimate from flow cytometry or other sources), expected heterozygosity, and ploidy. Differences in these characteristics may make genomes difficult to assemble, and may require us to use different strategies.

Today we are going to use our processed Illumina short-read data to count k-mers that can be used to explore some of these characteristics. K-mers are sequence sub-strings of length *k* (i.e., a DNA sequence of a specified length). Genome assembly algorithms typically use k-mers to connect sequences to one another. K-mer counting algorithms are typically based around the properties of short reads with consistent length and high quality, so these algorithms may perform poorly with long-read data that has more variable lengths and qualities. We'll be using a program called `Jellyfish` to count k-mers in our short-read data to assess characteristics of the genome we are assembling.

```
cd scripts/
```

Let's create a script for `Jellyfish` analysis with `nano jellyfish.sl`.

```
#!/bin/bash -e
#SBATCH --job-name=jellyfish
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=05:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=10 
#SBATCH --account=nesi02659

module purge
module load Jellyfish/2.3.0-gimkl-2020a

cd ~/obss_2022/genome_assembly/results/

# count 21-mers from read dataset
jellyfish count -C -m 21 -s 1G -t 10 -o kmer_21_illumina_reads.jf ~/obss_2022/genome_assembly/data/All_trimmed_illumina.fastq

# generate histogram of kmer counts
jellyfish histo -t 10 kmer_21_illumina_reads.jf  > jf_reads.histo
```

First, let's look at the SLURM resources. What can you tell about this job?

Now let's look at how we are navigating our directory structure. Where is the job being processed? Are any other paths included in this job?

Then let's look at the commands we'll be passing to Jellyfish. There are two parts to this job - in the first stage, Jellyfish counts the k-mers. In the second stage it computes the histogram of these counts. As you can see, there are a number of different parameters used, denoted by `-`.

A key part of bioinformatics is getting familiar with program manuals. These are often but not always hosted on GitHub, and provide information about program installation and usage. If we know that Jellyfish is used to count k-mers, can you use Google to find the manual?

{% comment %} https://github.com/gmarcais/Jellyfish - Usage link - discuss how manuals are essential to understand underlying algorithms, and the various parameters in the processes. {% endcomment %} 

Let's save our `jellyfish.sl` file, and run it.

```
sbatch jellyfish.sl
```

We can check where our job is in the queue using `squeue -u <nesi.id>`. Then when the job finishes, we can check our log files to see if there are any problems or errors. 

```
less jellyfish.*.out
less jellyfish.*.err
```

To visualise the results of k-mer analysis, we'll use [GenomeScope](http://qb.cshl.edu/genomescope/). Our input reads were produced from 150 bp Illumina sequencing.

What can you tell from these results, and how do your results compare to the examples at <https://www.nature.com/articles/s41467-020-14998-3/figures/1>?
