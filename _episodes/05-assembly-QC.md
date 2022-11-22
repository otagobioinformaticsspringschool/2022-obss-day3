---
title: "Assessing assembly quality"
teaching: 20
exercises: 30
questions:
- "How do we know whether a genome assembly is of sufficient quality to be biologically accurate?"
objectives:
- "Explain key metrics used for assessing genome assembly quality."
- "Generate metrics for the genome assembly produced."
keypoints:
- "All genome assemblies are imperfect, but some are more imperfect than others. The quality of an assembly must be assessed using multiple lines of evidence."
---

{% include links.md %}

# Assembly quality assessment

{% comment %} I think it's important to have this step here rather than last for a) timing reasons, but also to indicate that we may want to collect metrics at many steps in the process. {% endcomment %}

There are many factors to consider when assessing the quality of a genome assembly. Not only are we looking for an assembly that spans the expected genome size, but also we need to consider a range of other metrics.

In an ideal world, a genome assembly will span the expected genome size, be highly contiguous with one scaffold representing one chromosome and generally be comprised of long continuous scaffolds. The AT:GC ratio should be representative of that expected for the focal taxon. The assembly should contain the full set of expected gene orthologs (those genes common across species within a taxon). An ideal assembly will have consistent coverage across its length when the raw data is mapped gainst it, with no evidence of assembly errors. It will also contain only sequences from the focal organism, with no contaminating sequences from other species (e.g., no fungal sequences in a bird genome).

<img src="../fig/nat-dream-genome-space.png">

Even if an assembly is highly contiguous, or shows high completeness in terms of gene orthologs, there may be other issues that may not be identified using only those metrics. To investigate all of these aspects of an assembly requires a range of processes. 

## 01. Basic metrics

An initial step in the assembly QC process is to generate basic assembly metrics. Today we will be using the `assemblathon_stats.pl` script to extract these metrics. You will need to modify the name of the input file depending on the read set you used for the assembly.

```
assemblathon_stats.pl ~/obss_2022/genome_assembly/results/flye_raw_*/assembly.fasta > ~/obss_2022/genome_assembly/results/flye_raw_assembly_QC.txt
```

Now let's view the results.

```
less ~/obss_2022/genome_assembly/results/flye_raw_assembly_QC.txt
```

> ## Exercise
> Enter the metrics for your assembly into the shared spreadsheet. 
> Are everyone's results the same? If not, why may this be? How does your assembly compare to others? 
>> ## Solution
>> We expect to see different results between groups based on the input data. 
> {: .solution}
{: .challenge}

## 02. Gene orthologs

To assess assembly completeness in terms of the set of gene orthologs, we use a program called [BUSCO](https://busco.ezlab.org/). This identifies the presence of genes from a common database of gene common among the focal taxon (e.g., fungi) within the assembly. These genes are those that arose in a common ancestor and been passed on to descendant species. This means that selecting an appropriate database to make these comparisons is essential. For example, if I was assessing a bird genome, I would be best to select the 'aves' database rather than the 'vertebrate' database. Within an assembly, we expect these genes to be complete and present as a single copy.  

```
#!/bin/bash -e

#SBATCH --job-name=busco
#SBATCH --account=nesi02659
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:20:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16

module purge
module load BUSCO/5.1.3-gimkl-2020a

cd ~/obss_2022/genome_assembly/results/ 

# Don't forget to correct the name of the input assembly directory
busco  -i flye_raw_*/assembly.fasta -c 10 -o flye_raw_assembly_busco -m genome -l fungi_odb10
```

> ## Exercise
> These results act as a proxy to tell us how 'complete' the assembly is. 
> 
> Can you find the file containing a summary of the BUSCO results? Add the metrics for your assembly to the shared spreadsheet. 
> 
> How do the results for your assembly differ from others? Why might this be?
> 
> Are everyone's results the same? If not, why may this be? How does your assembly compare to that of other groups'? 
> 
>> ## Solution
>> We expect to see different results between groups based on the differences in input data. 
> {: .solution}
{: .challenge}


## 03. Plotting coverage

Another form of evidence of assembly quality can come through plotting the sequence data back against the assembly. Spikes in coverage may indicate that repetitive regions have been collapsed by the assembler. On the other hand, regions of low coverage may indicate assembly errors that need to be corrected.

To investigate coverage, we will take our input Illumina data and map this to the genome assembly using `bwa mem`. We will then use `bedtools` to calculate mean coverage across genomic 'windows' - regions with a pre-determined size. 

First, let's make a new output directory for the results of this process.

```
mkdir ~/obss_2022/genome_assembly/results/coverage/ 
cd ~/obss_2022/genome_assembly/results/coverage/
```

We have quite a lot of short-read data compared with the data we used for mapping yesterday, so we will use a SLURM script for this step. 

```
#!/bin/bash -e

#SBATCH --job-name=mapping
#SBATCH --account=nesi02659
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=00:45:00
#SBATCH --mem=3G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8

module purge
module load BWA/0.7.17-gimkl-2017a SAMtools/1.15.1-GCC-11.3.0

DATA=~/obss_2022/genome_assembly/data/
OUTDIR=~/obss_2022/genome_assembly/results/coverage/

cd ~/obss_2022/genome_assembly/results/flye_raw_*/

echo "making bwa index"
bwa index assembly.fasta

echo "mapping ONT reads"
bwa mem -t 8 assembly.fasta ${DATA}All_trimmed_illumina.fastq.gz | samtools view -S -b - > ${OUTDIR}mapped.bam

echo "sorting mapping output"
samtools sort -o ${OUTDIR}mapped.sorted.bam ${OUTDIR}mapped.bam

echo "getting depth"
samtools depth ${OUTDIR}mapped.sorted.bam > ${OUTDIR}coverage.out
```

While our data is mapping, we can prepare the genomic windows using `BEDTools`. First, let's load the modules for the programs we need: BEDTools and SAMtools.

```
module load BEDTools/2.30.0-GCC-11.3.0
module load SAMtools/1.15.1-GCC-11.3.0
```

Now we want to gather information about the length of the contigs in the genome assembly. To do this, we can use `SAMtools` to index the genome. The index file produced contains the information we need to pass to BEDTools. 

```
samtools faidx ~/obss_2022/genome_assembly/results/flye_raw_*/assembly.fasta
```

Now we will generate the genomic windows that we will calculate mean coverage across. We are going to set each window to 10 kb, making 5 kb steps before starting a new window. 

```
bedtools makewindows -g ~/obss_2022/genome_assembly/results/flye_raw_*/assembly.fasta.fai -w 10000 -s 5000 > genomic_10kb_intervals.bed
```

Once our script has finished running, we can then continue the process to get the mean coverage across the genome assembly. First, let's use head to see what information is in the final `coverage.out` file. 

We will use `awk` to extract data from the `coverage.out` file produced by `samtools depth` and make the format readable by the BEDtools. `awk` is a widely used command line language that is really useful for manipulating files and extracting useful information. In the following command, we are producing a new file that will be in .bed format, containing 4 columns: 
1. the contig name
2. the start position of the window within that contig
3. the end position of the window within that contig
4. the mean depth of coverage within that window

We tell `awk` to print the first column of `coverage.out`. It then inserts a tab space before printing the content of the second column. It inserts a tab again, before printing another column containing the sum of columns 1 and 2. Finally it adds another tab and then prints the third column of the coverage.out file. 

```
awk '{print $1"\t"$2"\t"$2+1"\t"$3}' coverage.out > coverage.bed
```

Use `head` to compare the format and content of `coverage.out` and `coverage.bed`. Do you see what we did using `awk`?

We can then use BEDTools to do a function that calculates the mean coverage within each of the predefined windows. This makes it more efficient to plot coverage, as instead of plotting a point for every single base in the assembly, we will have one point (the mean depth) for every 10 kb window.

```
bedtools map -a genomic_10kb_intervals.bed -b coverage.bed -c 4 -o mean -g ~/obss_2022/genome_assembly/results/flye_raw_*/assembly.fasta.fai > 10kb_window_5kb_step_coverage.txt
```

Use `head` to check the format of this new output file, `10kb_window_5kb_step_coverage.txt`.

Along with being a powerful tool for statistical analysis, R also is useful for data visualisation. We will then plot the mean coverage results file in RStudio. To do this, we need to open a new tab in Jupyter and select the RStudio launcher. Open a new Rscript file.

The first thing we need to do is to change our working directory to the directory where we have written all the coverage results. 

```
setwd("~/obss_2022/genome_assembly/results/coverage/")
```

Next we want to pass the `10kb_window_5kb_step_coverage.txt` file we just made in to R, ready to do some processing.

```
coverage <- read.delim("10kb_window_5kb_step_coverage.txt", header = FALSE)
```

We can use the R command `head()` in the same way that we use the bash command `head` to check the contents of the data file. 

```
head(coverage)
```

Because we used `header = FALSE` when we read in the data file, you can see that R has automatically named our columns V1-V4. Let's rename these to something more meaningful. To do this, we will use a package called `reshape`.

```
# load the reshape package
library(reshape)

# provide new column names
coverage<-rename(coverage,c(V1="contig", V2="startPos", V3="endPos",V4="meanDepth"))

# check the new headers
head(coverage)
```

To start with, we just want to look at the coverage across one contig. Each contig will have the same starting position (at 1), but the contigs are all different lengths. If we tried to plot the coverage for all the contigs in one plot, it would be messy and difficult to interpret. 

<figure>
<img src="../fig/example-coverage-all-contigs.PNG" style="width:auto;height:auto;max-width:65%;" >
<figcaption align = "center"><b>A coverage plot of all of the contigs in the genome assembly producing using Nanopore dataset D. Position on the contig is on the x-axis, and mean coverage depth on the y-axis. Each contig is assigned a different colour. As the contigs start at position 1 but are all different lengths, it is very difficult to interpret any patterns or identify any anomalies.</b></figcaption>
</figure>

Instead, we will do better to plot mean depth for each contig independently. Here we will just look at the first contig as an example, but you could choose the name of any contig in the genome assembly to investigate. To do this, we will use a package called `dplyr` that can help us filter our dataset.

First let's find out what the names of all the contigs are in the genome assembly. NOTE: Assemblies produced using different datasets may have different numbers of contigs, and the contigs may have different names. 

```
unique(coverage$contig)  
```

Now we know the contig name, we can extract the data for just that contig.

```
# load the dplyr package
library(dplyr)

# From the coverage dataset, we will make a new object called 'contig_10' that only contains the data for contig_10
contig_10 <- coverage %>% 
  filter(contig == "contig_10") 
```

Now we are ready to create a coverage plot. To make the most basic plot in R we can use the `plot` command.

```
plot(contig_10$startPos, contig_10$meanDepth, type = "l")
```

Can you identify any regions with very high or very low coverage compared to the overall pattern across the contig? What may have caused this?

## 04. Summary

So now we have an idea of our assembly characteristics. It can be useful to compare these metrics against those for other closely related taxa in the published literature. Like genome size, there may be some genome characteristics that are inherent to a taxon. For example, birds have highly conserved genomes, so we expect all bird genome assemblies to be around 1-1.5 Gb, and to see a particular skew to the AT:GC ratios. 

It is very useful to generate these metrics after each step in downstream genome assembly processing, so we can identify the extent of the improvements we are making to the genome. We will expect to see diminishing improvements as we progress.
