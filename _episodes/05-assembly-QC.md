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
~/obss_2022/genome_assembly/results/flye_raw_assembly_QC.txt
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

Another form of evidence of assembly quality can come through plotting the sequence data back against the assembly. Spikes in coverage may indicate that repetitive regions have been collapsed by the assembler. On the other hand, regions of low coverage may indicate assembly errors that need to be broken.


```
PLACEHOLDER MAPPING SCRIPT - TO BE TESTED & SLURMIFIED

# Create a genome file (this has two columns: chromosome name and chromosome length) 
samtools faidx assembly.fasta

awk '{print $1 "\t", $2} assembly.fasta.fai > genome.txt

# Make a bedfile with window size and step size defined  
bedtools makewindows -g genome.txt  -w 10000 -s 5000 > genomic_10kb_intervals.bed

# Run samtools depth to obtain read counts at every position
samtools depth  sorted.alignment.bam > coverage.out

# Convert coverage file to BED format
awk '{print $1"\t"$2"\t"$2+1"\t"$3} coverage.out > coverage.bed

# Calculate mean coverage within each genomic window
# Genome file is optional here, including it gives control over order of scaffolds in output.
bedtools map -a  genomic_10kb_intervals.bed -b coverage.bed-c 4 -o mean -g genome.txt > 10kb_window_5kb_step_coverage.txt

## THEN PLOT - with R.

```

So now we have an idea of our assembly characteristics. It can be useful to compare these metrics against those for other closely related taxa in the published literature. Like genome size, there may be some genome characteristics that are inherent to a taxon. For example, birds have highly conserved genomes, so we expect all bird genome assemblies to be around 1-1.5 Gb, and to see a particular skew to the AT:GC ratios. 

It is very useful to generate these metrics after each step in downstream genome assembly processing, so we can identify the extent of the improvements we are making to the genome. We will expect to see diminishing improvements as we progress.
