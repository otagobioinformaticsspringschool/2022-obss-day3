---
title: "Assessing assembly quality"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---

{% include links.md %}

# Assembly quality assessment

{% comment %} I think it's important to have this step here rather than last for a) timing reasons, but also to indicate that we may want to collect metrics at many steps in the process. {% endcomment %}

There are many factors to consider when assessing the quality of a genome assembly. Not only are we looking for an assembly that spans the expected genome size, but also we need to consider a range of other metrics.

In an ideal world, a genome assembly will span the expected genome size, be highly contiguous with one scaffold representing one chromosome and generally be comprised of long continuous scaffolds. The AT:GC ratio should be representative of that expected for the focal taxon. The assembly should contain the full set of expected gene orthologs (those genes common across species within a taxon). An ideal assembly will have consistent coverage across its length when the raw data is mapped gainst it, with no evidence of assembly errors. It will also contain only sequences from the focal organism, with no contaminating sequences from other species (e.g., no fungal reads in a bird genome).

[Nat dream genome space figure]

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

Enter the metrics for your assembly into the shared spreadsheet. Are everyone's results the same? If not, why may this be? How does your assembly compare to others? 

## 02. Gene orthologs

To assess assembly completeness in terms of the set of gene orthologs, we use a program called [BUSCO](https://busco.ezlab.org/). This identifies the presence of genes from a common database of gene common among the focal taxon (e.g., fungi) within the assembly. These genes are those that arose in a common ancestor and been passed on to descendant species.  

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

These results act as a proxy to tell us how 'complete' the assembly is. Can you find the file containing a summary of the BUSCO results? Add the metrics for your assembly to the shared spreadsheet. How do the results for your assembly differ from others? Why might this be?

## 03. Plotting coverage

Another form of evidence of assembly quality can come through plotting the sequence data back against the assembly. Spikes in coverage may indicate that repetitive regions have been collapsed by the assembler. On the other hand, regions of low coverage may indicate assembly errors that need to be broken.


[ADD BWA MAPPING & DEPTH SCRIPT}

So now we have an idea of our assembly characteristics. It can be useful to compare these metrics against those for other closely related taxa in the published literature. Like genome size, there may be some genome characteristics that are inherent to a taxon. For example, birds have highly conserved genomes, so we expect all bird genome assemblies to be around 1-1.5 Gb, and to see a particular skew to the AT:GC ratios. 
