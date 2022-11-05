---
title: "Draft genome assembly"
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

## 4.1 Genome assembly options

There are a huge range of genome assembly tools available. Many assemblers are designed for specific data types (e.g., short- vs long-reads). It is important to select an assembler that is appropriate to the data type(s) and characteristics of the focal genome. Reading assembler manuals and associated articles can be really helpful to understand the strengths and weaknesses of the underlying algorithms implemented by the assembly program.

[FIG TO INCLUDE]

## 4.2 Assembling a genome

Today we will be using [Flye](https://github.com/fenderglass/Flye), an assembler designed to work with long-read data. Flye has modes for Nanopore and PacBio long-read data, that take into account the characteristics specific to each data type. 

Although genome assembly algorithms are complex, using an assembler is typically very straightforward. Let's make the script `flye.sl`, containing the content below.

```
#!/bin/bash -e

#SBATCH --job-name=FLYE
#SBATCH --account=nesi02659
#SBATCH --output=%x.%j.flye.out
#SBATCH --error=%x.%j.flye.err
#SBATCH --time=2:00:00
#SBATCH --mem=70G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32 

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

# Change * to your group's dataset
flye --nano-raw ~/obss_2022/genome_assembly/data/all_trimmed_ont_*.fastq --out-dir ~/obss_2022/genome_assembly/results/flye_raw_* -t 24
```
