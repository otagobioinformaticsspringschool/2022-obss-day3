---
title: "Draft genome assembly"
teaching: 10
exercises: 10
questions:
- "What considerations should be taken into account when selecting a genome assembler?"
objectives:
- "Assemble a small fungal genome."
keypoints:
- "There are a wide variety of genome assembly programs available, each with specific strengths, weaknesses, and requirements. Building familiarity with a range of assemblers by reading program manuals and associated articles can be helpful."
---

{% include links.md %}

## 4.1 Genome assembly options

A wide (and continually expanding) range of genome assembly tools are available. Many assemblers are designed for specific data types (e.g., short- vs long-reads, Nanopore vs PacBio data). It is important to select an assembler that is appropriate to the data type(s) and characteristics of the focal genome. Becoming familiar with genome assembly manuals and associated research can be helpful in understanding the strengths and weaknesses of the underlying algorithms implemented by the various assembly program.

<img src="../fig/genome-asm-options.png">

## 4.2 Assembling a genome

Today we will be using [Flye](https://github.com/fenderglass/Flye), an assembler designed to work with long-read data. Flye has modes for Nanopore and PacBio long-read data, that take into account the characteristics specific to each data type. We will use the `--nano-raw` mode to pass just our Nanopore data (labelled 'ont'). 

Although genome assembly algorithms are complex, using an assembler is typically very straightforward. 

We will be making a series of scripts today to run different processes. Let's set up our directory structure so we can keep things tidy. We will store all our scripts in the `scripts` directory. We also need to make a `results` directory to put the outputs we will generate.

```
cd ~/obss_2022/genome_assembly/
mkdir results/
cd scripts/
```

Using `ls` we can see that there are already some files in this directory. Take a look in your `data` directory to see what input data you have been given. 

Now we will make our first script that will assemble a genome. To make a script, we will use a text editor called `nano`. Use the command `nano flye.sl` to open a new file called 'flye.sh' with nano. Now copy in the contents below: 

```
#!/bin/bash -e

#SBATCH --job-name=FLYE
#SBATCH --account=nesi02659
#SBATCH --output=%x.%j.out
#SBATCH --error=%x.%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=70G
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32 

module purge
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

# Change * to your group's dataset
flye --nano-raw ~/obss_2022/genome_assembly/data/all_trimmed_ont_*.fastq.gz --out-dir ~/obss_2022/genome_assembly/results/flye_raw_* -t ${SLURM_CPUS_PER_TASK}
```

You may have noticed that in our `data/` directory, there is more than one .fastq.gz file. We are going to be working in groups today. The letter you are assigned corresponds to one of these data files, so you need to replace the two wildcards in the script above with the letter your group is assigned.

Now let's save our new script and exit, using Ctrl+X. We have already named it `flye.sh` when we opened nano, so hit `Y`, and nano will save this and close. Check the contents of your `scripts` directory, and you should see this file.

As you can see from the SLURM resources, we expect this job to take around 2 hrs to run, so let's check back later.

## 4.3 Later...

<img width="55%" src="../fig/2hrs-later-dna.png">

Now that we've got more familiar with looking at scripts in Section 3, let's have a refresher on what the `flye.sh` script was doing. You should have an understanding of the SBATCH header now, but what about Flye itself? Let's load the Flye module, and take a look at the help documentation to learn more about the program.

```
module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

flye -h
```

Flye has really helpful help documentation! From this we can see that Flye can use both Nanopore and PacBio data to assemble a genome, but that we have to tell it what type of data we are giving it. Now let's investigate the logs produced from running Flye.  

> ## Exercise
> Look at your `.err` and `.out` logs for Flye. Were there any issues or errors? What can you tell about the process? What can you tell about the genome assembly? 
>> ## Solution
>> Take careful note of any ERROR messages. The logs describe the various processes implemented by Flye, and make some comments about the data along the way. At the end of the log there is a print-out of some brief assembly summary metrics, that will likely vary between groups. 
> {: .solution}
{: .challenge}


After confirming that Flye ran correctly without any errors, we can investigate the outputs.

```
ls -lt ~/obss_2022/genome_assembly/results/flye_raw_*/
```

You should see a directory that looks like this:
```
  `-flye.log
  `-assembly_info.txt
  `-assembly.fasta
  `-assembly_graph.gfa
  `-assembly_graph.gv
  `-params.json
  |-40-polishing/
  |-30-contigger/
  |-20-repeat/
  |-10-consensus/
  |-00-assembly/    
```

It is always important to check the log files for any other errors that may have occurred, or provide more details on resources used, or data characteristics. Flye logs contain very descriptive information about every step in the assembly process. This kind of transparent processing is valuable. You can see that a brief summary of these logs was printed to our `.err` file. If we jump to the end of the Flye log, we can see it prints the output directory contents, along with the brief summary of assembly statistics included in our `.err` file. 

Most importantly, we can see there is a file called `assembly.fasta`. This contains the genome assembly produced by Flye in FASTA format. Let's take a quick look at its contents.

```
head ~/obss_2022/genome_assembly/results/flye_raw_*/assembly.fasta
```

As expected for a FASTA file, we can see a header line denoted by `>` followed by a sequence line representing an assembled contig.

Congratulations, you've assembled a genome! But wait, there's more...
