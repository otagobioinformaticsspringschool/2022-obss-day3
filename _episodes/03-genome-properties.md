---
title: "Exploring genome properties with kmers"
teaching: 0
exercises: 0
questions:
- "Key question (FIXME)"
objectives:
- "First learning objective. (FIXME)"
keypoints:
- "First key point. Brief Answer to questions. (FIXME)"
---
FIXME

{% include links.md %}


### 2. 

Use Jellyfish to count kmers (21-mers) in filtered 150bp Illumina reads and produce a histogram that can be imported into Genomescope http://qb.cshl.edu/genomescope/


NB Genomescope website was refusing to import histograms for me today. I tried with various inputs including their own test, think website is temporarily down. There is an Rscript we can run on the command line, might want to explore that as a backup.


```
#!/bin/bash -e
#SBATCH --job-name=jellyfish
#SBATCH --output=GA_%j.jf.out
#SBATCH --error=GA_%j.jf.err
#SBATCH --time=05:00
#SBATCH --mem=4G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10 
#SBATCH --account=nesi02659

module purge
module load Jellyfish/2.3.0-gimkl-2020a

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs

#count 21-mers from read dataset
jellyfish count -C -m 21 -s 1G -t 10 -o kmer_21_illumina_reads.jf All_trimmed_illumina.fastq

#generate historgram of kmer counts
 jellyfish histo -t 10 kmer_21_illumina_reads.jf  > jf_reads.histo
    
```


```
#Example for subsetting ONT  reads (here by min length 1000) and then generating summary stats.

#!/bin/bash -e

#SBATCH --job-name=nanofilt_stat
#SBATCH --output=AW_%j.out
#SBATCH --error=AW_%j.err
#SBATCH --time=40:00
#SBATCH --mem=10G
#SBATCH --ntasks 1
#SBATCH --profile=task
#SBATCH --account=nesi02659


module purge

module load NanoStat nanofilt 

cd /nesi/project/nesi02659/GenomeAssembly/processed_sequence_data

NanoFilt -l 1000 all_q10_pass.pc.noCS.fastq > ../Test_outputs/all_q10_1kb_pass.pc.noCS.fastq

NanoStat --fastq ../Test_outputs/all_q10_1kb_pass.pc.noCS.fastq -t1 > ../Test_outputs/Nanostat_sup_after_PC_CS_Q10_10kb
nanolyse.sh (END)
```