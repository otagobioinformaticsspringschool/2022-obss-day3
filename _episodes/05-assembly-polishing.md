---
title: "Assembly polishing"
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



4.1 Long read polishing with medaka
4.2 Short read polishing with pilon
Map reads with BWA bwa index ../Test_outputs/medaka_polishing/consensus.fasta

Include check mapping stats with samtools flagstat command

```
#!/bin/bash -e

#SBATCH --job-name=medaka
#SBATCH --output=AW_%j.out
#SBATCH --error=AW_%j.err
#SBATCH --time=2:00:00
#SBATCH --mem=15G
#SBATCH --ntasks=1
#SBATCH --profile=task 
#SBATCH --account=nesi02659
#SBATCH --cpus-per-task=4
#SBATCH --gpus-per-node=1

module purge
module load medaka/1.6.0-Miniconda3-4.12.0

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs/

medaka_consensus -i  all_q10_1kb_pass.pc.noCS.fastq -d flye_raw_1kb_q10/assembly.fasta -o medaka_polishing  -t 4 -m r941_min_sup_g507
```

```
#!/bin/bash -e

#SBATCH --job-name=samtools
#SBATCH --output=AW_%j.out
#SBATCH --error=AW_%j.err
#SBATCH --time=1:00:00
#SBATCH --mem=10G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --profile=task
#SBATCH --account=nesi02659

module purge

module load BWA SAMtools

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs

#bwa mem medaka_polishing/consensus.fasta Na171_trimmed_R1.fastq.gz  Na171_trimmed_R2.fastq.gz > illumina_trimmed_mapped_to_medaka_consensus.sam

samtools sort  -@ 8 -o illumina_trimmed_mapped_to_medaka_consensus.sorted.bam illumina_trimmed_mapped_to_medaka_consensus.sam
samtools index illumina_trimmed_mapped_to_medaka_consensus.sorted.bam
```

```
#!/bin/bash -e
#SBATCH --job-name=FLYEq10
#SBATCH --output=AW_%j.flye.out
#SBATCH --error=AW_%j.flye.err
#SBATCH --time=4:00:00
#SBATCH --mem=48G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --profile=task 
#SBATCH --account=nesi02659

module purge


module load  Pilon/1.24-Java-15.0.2

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs

java -Xmx16G -jar $EBROOTPILON/pilon.jar \
    --genome consensus.fasta \
    --frags illumina_trimmed_mapped_to_medaka_consensus.sorted.bam \
    --changes --vcf --diploid --threads 16 \
    --output consensus_pilon_polished_rd1 \
    --minmq 30 \
    --minqual 30
```