---
title: "Draft assembly with FLYE"
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


```
#!/bin/bash -e
#SBATCH --job-name=FLYEq10
#SBATCH --output=AW_%j.flye.out
#SBATCH --error=AW_%j.flye.err
#SBATCH --time=4:00:00
#SBATCH --mem=50G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16 
#SBATCH --profile=task 
#SBATCH --account=nesi02659

module purge

module load Flye/2.9.1-gimkl-2022a-Python-3.10.5

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs
flye --nano-raw  /nesi/project/nesi02659/GenomeAssembly/processed_sequence_data/all_q10_1kb_pass.pc.noCS.fastq --out-dir flye_raw_1kb_q10  -t 16

```