---
title: "Assembly QC"
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


5.1 Genome assemby stats with QUAST
5.2 Completeness estimates with BUSCO
5.3 Coverage plots
could use samtools depth, but perhaps mosdepth?

5.4 Contaminant scans [demo ]

```
#!/bin/bash -e

#SBATCH --job-name=busco_example
#SBATCH --output=AW_%j.out
#SBATCH --error=AW_%j.err
#SBATCH --time=12:00:00
#SBATCH --mem=18G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --account=nesi02659

module purge

module load BUSCO/5.1.3-gimkl-2020a

cd /nesi/project/nesi02659/GenomeAssembly/Test_outputs 

busco  -i /nesi/project/nesi02659/GenomeAssembly/Test_outputs/medaka_polishing/consensus.fasta  -c 8 -o medaka_fungi_odb10_busco -m genome -l fungi_odb10
```