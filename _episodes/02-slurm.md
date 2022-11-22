---
title: "Using the slurm scheduler"
teaching: 10
exercises: 10
questions:
- "How to submit jobs to an HPC scheduler"
objectives:
- "Submit a job to the SLURM scheduler."
keypoints:
- "The SLURM scheduler faciliates the fair distribution of compute resource amongst users."
---


{% include links.md %}



Derived from [https://support.nesi.org.nz/hc/en-gb/articles/360000684396-Submitting-your-first-job](https://support.nesi.org.nz/hc/en-gb/articles/360000684396-Submitting-your-first-job)


## Shared resources and schedulers

The Mahuika HPC cluster provided by NeSI is a shared compute environment which means there needs to be a way of distributing resources in a fair way so that all users get the opportunity to use the cluster for their research needs. A complete list of resources can be found at
[Mahuika HPC Cluster resources](https://support.nesi.org.nz/hc/en-gb/articles/360000163575-Mahuika) but a summarised configuration is in the table below

Node type | n | Cores/node| RAM/node
---|---|---|---
Login | 2 | 36 | 512GB
Compute | 226 | 36 | 128GB

In order for the resources to be fairly managed, a scheduler is used to allocate resources. The login nodes allow for small interactive tasks to be done, but the bulk of the compute resource needs to be requested through the `slurm` scheduler which will take job resource requests, queue them, and distribute them throughout the compute nodes as suitable allocations become available.

On the Mahuika cluster, there are different 'partitions' or 'queues' depending on the resources being requested such as walltime, number of cpus, amount of ram, or gpu usage. This link has a list of the [Mahuika Slurm Partitions](https://support.nesi.org.nz/hc/en-gb/articles/360000204076-Mahuika-Slurm-Partitions).

Up until this point we have been using NeSI in an interactive manner, but we need to start requesting larger allocations of resources in order run our jobs in a timely manner.

### Creating a batch script

Jobs on Mahuika and Māui are submitted in the form of a batch script containing the code you want to run and a header of information needed by our job scheduler Slurm.


Create a new file and open it with `nano myjob.sl`

```
#!/bin/bash -e
#SBATCH --job-name=SerialJob # job name (shows up in the queue)
#SBATCH --time=00:01:00      # Walltime (HH:MM:SS)
#SBATCH --mem=512MB          # Memory in MB

pwd # Prints working directory
```

Copy in the above text and save and exit the text editor with 'ctrl + x'.

Note:`#!/bin/bash` is expected by Slurm

Note: if you are a member of multiple accounts you should add the line `#SBATCH --account=<projectcode>`

### Submitting a job

Jobs are submitted to the scheduler using:

```
sbatch myjob.sl
```

You should receive a similar output

Submitted batch job 1748836

sbatch can take command line arguments similar to those used in the shell script through SBATCH pragmas

You can find more details on its use on the [Slurm Documentation](https://slurm.schedmd.com/sbatch.html)

### Using more CPUs

Often, bioformatics software has the option to use more CPUs or threads to make the processing go faster. In the case of SLURM, CPUs are a resource that need to be requested of if you are trying to use more than one.

Depending on the software, this is usually done using `#SBATCH --cpus-per-task=N`, where N is the number of CPUs you wish to use. You then can reference this number in your actual command using the variable `SLURM_CPUS_PER_TASK`. It's a good idea to refer to this variable in your command so that if you increase or decrease the number of CPU cores from slurm this gets reflected in the command too. Usually bioinformatics programs that can make use of extra cpu cores will have a `--threads` or `-t` (or similar) argument.

~~~
#!/bin/bash -e
#SBATCH --job-name=SerialJob # job name (shows up in the queue)
#SBATCH --account=nesi02659  # account for the job
#SBATCH --time=00:01:00      # Walltime (HH:MM:SS)
#SBATCH --mem=512MB          # Memory in MB
#SBATCH --cpus-per-task=2    # Number of logical CPU cores being requested (2 logical cores = 1 physical core)

echo "Number of CPUs requested: $SLURM_CPUS_PER_TASK"
~~~
{: .bash}


### Job Queue

You can view all jobs that are currently in the job queue by using

```
squeue
```

If you want to see your jobs specifically you can use the `-u` flag and your username

```
squeue -u murray.cadzow
```

You can also filter to just your jobs using

```
squeue --me
```

You can find more details on its use on the Slurm Documentation

You can check all jobs submitted by you in the past day using:

```
sacct
```

Or since a specified date using:

```
sacct -S YYYY-MM-DD
```

Each job will show as multiple lines, one line for the parent job and then additional lines for each job step.

### Cancelling

`scancel <jobid>` will cancel the job described by `<jobid>`. You can obtain the job ID by using sacct or squeue.

Tips
`scancel -u [username]` Kill all jobs submitted by you.

`scancel {[n1]..[n2]}` Kill all jobs with an id between `[n1]` and `[n2]`



### Job Output

When the job completes, or in some cases earlier, two files will be added to the directory in which you were working when you submitted the job:

`slurm-[jobid].out` containing standard output.

`slurm-[jobid].err` containing standard error.


### Job efficiency

As you run your jobs, it's a good idea to keep efficiency in mind. You want to be requesting sufficient resources to complete your job, but not significantly more than you need "just in case" - the scheduler will assign all of the resource to you that you ask for regardless of if your task actually uses it all. The unused resource is unavailable to others to use so is wasted which is not ideal is a large, multi-user system like NeSI.

More is not always better. Requesting more CPU cores won't always make things go faster - many programs don't scale in a linear fashion, instead they have deminishing returns and will eventually hit a point where there is no increased performance. 

We can determine the efficiency of our resource request usage through the `nn_seff` command. This will report back the job based on the resources requested and what was used.

The command is run like so: `nn_seff <job_id>`