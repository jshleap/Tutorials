# Differential Expression Analysis: From basics to pipeline
This tutorial will walk you through the theory and praxis of running a 
differential expression analyses on unix (linux/mac) systems. It will start
by giving you a brief refresher on basic bash commands (for more complete
tutorial on bash basics, please see [here](
https://jshleap.github.io/programming/writing-jBasic_BASH/)), moving towards the
usage of [genpipes](https://genpipes.readthedocs.io/en/genpipes-v3.1.5/), a
popular HPC pipeline for multiple genomic analyses. I will explain the theory
behind each of the steps, as well as explain the options of some software to be
used. 
Here I will be focusing on the use in [Compute Canada](
https://docs.computecanada.ca/wiki/Getting_started) systems, but should be easily
extendable to other kinds of HPC systems (particularly if they use SLURM as a 
scheduler).


Table of Contents
=================

   * [Intro](#intro) 
        * [Unix generalities](#unix-generalities)
        * [Principles of RNA-seq](#principles-of-rna-seq)
        * [RNA-seq standard analysis](#rna-seq-standard-analysis)
    
Quality control (FASTQC) (30 mins)
Working with FASTQC 
Understanding the report
Generating a report for your files (assignment)

Trimming (Trimmomatic / cutadapt) (30 mins)
Introduction to Trimmomatic
Understanding Trimmomatic options
Working with Trimmomatic
Trimming your reads (assignment)

Alignment and junction discovery (STAR) (1 Hour - 2 hours)
Introduction to STAR
Understanding STAR options: Generating indices
Understanding STAR options: Mapping
Working with STAR
Generating your indices (Assignment 1)
Generating your mapping (Assignment 2)

Clean alignment (GATK / Picard) (30 mins - 45 mins)
Introduction to picard (only the relevant parts as this is a very big tool)
Understanding picard’s markduplicates
Understanding picard’s RNA metrics
Cleaning your data and generate metrics (Assignment)

Post-alignment quality control (30 mins)
Introduction to RNASeQC
Understanding RNASeQC options
Generating a report for your files (assignment)

Transcript assembly with Cufflinks (2 - 3 Hours)
Introduction to transcriptome assembly
Understanding Samtools hardclip
Understanding Cufflinks 
Understanding Cuffmerge
Understanding Cuffdiff
Understanding Cuffnorm
Running Cufflink bundle on your data (Assignment)

Differential expression analysis using DESeq2 (1 hour)
Differential expression analysis
Understanding R
Understanding the R package DESeq2
Doing DE analysis on your data (Assignment)

Putting all together and more: GenPIPEs (1 - 1.5 hours)
Understanding genpipes
Using GenPIPES on compute canada
Running GenPIPEs on your data (Assignment)



## Intro
### Unix generalities
### Principles of RNA-seq
### RNA-seq standard analysis