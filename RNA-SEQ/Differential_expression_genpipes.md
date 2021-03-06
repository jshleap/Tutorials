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

* [Quality control check with FASTQC](#quality-control-check-with-fastqc)
  * [Before we start: File formats](#before-we-start-file-formats)
  * [Working with FASTQC](#working-with-fastqc)
  * [Understanding the report](#understanding-the-report)
    
* [Trimming and adapter removal with Trimmomatic](#trimming-and-adapter-removal-with-trimmomatic)
  * [Introduction to Trimmomatic](#introduction-to-trimmomatic)
  * [Understanding Trimmomatic options](#understanding-trimmomatic-options)
  * [Working with Trimmomatic](#working-with-trimmomatic)
    
* [Alignment and junction discovery using STAR](#alignment-and-junction-discovery-using-star)
  * [Understanding STAR options: Generating indices](#understanding-star-options-generating-indices)
  * [Understanding STAR options: Mapping](#understanding-star-options-mapping)
  * [Running STAR mapping on Compute Canada](#running-star-genome-indexing-on-compute-canada-cluster)
    
* [Cleaning the alignment with Picard](#cleaning-the-alignment-with-picard)
  * [Introduction to picard (only the relevant parts as this is a very big tool)](#introduction-to-picard-only-the-relevant-parts-as-this-is-a-very-big-tool)
  * [Picard’s MergeSamFiles](#picards-mergesamfiles)
  * [Picard’s SortSAM](#picards-sortsam)
  * [Picard’s CollectRnaSeqMetrics](#picards-collectrnaseqmetrics)
  * [Cleaning your data and generate metrics](#cleaning-your-data-and-generate-metrics)
    
* [Quantifying Ribosomal RNA](#quantifying-ribosomal-rna)
  * [Using BWA](#using-bwa)

* [Running ALL in one pipeline! GENPIPES](#running-all-in-one-pipeline-genpipes)
    * [Setting up GENPIPES in Compute Canada](#setting-up-genpipes-in-compute-canada)
    * [Running GENPIPES in Compute Canada](#running-genpipes-in-compute-canada)
        * [Configuration files](#configuration-files)


<!--    
* [Post-alignment quality control](#post-alignment-quality-control)
  * [Introduction to RNASeQC](#introduction-to-rnaseqc)
  * [Understanding RNASeQC options](#understanding-rnaseqc-options)
* [Raw Counts](#raw-counts)
  * [Counting aligned reads with HTSeq](#counting-aligned-reads-with-htseq)
     * [Usage](#usage)
     * [Output](#output)
     * [Running HTSeq in Compute Canada](#running-htseq-in-compute-canada)

    

* [Transcript assembly with Cufflinks](#transcript-assembly-with-cufflinks)
  * [Introduction to transcriptome assembly](#Introduction-to-transcriptome-assembly)
  * [Understanding Samtools hardclip](#understanding-samtools-hardclip)
  * [Understanding Cufflinks](#understanding-cufflinks) 
  * [Understanding Cuffmerge](#understanding-cuffmerge)
  * [Understanding Cuffdiff](#understanding-cuffdiff)
  * [Understanding Cuffnorm](#understanding-cuffnorm)
  * [Running Cufflinks bundle on your data](#running-cufflinks-bundle-on-your-data)


Differential expression analysis using DESeq2 (1 hour)
Differential expression analysis
Understanding R
Understanding the R package DESeq2
Doing DE analysis on your data (Assignment)

Putting all together and more: GenPIPEs (1 - 1.5 hours)
Understanding genpipes
Using GenPIPES on compute canada
Running GenPIPEs on your data (Assignment)
-->

# Intro
## Unix generalities
This tutorial expect you to be more or less comfortable in the terminal. Let's just
touch up on a few things:
### Connecting to a remote server
On unix-like systems we can use the secure shell command `ssh` to connect to a 
remote server. In the terminal we can:
```commandline
ssh username@remotehost
```
username being your username in the remote host and remote host the IP or name 
of the remote server. For example, let's say that you have an account in the
[Graham](https://docs.computecanada.ca/wiki/Graham) cluster under the username
someuser, you could connect to graham through a unix terminal:
```commandline
ssh someuser@graham.computecanada.ca
```
This will connect you to one of the Compute Canada HPC called Graham. For Windows
users without the [linux subsystem](https://docs.microsoft.com/en-us/windows/wsl/about),
I recommend the applicatio [MobaXterm](https://mobaxterm.mobatek.net/) that allows
you to connect to remote servers and move files between your local computer and
a cluster for example. Similar to the example above, you need toc create an ssh
session:

![MobaxtermX11](https://raw.githubusercontent.com/jshleap/Tutorials/main/RNA-SEQ/images/800px-MobaXterm_X11.png)

### Moving files from and to a remote server
There are multiples ways to copy files from and to a remote server. One useful 
one is the rsync command, which allows you to also move only files and folder that
are not up to date in the destination. In unix like systems the general command
is:

```commandline
rsync source destination
```

To better understand this, let's assume that in our local computer (i.e. your 
laptop) you have a file called `afile.txt` in the path `/home/someuser/test`. 
Let's say that you want to move it to Graham supercomputer. You can so this by:

```commandline
rsync /home/someuser/test/afile.txt someuser@graham.computecanada.ca
```

This will put afile.txt in Graham's `/home/someuser`. Now let's say that you want
to put it in a different location, say that you have a folder within you home (
`/home/someuser`) called `testing`. You can do:

```commandline
rsync /home/someuser/test/afile.txt someuser@graham.computecanada.ca:/home/someuser/testing
```

In MobaXterm, you have a sidebar for sftp where you can drag and drop your files:

![sftp](https://mobaxterm.mobatek.net/img/moba/features/feature-sftp-browser.png)

### Creating folders, files, and moving around a unix server
Just a brief mention of how you can create a folder and move around a unix system.
Say you just connected to the cluster (i.e. Graham) and want to create a folder
in your home call test. You can do this by:

```commandline
mkdir test
```

Now you have created test in `/home/someuser` (which is the home of `someuser`), 
now to move into that folder we can use the change directory command `cd`:

```commandline
cd /home/someuser/test
```

You can check your current path with the command `pwd`. Paths can be absolute (
it gives you the full path from `/` to where you are), or relative to where you
are. For example, to get back to your home `/home/someuser`, you could:

```commandline
cd /home/someuser
```

as an absolute path or:

```commandline
cd ..
```

relative to the working directory that you were in `/home/someuser/test`.

Your home also has a special character to refer to without using an absolute path
with is the tilde (`~`). No matter where you are, you can always refer to you home
bu `~`. For example to move to the test folder you created earlier you can do:

```commandline
cd ~/test
```

Now, let's say that instead of a folder, you want to create a file within the test
folder (where we are). You can create the file in your own computer and move it as
explained in [Moving files from and to a remote server](#moving-files-from-and-to-a-remote-server)
or you can use the editors available in the cluster (i.e. nano, vin, emacs). To
use nano, you simply type `nano` in the terminal, and a blank space will show up:

![Nano](https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/Nano.png)

### Downloading files from the web
There are many unix commands to download files from the web. Two of the most
popular are [curl](https://github.com/curl/curl) and [wget](
https://www.gnu.org/software/wget/). I'll focus on the latter. Wget allows you
to get content from the web. **If the content is not a file, it will download the
full html page** so beware of that. For this tutorial let's download a file that 
we will use later on, located at this url: 
https://www.sharcnet.ca/~jshleap/tips_n_tricks/right.norm.fq. In your computer, 
you can just click on the link and will start downloading, but if you want it to 
download it directly on the terminal (or in a server's terminal), you can use 
`wget` as:

```commandline
wget https://www.sharcnet.ca/~jshleap/tips_n_tricks/right.norm.fq
```

That will start the download of this fastqfile. Wget has a lot of options that
can be seen using the manual:

```commandline
man wget
```
but their explanation is out of the scope of this tutorial.


#### Filesystem structure in Compute Canada
***If you are not using Compute Canada, you might skip this intermission.***

Most of Compute Canada resources (i.e. Graham, Beluga and Cedar) have 3 filesystem
spaces:
1. Home: Is a limited space where configuration files (like .bashrc) are stored
2. Project: This is a Lustre filesystem with bigger space (depends on the type 
   of account) and is permanent (not erased)
3. Scratch: Also a Lustre Filesystem with 20Tb of volatile space (gets purge 
   every 60 days). Is intended for computation and storage of intermediate files
   during a run
   
## Submitting jobs with SLURM scheduler
Compute Canada systems (along with many other HPC systems) have a [SLURM](
https://slurm.schedmd.com/documentation.html) scheduler. This is required since
HPC systems are clusters of interconnected computers and with many users a queue
is required to keep the system running smoothly. The heavy computation happens
in the compute nodes, and the instructions on how to run the job have to be sent
through SLURM commands such as `sbatch`. This command takes a batch script as 
input, which contains SBTACH directives along with your program's execution.
Let's say that you normally run a program called `myprogram`, and it takes a file
called `afile.txt` as input. Perhaps you can run this program in your own computer 
with 4 cpus like:

```commandline
myprogram -t 4 afile.txt
```

`-t` being the option for the number of cpus to use. Let's say that you'd like to
use this program in a HPC system using 32 cpus. Let's say that you know that with
the input `afile.txt` and those 32 cpus, you will be using 100Mb of memory, and
it will take roughly 2 hours. Then, in a nano window, you should write something
like this:

```commandline
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=0                 #<-- Amount of memory. In this case reserve all memory
myprogram -t 32 afile.txt       #<-- Your command goes after the directives
```
After that, you can exit nano by pressing the key `ctrl` and `X`. Nano will ask
if you'd like to save the modified file, and will ask you to name the file. Let's
assume we named it `submit.sh`. Now we have a submission script ready, and we 
can proceed for the submission using sbatch:

```commandline
sbatch submit.sh
```

You might be thinking, why `--mem=0`? This is because **in Compute Canada** most
nodes have around 125Gb of RAM, and is more efficient to reserve a full node than 
part of one (easier to schedule). So by asking for 32 CPUS (most Graham's nodes
have 32 cpus) and `mem=0` and  `--cpus-per-task=32` we reserved a full node in
this cluster. You will need to know the node characteristics of your systems to
tune this.

Another thing to point our is that the program `myprogram` needs to be on path
for this execution to work. Most HPC systems have a very vast software stack, and
you will have to get familiar with it. Compute Canada uses the [LMOD](
https://lmod.readthedocs.io/en/latest/) system to load programs as modules. It 
has a convenient search command called `module spider`. So let's say that instead
of `myprogram` you would like to run [fastqc](
https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), then your `submit.sh`
script should be changed to include the loading of fastqc. In Compute Canada systems:

```commandline
module spider fastqc
```
will render:
```angular2html
---------------------------------------------------------------------------------------------------------------------------
  fastqc:
---------------------------------------------------------------------------------------------------------------------------
    Description:
      FastQC is a quality control application for high throughput sequence data. It reads in sequence data in a variety of
      formats and can either provide an interactive application to review the results of several different QC checks, or
      create an HTML based report which can be integrated into a pipeline.

     Versions:
        fastqc/0.11.5
        fastqc/0.11.8
        fastqc/0.11.9

---------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "fastqc" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider fastqc/0.11.9
---------------------------------------------------------------------------------------------------------------------------
```

Which means that there are 3 different versions of the program. To see an specific 
version you can:

```commandline
module spider fastqc/0.11.9
```

rendering:
```angular2html
---------------------------------------------------------------------------------------------------------------------------
  fastqc: fastqc/0.11.9
---------------------------------------------------------------------------------------------------------------------------
    Description:
      FastQC is a quality control application for high throughput sequence data. It reads in sequence data in a variety of
      formats and can either provide an interactive application to review the results of several different QC checks, or
      create an HTML based report which can be integrated into a pipeline.

    Properties:
      Bioinformatic libraries/apps / Logiciels de bioinformatique

    You will need to load all module(s) on any one of the lines below before the "fastqc/0.11.9" module is available to load.

      StdEnv/2020
      nixpkgs/16.09
 
    Help:
      
      Description
      ===========
      FastQC is a quality control application for high throughput
      sequence data. It reads in sequence data in a variety of formats and can either
      provide an interactive application to review the results of several different
      QC checks, or create an HTML based report which can be integrated into a
      pipeline.
      
      
      More information
      ================
       - Homepage: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/
      


```

There it says that you need to load ***EITHER*** `StdEnv/2020` ***OR*** 
`nixpkgs/16.09` before you can load `fastqc/0.11.9`. Now, we can try to load it:

```commandline
module load StdEnv/2020 fastqc/0.11.9
```

To see what is loaded you can use the command `module list`, which will render:
```
Currently Loaded Modules:
  1) CCconfig             4) imkl/2020.1.217  (math)   7) libfabric/1.10.1      10) java/13.0.2   (t)
  2) gentoo/2020    (S)   5) intel/2020.1.217 (t)      8) openmpi/4.0.3    (m)  11) fastqc/0.11.9 (bio)
  3) gcccore/.9.3.0 (H)   6) ucx/1.8.0                 9) StdEnv/2020      (S)

  Where:
   S:     Module is Sticky, requires --force to unload or purge
   bio:   Bioinformatic libraries/apps / Logiciels de bioinformatique
   m:     MPI implementations / Implémentations MPI
   math:  Mathematical libraries / Bibliothèques mathématiques
   t:     Tools for development / Outils de développement
   H:                Hidden Module
```

Now that we know how to load it, let's add it to our submission script. Let's 
say that we want to process the fastq file that we downloaded earlier. Then our
new submission script will become:

```
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-0:05:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=1       #<-- Here is where you put the number of cpus you want
#SBATCH --mem=4G                #<-- Amount of memory.

module load StdEnv/2020 fastqc/0.11.9
fastq  ~/test/right.norm.fq     
```
In this case the program (`fastq`) does not require too much memory (hence the 
`mem=4G`) and time (hence the 5 minutes).

### Checking your job status with SLURM
Slurm has a command to tell you how busy the cluster is, and what is the status
of your job, is called `squeue`. If you cast squeue without any options it will
list all users that have submitted a job. However, you can narrow the search to
yourself by passing your username with the `-u` option. For example, if your
username is `someuser`, you can do:

```commandline
squeue -u someuser
```
A more real example would be my own. I have started an interactive job (yes 
there is a way) with my username `jshleap` and my own group `def-jshleap`, then:

```commandline
squeue -u jshleap
```

renders:

```
 JOBID     USER        ACCOUNT      NAME  ST  TIME_LEFT NODES CPUS TRES_PER_N MIN_MEM NODELIST (REASON) 
44129014  jshleap  def-jshleap_cpu   sh   R      59:45     1    1        N/A    256M   gra796   (None) 
```

Which gives you the jobid (the ID of the launched job) the user, the account used
the name you gave to the job (it defaults to the same of your submit script), 
the status (`st`), the requested resources (time, nodes, cpus, other resources, 
memory), the name of the nodes (if running), and the reason why it hasn't started 
(if it hasn't). The status can be Running (R), pending (PD), and clearing (C).

### Getting an interactive shell
Oftentimes you want to test or run a job that is shorter or that is not very 
resource intensive. For those cases is better to ask for a interactive shell. 
An interactive shell is basically an allocation of a compute node that you can 
access interactively. To schedule an interactive shell with SLURM, you can use
the `salloc` command. You will need to pass the account through the option `--accoun`
or `-A` for short. So in the example above, I would do either:

```commandline
salloc -A def-jshleap
```

or 

```commandline
salloc --account=def-jshleap
```

This will schedule an interactive job on the `def-jshleap` account with the
default memory (250Mb for Compute Canada clusters), one cpu in one node, and
one hour. The same sbatch directive we passed to in the `submit.sh` file can be
passed to the `salloc` command to modify the request. For example, if we wanted
8 cpus, 10Gb of memory, and two hours of computation, we can request an 
interactive shell like this:

```commandline
salloc -A def-jshleap --mem=10G --cpus-per-task=8 --time=00-02:00:00
```

### SLURM environmental variables
Unix systems already have a lot environmental variables that can be useful, you
can check a non-comprehensive list [here](https://www.guru99.com/linux-environment-variables.html).
Likewise, inside an slurm job there is a series of environmental variables that
you can use in your script:
* SLURM_CPUS_PER_TASK: Number of cpus requested per task. Only set if the 
  `--cpus-per-task` option is specified. This can be useful to assign threads
  to your programs, so you do not have to change the value in multiple parts of
  the script.
* SLURM_JOB_NAME: Name of the job. Can be useful for input/output names
* SLURM_MEM_PER_NODE: Same as `--mem`. Useful to assign maximum memory. It does
  not work when `mem` has been assigned to 0.
* SLURM_TMPDIR: Path to a local temporary directory. Usually is a SSD, which 
  makes IO much faster. However, is volatile  
  


For a more comprehensive list, check the slurm manual [here](https://slurm.schedmd.com/sbatch.html).

## Principles of RNA-Seq
RNA-Seq or RNA sequencing is a molecular technique focused in obtaining  a 
collection of RNA molecules (often refer to as library) from a set of samples.
Often the downtream bioinformatic process of analysing and comparing said libraries
is bundle within the term.

RNA-Seq is the main technique for transcriptomics analyses and provides information
about genes being expressed (profiling), identifying variants in gene transcription
(i.e alternative splicing), and -the subject of this tutorial- comparing the 
expression between two treatments of conditions.

<p align="center">
  <img src="https://github.com/jshleap/Tutorials/blob/main/images/RNA_seq.png?raw=true"><br>
  <sup>Image by <a href="http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004393, CC BY 2.5, https://commons.wikimedia.org/w/index.php?curid=53055894"> 
Malachi Griffith, Jason R. Walker, Nicholas C. Spies, Benjamin J. Ainscough, 
Obi L. Griffith</a></sup>
</p>
 
Most often the transcriptome sequencing is done on short fragment technologies
such Illumina, which produces hundreds of million of short reads of the cDNA that
was created in the laboratory. This tutorial assumes that you are already familiar
with the process from total RNA to cDNA library construction and we will focus
on what happes from sequencing to dowtream analyses.

### Illumina sequencing of cDNA
To generate the cDNA library, the most common protocol uses tagmentation. In brief,
you fragment you total (or rRNA depleted) RNA and reverse transcribed using primer
that contains a known tagging sequence in the 5' end, and a random hexamer sequence
in the 3' end. Once the reverse transcription have generated a di-tagged (tagged 
with at least two primers), single stranded cDNA, the fragment is purified and 
amplified, adding Illumina adapters and barcodes:

<p align="center">
  <img src="https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fnmeth.f.355/MediaObjects/41592_2012_Article_BFnmethf355_Fig1_HTML.jpg?as=webp"><br>
  <sup> Courtesy of <a href="https://doi.org/10.1038/nmeth.f.355"> https://doi.org/10.1038/nmeth.f.355</a></sup>
</p>

This is important for downstream analyses, since it gives us information about
the architechture of the constructs. We have to remove sources of noise such as:
1. Spurious amplification: Fragments without adapter sequence
2. Non-biological variation: Such as the one created by the barcodes, adaptor 
   and primer tag.
   
Once the library is generated, is loaded into a flowcell to be placed into the
sequencer. The flow cell is design with lines which in turn contain nanowells, 
and each nanowell contains oligonucleotide anchors that are complementary to
the Illumina adapters.

<p align="center">
  <img src="https://www.hackteria.org/wiki/images/a/a5/FlowCell.jpg"><br>
  <sup> taken from <a href="https://www.hackteria.org/wiki/HiSeq2000_-_Next_Level_Hacking"> www.hackteria.org/wiki/HiSeq2000_-_Next_Level_Hacking</a></sup>
</p>

The cDNA library you have created "flows" within the flow cell and the cDNA
fragments attach to the wells by affinity between the anchors and the adaptors
attached to the RNA fragments during the tagmentation process. This creates a 
bridge between both ends (reverse and forward adaptors) that allows the 
polymerase to generate the new fragment in both ways, essentially generating 
the forward (often called R1 in the resulting files) and reverse (often called 
R2 in the resulting files) reads of the target cDNA:

<p align="center">
  <img src="https://upload.wikimedia.org/wikipedia/commons/6/65/Cluster_Generation.png"><br>
  <sup>author: <a href="https://commons.wikimedia.org/w/index.php?title=User:DMLapato&action=edit&redlink=1"> DMLapato </a></sup>
</p>

After the sequencing process, you will receive a pair (or more) files with your
sequences. If you multiplexed your samples (pool samples together and added a 
barcode) you need to demultiplex them, i.e. split each sample from the run. This
is oftentimes done with the sequencing service, but is not always the case. 
In this tutorial I would skip the demultplexing, but in a nutshell demultiplexing
splits your data into your barcoded samples, often using the command `bcl2fastq`:

```commandline
bcl2fastq --run-folder-dir <FOLDER WITH ILLUMINA DATA> -p <THREADS/CPUS> --output-dir <OUTPUT DIRECTORY> --no-lane-splitting
```

replacing the terms encapsuled between `<>` with your data. Since this is not the 
most common scenario, I am not going to go into details. 

### RNA-seq standard analysis
<p align="center">
  <img src="https://rna-seqblog.com/wp-content/uploads/2016/02/typical.jpg"><br>
  <sup>from <a href="https://rna-seqblog.com/review-of-rna-seq-data-analysis-tools/"> rna-seqblog.com/review-of-rna-seq-data-analysis-tools</a></sup>
</p>

The steps in the figure above show the standard analysis for differential 
expression. That can be summarizeed in 4 main steps:



1. Preprocessing: This step is very important as it weeds out noise from true 
   signal. In this step we have to analyse and visualize the quality of the raw
   data and remove low quality reads that might affect our downstream analyses.
   We will touch on how to visualize and identify contaminants in the section
   [Quality control check with FASTQC](#quality-control-check-with-fastqc). Once 
   we already established contamination and the range of acceptable quality 
   thresholds, we need to trimm our reads. During this process we will remove all
   non-biological sequences (i.e. adapters, barcodes, etc) as well as low quality
   reads. We will cover this in detail in the section [Trimming and adapter removal with Trimmomatic](#trimming-and-adapter-removal-with-trimmomatic).
   ![report](https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/report.jpg)
   <sup>Image from http://cgga.org.cn:9091/gliomasdb/images/figure_1.jpg


2. Alignment and Assembly: Since we sheared our transcriptome, we now have little 
   pieces, and now we have the task to reconstruct full size transcriptomes. To 
   do this we need to map our reads to a reference genome (if we have one) or 
   do de novo assembly. The former is much more accurate if a suitable reference
   is available. In this tutorial we will focus on mapping to the human genome.
   Once we know where each reads go relative to the reference genome (mapping),
   we can piece together our transcriptome (assembly).
   ![mapping](https://home.cc.umanitoba.ca/~frist/PLNT7690/lec12/MapVsAssemble.png)
   <sup>Image from http://jura.wi.mit.edu/bio/education/hot_topics/RNAseq/RNA_Seq.pdf

   
3. Analysis: In this step we need to quantify and compare abundances of transcripts
   mapped to individual genes across treatments/conditions. We can estimate which
   genes have been upregulated (more transcripts produced) or downregulated (less
   transcripts produced) on your base or control condition vs your treatment or
   experimental condition.
   ![DE](https://hbctraining.github.io/DGE_workshop/img/de_theory.png)
   <sup>Image credit: Paul Pavlidis, UBC
   

4. PostProcessing (not in the figure): Visualizing your results and generating 
   figures. It is important to be able to explore the results of your pipeline, 
   and is easier to do it visualy. In the postprocessing step, you can generate
   figures and summaries of these results.
   ![figure](https://galaxyproject.github.io/training-material/topics/transcriptomics/images/rna-seq-viz-with-volcanoplot/volcanoplot.png)
   <sup>Image from https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/rna-seq-viz-with-volcanoplot/tutorial.html

# Quality control check with FASTQC
## Before we start: File formats
Before we start with the actual analysis, we need to understand the file formats
that we will be working with. On bioinformatics there are two main sequence
formats (Fastq and Fasta) and two main mapping formats (BAM and SAM).
### Simple sequence formats: Fasta and Fastq
Before the advent of next generation sequencing (NGS) technologies, most sequence
data was stored in a simple sequence file called a FAST**A** file. Fasta files
contain a header, with information about the sequence, and the sequence:

![fasta](https://www.researchgate.net/profile/Morteza_Hosseini17/publication/309134977/figure/fig1/AS:417452136648705@1476539753111/A-sample-of-the-Multi-FASTA-file_W640.jpg)
<br><sup>from DOI: https://doi.org/10.3390/info7040056

Fasta files can be multiline (as above), where the sequence is broken down into
lines, or single line. The multiline fasta files is by far the most common. The
sequence can be either protein, DNA or RNA. The header in a fasta file start 
with a > symbol, and is followed by the sequence ID (often the accession number
in some database). More information can be place there  and the format varies 
from source to source. 

Fasta files do not have quality information as the chain termination method 
sequencing (also known as Sanger sequencing) did not produce that information.
However, this changed with teh NGS technologies, where a probability of correctly
calling each base can be computed. A new file format was needed for such information,
and the FAST**Q** file was born. Since the earlier NGS technologies only produced
short reads, fastq files sequence and qualities are represented in a single line:

![fastq](https://www.researchgate.net/profile/Morteza-Hosseini-6/publication/309134977/figure/fig2/AS:417452136648711@1476539753452/A-sample-of-the-FASTQ-file_W640.jpg)
<br><sup>from DOI: https://doi.org/10.3390/info7040056

The fastq files' header starts with an @ symbol, followed by the sequence identifier
(often represents the techonology used, the lane, and other information). The next
line in a fastq file is the sequence, which is usually DNA only (RNA is 
retrotranscribed into cDNA). The third line, identified by a leading + sign, serves
as secondary information for each read, but is often empty. The last line is a 
series of alphanumeric characters that represent the quality of each base. Now 
you might be asking, wasn't the quality a probability? How can a probability be
a character? This is a good question, which brings us to a concept called encoding.

#### Encoding
Encoding is the "mapping" of values to some other value, often time more concise, 
the way the quality of the bases are written. There are many encodings, but 
the most popular are Sanger, Solexa, Ilumina 1.3+, Illumina 1.5+, and illumina 
1.8+. In summary, is a character that represents the confidence you have in a 
given base call:

![phred](https://github.com/CristescuLab/Tutorials/raw/master/NGS_QC/images/fastq_phread-base.png)
<br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding

As you can see, the point at which the quality starts (quality score of 0) is 
mapped differently in the [ASCII table](https://en.wikipedia.org/wiki/ASCII) 
depending on the encoding. For example, Illumina 1.3+ and 1.5+ are a PHRED-64
encoding which means that 0 aligned with the 64th character in the ASCII table
(@). However, Illumina 1.5+ actually starts at quality 3 (it does not report 
anything below that).

#### Phred Quality Score
But what does a quality score means? It s related to the probability of an error:

|Phred Quality Score |Probability of incorrect base call|Base call accuracy|
|--- |--- |--- |
|10|1 in 10|90%|
|20|1 in 100|99%|
|30|1 in 1000|99.9%|
|40|1 in 10,000|99.99%|
|50|1 in 100,000|99.999%|
|60|1 in 1,000,000|99.9999%|

As a rule of thumb a Phred score above 20 (99% chances to be right) is considered 
acceptable and above 30 (99.9% chances to be right) as good.

## Working with FASTQC
I have shown several examples working with `fastqc`. In this section I will walk
you through some options, starting with getting help: 

```commandline
module load fastqc/0.11.9
fastqc -h
```
This will render:

```angular2html
            FastQC - A high throughput sequence QC analysis tool

SYNOPSIS

	fastqc seqfile1 seqfile2 .. seqfileN

    fastqc [-o output dir] [--(no)extract] [-f fastq|bam|sam] 
           [-c contaminant file] seqfile1 .. seqfileN

DESCRIPTION

    FastQC reads a set of sequence files and produces from each one a quality
    control report consisting of a number of different modules, each one of 
    which will help to identify a different potential type of problem in your
    data.
    
    If no files to process are specified on the command line then the program
    will start as an interactive graphical application.  If files are provided
    on the command line then the program will run with no user interaction
    required.  In this mode it is suitable for inclusion into a standardised
    analysis pipeline.
    
    The options for the program as as follows:
    
    -h --help       Print this help file and exit
    
    -v --version    Print the version of the program and exit
    
    -o --outdir     Create all output files in the specified output directory.
                    Please note that this directory must exist as the program
                    will not create it.  If this option is not set then the 
                    output file for each sequence file is created in the same
                    directory as the sequence file which was processed.
                    
    --casava        Files come from raw casava output. Files in the same sample
                    group (differing only by the group number) will be analysed
                    as a set rather than individually. Sequences with the filter
                    flag set in the header will be excluded from the analysis.
                    Files must have the same names given to them by casava
                    (including being gzipped and ending with .gz) otherwise they
                    won't be grouped together correctly.
                    
    --nano          Files come from nanopore sequences and are in fast5 format. In
                    this mode you can pass in directories to process and the program
                    will take in all fast5 files within those directories and produce
                    a single output file from the sequences found in all files.                    
                    
    --nofilter      If running with --casava then don't remove read flagged by
                    casava as poor quality when performing the QC analysis.
                   
    --extract       If set then the zipped output file will be uncompressed in
                    the same directory after it has been created.  By default
                    this option will be set if fastqc is run in non-interactive
                    mode.
                    
    -j --java       Provides the full path to the java binary you want to use to
                    launch fastqc. If not supplied then java is assumed to be in
                    your path.
                   
    --noextract     Do not uncompress the output file after creating it.  You
                    should set this option if you do not wish to uncompress
                    the output when running in non-interactive mode.
                    
    --nogroup       Disable grouping of bases for reads >50bp. All reports will
                    show data for every base in the read.  WARNING: Using this
                    option will cause fastqc to crash and burn if you use it on
                    really long reads, and your plots may end up a ridiculous size.
                    You have been warned!
                    
    --min_length    Sets an artificial lower limit on the length of the sequence
                    to be shown in the report.  As long as you set this to a value
                    greater or equal to your longest read length then this will be
                    the sequence length used to create your read groups.  This can
                    be useful for making directly comaparable statistics from 
                    datasets with somewhat variable read lengths.
                    
    -f --format     Bypasses the normal sequence file format detection and
                    forces the program to use the specified format.  Valid
                    formats are bam,sam,bam_mapped,sam_mapped and fastq
                    
    -t --threads    Specifies the number of files which can be processed
                    simultaneously.  Each thread will be allocated 250MB of
                    memory so you shouldn't run more threads than your
                    available memory will cope with, and not more than
                    6 threads on a 32 bit machine
                  
    -c              Specifies a non-default file which contains the list of
    --contaminants  contaminants to screen overrepresented sequences against.
                    The file must contain sets of named contaminants in the
                    form name[tab]sequence.  Lines prefixed with a hash will
                    be ignored.

    -a              Specifies a non-default file which contains the list of
    --adapters      adapter sequences which will be explicity searched against
                    the library. The file must contain sets of named adapters
                    in the form name[tab]sequence.  Lines prefixed with a hash
                    will be ignored.
                    
    -l              Specifies a non-default file which contains a set of criteria
    --limits        which will be used to determine the warn/error limits for the
                    various modules.  This file can also be used to selectively 
                    remove some modules from the output all together.  The format
                    needs to mirror the default limits.txt file found in the
                    Configuration folder.
                    
   -k --kmers       Specifies the length of Kmer to look for in the Kmer content
                    module. Specified Kmer length must be between 2 and 10. Default
                    length is 7 if not specified.
                    
   -q --quiet       Supress all progress messages on stdout and only report errors.
   
   -d --dir         Selects a directory to be used for temporary files written when
                    generating report images. Defaults to system temp directory if
                    not specified.
                    
BUGS

    Any bugs in fastqc should be reported either to simon.andrews@babraham.ac.uk
    or in www.bioinformatics.babraham.ac.uk/bugzilla/
```
There are a few options I will not discuss (they are either too advanced or self
evident) like `--help`, `--quiet`, `--limits`, `--min_length`, `--nofilter`, 
`--nano`, and `--casava`. The last two are from a specific technologies. Nano 
refers to [nanopore](https://nanoporetech.com/), and casava to 
[Casava](https://gaow.github.io/genetic-analysis-software/c/casava/) 
software, both cases outside of the scope of this tutorial.

The `-o` or `--outdir` options allow you to modify the output directory. Say that
your current working directory is `/home/someuser`, but you want your outputs to
be in `/home/someuser/fastqc_results`, by passing the full path to `-o` it will
generate the files there:

```commandline
mkdir -p fastqc_results
fastqc -o fastqc_results file.fastq.gz
```
Note that the output directory must exist (i.e. you need to create it ahead of 
time). By default, the individual results will be created in a compressed folder. 
If you wish to have it uncompress, you can use the option `--extract`. On any 
non-interactive run of fastqc (like in Compute Canda) `--extract` is set by 
default, if you wish ***not*** to extract the archive, you can use the option
`--no-extract`.

Often, the java binary (the java program) is not set in your path. Fastqc allows
you give a specific path with the `--java` option (on Compute Canada systems is 
on path if you loaded fastqc correctly). Say for example that java is not in your
path (i.e. if you type `java --version` it gives you an error) but you installed
it in `/opt`, you can pass it to fastqc by:

```commandline
fastqc -j /opt/java file.fastq.gz
```

It is important to know that for efficiency and visualization purposes, FastQC
adjust computations. For example, most of the analyses are made on a subset of 
the data. Also, bases (i.e. quality, etc) are grouped into chunks. This behaviour
can be disable with the `--no-group` option, with the caveat that the plots can 
become very big, and the analyses might fail.

Another useful flag is when analysing reads that have been already mapped to a 
reference or your reads are in bam or sam formats. By providing a valid format 
(bam,sam,bam_mapped,sam_mapped and fastq) to the `--format` option, the report 
will be done on inputs with that format.

Like with many other tools, FastQC have the ability to use many threads. This
is important in HPC systems such as Compute Canada where you can use multiple
threads by requesting them to the scheduler. You can pass multiple threads to 
the `--threads` option. FastQC however, does parallelization differently. It 
provides a single file per each thread that you request, so is lost on single 
file reports.

A very important option that FastQC has, is the ability to provide a set of known
contaminants. Say that you are working with humans, but in your lab, there is 
significant work with mice. You can put known mice sequences to be detected in
the report so you can test for contaminants. The way to do this is passing a tab
delimited file like this:

```text
# Lines prefixed with a hash will be ignored.
# The format is Name [tab] SEQUENCE
miceSeq1    ACTGATGACCAGTAGCTGATGTTGGTAGT
miceSeq2    ACTGATGACCAAAAAAAAATGATGTTGGTAGT
miceSeq3    ACTGATGACTGTTGGTAGT
```

If you called your file `contaminants.txt`, then you can pass it to fastqc throuh
the `--contaminants` option:

```commandline
fastqc -c contaminants.txt
```

By default it will use a database of usual contaminants/overrepresented sequences.

Another important option is `--adapters`. This option allow you to provide a set
of non canonical adapters to be test if they are present. The file follows the
same structure as the contaminants file.

Oftentimes we want to detect [kmers](https://en.wikipedia.org/wiki/K-mer) or
repeated sequences of a particular lenght. FastQC defaults at kmers of size 7, 
but can be modified with the option `kmers`. This option will only take kmers 
from 2 to 10 bases.

Finally, the `--dir` option allows you to write temporary files in a different
location than the default `/tmp`. This is important when you know that the 
temporary directory is too small.

## Understanding the report
As a test, let's run FastQC on the fastq file `right.norm.fq` downloaded earlier.
Let's assume that we want to not extract the contents, use kmers of size 5, use
the SLURM temporary directory to write the temporary files and to not group bases,
and I want the results to be in a folder called `fastqc_results`:

```commandline
mkdir -p fastqc_results
fastqc --kmers 5 --dir ${SLURM_TMPDIR} --noextract --nogroup \
  --outdir fastqc_results right.norm.fq
```
Will generate an `html` file with the report, and a zip file with the same 
report, stand alone images and raw data of the report. Let's dig in a bit more 
into the fastq report and its contents. FastQC contains several modules that test
the quality of your reads. If the sequences pass the (mostly rules of thumb) 
statiscics, you will see a checkmark next to the module, otherwise an X (all
subsequent plots were obtained from https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):

### Basic Statistics
This just gives you some basic information about your reads, like name, encoding,
type of file, number of sequences, poor quality ones, lenght, and GC content:

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/Basic_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/Basic_bad.png">

### Per base sequence quality
Its name is self-explanatory: this module evaluates the quality at each base 
for a sample of the reads reads. FastQC gives you a box plot of the qualities, 
representing the inter-quartile range (25-75%) (yellow box), the extremes 10 
and 90th percentiles are represented by the whiskers, the median value by a red
line, and the mean quality by the blue line.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/quality_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/quality_bad.png">

From the documentation of this module:

>##### Warning
>A warning will be issued if the lower quartile for any base is less than 10, 
> or if the median for any base is less than 25.
>
>##### Failure
>This module will raise a failure if the lower quartile for any base is less 
> than 5 or if the median for any base is less than 20.

### Per tile sequence quality
This a feature that is exclusive to Illumina technologies. Their flow cells 
typically have 8 lanes,with 2 columns and 50 tiles:

![Flow cell pattern](https://github.com/jshleap/CristescuLab_misc/raw/master/Tutorials/NGS_QC/images/illumina_flowcell.png)
</br> <sup>Courtesy of http://zjuwhw.github.io/2016/08/13/Illumina_sequencer.html

When systematic error occur in a tile, it can indicate sequencing error such as
bubbles, smudges, or dirt. When the errors occur very sparsely and not too 
widespread, is often OK to overlook this error. When a full lane has a problem,
oftentimes is a sequencing error and this cannot be fixed with bioinformatics. 
The problem can occur as as well when the flowcell is overloaded.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/tiles_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/tiles_bad.png">

Not the best quality, but there is no systematic bias... we might be able to fix this with some quality trimming.

From FastQC documentation:
> ##### Warning
> This module will issue a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.
> ##### Failure
> This module will issue a warning if any tile shows a mean Phred score more than 5 less than the mean for that base across all tiles.

### Per sequence quality scores

This module allows you to explore if a significant portion of your reads are of
poor quality. Often times warnings occur when your sequence is shorter than 
your read length, and therefore the end of reads (or the end of the flowcell) 
is of poor quality.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/per_sequence_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/per_sequence_bad.png">


From FastQC documentation:
>##### Warning
>A warning is raised if the most frequently observed mean quality is below 27 -
> this equates to a 0.2% error rate.
>##### Failure
>An error is raised if the most frequently observed mean quality is below 20 - 
> this equates to a 1% error rate.


## Per base sequence content
This module shows the proportion of bases in each position. In an unbiased 
library, the proportion of A, T, C, G, should run parallel to each other. If 
there is a bias, this could imply that the primers or adaptors were not remove,
and therefore there would be a strong bias towards a certain composition. It 
could also mean that you have an over-fragmented library, creating 
over-represented k-mers, or a dataset that has been trimmed too aggressively. 
In amplicon sequencing, there tends to be biases in the composition of the 
given amplicon, especially when dealing with mitochondrial DNA.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/content_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/content_bad.png">


From FastQC documentation:
>#### Warning
>This module issues a warning if the difference between A and T, or G and C is greater than 10% in any position.
>#### Failure
>This module will fail if the difference between A and T, or G and C is greater than 20% in any position.

## Per sequence GC content
This module intends to show the proportion of GC content in the reads. The blue
line represents a theoretical distribution (Normal) of your observed data. 
Deviations from this theoretical distribution often implies contamination of 
some kind (adapter/primer dimers, multiple species in the run). FastQC assumes 
that you are analyzing a  single genome, and therefore will issue a warning in 
multispecies libraries.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/GC_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/GC_bad.png">

From FastQC documentation:
>#### Warning
>A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads.
>#### Failure
>This module will indicate a failure if the sum of the deviations from the normal distribution represents more than 30% of the reads.
>#### Common reasons for warnings
>Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.

## Per base N content
Some sequencer technologies would produce an N when it cannot define which of 
the four bases it has confidence on based on the phenogram. Illumina does not 
produce this, and therefore the plot should be flat and the module should 
always pass.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/N_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/N_bad.png">

From FastQC documentation:
>#### Warning
>This module raises a warning if any position shows an N content of >5%.
>#### Failure
>This module will raise an error if any position shows an N content of >20%.

Failure or warning in this module suggest that the sequencing should probably be
repeated since a significant portion of your reads have no information in them.

## Sequence Length Distribution
It self explanatory title describes well this module. It plots the distribution
of sequence length for your reads. Illumina produces the same length throughout
all the lanes, however, other sequencing platforms produce a distribution of 
them. This module can be safely ignore if you know that you are expecting a 
population of lengths in your reads. If you are using illumina and this module 
fails or gives you a warning, you should talk to the provider of the sequencing.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/lenght_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/lenght_bad.png">

## Sequence Duplication Levels
This module allows you to see the level of duplication of your library. Ideally, 
the blue (total sequences), and the red (deduplicated sequences) should match. 
This would mean that you had a diverse library, and that each sequence has been
sequenced in the proper depth. However, this assumes that you are working with 
a genome of single species, and therefore the warnings and failures of this 
module should only be worrysome then, since it will show a bias (i.e. PCR 
artefacts, resequencing parts of genome). In enriched libraries, you would 
expect some level of duplication, especially when this module only takes the 
first 50 bases and the first 100K sequences to run the tests. 

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/dup_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/dup_bad.png">

From FastQC docs:

>#### Warning
>This module will issue a warning if non-unique sequences make up more than 20%
> of the total.
>#### Failure
>This module will issue a error if non-unique sequences make up more than 50% 
> of the total.

## Overrepresented sequences
This cool module shows you sequences that are present in over 0.1% of your total
reads. The coolest thing about it is that it will run a search for common 
contaminants and report them. In a single species, diverse, uncontaminated 
library, you should expect not to have any overrepresented library.

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/overrep_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/overrep_bad.png">


## Adapter Content
Another self-explanatory module. Here, the most commonly used adapters are 
screened for. They are mostly illumina adapters (Universal, Small 3' RNA, Small
5' RNA, Nextera) and SOLiD small RNA adapter. 

Good Sequence            |  Bad Sequence
:-------------------------:|:-------------------------:
<img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/adapter_good.png" >  | <img src="https://github.com/jshleap/Tutorials/raw/main/RNA-SEQ/images/adapter_bad.png">


From the docs:
>Any library where a reasonable proportion of the insert sizes are shorter than
>the read length will trigger this module. This doesn't indicate a problem as 
>such - just that the sequences will need to be adapter trimmed before 
>proceeding with any downstream analysis.

### Generating a report for your files
Let's generate the report with your own files, and try to make sense of them
now that you know what each one means.

# Trimming and adapter removal with Trimmomatic
As we saw in the previous session, we often need to remove non-biological sequences,
such adapters, from our reads before analyses. Also, checking your reports, you 
might have seen that you required to remove sequences or to trim low quality edges.
This can be done with many programs such as [Trimmomatic](https://github.com/timflutre/trimmomatic)
or [CutAdapt](https://cutadapt.readthedocs.io/en/stable/guide.html). I prefer the
latter as is more versitile and allows you more tailoring. However, we will focus
on Trimmomatic as is very easy to use, is widely adopted, and most RNA-Seq 
pipelines use it.

## Introduction to Trimmomatic
This program does adaptive quality trimming, head and tail crop, and adaptor 
removal. You can check the documentation and download the program [here](
http://www.usadellab.org/cms/index.php?page=trimmomatic). One of the advantages
of trimmomatic is that it allows you to work with pair end sequences, retaining 
only matching pairs. Another advantage is that it allows partial and overlapping
matches for the searching of adapters. Before we run the program, let's check at
some of the options.

## Understanding Trimmomatic options
Here I am going to focus in pair-end reads, but it also works in single end:
### Efficiency and format flags
This flags go before the invocation of the output/input files:
 1. `-threads`: this flag modifies the number of cpu threads that trimmomatic 
    should use in the computations. A typical laptop computer have about 2 cores 
    which should amount to 4 available threads.
 2. `[-phred33 | -phred64]`: this flags tells trimmomatic the encoding of the 
    file (see [above](#encoding))
 
### Change encoding option
If you want to read your file in one encoding and output it in a different one,
this options are the ones you need to use:
- TOPHRED33: Convert quality scores to Phred-33
- TOPHRED64: Convert quality scores to Phred-64
This options (and all the following on CAPS) must go after the input/output 
  files.

### Cropping
Trimmomatic has several options that can be use simultaneously or not:
-   LEADING: Cut bases off the start of a read, if below a threshold quality
-   TRAILING: Cut bases off the end of a read, if below a threshold quality
-   CROP: Cut the read to a specified length
-   HEADCROP: Cut the specified number of bases from the start of the read

LEADING and TRAILING are adaptive cropping. That means that they will cut your 
read's head and/or tail if they fail the specified quality. This differs from 
CROP and HEADCROP, which would cut at an specified length or specified number 
of bases respectively. For the latter two, the program will perform the cropping
for all reads.

### Adaptive length filtering
Trimmomatic has the option MINLEN which will drop reads that fall under the 
specified length:
- MINLEN: Drop the read if it is below a specified length

Say that you have sequences from Illumina, all made up of 150 bp. After quality
trimming, you might end up with some being shorther than 70, others greater than
100, etc. You can set MINLEN, such that every sequence that falls under that
value threhold is dropped from your dataset.

### Adaptive quality trimming
The SLIDINGWINDOW option allows you to trimm reads based on their average 
quality in a window:

- SLIDINGWINDOW: Perform a sliding window trimming, cutting once the average 
  quality within the window falls below a threshold.

It takes two values like `SLIDINGWINDOW:4:15` which means  "Scan the read with 
a 4-base wide sliding window, cutting when the average quality per base drops 
below 15"

### Adapter trimming
Finally, trimmomatic will take a file with the sequences of your adapters and 
will trimm them out. It follows the following call: 
`ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
threshold>:<simple clip threshold>`. 

From Trimmomatic's documentation:
> - fastaWithAdaptersEtc: specifies the path to a fasta file containing all the adapters,
PCR sequences etc. The naming of the various sequences within this file determines
how they are used. See the section below or use one of the provided adapter files
> - seedMismatches: specifies the maximum mismatch count which will still allow a full
match to be performed.
>- palindromeClipThreshold: specifies how accurate the match between the two 'adapter
ligated' reads must be for PE palindrome read alignment.
>- simpleClipThreshold: specifies how accurate the match between any adapter etc.
sequence must be against a read.

<img src="https://raw.githubusercontent.com/CristescuLab/Tutorials/master/NGS_QC/images/trimmomatic_adapter.png" width="650" height="500">

As can be seen in the figure, there are 4 possible scenarios that trimmomatic 
cover:

<ol type="A">
<li>Technical sequence is completely covered by the read and therefore a simple
alignment will identify it.</li>
<li>Only a partial match between the technical sequence and the read, and 
therefore a short alignment is needed.
</li>
<li>Both pairs are tested at once, hence allowing for "is thus much more 
reliable than the short alignment in B, and allows adapter read-though to be 
detected even when only one base of the adapter has been sequenced."</li>
<li>Similar to C</li>
</ol>


The `palindrome clip
threshold` essentially tells how accurate the alignment of the adapters must be. This is the log10 probability against getting a match by random chance, and therefore values around 30 are recommended 

## Working with Trimmomatic
To run trimmomatic in any terminal, you have to:

```commandline
java -jar <path to trimmomatic.jar> PE [-threads <threads>] [-phred33 | -phred64] [-trimlog <logFile>] <input 1> <input 2> <paired output 1> <unpaired output 1> <paired output 2> <unpaired output 2> <OPTIONS>
```
As you can see you need java to be able to run this program.

### Trimmomatic in Compute Canada Clusters
If you are trying to run trimmomatic in Compute Canada systems, we have to first 
load an apropriate version through our modules. So first let's check the availability:

```commandline
module spider trimmomatic
```

as in February 2021, this will render:
```angular2html
---------------------------------------------------------------------------------------------------------------------------
  trimmomatic:
---------------------------------------------------------------------------------------------------------------------------
    Description:
      Trimmomatic performs a variety of useful trimming tasks for illumina paired-end and single ended data.The selection
      of trimming steps and their associated parameters are supplied on the command line. 

     Versions:
        trimmomatic/0.36
        trimmomatic/0.39

---------------------------------------------------------------------------------------------------------------------------
  For detailed information about a specific "trimmomatic" package (including how to load the modules) use the module's full name.
  Note that names that have a trailing (E) are extensions provided by other modules.
  For example:

     $ module spider trimmomatic/0.39
---------------------------------------------------------------------------------------------------------------------------
```

When we query the specific version, it tells us that the only requirement is 
`StdEnv/2020` (***NOTE: This is ONLY for Compute Canada users and may vary***)

When we load `trimmomatic/0.39` we obtain:

```commandline
module load trimmomatic/0.39
To execute Trimmomatic run: java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar
```

This command tells us how to execute our Trimmomatic run. So say that you have
a pair of fastq files named `sample1.R1.fastq.gz` and `sample1.R2.fastq.gz`, and 
you want your paired outputs sent to `sample1.trimmed.R1.fastq.gz` and 
`sample1.trimmed.R2.fastq.gz` and you unpaired outputs to 
`sample1.trimmed_unpaired.R1.fastq.gz` and `sample1.trimmed_unpaired.R2.fastq.gz`. 
You also want to run in on a SLURM cluster (like
the ones in Compute Canada) using 32 cpus. You also want to remove the last 10bp
of your 100bp reads and any trailing and leading base pairs that fall below a 
quality score of 25. Additionally, you want to scan your reads on a sliding window
of size 10 and cut is the average quality drops below 20. You can create your 
submission script like this:

```commandline
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-1:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=125G              #<-- Amount of memory.

module load StdEnv/2020 trimmomatic/0.39
java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE -threads ${SLURM_CPUS_PER_TASK} \
   sample1.R1.fastq.gz sample1.R2.fastq.gz sample1.trimmed.R1.fastq.gz sample1.trimmed_unpaired.R1.fastq.gz \
   sample1.trimmed.R2.fastq.gz sample1.trimmed_unpaired.R2.fastq.gz CROP:10 LEADING:25 \
   TRAILING:25 SLIDINGWINDOW:10:20
```

And voila!


# Alignment and junction discovery using STAR
Now that we have our reads clean, we need to reconstruct our transcripts. Since
the RNA-seq strategy with illumina involves shearing our DNA into smaller
sequentiable bits, now we need to stich them back toguether. There are two options
for doing this: De-novo and reference alignment. In here we will focus in the latter
since is more accurate. The first thing we want to do is to align our reads to 
a reference genome. There are many software available, but here we will be using 
[STAR](https://github.com/alexdobin/STAR) (Spliced Transcripts Alignment to a 
Reference). This software is specifically designed to work with RNA-seq data, 
taking into account spliced variants.

### STAR Strategy
STAR is faster and more accurate than other RNA-Seq alignment software, but it 
often requires a significant ammount of memory. 

STAR has a two-step process: Seeding and reconstructing.

#### Seeding
1. longest exact match or Maximal Mappable Prefixes (MMPs):
   
    ![alt](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step1.png)
    <br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding
   
2. Next MMPs or seed 2 from the unmapped part of the read:
   
    ![alt](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step2.png)
    <br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding

3. If sub-step 2 does not find an exact matching seed 1 gets extended:

    ![alt](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step3.png)
    <br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding
   
4. If extension does not align well it will be softclipped:

    ![alt](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step4.png)
    <br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding

#### Reconstructing
1. Seeds are clustered by proximity using anchor seeds (seeds that do not map
   to multiple sites)
2. The clustered seeds are "stiched" toguether based on the alignment and scoring

    ![alt](https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/img/alignment_STAR_step5.png)
    <br><sup>from https://en.wikipedia.org/wiki/FASTQ_format#Encoding

## Understanding STAR options: Generating indices
Indices are generated to make the accessing of particular regions of the reference
genome faster. Think about it as the index of a book: if you are looking for
something in particular within the book, it is faster to have a page with themes
and the page to where that theme can be found, that skimming the whole book until
you find it.
In the same way many mapping software first require to index the genome for 
efficiency purposes. 

In STAR, the basic options to generate genome indices are:

```
--runThreadN: number of threads
--limitGenomeGenerateRAM: maximum available RAM (bytes) for genome generation
--runMode: genomeGenerate mode
--genomeDir: /path/to/store/genome_indices
--genomeFastaFiles: /path/to/FASTA_file
--sjdbGTFfile: /path/to/GTF_file
--sjdbOverhang: readlength -1
```
##### runThreadN
As we have seen before, the number of threads or cpus is the number of processing
units you can use with your program. This means that your program can run in
parallel, and you are encouraged to use them. Just as a caution, do not request
more threads that you have (or requested if you are in a cluster), since this
will hinder the efficiency as the threads will be competing with each other.

##### limitGenomeGenerateRAM
Often you would like to limit the amount of memory you are using, especially in
HPC systems. By default (if you do not set anything in this option), the program
will use up to 31Gb of memory, so if you requested (or have) less than that, your
run could be killed due to the lack of memory. If you set this variable make sure
is in bytes.

##### runMode
For indexing this option has to be set to `genomeGenerate`. It will tell STAR
that instead of mapping, the current run is to index the fasta file found in
`genomeFastaFiles` and to place the indices into `genomeDir`. The other option
for this flag in STAR is `alignReads`, the default, which we will see in the 
aligning section.

##### genomeDir
This option tells STAR where to store the genome indices. As with 
[FastQC](#quality-control-check-with-fastqc), the folder must exist prior to the
run (i.e. will not create it for you), so make sure to use `mkdir` or to point
to an existing folder.
>This directory path will have
to be supplied at the mapping step to identify the reference genome.

##### genomeFastaFiles
This option must include the path to where the genome (the one you want to index)
fastas are located. The reference can be in a single fasta or multiple fastas 
(usually one per chromosome).
>The tabs are not allowed in chromosomes’ names, and spaces are not recommended

##### sjdbGTFfile
This optional flag will provide a path to the annoation file in 
[GTF format](https://useast.ensembl.org/info/website/upload/gff.html). This will
help STAR identify splice junctions and improve accuracy.
>While this is optional, and STAR can be run without annotations,
using annotations is highly recommended whenever they are available

##### sjdbOverhang
This option tells STAR the expected length around the annotated spliced junction
to construct the database. Ideally, `ReadLength - 1` (meaning the length of your
reads), but with reads of variable length, it should default max(ReadLength) -1.
> * In most cases, the default value of 100 will work as well as the ideal value
> * For Illumina 2x100b paired-end reads, the ideal value is 100-1=99

### Running STAR genome indexing on Compute Canada cluster
Show me how you would do it! 

<details open>
<summary>Peak at the solution after you have tried!</summary>
<br>
<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=0                 #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name STAR_index 	#<-- Job name
#SBATCH -o %j.out			    #<-- File to which standard out will be written
#SBATCH -e %j.err 		        #<-- File to which standard err will be written

module load StdEnv/2020 star/2.7.5a

mkdir -p hg38_index
STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
--runMode genomeGenerate \
--genomeDir hg38_index \    # This is the folder where the indices will be written
--genomeFastaFiles reference_data_ensembl38/*.gz \ #This is the folder with the genome's fastas
--sjdbGTFfile reference_data_ensembl38/Homo_sapiens.GRCh38.92.gtf \ # annotation file
--sjdbOverhang 99 # Assumes that your reads are 100
</code>
</pre>
</details>


## Understanding STAR options: Mapping
Now that we have indexed our reference genome, we can map our reads to it. This
is the default mode in star (`--runMode alignReads`), so the runMode option need
not to be passed. The basic options are:

```ignorelang
--runThreadN 
--genomeDir /path/to/generated/genome_indices
--readFilesIn /path/to/read1 [/path/to/read2]
--readFilesCommand UncompressionCommand
```

There are many advanced options, but those are out of the scope of this tutorial, 
so let's focus on these 4:

##### runThreadN
This is exactly as during genome indexing: the number of threads or cpus requested.
>Just as a caution, do not request more threads that you have (or requested if you 
> are in a cluster), since this will hinder the efficiency as the threads will 
> be competing with each other.

##### genomeDir
This should be the path were you generated the indices. If you do not have that
folder with the corresponding genome indices, go back to the previous section
and run STAR on `runMode genomeGenerate`

##### readFilesIn
In this option you pass the path to where your reads (the sequences to be mapped
) can be found.
>If using Illumina paired-end reads, the read1 and read2 files have to
be supplied. 

When using pair-end reads, you will have to pass both (forward and reverse) paths
to this option.

STAR allows you to pass multiple samples at the same time to be mapped in a single
job, using comma separated list of paths

>For single-end reads use a comma separated list (no spaces around commas), e.g.
> --readFilesIn sample1.fq,sample2.fq,sample3.fq. For paired-end reads, use comma 
> separated list for read1 /space/ comma separated list for read2, e.g.: 
> --readFilesIn sample1read1.fq,sample2read1.fq,sample3read1.fq sample1read2.fq,sample2read2.fq,sample3read2.fq.

You will need to include the paths (or execute the program where the reads are)

##### readFilesCommand
More often than not, your reads will be compressed and will look like read1.zip
or read.gz, or similar. Since there are many ways to compress files, you will
need to tell STAR how to decompress the reads and send them to the standard 
output (i.e. printing it to screen).

>files (*.gz) use --readFilesCommand zcat OR --readFilesCommand gunzip -c. For bzip2-
compressed files, use --readFilesCommand bunzip2 -c.
 

### Running STAR mapping on Compute Canada
Show me how you would do it! 

<details open>
<summary>Peak at the solution after you have tried!</summary>
<br>
This assumes that the reads are located in /scratch/someuser/sequences, that
they are pair end, and ran on 3 diffrent samples:
<ul>
<li>sample1_R1.gz and sample1_R2.gz</li>
<li>sample2_R1.gz and sample2_R2.gz</li>
<li>sample3_R1.gz and sample3_R2.gz</li>
</ul>
Note that all reads are compressed with gunzip (.gz). You can run this in Compute
Canada clusters like this:
<br>
<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=0                 #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name STAR_index 	#<-- Job name
#SBATCH -o %j.out               #<-- File to which standard out will be written
#SBATCH -e %j.err               #<-- File to which standard err will be written

module load StdEnv/2020 star/2.7.5a


STAR --runThreadN ${SLURM_CPUS_PER_TASK} \
--genomeDir hg38_index \    # This is the folder where the indices will be written
--readFilesIn /scratch/someuser/sample1_R1.gz,/scratch/someuser/sample2_R1.gz,/scratch/someuser/sample3_R1.fq.gz \
              /scratch/someuser/sample1_R2.gz,/scratch/someuser/sample2_R2.gz,/scratch/someuser/sample3_R2.fq.gz \
--readFilesCommand gunzip -c
</code>
</pre>
</details>

## Cleaning the alignment with Picard
### Introduction to picard (only the relevant parts as this is a very big tool)
Picard/GATK helps you manipulate high-throughput sequencing data in formats like
SAM/BAM/CRAM and VCF. 

#### Before we start, the SAM format
To understand what picard will be doing, we need to familiarized ourselves with
the SAM/BAM formats. Both formats are the same, with the exception that the BAM
format is the binary version of the SAM. SAM stands for Sequence Alignment Map
and is the main format of -- you guessed it -- mapping software. Similar to 
other sequence files lile fasta and fastq, SAM consist on header and alignemnt 
fields, the former always coming before the latter. 
<p align="center">
  <img src="https://raw.github.com/ecerami/samtools_primer/master/figs/sam_format_example.png"><br>
  <sup>Image from <a href="https://du-bii.github.io/module-5-Methodes-Outils/seance1/images/SAM_format.jpg"> https://du-bii.github.io/module-5-Methodes-Outils/seance1/images/SAM_format.jpg</a></sup>
</p>

As you can see, is a  it is a TAB-delimited text format. As with fastq, the 
header starts with the symbol `@`, but unlike fastq, both the alignment and the
header can have multiple lines. The header is optional and give extra information.

Each alignment section has 11 mandatory fields:

Column|	Field	|Type	|Brief Description
---   |---      |---    |---              
1	  |QNAME	|String	|Query template NAME
2	  |FLAG	    |Int	|[bitwise FLAG](https://en.wikipedia.org/wiki/SAM_(file_format)#Bitwise_flags)
3	  |RNAME	|String	|References sequence NAME
4	  |POS	    |Int	|1- based leftmost mapping POSition
5	  |MAPQ	    |Int	|[MAPping Quality](https://genome.sph.umich.edu/wiki/Mapping_Quality_Scores)
6	  |CIGAR	|String	|[CIGAR string](https://genome.sph.umich.edu/wiki/SAM#:~:text=The%20CIGAR%20string%20is%20a,are%20not%20in%20the%20reference)
7	  |RNEXT	|String	|Ref. name of the mate/next read
8	  |PNEXT	|Int	|Position of the mate/next read
9	  |TLEN	    |Int	|observed Template LENgth
10	  |SEQ	    |String	|segment SEQuence
11	  |QUAL	    |String	|[ASCII of Phred-scaled base QUALity+33](https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm)

The header is optional, but when present it has important metada tags, and marks
(like when you mark duplicates, etc). The options are too mant co cover them here
(for a more compreshensive explanation of the SAM file format check [here](
http://samtools.github.io/hts-specs/SAMv1.pdf), including headers), but some of
the main ones:
* `@HD`: File-level metadata. If present, there must be only one @HD line and it 
  must be the first line of the file. It gives information about the version of
  the file (VN), how is the file sorted (unknown (default), unsorted, queryname
  and coordinate), how the sequences are grouped (GO), etc.
  
* `@SQ`: Reference sequence dictionary. The order of @SQ lines defines the 
  alignment sorting order. It provides with name of the sequence (SN), refence
  sequence lenght (LN), alternate locus (AH), among many other options.
  
* `@RG`: Read group. Unordered multiple @RG lines are allowed. It gives 
  information about the groupiong of the reads such as read group identifier (ID)
  barcode sequence identifying the sample library (BC), library (LB), etc.
  

##### Command Syntax
Picard is also written in Java, and has to be invoked in a similar manner than
[Trimmomatic](#introduction-to-trimmomatic). In a regular computer, the general
syntax is:

```bash
java jvm-args -jar picard.jar PicardToolName \
	OPTION1=value1 \
	OPTION2=value2
```

where `jvm-args` are arguments passed to the Java engine, `PicardToolName` 
refers to the actual name of the tool that want to be use, and `OPTION1=value1`,
are the option-value pairs of the options for the tool.

In Compute Canada systems you can load picard with the LMOD command `module load`:

```bash
module load picard/2.23.3
To execute picard run: java -jar $EBROOTPICARD/picard.jar
```
As with trimmomatic, it prints the way to invoke the tool, then the general syntax
in Compute Canada clusters is:

```bash
java jvm-args -jar $EBROOTPICARD/picard.jar PicardToolName \
	OPTION1=value1 \
	OPTION2=value2
```
Note that you can modify Java's behavious in Compute Canada by passing options
(i.e. jvm-args)

After STAR, most RNA-seq pipeline would use picard to:
1. Merge resulting SAM files (if multiple)
2. Sort the merged SAM file
3. Mark duplicates
4. Generate RNA metrics

### Picard’s MergeSamFiles
Often we need to merge SAM/BAM files from parallel runs of the mapping software.
Like with most pipelines, there are multiple tools to do this, but Picard's 
MergeSamFiles is a popular choice. To run it, you have to provide multiple 
inputs (`I`) and an output (`O`). You can pass a single thread/cpu to help with
compressing and writting, but in general, this tools is not parallelizeable.
Assuming that you have two SAM files names `file1.sam` and `file2.sam`, and 
that you would like the output to be named `output.sam`, you can:

```bash
module load picard/2.23.3
java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
	I=file1.bam \
	I=file2.bam \
	O=output.bam \
	USE_THREADING=true
```
This will generate a single unsorted SAM file called output.sam. The same 
procedure can be perform to merge bamfiles. You can check more options [here](
http://broadinstitute.github.io/picard/command-line-overview.html#MergeSamFiles
).

For RNAseq, this should be done if you have more than one file PER SAMPLE. Do
not merge files from different samples as it will create problems down the road.

### Picard’s SortSAM
As mentioned in the previous section, once file are merged or mapped, they might
be unsorted. Sorted files are important for efficiency resons, and many programs
require the SAM/BAM files to be sorted. This can also be achieved through Picard's
functionality through the tool called SortSAM. This tool:
>This tool sorts the input SAM or BAM file by coordinate, queryname (QNAME), or
> some other property of the SAM record. The SortOrder of a SAM/BAM file is found 
> in the SAM file header tag @HD in the field labeled SO.
> 
> For a coordinate sorted SAM/BAM file, read alignments are sorted first by the
> reference sequence name (RNAME) field using the reference sequence dictionary
> (@SQ tag). Alignments within these subgroups are secondarily sorted using the
> left-most mapping position of the read (POS). Subsequent to this sorting 
> scheme, alignments are listed arbitrarily. 
> 
> For queryname-sorted alignments, all alignments are grouped using the 
>queryname field but the alignments are not necessarily sorted within these 
> groups. Reads having the same queryname are derived from the same template.

#### Example usage
```bash
module load picard/2.23.3
java -jar $EBROOTPICARD/picard.jar SortSam \
      I=input.bam \
      O=sorted.bam \
      SORT_ORDER=coordinate
```
The possible values of SORT_ORDER are:
1. unsorted: Leave the file as is 
2. queryname: Alphabetical sorting based on the query name 
3. coordinate: Numerical by chromosome and start position
4. duplicate: Sorts the reads so that duplicates reads are adjacent

### Picard’s MarkDuplicates
This tool check and mark (not remove) duplicates. Duplicate reads are defined as:

> originating from a single fragment of DNA. Duplicates can arise during sample
> preparation e.g. library construction using PCR. Duplicate reads can also 
> result from a single amplification cluster, incorrectly detected as multiple 
> clusters by the optical sensor of the sequencing instrument. These duplication
> artifacts are referred to as optical duplicates.
 
We want to know (and later softclip them) where the duplicates are, but since we
are working with RNA-SEQ, most reasercher choose not to remove them. This is 
because it has been shown that retaining some unatural duplicates does not cause 
significant artifacts if the library complexity, while removing them might also
remove natural duplicates in the transription data, distorting our analyses.

#### Usage example
```bash
module load picard/2.23.3
java -jar $EBROOTPICARD/picard.jar MarkDuplicates \
      I=input.bam \
      O=marked_duplicates.bam \
      M=marked_dup_metrics.txt
```
The metrics (M) option write down information about the duplication found.

If desired, you can remove optical duplicates with the option 
`REMOVE_SEQUENCING_DUPLICATES=true` and/or remove all duplicates with 
`REMOVE_DUPLICATES=true`. More options can be found on the tool's [documentation](
https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).

### Picard’s CollectRnaSeqMetrics
This tool produces important RNA alignment metrics for your files. It takes 
either a BAM or SAM of **aligned** reads (as we did with [STAR](#alignment-and-junction-discovery-using-star)).
This will describe the distribution of the bases within the transcripts (
nucleotides on each genomic region), the regions passing quality filters, etc.

Please see the [CollectRnaSeqMetrics definitions](http://broadinstitute.github.io/picard/picard-metric-definitions.html#RnaSeqMetrics) 
for details on how things are calculated.

Besides the aligned BAM/SAM file, this tool requires a REF_FLAT file, a 
tab-delimited file containing information about the location of RNA transcripts, 
exon start and stop sites, etc. For an example refFlat file for GRCh38, 
see refFlat.txt.gz at http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database. 

#### Usage example
```bash
module load picard/2.23.3
java -jar picard.jar CollectRnaSeqMetrics \
      I=input.bam \
      O=output.RNA_Metrics \
      REF_FLAT=ref_flat.txt \
      STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
      RIBOSOMAL_INTERVALS=ribosomal.interval_list
```

Where **I** and **O** are input and output as before, **REF_FLAT** is the path
to the flat reference file, and **RIBOSOMAL_INTERVALS** are the indices where ribosomal
regions are. The latter is important as ribosomal regions might inflate some of
the metrics. An example of that file can be downloaded from [here](https://gist.github.com/slowkow/b11c28796508f03cdf4b/raw/38d337698ff1e6915578dfa08826c73631c3e0b5/hg19.rRNA.interval_list).
On Compute Canada systems we have a set of intervarl at 
```
/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.${VERSION}/annotations/Homo_sapiens.${VERSION}.Ensembl${NSB_VER}.rrna.interval_list

```
Version being the specific reference genome version you are using (e.g GRCh38) and
NSV_VER is the ensemble version you want to use (e.g. 90).

The **STRAND** option refers to STRAND_SPECIFICITY:
>For strand-specific library prep. For unpaired reads, use FIRST_READ_TRANSCRIPTION_STRAND
> if the reads are expected to be on the transcription strand. Required. Possible values:
> {NONE, FIRST_READ_TRANSCRIPTION_STRAND, SECOND_READ_TRANSCRIPTION_STRAND}

### Cleaning your data and generate metrics
Now you have your aligned reads, and you want to run the picard tools to do some 
cleaning. That will look like this in slurm:

<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-6:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=0                 #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name Picard_steps	#<-- Job name
#SBATCH -o %j.out               #<-- File to which standard out will be written
#SBATCH -e %j.err               #<-- File to which standard err will be written

module load StdEnv/2020 picard/2.23.3

# First merge files PERSAMPLE if required
java -Xmx${SLURM_MEM_PER_NODE} $EBROOTPICARD/picard.jar \
    MergeSamFiles MergeSamFiles \
    I=SAMPLE1_1Aligned_out.bam \
    I=SAMPLE1_2Aligned_out.bam \
    I=SAMPLE1_3Aligned_out.bam \
    O=SAMPLE1.bam \
    USE_THREADING=true

# Sort the merged file (or your aligned file if only one)
java -Xmx${SLURM_MEM_PER_NODE} $EBROOTPICARD/picard.jar SortSam \
    I=SAMPLE1.bam \
    O=SAMPLE1_sorted.bam \
    SORT_ORDER=coordinate 

# Mark duplicates in the sorted file
java -Xmx${SLURM_MEM_PER_NODE} $EBROOTPICARD/picard.jar MarkDuplicates \
    I=SAMPLE1_sorted.bam \
    O=SAMPLE1_marked_dup.bam \
    M=marked_dup_metrics.txt 

# Get some metrics
java -Xmx${SLURM_MEM_PER_NODE} $EBROOTPICARD/picard.jar CollectRnaSeqMetrics \
    I==SAMPLE1_marked_dup.bam \
    O=SAMPLE1.RNA_Metrics \
    REF_FLAT=PATH/TO/REF/ref_flat.txt \
    STRAND=SECOND_READ_TRANSCRIPTION_STRAND \
    RIBOSOMAL_INTERVALS=/PATH/TO/INTERVAL/ribosomal.interval_list
</code>
</pre>

Now you will have your reads ready to proceed the pipeline!


# Quantifying Ribosomal RNA
Since ribosomal RNA is the most abundant of the RNAs, we want to quantify if the
depletion in the molecular process was enough. To do that we can use fast mapping
to known ribosomal sequences. In Compute Canada we have some versions of the human
genome and some annotation at:

```bash
$ ls /cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38
annotations  Homo_sapiens.GRCh38.Ensembl90.dbSNP150.ini
downloads    Homo_sapiens.GRCh38.ini
genome       log

$ cd /cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/
$ ls
Homo_sapiens.GRCh38.Ensembl77.rrna.fa
Homo_sapiens.GRCh38.Ensembl77.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl77.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl77.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl77.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl77.rrna.fa.sa
Homo_sapiens.GRCh38.Ensembl83.rrna.fa
Homo_sapiens.GRCh38.Ensembl83.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl83.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl83.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl83.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl83.rrna.fa.sa
Homo_sapiens.GRCh38.Ensembl85.rrna.fa
Homo_sapiens.GRCh38.Ensembl85.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl85.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl85.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl85.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl85.rrna.fa.sa
Homo_sapiens.GRCh38.Ensembl86.rrna.fa
Homo_sapiens.GRCh38.Ensembl86.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl86.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl86.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl86.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl86.rrna.fa.sa
Homo_sapiens.GRCh38.Ensembl87.rrna.fa
Homo_sapiens.GRCh38.Ensembl87.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl87.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl87.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl87.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl87.rrna.fa.sa
Homo_sapiens.GRCh38.Ensembl90.rrna.fa
Homo_sapiens.GRCh38.Ensembl90.rrna.fa.amb
Homo_sapiens.GRCh38.Ensembl90.rrna.fa.ann
Homo_sapiens.GRCh38.Ensembl90.rrna.fa.bwt
Homo_sapiens.GRCh38.Ensembl90.rrna.fa.pac
Homo_sapiens.GRCh38.Ensembl90.rrna.fa.sa
```

So we can use them to investigate the remnants of ribosomal RNA in our sample. To
do so, we can just use BWA to map our reads to the reference ribosomal sequences.

## Using BWA
As STAR, BWA is a genome aligner. As opposed to STAR, BWA is not geared towards
RNAseq and junction determination, but to do efficient mapping of reads to the 
reference. As STAR, BWA also have to modes `index` and `mem`. Fortunately, in the
path above we can see that index has been already called (files .amb, .ann, .bwt
and .pac), so we only need to run the mem option. Although bwa has multiple option
the only relevant for us here is the `-t` or threads, so the general usage is simple:

```bash
bwa mem -t <number of cpus> Path2reference Path2R1 PathtoR2 > outputsam
```
We can either "pipe" or transform the resulting sam into a bam using samtools:

```bash
bwa mem -t <number of cpus> Path2reference Path2R1 PathtoR2| samtools view -@ 32 -bS -o outputbam
```

Then we can use the flagstat tool of samtools to capture the alignemnt statistics.
The relevant number is the fifth column (successfully aligned reads), which should
be low. You can also pipe this last command (to avoid intermediate file):

```bash
bwa mem -t <number of cpus> Path2reference Path2R1 PathtoR2 | \
  samtools view -@ <number of cpus> -bS - | samtools flagstat -@ <number of cpus> - > rrna.stats
```

Now, try to create your submission script!!!

<details open>
<summary>Peak at the solution in Compute Canada after you have tried!</summary>
<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=32      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=0                 #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name rRNA_quant	#<-- Job name
#SBATCH -o %j.out               #<-- File to which standard out will be written
#SBATCH -e %j.err               #<-- File to which standard err will be written

module load StdEnv/2020 bwa/0.7.17 samtools
path2ref=/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/rrna_bwa_index/Homo_sapiens.GRCh38.Ensembl90.rrna.fa
bwa mem -t ${SLURM_CPUS_PER_TASK} ${path2ref} YOUR_R1_Goes_HERE YOUR_R2_Goes_HERE | &#92;
  samtools view -@ ${SLURM_CPUS_PER_TASK} -bS -| &#92;
  samtools flagstat -@ ${SLURM_CPUS_PER_TASK} - > rrna.stats
</code>
</pre>
Check out the output and tell me if is OK!
</details>

<!--
# Post-alignment quality control
We have done significant work so far, however, you do not know how good your 
alignments are, and if general what is the quality of your data. For that we need
a report similar to fastq, but tailored for RNAseq. The tool **RNASeQC**.

## Introduction to RNASeQC
RNASeQC is a handy tool specialized in getting RNAseq quality control metrics.
It makes use of curated annotations of transcripts, however (from the tool's page):
>This tool requires that the provided GTF be collapsed in such a way that there
> are no overlapping transcripts on the same strand and that each gene have a 
> single transcript whose id matches the parent gene id. This is not a 
> transcript-quantification method. Readcounts and coverage are made towards 
> exons and genes only if all aligned segments of a read fully align to exons of
> a gene, but keep in mind that coverage may be counted towards multiple transcripts
> (and its exons) if these criteria are met.

Therefore we need first a GTF annotation file, and make sure the file is flat.
The first part (obtaining the GTF file) is simple in Compute Canada as we have 
those files in the same path explored above:

```bash
/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl90.gtf
```

To make it flat we can use the [GTEx collapse annotation script](https://github.com/broadinstitute/gtex-pipeline/raw/master/gene_model/collapse_annotation.py).
To be able to use we it we will to create a python virtual environment:
```bash
# Load the required modules
module load scipy-stack/2020b
# create a virtual environment
virtualenv --no-download env
# Activate the virtual environment
source env/bin/activate
# install prerequisites
pip install bx-python
# download the required script
wget https://github.com/broadinstitute/gtex-pipeline/raw/master/gene_model/collapse_annotation.py
# faltten the reference genome:
REF=/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl90.gtf
python collapse_annotation.py ${REF} Homo_sapiens.GRCh38.Ensembl90_genes.gtf
```
That will generate the file `Homo_sapiens.GRCh38.Ensembl90_genes.gtf` which can
then be passed to RNA-SeQC

## Understanding RNA-SeQC options
The basic options for RNA-SeQC are:
```bash
rnaseqc [OPTIONS] gtf bam output
```
where `gtf` is the collapsed annotation file, `bam` is the input bam file resulting
from your mapping and cleaning, and the `output` is the name of the output for the
metrics.

### Options
There are a number of options that you can use to modify the default behaviour.
Below a **shortened** list, you can check the helper option for a more 
comprehensive list:

##### Helper options
You can pass the options `--help` (or `-h`, for short) to get help. You can also
verify the version of the program with the `--version` option.

#### Sample
If you want to add a specific name to your sample, you can pass it through the
`--sample` (`-s` for short). It defaults to the BAM filename. 

#### BED file
This option allow you to provide an optional file in [BED format](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)
of non-overlapping exons to be included in the fragment size calculation. If you
have such a fe, you can pass it through the option `--bed`.

#### Mapping quality
You can change the default (255 or missing value) for the threshold of [mapping
quality](http://maq.sourceforge.net/qual.shtml) that will be reported. To do so,
you can pass it through `--mapping-quality` (`-q` for short).

#### Base Mismatches 
To control how many mismatches are allowed between a read and the reference, you 
can use the option `--base-mismatch`. Reads with more than this number of 
mismatches are excluded from coverage metrics, ait defaults to 6.

#### Coverage statistics for all transcripts

You can ask RNASeQC to give you either a summary statistic of all transcripts
(default), or you can ask for coverage statistics of each transcript. If this is
the case, you can use the `--coverage` flag (no input needed).

### Running RNASeQC on Compute Canada
Compute Canada does not have this program installed, however, you can install it
in your home directory. To do this, you can:
```bash
# Download the source code
wget https://github.com/getzlab/rnaseqc/releases/download/v2.4.2/rnaseqc.v2.4.2.full_source.tar.gz
# uncompress and unarchive
tar -xzf rnaseqc.v2.4.2.full_source.tar.gz
# get into the directory
cd rnaseqc/rnaseqc
# Load required modules
module load StdEnv/2020  gcc/9.3.0 boost/1.72.0 bamtools/2.5.1
# Modify the make file
sed -i "2s|$| -I${EBROOTBAMTOOLS}/include -I${EBROOTBOOST}/include|" Makefile
# compile the program
make -j 4
# create a binary folder in your home
mkdir ${HOME}/bin/
# copy the binary into folder
cp rnaseqc ${HOME}/bin
# make it executable
chmod +x ${HOME}/bin/rnaseqc
```
You can access now the program by pointing to its path `${HOME}/bin/rnaseqc`. In
your sbatch script you can add something like:

<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=10      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=32G               #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name RNASeQC  	#<-- Job name
#SBATCH -o %j.out               #<-- File to which standard out will be written
#SBATCH -e %j.err               #<-- File to which standard err will be written
module load StdEnv/2020  gcc/9.3.0 boost/1.72.0 bamtools/2.5.1
GTF_FILE=/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl90.transcript_id.gtf
${HOME}/bin/rnaseqc -s SAMPLE1 ${GTF_FILE} SAMPLE1_marked_dup.bam SAMPLE1.RNASeQC
</code>
</pre>



# Raw Counts
Now that we have aligned and processed the reads, we need to get some data out 
of them. One especially important statistic are the raw counts, or how many reads
aligned with genes. Like with everything in bioinformatics there are plenty of 
options, but in this tutorial we will use the program HTSeq.

## Counting aligned reads with HTSeq
HTSeq is a python program that takes as input a splicing-aware alignment/mapping
(like  the one we just did with STAR) and a 'Features' file in the GFF format.
These features are the annotations of the reference genome to know where the 
exons are present.

HTSeq counts the reads with 3 possible overlapping modes:
1. The union of all the sets. This mode is recommended for most use cases: mode
   `union`
2. The intersection of all the sets: mode `intersection-strict`
3. The intersection of non-empty sets: mode `intersection-nonempt`

Special care must be taken to decide how to deal with reads that align to or 
overlap with more than one feature. This decision depends largely on the type
of experiment and data, as well as to the confidence the researcher has in the
data.

Apart from the mode, there is the non-unique modifier (as per docs):
>--nonunique none (default): the read (or read pair) is counted as ambiguous 
   and not counted for any features. Also, if the read (or read pair) aligns to
   more than one location in the reference, it is scored as alignment_not_unique.
> 
>--nonunique all: the read (or read pair) is counted as ambiguous and is also 
   counted in all features to which it was assigned. Also, if the read (or read
   pair) aligns to more than one location in the reference, it is scored as 
   alignment_not_unique and also separately for each location.
   Notice that when using --nonunique all the sum of all counts will not be equal 
   to the number of reads (or read pairs), because those with multiple alignments
   or overlaps get scored multiple times.

The following figure illustrates the effect of these three modes and the 
`--nonunique` option:
![alt](https://htseq.readthedocs.io/en/release_0.9.1/_images/count_modes.png)

### Usage
Once you have HTSeq installed, you can:

```bash
htseq-count [options] <alignment_files> <gff_file>
```

The (relavant) options being:
##### Format (-f \<format>, --format=\<format>)  
This option informs of the format of the input data, being sam or bam the possible
options. It defaults to sam. 
```bash
htseq-count -f bam <alignment_files> <gff_file>
```
#### Sorting order (-r \<order>, --order=\<order>)
This option allow you to describe how was the sam/bam file sorted (if not sorted
yes, do it with samtools or [picard](#picards-sortsam)). Use `name` (default), 
if your mapping/alignment was sorted by query name. You can also pass `pos` if
you sorted by coordinate.
```bash
htseq-count --order=pos <alignment_files> <gff_file>
```

#### Minumum quality (-a \<minaqual>, --a=\<minaqual>)
Set the quality treshold to keep in alignments. By default is 10.
```bash
htseq-count -a 20 <alignment_files> <gff_file>
```

#### Mode of overlap (-m \<mode>, --mode=\<mode>
Mode to handle reads overlapping more than one feature, as explained 
[previously](#counting-aligned-reads-with-htseq). The default is union:
```bash
htseq-count --mode=union <alignment_files> <gff_file>
```

#### Strandness (-s \<yes/no/reverse>, --stranded=\<yes/no/reverse>)
Tell the program if you data comes from a strand-specific assay, the default being
`yes`:

>For stranded=no, a read is considered overlapping with a feature regardless of
> whether it is mapped to the same or the opposite strand as the feature. For 
> stranded=yes and single-end reads, the read has to be mapped to the same strand 
> as the feature. For paired-end reads, the first read has to be on the same 
> strand and the second read on the opposite strand. For stranded=reverse, 
> these rules are reversed.

>If your RNA-Seq data has not been made with a strand-specific protocol, this 
> causes half of the reads to be lost. Hence, make sure to set the option 
> --stranded=no unless you have strand-specific data!
```bash
htseq-count -s no <alignment_files> <gff_file>
```

#### Handle of Non-uniques (--nonunique=\<nonunique mode>)
As explained in [Counting aligned reads with HTSeq](#counting-aligned-reads-with-htseq),
this is a modifier on how to handle multiple overlaps. <nonunique mode> are `none`
and `all` with the former by default. 
```bash
htseq-count --nonunique=none <alignment_files> <gff_file>
```

#### Write SAM output (-o \<samout>, --samout=\<samout>)
If instead of only a report file we want to introduce the count as a feature in
the sam file, we can pass the filename to the `-o` or `--samout` option.
```bash
htseq-count -o afile.sam <alignment_files> <gff_file>
```

#### Miscelaneaous
There are other options that you can see with the help option (`-h`, `--help`), 
which will also show a usage summary. You can also turn off the verbosity of 
the program with the quiet option (`-q`, `--quiet`), suppressing the  progress 
report and warnings.

### Output

The main output of HTSeq (see [usage](#usage)) is a table with counts of each 
feature, in our case genes/exons. It contains some special counters:
>__no_feature: reads (or read pairs) which could not be assigned to any feature 
   (set S as described above was empty).
> 
>__ambiguous: reads (or read pairs) which could have been assigned to more than
   one feature and hence were not counted for any of these, unless the 
> --nonunique all option was used (set S had more than one element).
> 
>__too_low_aQual: reads (or read pairs) which were skipped due to the -a option,
>
>__not_aligned: reads (or read pairs) in the SAM file without alignment
> 
> __alignment_not_unique: reads (or read pairs) with more than one reported alignment. 
 
### Running HTSeq in Compute Canada
As mentioned earlier, HTSeq is not written in a compiled language like we have use
so far, but in python. Compute Canada carries a series of wheels can can be installed.
You can check the [documentation](https://docs.computecanada.ca/wiki/Python) for more
information. In a nutshell:
1. Load required modules
   ```bash
   module load python/3.7
   ```
2. Create a virtual environment
   ```bash
   virtualenv --no-download ${HOME}/htseq_env
   ```
3. Activate the virtual environment:
   ```bash
   source ${HOME}/htseq_env/bin/activate
   ```
4. Install the wheel
   ```bash
   pip install --no-index htseq
   ```
Now you can run `htseq-count` directly. Remember to deactivate the environment
when you are done using the command `deactivate`. As per the
[documentation](https://docs.computecanada.ca/wiki/Python), you can create the 
virtual environment and do the installation inside a SLURM job:

<pre>
<code>
#!/bin/bash                       
#SBATCH --account=def-someuser  #<-- this is the account of the PI
#SBATCH --time=0-2:00:00        #<-- you need to provide the expected time (dd-hh:mm-ss)
#SBATCH --cpus-per-task=10      #<-- Here is where you put the number of cpus you want
#SBATCH --mem=32G               #<-- Amount of memory. In this case reserve all memory
#SBATCH --job-name HTSeq      	#<-- Job name
#SBATCH -o %j.out               #<-- File to which standard out will be written
#SBATCH -e %j.err               #<-- File to which standard err will be written
module load StdEnv/2020 python/3.7
PATH2ALIGNMENT=some/path/to/bam
PATH2GFF=/cvmfs/ref.mugqic/genomes/species/Homo_sapiens.GRCh38/annotations/Homo_sapiens.GRCh38.Ensembl90.gtf
virtualenv --no-download ${SLURM_TMPDIR}/htseq_env
source  ${SLURM_TMPDIR}/htseq_env/bin/activate
pip install --no-index htseq
htseq-count -s no ${PATH2ALIGNMENT} ${PATH2GFF} > myalignment.counts
</code>
</pre>

# Processing the Raw Counts
Raw counts produced by [HTSeq](#counting-aligned-reads-with-htseq) need to be 
read in a features (e.g. genes) by treatments, something like this:
<p align="center">
  <img src="https://4va.github.io/biodatasci/img/countdatacoldata.png"><br>
  <sup>Image from <a href="https://4va.github.io/biodatasci/img/countdatacoldata.png"> https://4va.github.io/biodatasci/img/countdatacoldata.png</a></sup>
</p>

To read it in the correct format, we could use multiple programming languages, 
but since the DGE will be done in R, we will use the same package throughout:  
DESeq2. But before we delve into learning a bit of R and the package let's talk
about differential expression analyses.

## Differential Expression Analyses (DGE)
Also known as Differential Gene Expression (hence the DGE acronym), it is an 
analysis on the counts of copies of transcripts of a set of genes. In the context
of RNAseq, most sequenced genes.
In short, DGE takes normalised read counts to perform an statistical analysis 
to quantify significant changes in expression levels between experimental groups. 
It tries to tease appart of the observed differences are greater than what would 
be expected just due to natural random variation. [DESeq2](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html) 
bases its test on a negative binomial (NB) distributions, while other resources 
(e.g  baySeq, EBSeq) use a Bayesian approach.

### Normalization
Since the counts varied for a lot of reasons, the first step is normalization.
This is done to make the comparison among samples valid. This step allows us to 
focus the analyses into the relevant difference instead of the spurious ones.
However, DESeq2 requires the un-normalized values as it makes some statistical
modeling of the errors, and corrects accordingly. That being said, there are a 
few ways to normalized your data that you should be aware of :

<table>
  <thead>
    <tr>
      <th>Normalization method</th>
      <th>Description</th>
      <th>Accounted factors</th>
      <th>Recommendations for use</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>CPM</strong> (counts per million)</td>
      <td>counts scaled by total number of reads</td>
      <td>sequencing depth</td>
      <td>gene count comparisons between replicates of the same samplegroup; <strong>NOT for within sample comparisons or DE analysis</strong></td>
    </tr>
    <tr>
      <td><strong>TPM</strong> (transcripts per kilobase million)</td>
      <td>counts per length of transcript (kb) per million reads mapped</td>
      <td>sequencing depth and gene length</td>
      <td>gene count comparisons within a sample or between samples of the same sample group; <strong>NOT for DE analysis</strong></td>
    </tr>
    <tr>
      <td><strong>RPKM/FPKM</strong> (reads/fragments per kilobase of exon per million reads/fragments mapped)</td>
      <td>similar to TPM</td>
      <td>sequencing depth and gene length</td>
      <td>gene count comparisons between genes within a sample; <strong>NOT for between sample comparisons or DE analysis</strong></td>
    </tr>
    <tr>
      <td>DESeq2’s <strong>median of ratios</strong> [<a href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-10-r106">1</a>]</td>
      <td>counts divided by sample-specific size factors determined by median ratio of gene counts relative to geometric mean per gene</td>
      <td>sequencing depth and RNA composition</td>
      <td>gene count comparisons between samples and for <strong>DE analysis</strong>; <strong>NOT for within sample comparisons</strong></td>
    </tr>
    <tr>
      <td>EdgeR’s <strong>trimmed mean of M values (TMM)</strong> [<a href="https://genomebiology.biomedcentral.com/articles/10.1186/gb-2010-11-3-r25">2</a>]</td>
      <td>uses a weighted trimmed mean of the log expression ratios between samples</td>
      <td>sequencing depth, RNA composition, and gene length</td>
      <td>gene count comparisons between and within samples and for <strong>DE analysis</strong></td>
    </tr>
  </tbody>
</table>

from [hbctraining](https://hbctraining.github.io/Training-modules/planning_successful_rnaseq/lessons/sample_level_QC.html)


# Transcript assembly with Cufflinks
## Introduction to transcriptome assembly

## Understanding Samtools hardclip
samtools view -bh -F 256 -f 81 input_bam > clipped1
samtools view -bh -F 256 -f 161 input_bam > clipped2
java -jar $EBROOTPICARD/picard.jar MergeSamFiles \
	I=clipped1 \
	I=clipped2 \
	O=output.bam \
	USE_THREADING=true
 
https://www.samformat.info/sam-format-flag
## Understanding Cufflinks 
## Understanding Cuffmerge
## Understanding Cuffdiff
## Understanding Cuffnorm
## Running Cufflinks bundle on your data
-->

# Running ALL in one pipeline! GENPIPES
`GENPIPES` is a bundle of multiple bioinformatics pipelines developed and maintain
by [C3G](https://www.computationalgenomics.ca/) and it is available in Compute 
Canada. According to [their documentation](https://genpipes.readthedocs.io/en/genpipes-v-3.4.0/):

>GenPipes is an open source Python-based Workflow Management System (WMS) for 
> Next Generation Sequencing (NGS) genomics pipeline development. 
> As part of its implementation, GenPipes includes a set of high-quality, 
> standardized analysis pipelines, designed for high-performance computing (HPC)
> resources and cloud environments. 
> GenPipes pipelines have been tested, continuously improved through industry 
> standard benchmarks, and used extensively over the past 4 years. By combining 
> both WMS and extensively validated end-to-end genomic analysis workflows, into 
> a simple tool, GenPipes not only enables bioinformatics professionals but also 
> students, researchers in need of bioinformatics tools to perform turnkey 
> analyses for a wide range of bioinformatics applications in the genomics field. 
> It also allows flexible and robust extensions for advanced genomics research.

`GENPIPES` have a battery of bioinformatic analyses, but here we will focused on
[RNASeq](https://genpipes.readthedocs.io/en/genpipes-v-3.4.0/user_guide/pipelines/gp_rnaseq.html?highlight=rnaseq).

Since `GENPIPES` is essentially a workflow manager, the usage is significantly
different of what we have done so far. For example, you will run the pipeline on
the login node, and the pipeline itself will submit multiple individual jobs
to the scheduler.

`GENPIPES` follow the schema below:

![alt](https://genpipes.readthedocs.io/en/genpipes-v-3.4.0/_images/rnaseq.png)

## Setting up GENPIPES in Compute Canada
As mentioned, `GENPIPES` is available in all general purpose clusters in Compute
Canada (Graham, Cedar, and Beluga). However, there are a few steps to make all
the tools visible to you. Specifically, you need to export some variables, and
tell the modules software how to use it:

```bash
## GenPipes/MUGQIC genomes and modules
# Set the HOME where programs and modules are
export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
# Tell the module program where to look for genpipes modules
module use $MUGQIC_INSTALL_HOME/modulefiles
# Load the particular python that genpipes runs in
module load mugqic/python/2.7.14
# Load the desired genpipes version
module load mugqic/genpipes/<latest_version>
# Export some variables with your info
export JOB_MAIL=<my.name@my.email.ca>
export RAP_ID=<my-rap-id>
```

You will need to replace `<my.name@my.email.ca>` with your own email,
`<my-rap-id>` with your lab's account (usually `def-PIuser` or `rrg=PIuser`), 
and `<latest_version>` with the genpipes version you want (current latest version
is 3.4.0).
You can add some or all of the lines in your `~/.bash_profile` file if you want
the changes to be permanent. Otherwise, you will have to do the loads and exports
every time you want to run `GENPIPES`.

You can check available pipline at `$MUGQIC_INSTALL_HOME/software/genpipes/genpipes-3.4.0/pipelines`"

```bash
ls $MUGQIC_INSTALL_HOME/software/genpipes/genpipes-3.4.0/pipelines

ampliconseq  container.ini  dnaseq_high_coverage     __init__.py  rnaseq
chipseq      covseq         hicseq                   methylseq    rnaseq_denovo_assembly
common.py    dnaseq         illumina_run_processing  nanopore     rnaseq_light
```
You can see that there are multiple pipelines that you can choose from. For this
course we will focus on `$MUGQIC_INSTALL_HOME/software/genpipes/genpipes-3.4.0/pipelines/rnaseq`.

## Running GENPIPES in Compute Canada
All `GENEPIPES` pipelines use a configuration file with the extension `.ini`. So
before we proceed, we need to understand it. 
### Configuration files
#### Why does GenPipes need configuration file?
As per their docs:
>An ini file is a file that contains parameters needed to run a pipeline. Our 
> genome alignment pipeline contains over 20 steps, each involving over 5 
> parameters per step. Imagine having to type all 100 parameters to run a 
> pipeline! For simplicity, all the parameters are stored in an “ini” file 
> (extention.ini) that accompanies the pipeline.

The same is true for their RNA seq pipeline. As we saw in each of the steps that
we have follow until now. Luckily for you, some pre-filled config files are available
in Compute Canada for you to modify with you own needs. 
From this point onward, I will assume that you'll be running the pipeline in the
`Graham` cluster. If youa re not, adjust accordingly. Since we are intending to
run a full RNAseq analysis, we can find the apropriate template for the configuration
file at `$MUGQIC_INSTALL_HOME/software/genpipes/genpipes-3.4.0/pipelines/rnaseq`. 
Make a copy of the `rnaseq.graham.ini` configuration file and explore it:

```bash
cd ~/scratch
cp $MUGQIC_INSTALL_HOME/software/genpipes/genpipes-3.4.0/pipelines/rnaseq/rnaseq.graham.ini .
```

You will see that the template looks like this (this is a cropped version showing 
the head and the tail):

```text
[DEFAULT]
# Cluster
cluster_submit_cmd=sbatch
cluster_submit_cmd_suffix= | grep "[0-9]" | cut -d\  -f4
cluster_walltime=--time=24:00:00
cluster_cpu= -n 1 -N 1
# IMPORTANT: update $RAP_ID with your own Resource Allocation Project ID or set it in your $HOME/.bash_profile!
cluster_other_arg=--mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID
cluster_queue=--mem=4G
cluster_work_dir_arg=-D
.
.
.
cluster_cpu=-N 1 -n 1

[ihec_metrics]
cluster_walltime=--time=5:00:0
cluster_cpu=-N 1 -n 4
cluster_queue=--mem=16G

[verify_bam_id]
cluster_cpu=-N 1 -n 4
cluster_queue=--mem=38G
```

Each step/tool section has a header between square brackets (e.g. [DEFAULT]), where 
you can adjust each parameter as you see fit (although these are mostly already optimal).
Let's walt to each one of them:

#### `DEFAULT` block
This is the block which adjust parameters regarding the run as a whole (not per tool),
like submission type, etc. The `rnaseq.graham.ini` files contains:

```text
[DEFAULT]
# Cluster
cluster_submit_cmd=sbatch
cluster_submit_cmd_suffix= | grep "[0-9]" | cut -d\  -f4
cluster_walltime=--time=24:00:00
cluster_cpu= -n 1 -N 1
# IMPORTANT: update $RAP_ID with your own Resource Allocation Project ID or set it in your $HOME/.bash_profile!
cluster_other_arg=--mail-type=END,FAIL --mail-user=$JOB_MAIL -A $RAP_ID
cluster_queue=--mem=4G
cluster_work_dir_arg=-D
cluster_output_dir_arg=-o
cluster_job_name_arg=-J
cluster_cmd_produces_job_id=true
cluster_dependency_arg=--depend=afterok:
cluster_dependency_sep=:
cluster_max_jobs=3000
tmp_dir=${SLURM_TMPDIR}

java_other_options=-XX:ParallelGCThreads=4 -Dsamjdk.buffer_size=1048576

assembly_dir=$MUGQIC_INSTALL_HOME/genomes/species/%(scientific_name)s.%(assembly)s
```

Given that the template is done especifically for `graham` then you do not need 
to modify this configuration file. Most of these options are relating to the 
scheduler. For example, `cluster_submit_cmd` is the command to submit to the 
scheduler, and as you already know is `sbatch`

#### `picard_sam_to_fastq` block
This block just gives information of resources needed. In the `rnaseq.graham.ini`
it is 1 node, 4 cpus and 16Gb of memory:

```text
[picard_sam_to_fastq]
cluster_cpu=-N 1 -n 4
cluster_queue=--mem=16G
```

#### `samtools_cram_output` block
The samtools block requests one node and 2 cpus but 48 hours of computation time:
```time
[samtools_cram_output]
cluster_cpu=-N 1 -n 2
cluster_walltime=--time=48:00:0
```
Note that this step bam to cram is high IO and therefore it does not benefit much
from a higher CPU count.

#### *trimmomatic* block
You are well acquainted with trimmomatic. In this block, one node and 6 CPUs are
requested during 24 hours, and with 24G of memory:

```text
[trimmomatic]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 6
cluster_queue=--mem=24G
```

#### *STAR align* block
In this case 16 CPUs are requested and set as `OMP_NUM_THREADS` with a lot of
ram (100G minimum), one day of computation:
```text
[star_align]
threads=16
ram=100G
cluster_cpu=-N 1 -n 16
cluster_walltime=--time=24:00:0
cluster_queue=--mem=128G
```
In this case two possible changes might be beneficial, change the `mem` to 0 and
the `n` to 32. This way you will reserve a full node. If in your previous runs
it took significantly less than 24 hours, you could also decrease the time. The 
latter is only important if you think it will run in 12 hours or less.

#### *STAR index* block
You also know this one very well. This block asks for 16 threads and 128G of ram, 
during 15 hours. The adjustments can be similar to `star align`:

```text
[star_index]
threads=16
ram=100G
cluster_cpu=-N 1 -n 16
cluster_walltime=--time=15:00:0
cluster_queue=--mem=128G
```

#### *STAR junction* block
This is a small job to identify joints in the star alignment. It only requires 
one cpu for 5 hours. It does not benefit from extra resources:

```text
[star_junction]
cluster_cpu=-N 1 -n 1
cluster_walltime=--time=5:00:0
```

#### *picard merge_sam files* block
Another well known block to you, uses picard to merged required sam/bamfiles. Here
`GENPIPES` will request 48G of memory, and 12 CPUs in one node for one day:

```text
[picard_merge_sam_files]
ram=40G
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G
```

#### *picard sort_sam* block
This block requests 48G and 12 CPUs for 24 hours:

```text
[picard_sort_sam]
ram=40G
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G
```

#### *picard_mark_duplicates* block
This is also a heavy read/write block, hence it takes loner but with less CPU/MEM
resources:

```text
[picard_mark_duplicates]
ram=14G
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 5
cluster_queue=--mem=20G
```

#### *rnaseqc* block
This module runs a lot of analyses of your RNA seq, hence it requires a long run
with 12 CPUs, and 48G of memory.

```text
[rnaseqc]
cluster_walltime=--time=72:00:0
cluster_cpu=-N 1 -n 12
ram=40G
cluster_queue=--mem=48G
```

#### *bed_graph* block
To process some of the bed files produce in intermediary files, `bed_graph` block
ask for 8 CPUs and 12 hours of computation.
```text
[bed_graph]
cluster_walltime=--time=12:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=38G
```

#### *wiggle* block
Wiggle produces a nice visual report in the form of tracks, to be visualize in 
any genome browser. To do that `genpipes` requests 12 hours of 12 CPUs and 48G
of memory:

```text
[wiggle]
cluster_walltime=--time=12:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G
```

#### *htseq_count* block
`htseq_count` produce the raw counts of copies in each of the genes in the 
reference genome. In this block one day with 6 CPUs and 24G of memory are requested.
```text
[htseq_count]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 6
cluster_queue=--mem=24G
```
#### *tuxedo_hard_clip* block
`tuxedo_hard_clip` removes all unmapped and uncomplete reads from the mapping. It
has an intermediate memory intensity and that is why 32G are requested along with
8 CPUs for a day.
```text
[tuxedo_hard_clip]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G
```

#### *stringtie* block
`stringtie` is a fast and highly efficient assembler of RNA-Seq alignments into
potential transcripts it requires the same ammount of memory, CPU and time than
the tuxedo block
```text
[stringtie]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G
```


By now I assume that you got the hang of it, so here are the rest of the options
without explanation:

```text
[stringtie_merge]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[stringtie_abund]
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[cufflinks]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[cuffmerge]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[cuffquant]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[cuffdiff]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[cuffcompare]
cluster_walltime=--time=2:00:0
cluster_cpu=-N 1 -n 1
[cuffnorm]
cluster_walltime=--time=48:00:0
cluster_cpu=-N 1 -n 8
cluster_queue=--mem=32G

[picard_collect_multiple_metrics]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
ram=40G
cluster_queue=--mem=48G

[picard_collect_rna_metrics]
ram=40G
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G

[picard_rna_metrics]
ram=40G
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G

[estimate_ribosomal_rna]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 12
cluster_queue=--mem=48G
[bwa_mem_rRNA]
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 16
cluster_queue=--mem=64G

[picard_sort_sam_rrna]
ram=7G
cluster_cpu=-N 1 -n 2
cluster_queue=--mem=8G
java_other_options=-XX:ParallelGCThreads=1  -Dsamjdk.buffer_size=1048576

[metrics]
cluster_walltime=--time=5:00:0
cluster_cpu=-N 1 -n 1

[rpkm_saturation]
threads=15
other_options=1
cluster_walltime=--time=24:00:0
cluster_cpu=-N 1 -n 16
cluster_queue=--mem=128G

[differential_expression]
cluster_walltime=--time=10:00:0
cluster_cpu=-N 1 -n 1

[differential_expression_goseq]
cluster_walltime=--time=10:00:0
cluster_cpu=-N 1 -n 1

[gq_seq_utils_exploratory_analysis_rnaseq]
cluster_walltime=--time=00:30:0
cluster_cpu=-N 1 -n 1

[ihec_metrics]
cluster_walltime=--time=5:00:0
cluster_cpu=-N 1 -n 4
cluster_queue=--mem=16G

[verify_bam_id]
cluster_cpu=-N 1 -n 4
cluster_queue=--mem=38G
```

Notice that there are some cases in which the word cluster is not prepended. That
means that some extra variables are being defined. For example, when you see `threads`
this will assign the environmental variable `OMP_NUM_THREADS` that is read by some 
software.

### Design file
As per `GENPIPES` documentation:
>In addition to the configuration files and the input readset file, certain 
> pipelines such as ChIP-Seq and RNA sequencing (RNA-Seq), require a design 
> file that describes each contrast. Custom sample groupings can be defined 
> in the design file

This is beacause you need to tell the pipeline which one is your control and 
which ones your treatment, along with other information. Let's take a look at
the design file:

```text
Sample Contrast_AB Contrast_AC
sampleA 1 1
sampleB 2 0
sampleC 0 2
sampleD 0 0
```

As you can see, the design file is a tab-separated plain text file with one 
line per sample. The columns are `Sample` and as many contrasts as you need. 
In the example above we have two contrasts: `sampleA` vs `sampleB`, and 
`sampleA` vs `sampleC`. note that `0` means that the sample does not belong to 
any group, `1` means that the sample belongs to the control group, and `2`
that the sample belongs to the treatment test case group. In the example above
we state that sample A is the control in both contrasts

### Read Set File
The other file that you must provide to the pipeline is the read set file.
The readset file is a tab-separated file that contains the following information:

Field |Contents| 
---|---|
Sample:|Sample must contain letters A-Z, numbers 0-9, hyphens (-) or underscores (_) only; BAM files will be merged into a file named after this value; **mandatory**. |
Readset:|a unique readset name with the same allowed characters as above; **mandatory**.|
Library:|Information about the library. **Optional**.|
RunType:|PAIRED_END or SINGLE_END; **mandatory**.|
Run:|1 to use in this run 0 not to use it; **mandatory**.|
Lane:|Lane used; **mandatory**.|
Adapter1:| sequence of the forward trimming adapter; **mandatory**.|
Adapter2:|sequence of the reverse trimming adapter; **mandatory**.|
QualityOffset:|quality score offset integer used for trimming; **optional**.
BED:|relative or absolute path to BED file; **optional**.|
FASTQ1:|relative or absolute path to first FASTQ file for paired-end readset or single FASTQ file for single-end readset; **mandatory if BAM value is missing**.|
FASTQ2:|relative or absolute path to second FASTQ file for paired-end readset; **mandatory if RunType value is “PAIRED_END”**.|
BAM:|relative or absolute path to BAM file which will be converted into FASTQ files if they are  not available; **mandatory if FASTQ1 value is missing, ignored otherwise**.|

If you sequenced your reads in Genome Quebec and still have access to nanuq, 
there is a helper script `$MUGQIC_PIPELINES_HOME/utils/csvToreadset.R` that 
takes a csv file downloadable from nanuq and creates the Readset file. 
Otherwise you need to do it by hand.

### Finding genpipes in Compute Canada
To use GenPipes in Compute Canada clusters you will need to load it. However,
genpipes is manage by a third party, and therefore they have their own stack.
Luckily, the C3G folks have been kind enough to make it work in our systems 
with a few steps:
1. Export the base path where C3G software stack is located:
   `export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6`
2. Tell LMOD (the modules platform) to use the C3G modules:
   `module use $MUGQIC_INSTALL_HOME/modulefiles`
3. In the case of genpipes, they worked with a very specific version of python,
   so let's load it:
   `module load mugqic/python/2.7.14`
4. Load the desired version of genpipes (you can also use spider to check):
   `module load mugqic/genpipes/<latest_version>`
5. Genpipes looks for your email and account, so let's set it:
   ```
   export JOB_MAIL=<my.name@my.email.ca>
   export RAP_ID=<my-rap-id>
   ```
The RAP_ID is the account that you can submit jobs to our clusters, usually 
def-someuser, someuser being the PI of the group.

Once that is set, you can run the pipline creator by:

```
rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
$MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.graham.ini -r readset.rnaseq.txt \
-d design.rnaseq.txt > rnaseqCommands.sh
```

This will generate a bash file with all the submissions of each step to our 
queue just by doing:

```
bash rnaseqCommands.sh
```

then you will notice how the jobs start being submitted!!! now let's take a 
closer look to the options that the RNAseq pipeline builder `rnaseq.py` from
genpipes has.

#### CONFIG (-c)
The config option allows you to pass a list of config files with `.ini` 
extension.  We cover part of their contents earlier. It contains all the info 
for the resources and variables needed at each step. You can pass as many files
as you like and it will be overwritten in the order you pass them. For example,
in the example we saw earlier, we passed two config files, the `rnaseq.base.ini` 
and `rnaseq.graham.ini`. Anything variable that exist in both side will only
retain the value in `rnaseq.graham.ini` as is the last one passed.
For all intends and purposes you can use the ones provided by genepipes located at 
`$MUGQIC_PIPELINES_HOME/pipelines/rnaseq`. The `rnaseq.base.ini` contains 
general information, while the `rnaseq.graham.ini` conatains information 
pertaining to the Graham cluster. You can see other clusters in the same path.

#### STEPS (-s, --steps)
GenPipes has the ability to only run a particular range of interest in the 
steps of the pipeline:

![alt](https://bitbucket.org/mugqic/genpipes/raw/master/resources/workflows/GenPipes_rnaseq.png)

As you can see there are about 25 steps (some of which you'll need to skip, for
example step 1). Say you want to only run steps 2 to 6, then your'll pass 
`-s 2-6`.

#### OUTPUT_DIR (-o, --output-dir)
If you want the output in anything different than the current working directory,
pass it throgh the `-o` option.

#### Scheduler (-j,--job-scheduler)
The type of scheduler you are using (e.g  pbs, batch, daemon, slurm). In 
Compute Canada we use the default SLURM

#### Force Job re-run (-f, --force)
By default, genpipes will not generate jobs that have already been done. You
can force it to be redone, by using the `-f` flag

#### No JSON (--no-json)
By default, the pipeline creates one JSON (a type of file) per each sample. You
can disable this with the `--no-json` flag.

#### Create a Report (--report)
With this flag, `rnaseq.py` will generate a report by merging all job markdown
report files in the given step range into HTML, if they exist.
If `--report` is set, `--job-scheduler`, `--force`, `--clean` options and job 
up-to-date status are ignored

####  Clean up (--clean)
With this flag a job script will be created, with the 'rm' commands for all job
removable files in the given step range, if they exist. If --clean is
set, `--job-scheduler`, `--force` options and job up-to-date status are ignored. 

#### Level of information in Logs (-l, --log)
You can choose the level of information available in the logs, from debug,
info, warning, error, and critical. By default is set to `info`.

#### Sanity check (--sanity-check) 
Run the pipeline in `sanity check mode` to verify that all the input files 
needed for the pipeline to run are available on the system, a.k.a dry-run

#### Run in a container (--container)
You can use this option if you have a properly setup container as `wrapper` or
`singularity`. You need to provide valid singularity image path.

#### Design file (-d, --design)
The design of your RNAseq experiment (as explained above)

#### Assembly type (-t, --type)
Using wthis option you can select wether to use `cufflinks` or `stringtie` as 
RNA-seq assembly method, the latter being the default.

#### Read Set file (-r, --readsets)
The readseq file as explained above.

### Running an example before you use your reads
Let's try this out:
1. First, let's download some example data:
    ```
    wget https://www.computationalgenomics.ca/tutorial/c3g_analysis_workshop/C3GAW_RNA_TestData_Aug2018.zip
    unzip C3GAW_RNA_TestData_Aug2018
    cd C3GAW_RNA_TestData_Aug2018
    ```
2. Now follow the steps explained above to load genpipes (remember to change
   the variables between the `<` and `>` symbols)
   ```
   ## GenPipes/MUGQIC genomes and modules
   export MUGQIC_INSTALL_HOME=/cvmfs/soft.mugqic/CentOS6
    module use $MUGQIC_INSTALL_HOME/modulefiles
    module load mugqic/python/2.7.14
    module load mugqic/genpipes/<latest_version>
    export JOB_MAIL=<my.name@my.email.ca>
    export RAP_ID=<my-rap-id>
   ```
3. Explore the sample readset file:
   ```
   less readset.rnaseq.txt
   ```
4. Explore the design file:
   ```
   less design.rnaseq.txt
   ```
5. Explore the config files that we will be using:
   ```
   less $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini
   less $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.graham.ini
   ```
6. Run a sanity check to make sure all is good:
   ```bash
   rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
   $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.graham.ini -r readset.rnaseq.txt \
   -d design.rnaseq.txt --sanity-check
   ```
7. Run the pipeline builder:
   ```bash
   rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
   $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.cedar.ini -r readset.rnaseq.txt \
   -d design.rnaseq.txt > rnaseqCommands.sh
   ```
8. Check the commands created:
   ```bash
   less rnaseqCommands.sh
   ```
9. Launch all the required jobs (**Warning: this will send you a lot of emails**):
   ```bash
   bash rnaseqCommands.sh
   ```
10. When you stop receiving emails, check the results.
11. Generate a report in html:
    ```bash
    rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
    $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.graham.ini -r readset.rnaseq.txt \
    -d design.rnaseq.txt --report > report.sh
    bash report.sh
    ```
12. Move the report folder to your own computer (or double click in it if using
    MobaXterm) and visualize it in your browser 
13. Clean up:
    ```bash
    rnaseq.py -c $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.base.ini \
    $MUGQIC_PIPELINES_HOME/pipelines/rnaseq/rnaseq.graham.ini -r readset.rnaseq.txt \
    -d design.rnaseq.txt --clean > cleanup.sh
    bash cleanup.sh
    ```

And voila! now do it for your readsets.

Good luck!
