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
 
<!---
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

![MobaxtermX11](https://github.com/jshleap/Tutorials/blob/main/images/800px-MobaXterm_X11.png?raw=true)

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

![Nano](https://github.com/jshleap/Tutorials/blob/main/images/Nano.png?raw=true)

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

![fastq](https://www.researchgate.net/profile/Reinhard_Schneider2/publication/256095540/figure/fig7/AS:298005665206280@1448061495472/FASTQ-file-1st-line-always-starts-with-the-symbol-followed-by-the-sequence.png)
<br><sup>from https://www.researchgate.net/publication/256095540_1756-0381-6-13

The fastq files' header starts with an @ symbol, followed by the sequence identifyer
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

mkdir -p hg38_index
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
  <img src="https://du-bii.github.io/module-5-Methodes-Outils/seance1/images/SAM_format.jpg"><br>
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
	I=file1.sam \
	I=file1.sam \
	O=output.sam \
	USE_THREADING=true
```
This will generate a single unsorted SAM file called output.sam. The same procedure
can be perfoem to merge bamfiles.

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

### Picard’s MarkDuplicates
### Picard’s CollectRnaSeqMetrics
### Cleaning your data and generate metrics
