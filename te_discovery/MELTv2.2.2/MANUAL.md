Table of Contents
=================

1. [Introduction to MELT](#_Introduction_to_MELT)
	1. [Readme Conventions](#_Readme_Conventions)  
2. [MELT Quick-Start Guide](#_MELT_Quick-Start_Guide)  
	1. [Quick Install MELT](#_Install_MELT)  
	2. [Preparing for MEI Discovery](#_Preparing_for_MEI_Discovery)  
	3. [Running MELT Quick](#_Running_MELT_Quick)  
3. [MELT System Requirements](#_MELT_System_Requirements)  
	2. [Software Testing](#_Software_Testing)  
	3. [Hardware Requirements](#_Hardware_Requirements)  
4. [Installing MELT](#_Installing_MELT)  
	1. [MELT Installation](#_MELT_Installation)  
	2. [Installing MELT Using SGE](#_Installing_MELT_Using_SGE)  
	3. [Required Software Dependencies](#_Required_Software_Dependencies)  
	4. [Required External File Dependencies](#_Required_External_File_Dependencies)  
	5. [Optional External Files](#_Optional_External_Files)  
5. [Transposon Reference Files](#_Transposon_Reference_Files)  
	1. [General Description](#_General_Description)  
	2. [Creating Custom MEI.zip Files](#_Creating_Custom_MEI.zip_Files)  
		1. [The MELT .zip File Structure](#_The_MELT_.zip_File_Structure)  
		2. [Deciding on an MEI Consensus](#_Deciding_on_an_MEI_Consensus)  
		3. [MELT-BuildTransposonZIP](#_MELT-BuildTransposonZIP)  
6. [MELT Options](#_MELT_Options)  
7. [Running MELT](#_Running_MELT)  
	1. [Preprocessing .bam Files for MELT](#_Preprocessing_.bam_Files_for_MELT)  
	2. [Which Version of MELT Should I Use?](#_Which_Version_of_MELT_Should_I_Use)  
	3. [Running MELT Using MELT-SGE](#_Running_MELT_Using_MELT-SGE)  
		1. [Description](#_Description_SGE)  
		2. [Options](#_Options_SGE)  
	4. [Running MELT Using MELT-SPLIT](#_Running_MELT_Using_MELT-SPLIT)  
		1. [Description](#_Description_Split)  
		2. [IndivAnalysis](#_IndivAnalysis)  
		3. [GroupAnalysis](#_GroupAnalysis)  
		4. [Genotype](#_Genotype)  
		5. [MakeVCF](#_MakeVCF)  
		6. [Example MELT-SPLIT Workflow](#_Example_MELT-SPLIT_Workflow)  
	5. [Running MELT Using MELT-SINGLE](#_Running_MELT_Using_MELT-SINGLE)  
		1. [Description](#_Description_Single)  
		2. [Options](#_Options_Single)  
8. [MELT Output](#_MELT_Output)  
9. [Other MELT Tools](#_Other_MELT_Tools)  
	1. [MEI Species Classification](#_MEI_Species_Classification)  
		1. [CAlu](#_CAlu)  
		2. [LineU](#_LineU)  
	2. [MELT-Deletion](#_MELT-Deletion)  
		1. [MELT-Deletion Output](#_MELT-Deletion_Output)  
	3. [Transduction](#_Transduction)  
		1. [Source](#_Source)  
		2. [TransductionFind](#_TransductionFind)  
		3. [TransductionMerge](#_TransductionMerge)  
10. [What MELT is for](#_What_MELT_is_for)  
11. [What MELT is not for](#_What_MELT_is_not_for)  
12. [Reporting Issues With MELT](#_Reporting_Issues_With_MELT)  
13. [How to Cite MELT](#_How_to_Cite_MELT)  
14. [References](#_References)  

<span id="_Introduction_to_MELT" class="anchor"></span>1. Introduction to MELT
=======================

The Mobile Element Locator Tool (MELT) is a software package, written in
Java, that discovers, annotates, and genotypes non-reference Mobile
Element Insertions (MEIs) in Illumina DNA paired-end whole genome
sequencing (WGS) data. MELT was first conceived as a tool to identify
non-reference MEIs in large genome sequencing projects, specifically as
part of the 1000 Genomes Project, and has been further optimized to run
on a wide range of data types. MELT is optimized for performing
discovery in a large number of individual sequencing samples using the
Sun Grid Engine (SGE). MELT also has two additional workflows: analysis
without SGE (for adaptability to other parallel computing platforms) and
single genome analysis. MELT is highly scalable for many different types
of data, and the different workflows are outlined and detailed in this
documentation. MELT was also programmed to be user friendly, with only a
small number of common and easily installed bioinformatic software
dependencies (see [Required software
dependencies](#_Required_Software_Dependencies)). Input is standard
BWA (MEM or ALN) WGS alignment(1), while output is in the Variant Call
File (VCF) 4.2 format(2). This allows users to integrate MELT into
pre-existing pipelines that utilize BWA alignments and accept VCF as
input. Additionally, In keeping with the user-friendly design ideology,
MELT operates with several safeguards to ensure uninterrupted runtime
from beginning to end, with minimal user interaction.

MELT uses an approach that is fairly standard for many forms of
Structural Variation (SV) discovery, but has been highly specialized for
MEI discovery. In short, MELT first collects all discordant pairs from a
WGS alignment, and aligns them to provided MEI references using the
Bowtie2(3) sequencing aligner. It then ’walks’ across the reference
genome and creates a list of putative MEIs based on total read support
at each putative site. MELT will then merge the initial MEI calls across
the datasets provided, and determine specific information and exact
breakpoints for each of the putative MEIs. Finally, all sites are
genotyped according to a modified (for SV genotyping) version of Heng
Li’s genotype likelihood equation(4), and subsequently filtered based on
true positive calls. Output is then generated in VCF format.

<span id="_Readme_Conventions" class="anchor"></span>1.1. Readme Conventions
-----------------------

This readme uses the following conventions:

-   System commands you should input into your terminal are given as:

\$ cmd options for command

-   Specific MELT Runtimes (Programs within the MELT.jar) are
    highlighted in **bold red** when discussed in the text (They will
    not be highlighted in command line examples)**.**
    
-   Commands for running MELT are given as:

		java -jar MELT.jar MELTRuntime <options>
		
-   Important information is marked by ***Note:***

-   [References](#_References) are denoted by a number in parentheses: (1)

<span id="_MELT_Quick-Start_Guide" class="anchor"></span>2. MELT Quick-Start Guide
=========================

This guide will take you through the complete detection of MEIs using
MELT on a single human genome. For further instructions on how to run
multiple genomes at a time, or use any other MELT functionality, please
see the appropriate sections as outlined in the table of contents.

<span id="_Install_MELT" class="anchor"></span>2.1 Install MELT
----------------

***Note:*** Further installation instructions can be found in section
‘[4](#_Installing_MELT)’

Simply unpack the .tar.gz file using a command like:

\$ tar zxf MELTvX.X.tar.gz

This will create a directory (./MELTvX.X/) in your current folder. In
this directory will be the MELT.jar and several other folders and files
for MEI discovery.

To run MELT, java1.8 (or greater) and Bowtie2 are required to be
installed on your local system. To check if java1.8 is installed, simply
type:

\$ java –version

If the command returns ‘java version “1.8.\#\_\#” (where the ‘8’ can
also be any number greater than 8, and \# can be any integer), java is
properly installed. If not, please refer to the [Java website](http://docs.oracle.com/javase/8/docs/technotes/guides/install/install_overview.html) for
instructions for java installation for your respective platform.

To check if bowtie2 is installed, simply type:

\$ bowtie2 --version

If the command returns information from bowtie2, then bowtie2 is
properly installed.

<span id="_Preparing_for_MEI_Discovery" class="anchor"></span>2.2 Preparing for MEI Discovery
-------------------------------

To perform MEI discovery, MELT needs the reference used to perform
alignment of sequencing traces (from BWA mem or aln) and several MEI
references. In the folder ‘./MELTvX.X/me\_refs/1KGP\_Hg19/’ there are 3
such references for version Hg19 of the human reference. If using Hg19,
Create a list of these files (with full path) in a text file like so:

\$ ls /full/path/to/MELTvX.X/me\_refs/1KGP\_Hg19/\*.zip &gt;
mei\_list.txt

If not using Hg19, please refer to the section on creating MEI reference
files: [5. Transposon Reference Files](#_Transposon_Reference_Files).

<span id="_Running_MELT_Quick" class="anchor"></span>2.3 Running MELT Quick
----------------

To run MELT, we are going to use a ‘BAM’ file sequenced as part of the 1KGP. This BAM was aligned using BWA aln to the Hg19 version of the human reference. You can download the BAM file and index using either a command like wget:

***Note:*** Downloading with wget may take upwards of an hour, if possible it is recommended to use aspera or some other alternative available to you.

\$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low\_coverage.20121211.bam
\$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low\_coverage.20121211.bam.bai

or, if you have aspera, like:

\$ ascp -i \</Path/to/\>asperaweb\_id\_dsa.putty -T -l300M -Q -k1
fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low\_coverage.20121211.bam .

\$ ascp -i \</Path/to/\>asperaweb\_id\_dsa.putty -T -l300M -Q -k1 fasp-g1k@fasp.1000genomes.ebi.ac.uk:vol1/ftp/phase3/data/NA12878/alignment/NA12878.mapped.ILLUMINA.bwa.CEU.low\_coverage.20121211.bam.bai .

We also need the reference genome. Please download the reference file like so:

$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2\_reference\_assembly\_sequence/hs37d5.fa.gz

$ wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2\_reference\_assembly\_sequence/hs37d5.fa.gz.fai

unzip the main reference file to an easy to access location using:

$ gunzip hs37d5.fa.gz

and rename the index file like:

$ mv hs37d5.fa.gz.fai hs37d5.fa.fai

The following is the basic MELT command line for running a single
genome. (For the arguments to -h and -bamfile, set to the location where you downloaded the files from above. Please see [6.5. Running MELT Using
MELT-SINGLE](#_Running_MELT_Using_MELT-SINGLE) for complete instructions):

	java –jar ./MELTvX.X/MELT.jar Single \
		-a \
		–b hs37d5/NC_007605 \
		–c 8 \
		–h /path/to/hs37d5.fa \
		–bamfile /path/to/NA12878.mapped.ILLUMINA.bwa.CEU.low_coverage.20121211.bam \
		–n ./MELTvX.X/add_bed_files/1KGP_Hg19/hg19.genes.bed \
		–t mei_list.txt \
		–w /path/to/working\_dir/

This will perform MELT discovery on the test bam file, and should take
approximately 10-30 minutes depending on system resources. Output will be in the top level of the directory provided to '-w' in the form of three vcf files: LINE1.final\_comp.vcf, ALU.final\_comp.vcf, and SVA.final\_comp.vcf - one file for each MEI.

<span id="_MELT_System_Requirements" class="anchor"></span>3. MELT System Requirements
=========================================================================================================================

<span id="_Development" class="anchor"></span>3.1. Development
==============================================================

MELT is coded in Java, specifically Java SE Runtime Environment build
1.8.0_40-b26 (Java 1.8). MELT was developed using a machine with
the Red Hat Enterprise Linux Server release 5.9 OS with a 64-bit
(x86\_64) processor architecture with the following system specs:

-   8-Core Intel Xeon CPU with 3.00GHz/core

-   16GB RAM

***Note:*** Previous version of MELT used java1.7 or later, please update to java1.8 to continue to use MELT.

<span id="_Software_Testing" class="anchor"></span>3.2. Software Testing
-------------------------------------------------------------------------------------------------------------------

All three MELT workflows have been tested on the following linux setup:

-   Local Machine: Red Hat Enterprise Linux Server release 5.9

-   Java VM: Java HotSpot 64-Bit Server VM (build 20.6-b01, mixed mode)

-   Grid Architecture: Sun Grid Engine (GE vers. 6.2u5)

MELT-SPLIT and MELT-SINGLE workflows have been tested on the following
Macintosh setup:

-   Local Machine: Macintosh OSX Yosemite version 10.10.2

-   Java VM: Java HotSpot 64-Bit Server VM (build 25.0-b70, mixed mode)

MELT was coded in Java for its high level of portability. It is expected
that as long as Java 1.8 is provided as the JRE, MELT should function on
any standard \*nix machine.

***Note:*** Please do not attempt to run MELT on Windows PCs, as the
file structure and environment assumptions are much different then those
on \*nix based systems.

<span id="_Hardware_Requirements" class="anchor"></span>3.3. Hardware Requirements
--------------------------

If running locally, MELT minimally requires 6GB of system memory, and at
least a single x64 3.00 GHz CPU (This is simply for MELT, and does not
include system overhead). In addition, MELT requires \~20GB of TEMPORARY
storage space per 60X genome in order to run. This does not represent
final storage size for MELT data.

When running MELT using **SGE**, additional pieces of information should
be taken into account:

1.  Number of nodes.

2.  Number of memory on each node.

As an example, each MELT job submitted to the grid uses a maximum of 6GB
of memory. If the grid is made up of 4 nodes with 12GB of memory each,
then it is advised to limit the total number of jobs to 7 jobs at a time
(12 x 4 = 48, 32 GB of total memory, 7 jobs = 42 GB used by MELT with
room for overhead). MELT will communicate to the grid the amount of
memory needed for each job and throttle jobs based on provided user
parameters. SGE additionally accounts for this in job submission, but it
is advised to run MELT with an eye towards grid submission stability.
And as always, please take into account your fellow grid users!

<span id="_Installing_MELT" class="anchor"></span>4. Installing MELT
======================================================================================================================

<span id="_MELT_Installation" class="anchor"></span>4.1. MELT Installation
----------------------

MELT installation is fairly straightforward. First, extract the MELT
tarball:

\$ tar zxf MELTvX.X.tar.gz

This will create a MELTvX.X directory in your current directory
(./MELTvX.X/). This folder will contain the MELT.jar jar file. MELT runs
from a single .jar file, with the first argument representing the
specific runtime that you want to execute. For example:

	java –jar MELT.jar SGE <options>;

Would run the MELT-SGE runtime. Possible MELT Runtimes are described in
short below (or by providing --help/-help/-h to MELT.jar), and in more
detail further on in the README:

**Table 1: MELT Runtimes**

| Runtime Name             | Purpose                                            |
|:------------------------ |:-------------------------------------------------- |
| **BuildTransposonZIP**   | Build custom TE discovery repositories             |
| **Preprocess**           |Processes BAM Alignments for input into MELT        |
|  **Single**              | Run MELT analysis on a Single .bam file            |
|  **SGE**                 | SGE Executable                                     |
|  **IndivAnalysis**       | Step 1 of Non-SGE MELT Pipeline                    |
|  **GroupAnalysis**       | Step 2 of Non-SGE MELT Pipeline                    |
|  **Genotype**            | Step 3 of Non-SGE MELT Pipeline                    |
|  **MakeVCF**             | Step 4 of Non-SGE MELT Pipeline                    |
|  **Source**              | Create an input file for transduction discovery    |
|  **TransductionFind**    | Process .bam files for transductions               |
|  **TransductionMerge**   | Merge transductions discovered by TransductionFind |
|  **Deletion-Genotype**   | MELT Reference MEI Deletion Genotyper              |
|  **Deletion-Merge**      | MELT Reference MEI Deletion VCF writer             |
|  **CALU**                | CAlu Alu classification software                   |
|  **LINEU**               | LineU L1 classification software                   |

You can also get help for individual runtimes by doing:

	java –jar MELT.jar Runtime --help/-help/-h

There are also several other folders included in the ./MELTvX.X/
directory that will be explained at length later in this tutorial.

Next, decide whether or not you will be running MELT using SGE or not.
If running with SGE, please read the section [4.2. Installing MELT Using
SGE](#_Installing_MELT_Using_SGE), otherwise, you can skip ahead to
[4.3. Required Software dependencies](#_Required_Software_Dependencies).

<span id="_Installing_MELT_Using_SGE" class="anchor"><span id="_Installing_MELT_Using_1" class="anchor"><span id="_Toc314153960" class="anchor"></span></span></span>4.2. Installing MELT Using SGE
-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

To run on SGE, MELT requires the SGE library dir for your
grid system. You should be able to find this using the following
command:

\$ echo \$SGE\_ROOT

This will list a directory, and the lib directory will be within this.
As an example, \$ echo SGE\_ROOT could return
/usr/local/packages/sge-root/ and the required lib dir for MELT will be
/usr/local/packages/sge-root/lib/lx24-amd64/. There are typically
several folders in the lib directory – please check with your System
Admin to determine which one is correct (lx24-amd64 is for the 64-bit
AMD architecture present at University of Maryland, Baltimore).

<span id="_Required_Software_Dependencies" class="anchor"></span>4.3. Required Software Dependencies
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

The only external dependency of MELT is Bowtie2 (MELT was tested with
version 2.0.5 & 2.2.4. No differences were seen between versions during
tests). Bowtie must be included on your system path (and Grid path if
using SGE). To test this, you can simply try typing:

\$ bowtie2 --version

This should show the version number and options for bowtie2 in your
terminal window.

<span id="_Required_External_File_Dependencies" class="anchor"></span>4.4. Required External File Dependencies
--------------------------------------------------------------------------------------------------------------

Several files are required for MELT to run:

1.  Reference in FASTA format with .fai index file (users can make .fai from
    fasta using samtools faidx). This can be recovered from the BAM file
    header by using:

 \$ samtools view –H bam\_to\_be\_analyzed.sorted.bam

 It is listed, typically, on each @SQ line in the SAM header under the
 flag AS (sequence name) or UR (URL download location)

 ***Note:*** Must be the same reference as used to align reads! **This
 is extremely important!** Not following this rule can result in any
 number of unexpected errors! These may include increased number of
 false positive/false negative sites, as well as unrecoverable runtime
 errors.

2.  Gene Annotation in BED format. The human and chimpanzee versions of
    this file are included in the MELTvX.X.tar.gz (in /add\_bed\_files/), and must be 	sorted according to coordinate. This file can be downloaded
    from the UCSC Genome Browser website(5,6). Simple download the GENE
    track from the TABLES section of the website. The following table
    outlines the format of this file:

**Table 2: Gene File Format**

| Column       | Description |
|:------------ |:----------------------------------------------------------------------- |
| 1            | chromosome (in either chrN or N format, where N is the chromosome ID)   |
| 2            | Gene start position |
| 3            | Gene stop position |
| 4            | Gene ID (Any format) |
| 5            | Quality Score (default 0, no need to set, MELT does not use) |
| 6            | Strand |
| 7            | Coding start position |
| 8            | Coding stop position |
| 9            | Quality Score (default 0, no need to set, MELT does not use) |
| 10           | Number of exons |
| 11           | Exon lengths (comma delimited) |
| 12           | Exon starts (comma delimited) |

 ***Note:*** Please ensure the chromosome names match those in the
 reference .fasta file you plan to use for BWA alignment and MELT
 analysis. A common example is the ‘chr’ prefix in front of chromosome
 names.

3.  MEI .zip files which include the following files (A complete
    description is below in the section [5. Transposon Reference
    Files](#_Transposon_Reference_Files)):

	a.  FASTA sequence for the ME
	
	b.  Reference insertions in BED format (for filtering known hits)

<span id="_Optional_External_Files" class="anchor"></span>4.5. Optional External Files
----------------------------

MELT can also use a collection of MEI sites, or priors, from other
studies. This is supplied to MELT in the form of a standard site VCF
file. All columns past 9 are not used, and can be excluded when
providing priors. It is recommended to use MELT discovered sites for
prior input, as several pieces of information in the MELT VCF output
(see [7. MELT Output](#_MELT_Output) for more information) can help with
improved genotyping. MELT uses priors in two ways:

1.  Allow for genotyping of sites not discovered in the current
    data set.

2.  Provide more accurate breakpoints for use during genotyping. This is
    only possible if MELT formatted MEIs are included as priors.

As is outlined above, priors are used mainly to increase total number of
genotyped sites, and increase accuracy of that genotyping. Input of
priors depends on the MELT workflow used:

-   If using MELT-SGE or MELT-SINGLE, prior files are provided in the –t
    input file as a second column. Please see the respective workflow
    for a working example.

-   If using MELT-SPLIT, use the –v option while running
    **GroupAnalysis**.

***Note:*** Adding priors will likely increase the runtime of MELT, as
all new sites must be genotyped in new samples.

After initial discovery, discovered sites are checked against the list
of priors to see if any of the prior sites have been rediscovered.
Matches are determined based on a +/- 100 bp window around the prior
site. If a prior is matched as a hit, MELT then checks whether or not
the provided prior site is a MELT formatted VCF. If using a MELT
formatted VCF as prior input, the ASSESS score of the discovered site is
checked against the ASSESS score in the prior file. If the ASSESS score
of the prior is greater than the ASSESS score of the discovered site,
the position of the prior is used for further analysis. Otherwise the
discovered site position is used. If not using a MELT formatted VCF, the
name (Column 3) of the prior is simply cataloged for record keeping.

A collection of VCF priors for L1 (3,048 sites), *Alu* (12,748 sites),
and SVA (835 sites) are provided in the folder ./MELTvX.X/prior\_files/.
These sites were discovered as part of the 1000 Genomes Project, Phase
III. Please see the [How to Cite MELT](#_How_to_Cite_MELT) section for more
information on these sites, and how to cite when using them.

Please see the section [7. MELT Output](#_MELT_Output) for information on
how priors are reported.

<span id="_Transposon_Reference_Files" class="anchor"></span>5. Transposon Reference Files
===============================================================================================

<span id="_General_Description" class="anchor"></span>5.1. General Description
-------------------------

MELT requires transposon .zip files to direct MEI discovery. Such files
are included within the MELT tar-ball for the 1KGP versions of Hg19 and
Hg38, and are located in the ./MELTvX.X/me\_refs/ folder after
extraction. The provided .zip files listed above do not require further
modification in order to be used for MEI discovery. So, if you only want
to discover human transposons in the Hg19 or Hg38 reference
sequences or chimp transposons in the panTro4 reference, you can stop reading here and jump to the [next section](#_Running_MELT). If you want to discover non-reference
transposons in another species or reference, please read on.

***Note:*** If not using the aforementioned references **YOU MUST**
rebuild the transposon reference files. Chromosome names can vary from
reference to reference, and can cause errors when running MELT.

<span id="_Creating_Custom_MEI.zip_Files" class="anchor"></span>5.2. Creating Custom MEI.zip Files
----------------------------------

Advanced users can also create their own MEI reference files for use
with MELT (Tested MEI refs that will run out of the box are included in
the MELT tarball in the ./MELTvX.X/me\_refs/ folder). A rundown on how
to make your own .zip using the **BuildTransposonZIP** runtime is given
step by step below:

### <span id="_The_MELT_.zip_File_Structure" class="anchor"></span>5.2.1. The MELT .zip File Structure

Each MEI\_MELT.zip repository is simply a collection of 5 different
files/indices. This section is simply for reference, MELT handles the
creation of files other than the .fa consensus and the .bed mask. Please
note – file suffixes are crucial to creating the .zip file, do not
change them (i.e. .fa versus .fasta):

**Table 3: mei.zip file structure**

|**File**                            |**Short Description** |
|:---------------------------------- |:------------------------------------------------------------ |
|  Fasta sequence (.fa)              | Sequence of transposon of interest in fasta format |
|  Fasta index (.fa.fai)             | .fai index of the Fasta sequence created by samtools faidx |
|  Reference Transposon Mask (.bed)  | Bed file of all reference insertions in the genome of interest. |
|  Bowtie2 index (.bt2)              | Bowtie2 index of .fa file. Created using bowtie2-build |
|  Info file (.info)                 | Allows MELT to direct specificity and naming of transposon discovery |

### <span id="_Deciding_on_an_MEI_Consensus" class="anchor"></span>5.2.2. Deciding on an MEI Consensus

Creating a consensus sequence for the transposon you want to discover is
very important to accurate MEI discovery. The simplest way to get an
accurate consensus sequence is to simply look at the literature for the
element of interest. If information on the element is scarce on
non-existent, then the following steps may help to generate an accurate
consensus for MEI discovery.

Simply using the repbase(7) provided sequence may not be enough to
discover that particular transposon if it is highly variable from
insertion to insertion. This is typically true of the smaller SINE
elements from genome to genome, but typically not true for the larger
LINE elements. Using a combination of the MELT-Deletion caller and
RepeatMasker, one can evaluate the genome of interest for transposons
that are polymorphic in the reference genome (i.e. when sequencing
another individual, these transposons are ‘deleted’).

1.  Run the [MELT-Deletion](#_MELT-Deletion) caller on the
    pre-existing RepeatMasker masked sites in the genome.

2.  Extract the fasta sequence of around 100 polymorphic sites

3.  Perform a Multiple Sequence Alignment (MSA) using one of several
    available software packages(8,9)

4.  Generate a consensus from this MSA using a tool such as cons from
    the EMBOSS package(10).

5.  Remove any Ns from the consensus (these typically represent InDels
    in some, but not all, extracted elements), as they will hamper MELT
    MEI discovery.

6.  Use this generated consensus when running
    [BuildTransposonZIP](#_MELT-BuildTransposonZIP).

### <span id="_MELT-BuildTransposonZIP" class="anchor"></span>5.2.3. MELT-BuildTransposonZIP

General Usage:

	java –Xmx1G –jar MELT.jar BuildTransposonZIP \
		Transposon_Sequence.fa|Transposon_Sequence.fasta \
		Transposon_Mask.bed \
		NAME[A-Z,a-z,0-9] \
		ERROR[0-9+]

File suffixes are required by **MakeTransposonZIP** to match the
requirements outlined in the general usage above. Naming requirements
and a more detailed explanation of each element of the command line is
explained below:

#### *Fasta* File

Create a fasta file of the transposon reference sequence. The first line
of the file can be anything as long as it matches
‘&gt;\[Non-Whitespace\]’. No spaces allowed. This name will NOT be used
at all, and will be replaced by the NAME input. Following this, either
multiple or single lines of sequence can be placed below. For example,
Alu is 283 bp long; you can either place these full 283 bp on one line,
or split them so there are 60 (fairly standard fasta convention)
characters on each line. MELT does not accept ‘N’ for alignment (ATCG
**only**). ‘N’ will confuse the alignment algorithms during discovery.
Characters can be either upper OR lower case. Only include one sequence
in the fasta file you want to analyze, otherwise **MakeTransposonZIP**
will fail.

#### *BED File*

Create a repeat masked bed file for your element of interest. Download
the repeat masker(11) track (in BED format) from the UCSC Genome Browser
webpage(5,6):

1.  From the homepage go to Tables section (top of webpage)

2.  Choose your organism

3.  Select Repeats from the ‘groups’ dropdown box

4.  Select Repeatmasker from the ‘track’ dropdown box

5.  Choose BED (browser extensible data) from the ‘output format’
    dropdown box

6.  Name the file and click ‘get output’

***Note:*** MELT currently does not support discovery of Transposons
that have not been masked in their respective organism.

***Note:*** Please ensure the chromosome names match those in the
Reference .fasta file you plan to use for BWA alignment and MELT
analysis. A common example is the ‘chr’ prefix in front of chromosome
names.

Next, extract all of the reference coordinates of the type of element
that you will be looking for (an Alu Human example):

\$ grep Alu hg19.repeatmasker.bed &gt; ALU.bed

This will grab ALL reference coordinates for your MEI of interest. You
can also create a custom mask (i.e. masking coordinates you do not wish
to discover in) if you like, just ensure that the first three columns of
the file are chr, start, end in that order. This will prevent discovery
in and around the bed coordinates provided.

#### *Name*

Decide on a name for your element. This name will be a standard feature
of each of the files that **BuildTransposonZIP** will create, and will
be kept across all file types.

***Note:*** This will affect the ability of MELT to classify Family and
Subfamily for Alu and L1 (see [9.1. MEI Species
Classification](#_MEI_Species_Classification)), and twin priming for L1 if not
named properly. If you want to discover primate lineage L1 or Alu and
get all the features of MELT, pleasure ensure ALU or LINE1 is used for
name.

#### *Error*

ERROR tells MELT how many mismatches to allow per 100 bases of the MEI
reference during alignment. This is a specificity parameter, and can
only bet set here (not during MELT MEI Discovery). This number should
reflect the biology of the transposon. If we were performing this for
L1, we would set this number at 3, since L1 elements change at a much
slower rate. *Alu* elements, on the other hand, mutate at a high rate,
will be more divergent from element to element and we set ERROR to 10.
Higher numbers are less specific, lower are more. The ERROR number in
the .zip files included with MELT was empirically tested using the 1000
Genomes Phase III data.

**Table 4: ERROR input for human transposons**

|**Transposon**    |**Error Rate** |
|:---------------- |:---------------- |
| LINE-1           | 3 |
|  *Alu*           | 10 |
|  SVA             | 10 |

***Note:*** **BuildTransposonZIP** will rename the fasta and bed file to
\*.fa.old or \*.bed.old if they share the same name as **BuildTransposonZIP** intermediary files.

<span id="_MELT_Options" class="anchor"></span>6. MELT Options
=============================================================================================================

To ensure ease of use, MELT uses a standard set of options for all runtimes shown in Figure 1. Please refer to the below table for all MELT options unless otherwise specified. Not all options are required for each MELT runtime, this is simply a table of all possible options. For the options specific to a certain MELT runtime, please run:

	java –jar MELT.jar <runtime>

Please see the specific MELT runtime you are interested in using for further information/examples for these options. Required options have no default, and ***must*** be set by the user prior to running MELT. Discretionary options have a default value, and is shown like: **\[false\]** 

-   The table below uses the following format, and contains options for all MELT
    **Runtimes**:

| Option   | Required?    | Boolean   | Description                     |
|:-------- |:-----------  |:--------- |:---------------------------     |
| a        | YES          | NO        | A Required option               |
| b        | NO           | YES       | A boolean optional option       |
| c        | NO           | NO        | A non-boolean optional option   |

**Note**: Do not include any parameters for Boolean options (i.e. do
not do ‘–a false’ if reads were aligned with bwa-mem, just leave –a out
of the command line; Applicable for options without &lt;arg&gt; after
them). Only include in the command line if they are true.

<span id="_Table5" class="anchor"></span>**Table 5: MELT Options**

| Option | Required? | Boolean? | Description | 
|:----------- |:----------- |:----------- |:--------------------------------------------------------- |
| bamfile | YES | NO | .bam alignment to perform MEI discovery on. This is the .sorted.bam or .bam file, NOT the .disc, .disc.bai, or .fq file from MELT-Preprocess.jar. PreProcessing is assumed to have been run already. |
| bamlist | YES | NO | List of .bam alignments with one file per line. This is the .sorted.bam or .bam file, NOT the .disc, .disc.bai, or .fq file from MELT-Preprocess.jar. Also provide an additional column for coverage in the format “file&lt;TAB&gt;coverage”. Coverage is required. |
| bed | YES | NO | A bed file of reference MEs. |
| discoverydir | YES | NO | Directory used for performing initial MEI discovery. Please see [7.4. Running MELT Using MELT-SPLIT](#_Running_MELT_Using_MELT-SPLIT) for further information. |
| exec | YES | NO | Location of the MELT.jar file. |
| genotypingdir | YES | NO | Directory containing samples genotyped by MELT. These files are located in the top level of the dir provided by –w to MELT-Genotype-jar, and are suffixed with &lt;MEI\_EXAMINED&gt;.tsv. Please see [7.4. Running MELT Using MELT-SPLIT](#_Running_MELT_Using_MELT-SPLIT) for more information. |
| h | YES | NO | FULL PATH to the reference fasta. Has to be the same as the reference used to align samples being used for MEI detection, and must have a .fai index at the same location. |
| i | YES | NO | The sge directory, see [4.2 Installing MELT Using SGE](#_Installing_MELT_Using_SGE) for further explanation on what this is. |
| mergelist | YES | NO | List of output files to merge into the MELT-Deletion VCF. |
| n | YES | NO | The gene annotation for the FASTA reference (-h). Must be sorted by coordinate. An annotation for Hg19 is included in MELTvX.X.tar.gz (hg19.genes.bed). See [4.4. Required External File Dependencies](#_Required_External_File_Dependencies) for further information. |
| p | YES | NO | Full path to working directory from GroupAnalysis (Step 2). |
| source | YES | NO | A bed file of all possible source elements to examine for source-offspring relationships. See [MELT Transduction](#_Transduction) for more information. |
| t | YES | NO | List of MEI.zip files with one file per line. See above section on [Transposon Reference Files](#_Transposon_Reference_Files) for more information. |
| u | YES | NO | Additional options required by the SGE runtime at your institution. For example, run MELT using –u ‘–q &lt;que name&gt; -P &lt;project-name&gt; -V’. This should have the same syntax as when you run a standard grid job using qsub (exclude the qsub part in the MELT command!). **Always** provide the qsub flag ‘–V’, and **do not** include anything involving the qsub ‘-b’ option. |
| vcf | YES | NO | The LINE-1 VCF file produced by all MELT for determining source-offspring relationships. See [MELT Transduction](#_Transduction) for more information. |
| w | YES | NO | The working dir for this project. |
| a | NO | YES | Boolean flag to tell MELT that alignments were performed with BWA-ALN rather than BWA-MEM. **\[false\]** |
| ac | NO | YES | Filter out allele count 0 sites when making final MELT VCF. **\[false\]** |
| b | NO | NO | Exclusion list for chromosomes. By default MELT performs MEI discovery on all chromosomes greater than 1Mb. To exclude chromosomes greater than 1Mb, simply list them, delimited by a forward slash. For example, to exclude chromosomes 1,2, and 4, use –b 1/2/4. **\[null\]** |
| bowtie | NO | NO | Path to the bowtie2 execultable if not in PATH. **\[null\]** |
| c | NO | NO | Coverage level of supplied bam file. Can either be set using -c or will be calculated by MELT using the first contig in the reference .fai file starting at the position provided to -d (default position 1,000,000). **Note:** Coverage is important for proper discovery. If using a tool other than MELT, please determine approximate (+/- 2x) coverage prior to analysis. **\[null\]** |
| samtools | NO | NO | Path to the samtools execultable if not in PATH. **\[null\]** ***Note:*** Required if using CRAM format files. |
| cov | NO | NO | Standard Deviation cutoff for distance between left and right sides of MEI evidence. **\[35\]** ***Note:*** The default for –cov 35 was used for all MEI detection in the main MELT paper except for the ancient (Denisovan & Neanderthal) discovery, where it was set to 10. The default should be appropriate for all standard sequencing library preparations. |
| d | NO | NO | Minimum length of chromosome to analyze. **\[1000000\]** |
| e | NO | NO | Fragment length of the samples. **\[500\]** |
| exome | NO | YES | Is analysis being performed on exome samples? **\[false\]** |
| j | NO | NO | No call filter cutoff. Please see [8. MELT Output](#_MELT_Output) for a further explanation of the VCF output provided by MELT. **\[25\]** |
| k | NO | YES | .bam file has already ben preprocessed. See above for description of what this means. **\[false\]** |
| m | NO | NO | Number of jobs to be submitted to the grid at a time. **\[100\]** |
| mcmq | NO | NO | Allow MELT to use MC/MQ tags at greater than X percentage prevalance in the original bam file. **\[95.0\]** |
| minlen | NO | NO | Minimum length of an L1 to classify as full length. **\[5900\]** |
| nocleanup | NO | YES | Do not cleanup intermediate files. Can lead to a sizeable increase on your file-system, but useful for debugging. **\[false\]** |
| o | NO | NO | Output directory for VCF Files. **\[./\]** |
| q | NO | YES | Alignments are phred+64 quality encoding. **\[false\]** |
| r | NO | NO | Approximate read length of the supplied bam files. **\[100\]** |
| s | NO | NO | STDEV cutoff. Please see [8. MELT Output](#_MELT_Output) for a further explanation of the VCF output provided by MELT. **\[2.0\]** |
| sr | NO | NO | Sites with < X total split reads will be filtered from all downstream MELT analysis. Default X, 0, is to not filter any sites. **\[0\]** |
| v | NO | NO | Priors VCF for current MEI. **\[null\]** |
| verbose | NO | YES | Report all SGE output and error. Will be deposited in –w/output/ and –w/error/, respectively. This will result in a number of small text files equal to the number of .bam files provided in –bamlist \* \~6. **\[false\]** |
| x | NO | NO | Number of time for any job to run on SGE, in minutes. Set higher for bigger bam files. **\[240\]** |
| z | NO | NO | Maximum reads to load into memory at any given time. Best to leave alone, but if looking at a .bam file larger than \~200Gb, set moderately higher. **\[5000\]** |


<span id="_Running_MELT" class="anchor"></span>7. Running MELT
=============================================================================================================

As stated before, MELT can be run in one of three ways depending on the
number of .bam files you wish to analyze. If running MELT on multiple
.bam files (Figure 1 – in Red), it is possible to run with SGE
(MELT-SGE) or without SGE (MELT-SPLIT). If discovery is to be performed
on a single .bam file (Figure 1 – In Green), MELT-SINGLE is the preferred
method. An additional method is availble for genotyping reference MEs (Figure 1 - In Blue), and is discussed in more detail in [9.2 MELT-Deletion](#_MELT-Deletion).

**Figure 1: MELT MEI Discovery Workflows**

![MELT Workflows](https://github.com/eugenegardner/MELT/blob/master/SuppFigure2.png)

Each Runtime works with a slightly different set of options, but all
option letters carry over between .jar files (i.e. –h is
always for the reference fasta). Please see [6. MELT Options](#_MELT_Options) for more details, and all options discussed below are also available, albeit less verbosely, via command line:

\$ ./MELT.jar &lt;Runtime&gt; --help/-help/-h

MELT will warn you when required options are not properly provided.

<span id="_Preprocessing_.bam_Files_for_MELT" class="anchor"></span>7.1. Preprocessing .bam Files for MELT
--------------------------------------

Beyond the .jar files discussed below, you can additionally
‘pre-process’ bam files to speed up MELT’s runtime. This is a required
step for MELT-SPLIT and MELT-SGE, but is optional for MELT-SINGLE (you
can use the flag ‘-k’ to tell MELT that .bam files have already been
processed). To process a .bam file using **Preprocess**, simply use the
following syntax:

	java –Xmx2G –jar MELT.jar Preprocess -bamfile name.sorted.bam -h reference.fa;

.bam files are expected to be both sorted and indexed, with at least the
.bam suffix (.sorted.bam and .bam will both work). The .bai suffix is
expected to be AFTER the .bam suffix, not in place of it (i.e.
.sorted.bam.bai, NOT .sorted.bai). MELT will not run if this convention
is not used. **Preprocess** generates three files, all as suffixes to
the current bam file:

1.  sorted.bam.disc – discordant pairs from the current BAM file

2.  sorted.bam.disc.bai – index of these discordant pairs

3.  sorted.bam.disc.fq – fastq version of all discordant pairs for MEI
    alignment

These files will be in the same directory as the provided .bam file, and
can be used as input for further MEI discovery using MELT.

MELT also can make use of additional tags provided by `samtools fixmate` 
or similar postprocessing algorithms to substantially improve runtime, 
in particular when providing .cram format files as input. These tags are: 

* MC - Mate cigar string
* MQ - Mate mapping quality 

The fixmate postprocessing step is often part of the workflows of large
sequencing centers and comes as standard for bams produced by several
consortia (1KGP high coverage samples included). It is highly recommended
that, especially if using cram, you take advantage of these improvements 
as this can improve MELT IndivAnlysis runtime by a factor of 10x. The percentage
of sites which have MC/MQ tags required to trigger this option (default 95%) 
can be adjusted using the `-mcmq` flag. This should not be required but is
included for rare cases where it might be necessary.


<span id="_Which_Version_of_MELT_Should_I_Use" class="anchor"></span>7.2. Which Version of MELT Should I Use?
----------------------------------------

This depends largely on whether SGE is installed at your home
institution and the number of .bam files you wish to analyze. It also
depends on how much interactivity/input you would like to have running
MELT. MELT-SGE is simply a SGE scheduling wrapper, and MELT-SGE,
MELT-SPLIT, and MELT-SINGLE are identical programmatically. The only
major differences between them being:

1.  MELT-SGE handles running each of the individual steps

2.  MELT-SGE & MELT-SPLIT can run multiple MEIs at once

3.  MELT-SGE throttles the number of concurrent jobs on a computational
    grid using input from the user

4.  MELT-SGE is less interactive, and the only input from the user is at
    the beginning of a run, while MELT-SGE requires human input at
    each step.

5.  MELT-SGE requires SGE, while MELT-SPLIT can be run on any
    grid architecture.

6.  MELT-SINGLE can only run on a single .bam file.

You can also use MELT-SPLIT and MELT-SINGLE if you have SGE, you will
just have to submit and monitor the jobs manually. That means ensuring
that each step completes properly, and that you are not overloading
whatever system architecture/computer you are currently using. To use
MELT-SGE, please read the section [Running MELT using
MELT-SGE](#_Running_MELT_Using_MELT-SGE). To use MELT-SPLIT, please read the
section [Running MELT using MELT-SPLIT](#_Running_MELT_Using_MELT-SPLIT). To use
MELT-SINGLE, please read the section [Running MELT using
MELT-SINGLE](#_Running_MELT_Using_MELT-SINGLE).

<span id="_Running_MELT_Using_MELT-SGE" class="anchor"></span>7.3. Running MELT Using MELT-SGE
------------------------------------------------------------------------------------------------

### <span id="_Description_SGE" class="anchor"></span>7.3.1. Description

Perform MEI discovery using Sun Grid Engine (SGE). This will look in
multiple genomes (-bamlist) for multiple transposons (-t).

General usage:

	java –Xmx1g -jar MELT.jar SGE <options>;

### <span id="_Options_SGE" class="anchor"></span>7.3.2. Options

***Note:*** Due to the nature of SGE, it is recommended to provide FULL
PATH for all options.

***Note:*** While –b is not a required option, it is crucial to make
sure the genome you are using does not have any ‘decoy’ or ‘fake’
chromosomes. An example is the hs37d5 ‘decoy’ chromosome attached to the
1000 Genomes human genome fasta file. This ‘chromosome’ works to help
BWA to align reads properly to the reference(12), and if included in
MELT analysis will result in increased runtime and discovery of MEIs on
this fake chromosome.

An Example file\_list.txt file for –bamlist (with column 1 being full path to
the .bam file):

/full/path/to/alignment1.bam 
/full/path/to/alignment2.bam
.  
.  
.  
/full/path/to/alignmentN.bam

***Note:*** Optionally, coverage can be provided in this file, but can be if you want to use coverage calculated via an external source, like below:

/full/path/to/alignment1.bam\<TAB\>15  
/full/path/to/alignment2.bam\<TAB\>32  
.  
.  
.  
/full/path/to/alignmentN.bam\<TAB\>17

An example transposon\_file\_list.txt file for –t:

/full/path/to/LINE\_MELT.zip  
/full/path/to/ALU\_MELT.zip  
/full/path/to/SVA\_MELT.zip  

OR if utilizing priors with –t:

/full/path/to/LINE\_MELT.zip\<TAB\>/full/path/to/LINE1\_priors.VCF  
/full/path/to/ALU\_MELT.zip\<TAB\>/full/path/to/ALU\_priors.VCF  
/full/path/to/SVA\_MELT.zip\<TAB\>/full/path/to/SVA\_priors.VCF  

***Note:*** Priors do not have to be provided for each MEI class, but
please ensure that the priors provided are only specific to the MEI type
being discovered (i.e. only provide L1 priors on the LINE\_MELT.zip
line).

For detailed option descriptions, please see [6. MELT Options](#_MELT_Options).

Output will be found in the directory –w/ in the form of a VCF file for each MEI. For example, -w/ALU.final\_comp.vcf if performing ALU discovery

<span id="_Running_MELT_Using_MELT-SPLIT" class="anchor"></span>7.4. Running MELT Using MELT-SPLIT
----------------------------------

### <span id="_Description_Split" class="anchor"></span>7.4.1. Description

Perform MEI discovery without using Sun Grid Engine. Rather than one
runtime, MELT-SPLIT uses 4 different runtimes to perform analysis. This
is to ensure that the user gets the same population level MEI discovery
as MELT-SGE, but without having to use SGE as a parallel computing
platform. The pipeline follows the following 4 general steps:

I.  **IndivAnalysis** – MEI discovery in individual samples

II. **GroupAnalysis** – Merge discovery information across all genomes
    in project to determine accurate breakpoint information. Also
    determines several MEI specific pieces of information such as
    length, strand, subfamily, target site duplication (TSD),
    among others.

III. **Genotype** – Genotype all samples using merged MEI discovery
    information**.**

IV. **MakeVCF** – Perform final filtering and merging of individual
    samples into final VCF.

For more information on how each individual step is performed, please
see our publication (see [Citing MELT](#_How_to_Cite_MELT) below). For information on each option, please see For detailed options, please see [6. MELT Options](#_MELT_Options). For actual information on running each individual step, please continue following
the README.

### 7.4.2. <span id="_IndivAnalysis" class="anchor"></span>IndivAnalysis

General Usage:

	java –Xmx6G –jar MELT.jar IndivAnalysis <options>;

If performing an analysis across multiple .bam files, it is recommended
to use the same working directory for each step. For example:

	java –Xmx2G –jar MELT.jar IndivAnalysis –w /path/to/MEI_NAME/ -bamfile /path/to/file1.bam <other_options>;

	java –Xmx2G –jar MELT.jar IndivAnalysis –w /path/to/MEI_NAME/ -bamfile /path/to/file2.bam <other_options>;

	java –Xmx2G –jar MELT.jar IndivAnalysis –w /path/to/MEI_NAME/ -bamfile /path/to/file3.bam <other_options>;

Where –w is the working directory, and –bamfile is the .bam file of interest (described further in [Table 5](#_Table5)). This is due to the adjusted workflow of MELT-Split:

**IndivAnalysis** creates several files as output, and **GroupAnalysis**
looks for them all in one folder on initialization.

***Note:*** **IndivAnalysis** assumes [MELT Preprocess](#_Preprocessing_.bam_Files_for_MELT) has been run on the supplied .bam file (-bamfile).

***Note:*** It is not recommended to use the same working directory for multiple
transposons. i.e. for LINE1 use LINE1, for Alu use ALU. This is to
ensure GroupAnalysis can retrieve the proper files.

***Note:*** While –b is not a required option, it is crucial to make
sure the genome you are using does not have any ‘decoy’ or ‘fake’
chromosomes. An example is the hs37d5 ‘decoy’ chromosome attached to the
1000 Genomes hg19 fasta file. This ‘chromosome’ works to help
BWA to align reads properly to the reference(12), and if included in
MELT analysis will result in increased runtime and discovery of MEIs on
this fake chromosome.

### 7.4.3. <span id="_GroupAnalysis" class="anchor"></span>GroupAnalysis

General Usage:

	java –Xmx4g –jar MELT.jar GroupAnalysis <options>;

The key output from this process is an
&lt;MEI\_EXAMINED&gt;.pre\_geno.tsv file. This file is used in the input
of the final two steps.

The same rules used for **IndivAnalysis** apply here. i.e. Ensure you
only analyze one transposon (-t) at a time, and use the same working
directory (-w) for all samples.

### 7.4.4. <span id="_Genotype" class="anchor"></span>Genotype

General Usage:

	java –Xmx2G –jar MELT.jar Genotype <options>;

Using the pre\_geno.tsv file from the GroupAnalysis step, perform genotyping across each sample of interest.

### 7.4.5. <span id="_MakeVCF" class="anchor"></span>MakeVCF

General Usage:

	java –Xmx2G –jar MELT.jar MakeVCF <options>;

This step, beyond merging all MEI calls, also performs final filtering.
Please see [MELT Output](#_MELT_Output) for information on filtering
parameters.

### 7.4.6. <span id="_Example_MELT-SPLIT_Workflow" class="anchor"></span>Example MELT-SPLIT Workflow

As a simulated example to demonstrate MELT-SPLIT based MEI discovery,
suppose we have sequenced four human individuals (Person1.sorted.bam,
Person2.sorted.bam, Person3.sorted.bam, Person4.sorted.bam) to
approximately 30x coverage. We then used bwa-mem to align them to the
Hg19 human reference sequence. We now want to know if they have any
LINE1 non-reference Mobile Element Insertions. We first want to run
Preprocess on all four samples to ensure MELT has the proper information
to discover MEIs:

	java –Xmx2G –jar MELT.jar Preprocess /path/to/Person1.sorted.bam

	java –Xmx2G –jar MELT.jar Preprocess /path/to/Person2.sorted.bam

	java –Xmx2G –jar MELT.jar Preprocess /path/to/Person3.sorted.bam

	java –Xmx2G –jar MELT.jar Preprocess /path/to/Person4.sorted.bam

Next, we want to do initial MEI site discovery in all 4 genomes (only
required options shown):

	java –Xmx6G –jar MELT.jar IndivAnalysis \
		–bamfile /path/to/Person1.sorted.bam \
		–w ./LINE1DISCOVERY/ \
		-t /path/to/LINE_MELT.zip \
		–h /path/to/human_reference.fa;

	java –Xmx6G –jar MELT.jar IndivAnalysis \
		-bamfile /path/to/Person2.sorted.bam \
		–w ./LINE1DISCOVERY/ \
		-t /path/to/LINE_MELT.zip \
		–h /path/to/human_reference.fa;
		
	java –Xmx6G –jar MELT.jar IndivAnalysis \
		-bamfile /path/to/Person3.sorted.bam \
		–w ./LINE1DISCOVERY/ \
		-t /path/to/LINE_MELT.zip \
		–h /path/to/human_reference.fa;

	java –Xmx6G –jar MELT.jar IndivAnalysis \
		-bamfile /path/to/Person4.sorted.bam \
		–w ./LINE1DISCOVERY/ \
		-t /path/to/LINE_MELT.zip \
		–h /path/to/human_reference.fa;

Next, we want to combine the initial discovery across these four genomes
into a single piece of information to aid genotyping and filtering of
false positive hits:

	java –Xmx4G –jar MELT.jar GroupAnalysis \
		-discoverydir ./LINE1DISCOVERY/ \
		-w ./LINE1DISCOVERY/ \
		-t /path/to/LINE_MELT.zip \
		-h /path/to/human_reference.fa \
		-n /path/to/human_annotation.bed;

Next, we will genotype each of the samples using the merged MEI
information file:

	java –Xmx2G –jar MELT.jar Genotype \
		-bamfile /path/to/Person1.sorted.bam \
		-t /path/to/LINE_MELT.zip \
		-h /path/to/human_reference.fa \
		-w ./LINE1DISCOVERY/ \
		-p ./LINE1DISCOVERY/;

	java –Xmx2G –jar MELT.jar Genotype \
		-bamfile /path/to/Person2.sorted.bam \
		-t /path/to/LINE_MELT.zip \
		-h /path/to/human_reference.fa \
		-w ./LINE1DISCOVERY/ \
		-p ./LINE1DISCOVERY/;

	java –Xmx2G –jar MELT.jar Genotype \
		-bamfile /path/to/Person3.sorted.bam \
		-t /path/to/LINE_MELT.zip \
		-h /path/to/human_reference.fa \
		-w ./LINE1DISCOVERY/ \
		-p ./LINE1DISCOVERY/;

	java –Xmx2G –jar MELT.jar Genotype \
		-bamfile /path/to/Person4.sorted.bam \
		-t /path/to/LINE_MELT.zip \
		-h /path/to/human\_reference.fa \
		-w ./LINE1DISCOVERY/ \
		-p ./LINE1DISCOVERY/;

Finally, we provide the genotyping directory from Genotype (as -genotypingdir) to MakeVCF to generate the final VCF file:

	java –Xmx2G –jar MELT.jar MakeVCF \
		-genotypingdir ./LINE1DISCOVERY/ \
		-h /path/to/human_reference.fa \
		-t /path/to/LINE_MELT.zip \
		-w ./LINE1DISCOVERY/ \
		-p ./LINE1DISCOVERY/LINE1.pre_geno.tsv;

This will result in a final VCF file in the directory provided to –o, or
by default in your current directory: ./LINE1.final\_comp.vcf.

<span id="_Running_MELT_Using_MELT-SINGLE" class="anchor"></span>7.5. Running MELT Using MELT-SINGLE
---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

### 7.5.1. <span id="_Description_Single" class="anchor"></span>Description

MELT-SINGLE is a method to call MEIs on a single genome without having
to manually perform each step individually using MELT-SPLIT. This method
is advantageous for users who only have to analyze one genome for MEIs.
While MELT runs fine on low coverage data, it is recommended to use a set
of priors (included in -t; see below) to decrease the false negative rate associated with
low coverage sequencing.
 
General Usage:

	java –Xmx6G –jar MELT.jar Single <options>;

### 7.5.2. <span id="_Options_Single" class="anchor"></span>Options

***Note:*** MELT-SINGLE makes use of the –k flag. If you have already
run Preprocess on your .bam alignment, ensure you use –k or MELT-SINGLE
will re-perform this step.

***Note:*** While –b is not a required option, it is crucial to make
sure the genome you are using does not have any ‘decoy’ or ‘fake’
chromosomes. An example is the hs37d5 ‘decoy’ chromosome attached to the
1000 Genomes human genome fasta file. This ‘chromosome’ works to help
BWA to align reads properly to the reference(12), and if included in
MELT analysis will result in increased runtime and discovery of MEIs on
this fake chromosome.

An example transposon\_file\_list.txt file for –t:

/full/path/to/LINE\_MELT.zip  
/full/path/to/ALU\_MELT.zip  
/full/path/to/SVA\_MELT.zip  

OR if utilizing priors with –t:

/full/path/to/LINE\_MELT.zip\<TAB\>/full/path/to/LINE1\_priors.VCF  
/full/path/to/ALU\_MELT.zip\<TAB\>/full/path/to/ALU\_priors.VCF  
/full/path/to/SVA\_MELT.zip\<TAB\>/full/path/to/SVA\_priors.VCF  

For detailed option descriptions, please see [6. MELT Options](#_MELT_Options).

***Note:*** Priors do not have to be provided for each MEI class, but
please ensure that the priors provided are only specific to the MEI type
being discovered (i.e. only provide L1 priors on the LINE\_MELT.zip
line).

<span id="_MELT_Output" class="anchor"></span>8. MELT Output
==============================================================================================================================================================================================================

MELT provides one output VCF file for each MEI analyzed, and when run
with MELT-SGE, a BED file for each MEI for each sample. See VCF format
specifications for further details on default columns (columns 1-6, 9+).
MELT, by default, does not provide any information for the NAME (Column
3) column. However, if a prior file is provided, the NAME column will be
the prior name.

In column 7, MELT provides one of six FILTERS – PASS, rSD, s25, hDP, lc, or ac0:

1.  (-j - MakeVCF) No-call filter \[s25\]. If sites are no-call (./.) for greater
    than, default, 25% of individuals, then the call will be flagged
    with the r25 flag in the FILTER column.

2.  (-s - MakeVCF) 5’ and 3’ evidence filter \[rSD\]. If sites have too much
    evidence on one side (only evidence on 3’ or 5’ of predicted
    insertion site), default 2.0 standard deviations outside of the mean
    for this discovery project, then the call will be flagged with the
    rSD flag in the FILTER column. Will not filter a site with fewer
    than 10 reads of combined evidence. For a more detailed explanation
    of this, please see our paper (see [Citing MELT](#_How_to_Cite_MELT) below).

3.  Split discordant filter \[hDP\]. Flagged when the (total number of split discordant pairs / total number of standard discordant pairs) * 100 >= 55.

4.  Low complexity filter \[lc\]. Filters sites that are within +/- 25bps of a low complexity region. Currently in testing.

5.  Allele count 0 filter \[ac0\]. Filters sites without a genotyped allele.

If the site passes all of these filters PASS will be listed in column
7.

MELT also provides several custom fields in the VCF INFO column (column
8) to provide additional information to the user.

**Table 6: MELT INFO Fields**

| INFO       | Description |
|:---------- |:--------------- |
| TSD        | Provides the sequence of the Target Site Duplication (TSD) for each insertion. Will be ‘null’ if was not able to be determined. |
| ASSESS     | Provides an accuracy assessment, and the information used, to determine the breakpoint of each insertion. Please read our paper to learn more above individual scores (see [Citing MELT](#_How_to_Cite_MELT)). |
| INTERNAL   | Provides information on if the insertion is internal to a gene in the provided reference annotation file. The first piece of information is the gene name; the second is the location within that gene (either exon, intron, 5\_UTR, or 3\_UTR). |
| SVTYPE     | Gives basic information about what type of MEI is represented in the record. |
| SVLEN      | Length of this MEI. |
| MEINFO     | Provides estimate information on the start, stop, and length of the MEI against its own reference. If polarity is unknown, will be ‘null’. If either start or stop is not know, will be -1. |
| DIFF       | Provides differences for this element compared to the MEI reference provided. Please see [MEI Species Classification](#_MEI_Species_Classification) below for more information. |
| LP         | Number of read pairs supporting the left, 5’, side of this insertion. |
| RP         | Number of read pairs supporting the right, 3’, side of this insertion. |
| RA         | Ratio of left and right read pairs. |
| ISTP       | If the SVTYPE info field is L1 or LINE1, the site of twin priming in BP. If insertion is not twin-primed, will be 0. If no measurement was made, will be -1. |
| PRIOR      | If a site was found in priors file, will be true. If site was not found in priors file, or no priors file provided, will be false. |

An additional description if each field is provided in the MELT VCF header.

<span id="_Other_MELT_Tools" class="anchor"></span>9. Other MELT Tools
===================

<span id="_MEI_Species_Classification" class="anchor"></span>9.1. MEI Species Classification
----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

Built into MELT is a sequence assembly algorithm that rebuilds the
sequence of the MEI from discordant read pairs discovered by MELT, and
then classifies the insertion according to known family relationships.
This algorithm currently contains three parts:

1.  Sequence assembler – Assembles discordant pairs from WGS alignments
    into a complete sequence for an MEI.

2.  Difference Discovery – Looks for differences in assembled sequence
    according to the provided MEI reference.

3.  Classifier – Determines family and sub-family information about
    specific transposon types.

Part one only runs during MELT analysis, and currently is not
stand-alone.

For Parts two and three, analysis depends on the transposon type. MELT
can currently only classify Alu and LINE-1 insertions into family and
subfamily. If, during a MELT Run, the name provided in the .info file
matches either ‘Alu’, ‘LINE1’, or ‘L1’ (with any combination of lower
and uppercase) family classification will be performed. If not, MELT
will stop after part two as outlined above and provide a list of differences
from the MEI .fa file.

### 9.1.1 <span id="_CAlu" class="anchor"></span>CAlu

CAlu is a standalone tool that, given an unknown Alu sequence, will
classify it according to family and mutations. Sites are classified according
to several publications, and are included in the [MELT Manuscript](#_How_to_Cite_MELT).

General Usage (When using multi-fasta):

	java –Xmx2g –jar MELT.jar CALU –f /path/to/sequences.fasta

General Usage (When providing sequence from the command line):

	java –Xmx2g –jar MELT.jar CALU –s <Input sequence [ATCG]>;

-o can additionally be specified to direct output to a specific file.
Otherwise, will print to standard out.

CAlu output is a tab-delimited and contains the following information:

**Table 7: CAlu/LineU Output**

| Column       | Description |
|:------------ |:---------------------------------------------------------------------------------- |
| 1            | Name of input sequence |
| 2            | Family |
| 3            | Subfamily |
| 4            | Start |
| 5            | Stop|
| 6            | Number of diagnostic changes (used to determine classification in columns 2 & 3) |
| 7            | Number of overall changes |
| 8            | Changes, delimited by ‘\|’ |
| 9            | Original input sequence |

### 9.1.2 <span id="_LineU" class="anchor"></span>LineU

LineU is a standalone tool for LINE-1 subfamily analysis that will
perform a similar analysis to CalU. I/O is identical for both.

<span id="_MELT-Deletion" class="anchor"></span>9.2. MELT-Deletion
----------------------------------------------------------------------------------------------------------------------

MELT-Deletion is a tool to determine presence/absence of reference
mobile elements (those included in the reference release for a given
species). It is included in the ./MELTvX.X/execs/ folder.

MELT-Deletion uses the repeat masker track (the same file used as a mask
for MELT MEI discovery analysis) to determine if a reference MEI is
present or absent. Please refer to the detailed description of .zip
creation in the [MELT-BuildTransposonZIP](#_MELT-BuildTransposonZIP) section for specifics on how to generate this file. MELT-Deletion uses two methods depending on the length of the reference MEI. If the insertion is longer than the fragment read length,
MELT-Deletion uses a joint discordant and split read pair analysis to
estimate the genotype. If the insertion is shorter than the fragment
read length, MELT-Deletion only uses a split read pair genotyping
approach.

MELT-Deletion is run in two steps:

1.  **Genotype** – Analyze an individual .bam file for reference
    MEI deletions.

2.  **Merge** – Merge and filter individual deletion calls from multiple
    .bam files.

Each of these two processes is included in the MELT.jar file, and is
specified with either Deletion-Genotype or Deletion-Merge. MELT-Deletion
is first run on multiple .bam files, and then the resulting output files
are merged. For specific options, please see the MELTOptions Section above. For example:

	java –Xmx2g –jar MELT.jar Deletion-Genotype -bamfile file1.sorted.bam –w /path/to/working/dir/ -bed /path/to/coords.bed –h /path/to/ref.fasta <other options>;

	java –Xmx2g –jar MELT.jar Deletion-Genotype -bamfile file2.sorted.bam –w /path/to/working/dir/ -bed /path/to/coords.bed –h /path/to/ref.fasta <other options>;

	java –Xmx2g –jar MELT.jar Deletion-Genotype -bamfile file3.sorted.bam –w /path/to/working/dir/ -bed /path/to/coords.bed –h /path/to/ref.fasta <other options>;

\$ ls \*.del.tsv > list.of.outputs.txt

	java -jar -Xmx2G Deletion-Merge –mergelist list.of.outputs.txt -bed /path/to/coords.bed –h /path/to/ref.fasta <other options>;

Deletion-Genotype will create a single output tsv file, using the name
of the bam file. The input of file1.sorted.bam would result in the
output file1.del.tsv being placed in the working directory specified by
-w.

***Note:*** Do not change the suffix (.del.tsv), as this is required by
Deletion-Merge.

***Note:*** The default minimum chromosome size to analyze is 1Mb
(1000000) but can be set lower if necessary. However, it is **not**
recommended to set it below about 100Kb (100000), as this may cause
errors during merging.

### 9.2.1. <span id="_MELT_Deletion_Output" class="anchor"></span>MELT-Deletion Output

MELT-Deletion outputs in VCF 4.1 format using similar conventions to
regular MELT MEI discovery analysis. The primary difference lies in the
INFO column.

**Table 8: MELT-Deletion INFO Fields**

| INFO       | Description |
|:---------- |:------------ |
| SVTYPE     | Transposon subfamily according to Repeat Masker. |
| SVLEN      | Length of this MEI. |
| ADJLEFT    | If breakpoint determined by MELT-Deletion is different from Repeat Masker on the 5’ end, list adjusted breakpoint here. |
| ADJRIGHT   | If breakpoint determined by MELT-Deletion is different from Repeat Masker on the 3’ end, list adjusted breakpoint here. |

<span id="_Transduction" class="anchor"></span>9.3. Transduction
----------------------------------------------------------------------------------------------------------------------

Transduction is a tool for the discovery of 3' transductions in MELT discovery MEI's. These transductions can be used to definitively determine the source element from which an offspring LINE-1 is derived(13). MELT performs this analysis in 3 steps:

1. Generate a source list using **Source**.
2. Preprocess .bam files using this list of known possible source elements using **TransductionFind**.
3. Merge predictions, filter, and add additional information into the MELT VCF file for elements with a known source using **TransductionMerge**.

Please see [MELT Options](#_MELT_Options) for further description of options provided to the steps in this pipeline.

### 9.3.1. <span id="_Source" class="anchor"></span>Source

**Source** takes two inputs, a LINE1.final\_comp.vcf file generated by MELT (-vcf), and a bed file of *known* possible source elements (-bed). The file provided to -vcf can be any final LINE1 discovery output from any of the MELT processes listed above ([MELT-SGE](#_Running_MELT_Using_MELT-SGE), [MELT-SPLIT](#_Running_MELT_Using_MELT-SPLIT), or [MELT-SINGLE](#_Running_MELT_Using_MELT-SINGLE)). The file for -bed is a list of known Full-Length LINE1 elements in the reference genome. There is a sample file of 298 such elements from Hg19 in the MELT repository located at ./MELTvX.X/add_bed_files/1KGP_Hg19/Hg19.FL.bed.

These two files allow for discovery of transductions sourced from both non-reference (-vcf) and reference (-bed) elements. For a more detailed description and example of information provided by this kind of study, please see the [MELT Publication](#_How_to_Cite_MELT). The output list of compiled 'source' elements is provided as a supplement to the -bed input. For example:

	java -jar MELT.jar Source -vcf LINE1.final_comp.vcf -bed Hg19.FL.bed
	
Would create a source file called Hg19.FL.bed.source in the same directory as Hg19.FL.bed. This final is then provided as input to the next steps in the transduction discovery pipeline.

***Note:*** The default for -minlen is 5900, which is a rough approximation of the minimum nucleotide length for a full-length LINE-1 element. This may not be the case in other species, and may be below the your desired estimated accuracy. Please use this parameter to adjust the elements you wish to extract from both -bed and -vcf.

### 9.3.2. <span id="_TransductionFind" class="anchor"></span>TransductionFind

TransductionFind preprocesses bam files using the [Source File](#_Source) generated in the previous step. Simply run as follows (using the file generated by **Source** as input):

	java -Xmx2G -jar MELT.jar TransductionFind -bamfile <path/to/file.sorted.bam> -h <path/to/reference.fa> -source <path/to/Hg19.FL.bed.source>

This will generate a file with the suffix '.trans' at the location of the bam file given to -bamfile (i.e. running TransductionFind with file.sorted.bam as input would result in the file file.sorted.bam.trans).

***Note:*** You must run TransductionFind on all samples that are included in the vcf file that you wish to perform transduction discovery on in the [TransductionMerge](#_TransductionMerge) step. MELT will provide an error if this is not the case.

### 9.3.3. <span id="_TransductionMerge" class="anchor"></span>TransductionMerge

TransductionMerge takes as it's input a list of bam files (-bamlist) found in a given vcf file (-vcf). The format for the -bamlist file is identical to that used for MELT-SGE. Please see [MELT Options](#_MELT_Options) for a more detailed description. Do not provide a list of the *.trans files generated by TransductionFind, this will result in MELT failing. Simply provide a tab delimited list of bam files and the coverage for that bam file:

/path/to/file1.sorted.bam\<TAB\>10  
/path/to/file2.sorted.bam\<TAB\>9  
/path/to/file3.sorted.bam\<TAB\>11  

Where each file (*.sorted.bam) has been preprocessed by TransductionFind and has a .trans file at the same location. Then, provide this file as input for TransductionMerge:

	java -Xmx2G -jar MELT.jar TransductionMerge -bamlist <path/to/filelist.txt> -h <path/to/reference.fa> -source <path/to/Hg19.FL.bed.source> -vcf <path/to/LINE1.final_comp.vcf>
	
This will generate a new VCF file with the suffix .trans (e.g. LINE1.final_comp.trans.vcf). This will not overwrite your previously generated VCF file. Included in this file are two additional INFO fields for transduction information:

**Table 9: Transduction INFO Fields**

| INFO       | Description |
|:---------- |:------------ |
| METRANS    | Mobile element transduction info of the form CHR,START,END,POLARITY. *null* if no transduction. |
| MESOURCE   | Source element information if METRANS not equal to null in the form of SOURCECHR,SOURCEPOS,SOURCETYPE,GENOTYPE CONCORDANCE,SUPPORTING READS,TRANSDUCTION SIZE. *null* if no transduction. |

These two fields provide different pieces of information:

- METRANS provides the MELT ascertained transduction coordinates adjacent to the source element from which the offspring LINE-1 was derived.
- MESOURCE provides information on the source element, and the amount of evidence supporting the confidence of the transduction, and therby, the source-offspring relationship:
	- SOURCECHR and SOURCEPOS are the coordinates of the source element that gave rise to the offspring insertion
	- SOURCETYPE is either REF (for a source element found in the reference) or NONREF (for a source element found in the vcf file provided by -vcf).
	- GENOTYPE CONCORDANCE is given by the formula:  
	(% of individuals with element that have transduction support) * (% of individuals with transduction support that have element)
	Only sites with greater than 0.25 concordance are reported.
	- SUPPORTING READS is the total number of reads supporting a given transduction. Only transductions with greater than 5 reads supporting are reported.
	- TRANSDUCTION SIZE is the *adjusted* size of the transduction based on read support and location of the reference element. For example, read support for a transduction may begin internally to a source element, and then continue into uniquely mappable sequence. The TRANSDUCTION SIZE reflects this and is essentially the length of the transduction that *does not* contain LINE-1 sequence.

For additional information on the ascertainment of source-offspring relationships, please see the [MELT Publication](#_How_to_Cite_MELT).

	
<span id="_What_MELT_is_for" class="anchor"></span>10. What MELT is for
===================

MELT is an accurate MEI caller with FDR established by a blind study
performed by a very subjective group separate from our own. FDR
estimates for the MEI types already discovered and studied by the 1000
Genomes project are &lt;5%, with MELT designed to be just as sensitive
at discovering low allele frequency insertions, as high allele frequency
insertions. The sequencing depth of the samples that are provided are
the only limits to MELT. MELT is also highly optimized for ease of use,
and due to the portability of Java, users should not experience issues
when moving to other platforms.

MELT is also a suite of tools for performing various analyses of ascertained MEIs.
It can determing subfamily of human elements, source-offspring relationships, 5'
inversions of LINE-1 elements, and much more.

<span id="_What_MELT_is_not_for" class="anchor"></span>11. What MELT is not for
========================

1.  MELT is NOT repeat masker! It will not discover MEIs in the
    reference sequence from your recently sequenced genome. It either
    discovers polymorphic MEIs in large populations or genotypes
    reference MEIs as compared to the reference.

2.  MELT does very well on most transposons. However, MELT has only been
    tested for non-LTR transposons (*Alu*, L1, and SVA in humans) and
    not LTR transposons such as ERVs or DNA transposons such as Mariner.
    Due to the nature of HERV-K’s LTRs having high sequence homology
    with SVA, MELT can sometimes miss-call an SVA insertion as a
    HERV-K insertion. Software development is underway to ensure that
    these transposons reach the high standard that we have already set
    for non-LTR transposons. Feel free to try these types of MEIs, but
    results may vary. A HERV-K.zip file has been included in the MELT
    release, but a FDR for discovery has not been established.

3.  MELT has currently been tested only in human (hg19), chimpanzee
    (panTro4), and dog (canFam3.1). MELT should do an excellent job at
    discovering transposons in any sequenced organism, but tests to
    determine this have yet to be performed.

4.  MELT has only been tested on BWA-MEM and BWA-ALN alignments. MELT
    does not discriminate between these types of alignments, though
    performing BWA-MEM alignments is typically easier and
    less time-intensive. It is recommended that you use these types of
    alignments, as MELT has not been tested on other types
    of alignments. More information on performing BWA alignments is
    provided at the BWA website (<http://bio-bwa.sourceforge.net/>), and
    in the BWA publication (1).
    
5.  MELT is designed only for Illumina paired-end data. Please do not use it on 454, Sanger traces, Pacific Biosciences SMRT sequencing, or other non-paired end sequencing.

<span id="_Reporting_Issues_With_MELT" class="anchor"></span>12. Reporting Issues With MELT
==============================

Please first ensure you are using the most recent release of MELT. If issues still occur, please refer to the MELT website (<http://melt.igs.umaryland.edu>) and use the provided link for assistance. This is the prefered method of contact, and some issues may have already been addressed by other users.

<span id="_How_to_Cite_MELT" class="anchor"></span>13. How to Cite MELT
=========================================================================

If you use MELT to discovery MEIs, or use the VCF discovery sets provided as part of the MELT paper, please cite:

Gardner, E. J., Lam, V. K., Harris, D. N., Chuang, N. T., Scott, E. C., Mills, R. E., Pittard, W. S., 1000 Genomes Project Consortium & Devine, S. E. The Mobile Element Locator Tool (MELT): Population-scale mobile element discovery and biology. *Genome Research*, 2017. **27**(11): p. 1916-1929.

If you use the included 1KGP phase III MEI priors’ list, please cite:

1000 Genomes Project Consortium, A global reference for human genetic variation. *Nature*, 2015. **526**(7571): p. 68-74.

Sudmant, P. H., Rausch, T., Gardner, E. J., *et al.*, An integrated map of structural variation in 2,504 human genomes. *Nature*, 2015. **526**(7571): p. 75-81.

<span id="_References" class="anchor"></span>14. References
==============

(1) Li, H. & Durbin, R. Fast and accurate long-read alignment with
Burrows-Wheeler transform. *Bioinformatics (Oxford, England)* **26**,
589-595, doi:10.1093/bioinformatics/btp698 (2010).

(2) Danecek, P. *et al.* The variant call format and VCFtools.
*Bioinformatics (Oxford, England)* **27**, 2156-2158,
doi:10.1093/bioinformatics/btr330 (2011).

(3) Langmead, B. & Salzberg, S. L. Fast gapped-read alignment with Bowtie
2. *Nature methods* **9**, 357-359, doi:10.1038/nmeth.1923 (2012).

(4) Li, H. A statistical framework for SNP calling, mutation discovery,
association mapping and population genetical parameter estimation from
sequencing data. *Bioinformatics (Oxford, England)* **27**, 2987-2993,
doi:10.1093/bioinformatics/btr509 (2011).

(5) Kent, W. J. *et al.* The human genome browser at UCSC. *Genome
research* **12**, 996-1006, doi:10.1101/gr.229102. Article published
online before print in May 2002 (2002).

(6) Karolchik, D. *et al.* The UCSC Table Browser data retrieval tool.
*Nucleic acids research* **32**, D493-496, doi:10.1093/nar/gkh103
(2004).

(7) Jurka, J. *et al.* Repbase Update, a database of eukaryotic repetitive
elements. *Cytogenetic and genome research* **110**, 462-467,
doi:10.1159/000084979 (2005).

(8) Larkin, M. A. *et al.* Clustal W and Clustal X version 2.0.
*Bioinformatics (Oxford, England)* **23**, 2947-2948,
doi:10.1093/bioinformatics/btm404 (2007).

(9) Katoh, K., Kuma, K., Toh, H. & Miyata, T. MAFFT version 5: improvement
in accuracy of multiple sequence alignment. *Nucleic acids research*
**33**, 511-518, doi:10.1093/nar/gki198 (2005).

(10) Rice, P., Longden, I. & Bleasby, A. EMBOSS: the European Molecular
Biology Open Software Suite. *Trends in genetics : TIG* **16**, 276-277
(2000).

(11) Smit AFA, H. R., Green P. RepeatMasker Open-3.0. (1996-2010).

(12) Genovese, G., Handsaker, R. E., Li, H., Kenny, E. E. & McCarroll, S.
A. Mapping the human reference genome's missing sequence by three-way
admixture in Latino genomes. *American journal of human genetics*
**93**, 411-421, doi:10.1016/j.ajhg.2013.07.002 (2013).

(13) Moran, J. V., DeBerardinis, R. J. & Kazazian, H. H., Jr. Exon shuffling by L1 retrotransposition. *Science* **283**, 1530-1534 (1999).
