# Wonderful world of Structural Variant Calling! 

Lead: Fritz J Sedlazeck

Human Genome Sequencing Center,
Baylor College of Medicine

Computer Science Department,
Rice University 
***

## Day 3

### Goals of this module
The goal of this module is to get you familiarized with Structural Variations (SVs) identification using long-read based mapping. 
For SVs detection in long-read data we will be using Sniffles2. Other methods exists that use both ONT and PacBio data, such as cuteSV.

### IMPORTANT NOTES
All the instructions below requires you to think along. You are responsible for a clear data structure on your workspace and to find the individual files that you need. Read alignments will be provided, and in some cases, intermediary files that are produced by Sniffles2, with the aim of speeding the process. We are happy to assist but please take a moment to look around or think where a certain file could be. 

### The main steps in this Module are:
1. SV calling using [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
2. SV annotation using [AnnotSV](https://github.com/lgmgeo/AnnotSV)
3. Population level SV calling using Sniffles
4. allele frequency SV annotation using [SVAFotate](https://github.com/fakedrtom/SVAFotate)
5. Low-frequency SV calling using Sniffles

## Organism
We will utilize data from human for our exercise. Specifically, we are going to use a couple of reference samples from the Genome in a Bottle Consortium (GIAB)
who has the aim to provide the characterization of human genomes for use in benchmarking, including analytical validation and technology development, optimization, and demonstration.
For the first excercise we will use chromosome 22 of the human sample HG002, who is a male with Ashkenazi Jewish ancestry. For excecise three we will also use human samples HG003 and HG004, who are the father and mother of HG002 (respectively). Finally, for excercise five we will use a mix of human sample HG002 and H00733 (1000 genomes project) 
You can read more about [GIAB here](https://www.nist.gov/programs-projects/genome-bottle)


## Part 1: Long-read based SV detection: Sniffles 
As discussed in the lecture, assembly based SV detection is quite comprehensive. As discussed in the previous lectures, before we start calling SVs, we need to align the reads to a reference genome, in this case the human genome. Different versions of the [human genome exists](https://www.ncbi.nlm.nih.gov/genome/guide/human/). Currently the most used are, GRCh37 (released in ) and GRCh38 (released in ). The latest version, known as CHM13 or T2T, is the most complete to datem however it lacks the annotation resources of the former ones.

### 1.1 Get the data
For this excercise we will use reference genome GRCh37. To begin to call variants we need the reads to be aligned to a reference genome. In the following path you will find the read alignment for Chromosome 22 of human sample HG002.

```bash
TODO
```

### 1.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'SNF2' using the following command

```bash
TODO
```

We will use Sniffles to detect the signals we can obtain from abnormally mapped long-reads. To run Sniffles2, execute the following command:
```bash
TODO
```

#### 1.2.1 Sniffles2 command explained
TODO: explain each parameter from the command

#### 1.2.2 Sniffles2 help
If you want to see the help menu from Sniffles you can execute the following command:
```bash
TODO
```

#### 1.2.3 Sniffles2 parameters detailed

You can see the full list of parameters in detail from Sniffles you can execute the following command:
```bash
TODO
```

#### 1.2.4 Number of variants detected
Now that Sniffles has finalized, let us count how many SV we could identify: 
```
bcftools view --no-header sniffles2_hg002.vcf.gz
```


## Part 2: SV Annotation with AnnotSV
In this part we will use AnnotSV to identify potential implications of the Stuctural Variants that we detected with Sniffles2. Here, the genome version that was discussed in earlier points is important. We need to use the same genome version in both the read-alignment/SV calling and annotation.

### 2.1 Running AnnotSV

First, we need to load the conda environment for AnnotSV to properly work.

```bash
TODO
```

Once loaded, let us run AnnotSV

```bash
TODO
```

#### 2.1.1 AnnotSV command explained
TODO: explain each parameter from the command

#### 2.1.2 AnnotSV help
If you want to see the help menu from Sniffles you can execute the following command:
```bash
```

### 2.2 Filtering interesting hits
TODO:


## Part 3: Population level SV calling with Sniffles
Now, we are going to investigate SV not in a single sample, but in a collection of samples. The samples can be related (family) or not (population). Sniffles2 has a built in way to merge samples in order to get a **population level VCF.**

### 3.1 Get the data
For this excercise we will use two related samples, HG002 and HG003. The family relationship is that HG003 is the father of HG002. This family relationship is useful when exploring, for example, inheritance, as some of the variants from HG002 come from HG003.

For this excercise, we have already call SVs with Sniffles2 and produced a special intermediary file (.snf) that contains the needed information fo perform the population level SV calling. 

```bash
TODO
```


### 3.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'SNF2' using the following command

```bash
TODO
```

To run Sniffles2 population calling, execute the following command:
```bash
TODO
```

#### 3.2.1 Sniffles2 command explained
TODO: explain each parameter from the command



## Part 4: Population level allele frequency SV annotation with SVAFotate

TODO: 

## Part 5: Low-frequency SV detection with Sniffles
In part one, we investigated the SV in sample HG002. The varians we called there are **germline variants** which means that they are present in all the cells of the individual from which we collected the sample. However, there are cases in which the variants are present only in a single tissue, or in some cells. This SVs are **low-frequency** because they are not present in all the cells, TODO: **mosaic** or **somatic**


We are going to use Sniffles2 to investigate these types of SVs

### 5.1 Get the data
For this excercise we will use the same file as in part one.

For this excercise, we have already call SVs with Sniffles2 and produced a special intermediary file (.snf) that contains the needed information fo perform the population level SV calling. 

```bash
TODO
```


### 5.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'SNF2' using the following command

```bash
TODO
```

To run Sniffles2, execute the following command:
```bash
TODO
```

#### 5.2.1 Sniffles2 command explained
TODO: explain each parameter from the command


#### 5.2.2 Number of variants detected
Now that Sniffles has finalized, let us count how many SV we could identify: 
```
bcftools view --no-header sniffles2_hg002_mosaic.vcf.gz
```
