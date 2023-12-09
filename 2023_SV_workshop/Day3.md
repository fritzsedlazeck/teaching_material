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
For SVs detection in long-read data we will be using Sniffles2. Other methods exist that use both ONT and PacBio data, such as cuteSV.


### IMPORTANT NOTES
All the instructions below require you to think along. You are responsible for a clear data structure on your workspace and to find the individual files that you need. Read alignments will be provided, and in some cases, intermediary files that are produced by Sniffles2, with the aim of speeding the process. We are happy to assist but please take a moment to look around or think where a certain file could be.


### The main steps in this Module are:
1. SV calling using [Sniffles](https://github.com/fritzsedlazeck/Sniffles)
2. SV annotation using [AnnotSV](https://github.com/lgmgeo/AnnotSV)
3. Population level SV calling using Sniffles
4. allele frequency SV annotation using [SVAFotate](https://github.com/fakedrtom/SVAFotate)
5. Low-frequency SV calling using Sniffles

## Organism
We will utilize data from humans for our exercise. Specifically, we are going to use a couple of reference samples from the Genome in a Bottle Consortium (GIAB)
who has the aim to provide the characterization of human genomes for use in benchmarking, including analytical validation and technology development, optimization, and demonstration.
For the first exercise we will use chromosome 22 of the human sample HG002, who is a male with Ashkenazi Jewish ancestry. For exercise three we will also use human samples HG003 and HG004, who are the father and mother of HG002 (respectively). Finally, for exercise five we will use a mix of human sample HG002 and H00733 (1000 genomes project)
You can read more about [GIAB here](https://www.nist.gov/programs-projects/genome-bottle)


## Conda 
```bash
source Share/conda_init.sh
```

## Part 1: Long-read based SV detection: Sniffles
As discussed in the lecture, assembly based SV detection is quite comprehensive. As discussed in the previous lectures, before we start calling SVs, we need to align the reads to a reference genome, in this case the human genome. Different versions of the [human genome exists](https://www.ncbi.nlm.nih.gov/genome/guide/human/). Currently the most used are, GRCh37 (released in ) and GRCh38 (released in ). The latest version, known as CHM13 or T2T, is the most complete to date however it lacks the annotation resources of the former ones.


### 1.1 Get the data
For this exercise we will use the reference genome GRCh37. To begin to call variants we need the reads to be aligned to a reference genome. In the following path you will find the read alignment for Chromosome 14 of human sample HG002.


```bash
mkdir 01_sv_germline
cd 01_sv_germline
# input
ln -T ~/Share/data/day3/01/hg002.bam -s hg002.bam
ln -T ~/Share/data/day3/01/hg002.bam.bai -s hg002.bam.bai
# reference
ln -T ~/Share/data/ref/grch37.fasta.gz -s grch37.fasta.gz
ln -T ~/Share/data/ref/grch37.fasta.gz.fai -s grch37.fasta.gz.fai
ln -T ~/Share/data/ref/grch37.fasta.gz.gzi -s grch37.fasta.gz.gzi
```

### 1.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'long' using the following command

```bash
conda activate long
```

We will use Sniffles to detect the signals we can obtain from abnormally mapped long-reads. To run Sniffles2, execute the following command:
```bash
sniffles \
   --input hg002.bam \
   --vcf hg002_germline.vcf.gz \
   --snf hg002_germline.snf \
   --reference grch37.fasta.gz \
   --threads 1
```

#### 1.2.1 Sniffles2 command explained
Long read alignment file in bam format used as input:
`--input hg002.bam`

Output file in 'Variant Call Format' used to store the Structural Variants that we detect:
`--vcf hg002_germline.vcf.gz`

Reference genome used during the alignment in order to have deletion sequences in the results:
`--reference grch37.fasta.gz`

Number of logical threads used during the analysis (parallelization):
`--threads 1`

You can read more about the specifications [here](https://samtools.github.io/hts-specs/)

#### 1.2.2 Sniffles2 help
If you want to see the help menu from Sniffles you can execute the following command:

```bash
sniffles
```

#### 1.2.3 Sniffles2 parameters detailed


You can see the full list of parameters in detail from Sniffles you can execute the following command:

```bash
sniffles --help
```

#### 1.2.4 Number of variants detected
Now that Sniffles has finalized, let us count how many SV we could identify:
```bash
bcftools view --no-header hg002_germline.vcf.gz | wc -l

bcftools view --no-header hg002_germline.vcf.gz | less -S
```


## Part 2: SV Annotation with AnnotSV
In this part we will use AnnotSV to identify potential implications of the Structural Variants that we detected with Sniffles2. Here, the genome version that was discussed in earlier points is important. We need to use the same genome version in both the read-alignment/SV calling and annotation.


### 2.1 Running AnnotSV

For this exercise we will use the same conda environment as before.

```bash
conda activate long
```

Once loaded, let us run AnnotSV

```bash
~/Share/AnnotSV/bin/AnnotSV \
  -SVinputFile hg002_germline.vcf.gz \
  -outputFile  hg002_germline_annotsv.tsv \
  -genomeBuild GRCh37 \
  -annotationMode full \
  -outputDir .
```

#### 2.1.1 AnnotSV command explained
Structural Variants detected by Sniffles2:
`-SVinputFile hg002_germline.vcf.gz` 

Output file that will contain the annotation:
`-outputFile  hg002_germline_annotsv.tsv`

Genome version to ensure the genomic coordinates of the SVs match the annotation:
`-genomeBuild GRCh37`

Two types of lines are produced by AnnotSV: `-annotationMode split`
* '**full**' length of the SV, every SV are reported, even those not covering a gene
* 'split' by gene, which gives focuses on each gene overlapped by the SV

Output directory: 
`-outputDir .`

#### 2.1.2 AnnotSV help
If you want to see the help menu from AnnotSV you can execute the following command:
```bash
~/Share/AnnotSV/bin/AnnotSV  --help
```

### 2.2 Filtering interesting hits
According to AnnoSV [manual](https://lbgi.fr/AnnotSV/Documentation/README.AnnotSV_latest.pdf)
AnnotSV_ranking_score (field 116): SV ranking score following the 2019 joint consensus
recommendation of ACMG and ClinGen. Scoring: 
* pathogenic if ≥0.99
* likely pathogenic [0.90;0.98]
* variant of uncertain significance [0.89;-0.89]
* likely benign [-0.90;-0.98]
* benign if ≤-0.99

```bash
cat hg002_germline_annotsv.tsv | cut -f 116 | sort | uniq -c

cat hg002_germline_annotsv.tsv | awk -F '\t' '{if($117 >= 0.99) print $_ }' | less -S
```

Based on the American College of Medical Genetics (ACMG_class, field 118):
* class 1 (benign)
* class 2 (likely benign)
* class 3 (variant of unknown significance)
* class 4 (likely pathogenic)
* class 5 (pathogenic)
* class NA (Non Attributed)
```bash
cat hg002_germline_annotsv.tsv | cut -f 118 | sort | uniq -c

cat hg002_germline_annotsv.tsv | awk -F '\t' '{if($118 == 1) print $_ }' | less -S
```


## Part 3: Population level SV calling with Sniffles
Now, we are going to investigate SV not in a single sample, but in a collection of samples. The samples can be related (family) or not (population). Sniffles2 has a built in way to merge samples in order to get a **population level VCF.**

### 3.1 Get the data
For this exercise we will use two related samples, HG002 and HG003. The family relationship is that HG003 is the father of HG002. This family relationship is useful when exploring, for example, inheritance, as some of the variants from HG002 come from HG003.

For this exercise, we have already called SVs with Sniffles2 and produced a special intermediary file (.snf) that contains the needed information to perform the population level SV calling.

```bash
mkdir 03_population
cd 03_population
# input
ln -T ~/Share/data/day3/03/hg002.snf -s hg002.snf
ln -T ~/Share/data/day3/03/hg003.snf -s hg003.snf
# reference
ln -T ~/Share/data/ref/grch37.fasta.gz -s grch37.fasta.gz
ln -T ~/Share/data/ref/grch37.fasta.gz.fai -s grch37.fasta.gz.fai
ln -T ~/Share/data/ref/grch37.fasta.gz.gzi -s grch37.fasta.gz.gzi
```

### 3.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'long' using the following command

```bash
conda activate long
```

To run Sniffles2 population calling, execute the following command:

```bash
sniffles \
   --input hg002.snf hg003.snf \
   --vcf merge.vcf.gz \
   --reference grch37.fasta.gz \
   --threads 1
```

#### 3.2.1 Sniffles2 command explained
Two intermediary files (produced with sniffles, as in the previous example) that contain all the Structural Variant candidates. We will use this file to call SV in both samples at the same time.
`--input hg002.snf hg003.snf`

Output file in 'Variant Call Format' used to store the Structural Variants that we detect:
`--vcf merge.vcf.gz`

Reference genome used during the alignment in order to have deletion sequences in the results:
`--reference grch37.fasta.gz`

Number of logical threads used during the analysis (parallelization):
`--threads 1`

### 3.3 How many variants we found
```bash
bcftools view --no-header merge.vcf.gz | wc -l

bcftools view --no-header merge.vcf.gz | less -S
```


## Part 4: Population level allele frequency SV annotation with SVAFotate
SVAFotate is a tool for annotating structural variant VCFs with population level allele frequency information. This is to be able to distinguish whether a mutation or variant is common or rare in the general population. 

### 4.1 Get the data
```bash
mkdir 04_population_af
cd 04_population_af
# input
ln -T ~/Share/data/day3/04/popAFs_GRCh37.bed.gz -s population_af.bed.gz
```


### 4.2 Running SVAFotate
To run SVAFotate, we will load the conda environment 'svafotate' using the following command

```bash
conda activate svafotate
```

Next lets use the following command to annotate the VCF from exercise 3:

```bash
svafotate annotate \
   --vcf ../03_population/merge.vcf.gz \
   --out merge_af_annotated.vcf.gz \
   --bed population_af.bed.gz
```

#### 4.2.1 SVAFotate command explained
Sub-command from SVAFotate to annotate VCFs: `annotate`

Input SV VCF from a previous analysis: `--vcf ../03_population/merge.vcf.gz`

Output VCF that will contain population AF annotations: `--out merge_af_annotated.vcf.gz`

Population AF: `--bed population_af.bed`

#### 4.2.2 Let's check a couple of SVs
```bash
bcftools view --no-header merge_af_annotated.vcf.gz | less -S


#### 4.3 SVAFotate help
To get the help menu from svafotate we can type:
```bash
svafotate -h
```


## Part 5: Low-frequency SV detection with Sniffles
In part one, we investigated the SV in sample HG002. The variants we called there are **germline variants** which means that they are present in all the cells of the individual from which we collected the sample. However, there are cases in which the variants are present only in a single tissue, or in some cells. This SVs are **low-frequency** because they are not present in all the cells. We also call this low-grequency SV **mosaic** or **somatic** depending on the context.

For htis exercise, we are going to use Sniffles2 to investigate low-frequency types of SVs

### 5.1 Get the data
For this exercise we will use the same file as in part one.

```bash
mkdir 05_sv_mosaic
cd 05_sv_mosaic
# input
ln -T ~/Share/data/day3/05/hg002.bam -s hg002.bam
ln -T ~/Share/data/day3/05/hg002.bam.bai -s hg002.bam.bai
# reference
ln -T ~/Share/data/ref/grch37.fasta.gz -s grch37.fasta.gz
ln -T ~/Share/data/ref/grch37.fasta.gz.fai -s grch37.fasta.gz.fai
ln -T ~/Share/data/ref/grch37.fasta.gz.gzi -s grch37.fasta.gz.gzi
```

### 5.2 Running Sniffles2
To run Sniffles2, we will load the conda environment 'long' using the following command

```bash
conda activate long
```

To run Sniffles2, execute the following command:
```bash
sniffles \
   --input hg002.bam \
   --vcf hg002_mosaic.vcf.gz \
   --mosaic \
   --reference grch37.fasta.gz \
   --threads 1
```

#### 5.2.1 Sniffles2 command explained
Long read alignment file in bam format used as input:
`--input hg002.bam`

Output file in 'Variant Call Format' used to store the Structural Variants that we detect:
`--vcf hg002_mosaic.vcf.gz`

Flag parameter that indicates that we will be calling low frequency Structural VAriants (.05 >= AF >= 0.2):
`--mosaic`

Reference genome used during the alignment in order to have deletion sequences in the results:
`--reference grch37.fasta.gz`

Number of logical threads used during the analysis (parallelization):
`--threads 1`


#### 5.2.2 Number of variants detected
Now that Sniffles has finalized, let us count how many SV we could identify:
```bash
bcftools view --no-header hg002_mosaic.vcf.gz | wc -l

bcftools view --no-header hg002_mosaic.vcf.gz | less -S
```

