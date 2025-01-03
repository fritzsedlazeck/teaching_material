# Wonderful world of Structural Variant Calling! 

***
Lead: Fritz J Sedlazeck

Human Genome Sequencing Center,
Baylor College of Medicine

Computer Science Department,
Rice University 
***

## Physalia-Courses
This is part of the Physalia-Courses lectured in Dec 2024
https://www.physalia-courses.org/


## Goals of this module
The goal of this module is to get you familiarized with variant identification (especially Structural Variations) across assembly, short read based mapping and long read based mapping. 
For Structural Variations (SVs) detection we will be using multiple methods (Assemblytics, Dipcall, Manta, Sniffles, etc) and then later compare them across to obtain more insights about advantages and disadvantages across the differnet approaches.

### IMPORTANT NOTES
All the instructions below requires you to think along! You are responsible for a clear data structure on your workspace and to find the individual files that you need. I am happy to assist but please take a moment to look around or think where a certain file could be. 

### The main steps in this Module are:
1. Assembly based SV detection using [Assemblytics](http://assemblytics.com/)
2. Assembly based SV detection using [Dipcall](https://github.com/lh3/dipcall)

Dipcall is already installed!


## Organism
We will utilize data from Cryptosporidium for our exercise. Cryptosporidium is a microscopic parasite that causes the diarrheal disease cryptosporidiosis. Both the parasite and the disease are commonly known as “Crypto.” There are many species of Cryptosporidium that infect animals, some of which also infect humans.
The data is from a recent publication in Gigascience: [Fully resolved assembly of Cryptosporidium parvum](https://doi.org/10.1093/gigascience/giac010). Here we produced short read Illumina data together with Oxford Nanopore long read data to obtain a fully assembled (T2T) assembly and thus improve the reference overall. 
We will study this sample of Cryptosporidium compared to the previously available reference sequence (GCF_000165345.1.). 


## Conda
Keep in mind that we (ie. Luis) has already organized the installation of all the needed programs. For this we are using different enviroments that you can list and activate. For installation of additional software for your home cluster use https://bioconda.github.io/
Here is a link to give you more insights into conda enviroments https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#activating-an-environment 


Please activate your coda enviroment using: 
```bash
source /home/ubuntu/share/conda_init.sh
```
You will need to do this everytime you open a new terminal!

If you don't see (base) next yo your username then type
```bash
conda activate base
```


## Part 1: Assembly based SV detection (Assemblytics)
As discussed in the lecture, assembly based SV detection is quite comprehensive. Nevertheless, before we start, we need to align a new assembly to an existing e.g. reference genome to identify structural differences. For this we will be using MUMmer (nucmer). Alternatively, you could also utilize Dipcall as a program especially if you have a phased assembly at hand. MUMmer is a commonly used package that allows you to rapidly compare two sequences together and includes multiple packages for summary reports over the aligned sequences. This includes variant reporting, coordinate reporting or even the generation of dot plots. As shown in the lecture dot plot is very helpful and intuitive method to compare sequences for us.

To begin to call variants we need to first compare the reference to our assembly. 
Next, we take assembly and reference there: (NOTE! the workshop organizers might have organized the data already for you)
```
# day 1 data
mkdir day1
cd day1
ln -T /home/ubuntu/share/data/day1 -s day1_data

```

Now we need to initiate the conda environment `sr` (short reads)
```bash
conda activate sr
```

Now we can initiate the alignment:
```
nucmer -maxmatch -l 100 -c 500 day1_data/GCF_000165345.1.fa day1_data/assembly/crypto_BCM2021_v2_min100k.fasta --prefix mummer_out
```

We can now switch over to Assemblytics. For simplicity had created a session for us, moreover the web application recently was sutdown, but we have saved a copy for you to see. The command line tool still works. 
```
# assemblytics <delta_file> <output_prefix> <unique_anchor_length> <min_variant_size> <max_variant_size>

assemblytics mummer_out.delta  output_assemblytics 10000  50  10000
```

To convert the Assemblytics file for SV we will need SURVIVOR: We want to discard all SV impacting less than 50 bp. 
```
SURVIVOR convertAssemblytics output_assemblytics.Assemblytics_structural_variants.bed 50 assemblytics.vcf
```

Lets count how many SV we could identify: 
```
grep -vc '#' assemblytics.vcf
```

You can look at all the SVs:
```
grep -v '#'  assemblytics.vcf
```


In case you were not able to run Assemblytics, you can also see that VCF file from here if you had any difficulties:
```
SURVIVOR convertAssemblytics day1_data/assembly/crypto.Assemblytics_structural_variants.bed 50 assemblytics.vcf

grep -vc '#' day1_data/assembly/assemblytics.vcf

grep -v '#'  day1_data/assembly/assemblytics.vcf
```

Congratulations we now have SV calls from Assemblytics. Lets discuss these results together!


## Part 2: Assembly based SV detection (Dipcall)

After Assemblytics we also want to use another method to call SV based on our assembly. Dipcall is quite state of the art method to identify variants between an assembled genome (mainly phased) and a reference genome. Dipcall was implemented by Heng Li and can be found here: [GitHub Dipcall](https://github.com/lh3/dipcall).


Before using Dipcall make yourself familiar with the command line options! Here Dipcall asks us to provide two input file one for each haplotype. Our toy example assembly from Crypto doesnt have two haplotypes resolved so we will just use it for both haplotypes. 

First, we will need to initate our alignment and variant calling step. Note that here Dipcall only writes the commands in one file given your input parameters. 
```
conda activate lr
run-dipcall  -t 2 -a outputdipcall day1_data/GCF_000165345.1.fa day1_data/assembly/crypto_BCM2021_v2_min100k.fasta day1_data/assembly/crypto_BCM2021_v2_min100k.fasta > prefix.mak
```

Thus, next we need to execute this commad to actually run the comparson and obtain the variant calls:
```
make -j2 -f prefix.mak

```




