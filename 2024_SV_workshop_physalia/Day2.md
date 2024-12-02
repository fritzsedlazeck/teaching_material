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
For Structural Variations (SVs) detection we will be using multiple methods (Assemblytics, Manta, Sniffles) and then later compare them across to obtain more insights about advantages and disadvantages across the different the differnet approaches.

### IMPORTANT NOTES
All the instructions below requires you to think along. You are responsible for a clear data structure on your workspace and to find the individual files that you need. I am happy to assist but please take a moment to look around or think where a certain file could be. 

### The main steps in this Module are:
1. Short read mapping based SV detection (using [Manta](https://github.com/Illumina/manta))
2. Long read based mapping based SV detection (using [Sniffles](https://github.com/fritzsedlazeck/Sniffles))
3. SV comparison (using [SURVIVOR](https://github.com/fritzsedlazeck/SURVIVOR))

## Organism
We will utilize data from Cryptosporidium for our exercise. Cryptosporidium is a microscopic parasite that causes the diarrheal disease cryptosporidiosis. Both the parasite and the disease are commonly known as “Crypto.” There are many species of Cryptosporidium that infect animals, some of which also infect humans.
The data is from a recent publication in Gigascience: [Fully resolved assembly of Cryptosporidium parvum](https://doi.org/10.1093/gigascience/giac010). Here we produced short read Illumina data together with Oxford Nanopore long read data to obtain a fully assembled (T2T) assembly and thus improve the reference overall. 
We will study this sample of Cryptosporidium compared to the previously available reference sequence (GCF_000165345.1.). 



## Conda 
```bash
source /home/ubuntu/share/conda_init.sh
conda activate sr
```


## Part 2: Illumina based SV detection 
In this part we will use Manta to interpret the signals we can obtain from abnormally mapped paired end Illumina reads. 

Again you might get the file already on your account or you can download from here:
```
# make a directory
mkdir day2
cd day2
ln -T /home/ubuntu/share/data/day2 -s day2_data

```

It contains both the reads and the reference genome


### 1. Initiate the run:
Next we initiate the Manta run:
```
configManta.py --bam=day2_data/short-reads/our_refCrypto_reads.sort.bam --referenceFasta=day2_data/GCF_000165345.1.fa  --runDir=Out_Manta
```
This should just take seconds as it initiates the folder structure and specifies for the subsequent process to use our mapped reads and our reference file. In addition, we specify the output to be written in `Out_Manta`

### 2. Run the analysis:
```
python2 Out_Manta/runWorkflow.py -j 2 -m local -g 6
```

This will launch the Manta pipeline that we previously configured. `-j` specifies the number of CPU threads (2 in our case), `-m` local indicates that it should not try to run things on different nodes or instances and `-g 30` specifies the available memory for the process in GB.

Manta now searches for abnormal paired-end reads and split-reads across our mapped reads. These will be analyzed together and clustered to identify SV in these samples. After ~2-3 minutes you should see that the program has finished.

Our SV calling results can be found here (Out_Manta/results/variants/). Let us open quickly the output of the highest quality SV files:
```
ls  Out_Manta/results/variants/
```
As you can see, we have multiple VCF files. These represent the different stages of Manta and the confidence level for the SV calls. We typically use the `diploidSV.vcf.gz` file.


Let us open quickly the output of the highest quality SV files:
```
gzip -dc Out_Manta/results/variants/diploidSV.vcf.gz > illumina.vcf
less illumina.vcf
```

You can get the file here if you had difficulties:
```
ln -T day2_data/short-reads/illumina.vcf -s illumina_results.vcf
```

Lets count how many SV we could identify: 
```
grep -vc '#' illumina.vcf

# or if you linked the results file
grep -vc '#' illumina_results.vcf
```

Let us all discuss the different variant types and how they are reported! 


## Part 3: Long read based SV detection 
Finally we are ready for the Oxford Nanopore detection using sniffles. For this use the "ont_mapped.sort.bam" file that I have previously mapped using minimap2. 
You might have that file on your account, but if not you can download it here:
```
conda activate lr
```

Using Sniffles v2 this should be a simple command like:

```
sniffles -i day2_data/ont-reads/ont_prev.sort.bam -v sniffles.vcf
```

To check the number of SVs detected by Sniffles:
```
grep -vc '#' sniffles.vcf
```

How many SV did you detect? 


You can also check the file from here if you had issues:
```
less -S day2_data/ont-reads/sniffles.vcf
```

## Part 4: Structural Variant comparison

Now that we generated Assembly, Illumina  and ONT based SV calls it is time to compare them. One tool that you can use for this very easily is SURVIVOR. SURVIVOR is a very simple method to compare SV but also includes multiple other methods that might be useful. Also feel free to check out Truvari, which is one of our newest methods to compare SV. When comparing SV it is quite important to take multiple things into considerations. Feel free to ask if you want to discuss certain points more. 

For SURVIVOR we want to use the merge option. Before doing this, the merge option requires a file including all paths and VCF files that you want to compare. Thus, we generate the file like this:
```
ls sniffles.vcf > vcf_files
ls illumina.vcf >> vcf_files
ls day2_data/assembly/assemblytics.vcf >> vcf_files
```

Next we can initiate the compare with 100bp wobble and requiring that we are only merging with SV type agreement. Furthermore, we will only take variants into account with 50bp+. 
```
SURVIVOR merge vcf_files 100 1 1 0 0 50 sample_merged.vcf
```
Lets check how good/bad the overlap is:
```
perl -ne 'print "$1\n" if /SUPP_VEC=([^,;]+)/'  sample_merged.vcf | sort | uniq -c 

```

As you can see you will get the pattern and the number of times the pattern occurs. The first number is the number of times it can be observed in the VCF file. The 2nd number is the pattern (0 or 1 depending if it was observed or not) 

