# Wonderful world of Structural Variant Calling! 

***
Lead: Fritz J Sedlazeck

Human Genome Sequencing Center,
Baylor College of Medicine

Computer Science Department,
Rice University 
***

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
source Share/conda_init.sh
```



## Part 2: Illumina based SV detection 
In this part we will use Manta to interpret the signals we can obtain from abnormally mapped paired end Illumina reads. 

Again you might get the file already on your account or you can download from here:
```
wget https://www.dropbox.com/s/1mjb0tjbtmrcugj/our_refCrypto_reads.sort.bam?dl=0
mv our_refCrypto_reads.sort.bam?dl=0 our_refCrypto_reads.sort.bam
```

You might need to index the reference and the bam file that we provided to you:
```
samtools faidx GCF_000165345.1.fa
samtools index our_refCrypto_reads.sort.bam
```

### 1. Initiate the run:
Next we initiate the Manta run:
```
configManta.py --bam=our_refCrypto_reads.sort.bam --referenceFasta=GCF_000165345.1.fa  --runDir=Out_Manta
```
This should just take seconds as it initiates the folder structure and specifies for the subsequent process to use our mapped reads and our reference file. In addition, we specify the output to be written in `Out_Manta`

### 2. Run the analysis:
```
python2 Out_Manta/runWorkflow.py -j 2 -m local -g 10
```

This will launch the Manta pipeline that we previously configured. `-j` specifies the number of CPU threads (2 in our case), `-m` local indicates that it should not try to run things on different nodes or instances and `-g 30` specifies the available memory for the process in GB.

Manta now searches for abnormal paired-end reads and split-reads across our mapped reads. These will be analyzed together and clustered to identify SV in these samples. After ~2-3 minutes you should see that the program has finished.

Our SV calling results can be found here (Out_Manta/results/variants/). Let us open quickly the output of the highest quality SV files:
```
cd  Out_Manta/results/variants/
ls 
```
As you can see, we have multiple VCF files. These represent the different stages of Manta and the confidence level for the SV calls. We typically use the `diploidSV.vcf.gz` file.


Let us open quickly the output of the highest quality SV files:
```
gunzip diploidSV.vcf.gz
mv  diploidSV.vcf illumina.vcf
less illumina.vcf
```

Please dont forget to switch back and copy that file. 
```
cd ../../../
cp Out_Manta/results/variants/illumina.vcf .
```

You can get the file here if you had difficulties:
```
wget https://www.dropbox.com/s/wcro3nzkird5nx8/illumina.vcf?dl=0
mv illumina.vcf?dl=0 illumina.vcf
```

Lets count how many SV we could identify: 
```
grep -vc '#' illumina.vcf
```

Let us all discuss the different variant types and how they are reported! 


## Part 3: Long read based SV detection 
Finally we are ready for the Oxford Nanopore detection using sniffles. For this use the "ont_mapped.sort.bam" file that I have previously mapped using minimap2. 
You might have that file on your account, but if not you can download it here:
```
wget https://www.dropbox.com/s/ttafrqaikst8xea/ont_prev.sort.bam?dl=0 
mv ont_prev.sort.bam?dl=0 ont_prev.sort.bam
```

Using Sniffles v2 this should be a simple command like:

```
samtools index ont_prev.sort.bam
sniffles -i ont_prev.sort.bam -v sniffles.vcf
```

You can also download the file from here if you had issues:
```
wget https://www.dropbox.com/s/7fpgnoq818mxsnk/sniffles.vcf?dl=0
mv sniffles.vcf?dl=0 sniffles.vcf
```

Next we can inspect the file with e.g.:
```
less -S sniffles.vcf
```
or 
```
grep -vc '#' sniffles.vcf
```

How many SV did you detect? 

## Part 4: Structural Variant comparison

Now that we generated Assembly, Illumina  and ONT based SV calls it is time to compare them. One tool that you can use for this very easily is SURVIVOR. SURVIVOR is a very simple method to compare SV but also includes multiple other methods that might be useful. Also feel free to check out Truvari, which is one of our newest methods to compare SV. When comparing SV it is quite important to take multiple things into considerations. Feel free to ask if you want to discuss certain points more. 

For SURVIVOR we want to use the merge option. Before doing this, the merge option requires a file including all paths and VCF files that you want to compare. Thus, we generate the file like this:
```
ls sniffles.vcf > vcf_files
ls assemblytics.vcf >> vcf_files
ls illumina.vcf >> vcf_files
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

