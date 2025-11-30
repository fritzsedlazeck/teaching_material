# Wonderful world of Structural Variant Calling!

Lead: Fritz J Sedlazeck

Human Genome Sequencing Center,
Baylor College of Medicine

Computer Science Department,
Rice University
***

## SAMTools and BCFTools

### SAMTools
[SAMtools](https://samtools.github.io/) is a set of utilities for interacting with and post-processing short DNA sequence read alignments in the [SAM](https://www.htslib.org/doc/sam.html) (Sequence Alignment/Map), BAM (Binary Alignment/Map) and CRAM formats.

#### Useful commands

View the sequences of a BAM file
```bash
samtools view file.bam | less -S
```

View the header of a BAM file
```bash
samtools view -H file.bam | less -S
```

View the sequences of a BAM file for a specific chromosome
```bash
samtools view file.bam chromosome_name | less -S
```

View the sequences of a BAM file for a specific chromosomal region
```bash
samtools view file.bam chromosome:start-end | less -S
```

Count reads
```bash
samtools view -c file.bam 
```

---

### BCFTools
[BCFtools](https://samtools.github.io/bcftools/bcftools.html) is a set of utilities that manipulate variant calls in the Variant Call Format (VCF) and its binary counterpart BCF. All commands work transparently with both VCFs and BCFs, both uncompressed and BGZF-compressed.

[Here](https://samtools.github.io/bcftools/howtos/index.html)are a couple of "How to" recipes that are useful

#### Useful commands

**Note: The VCF need to be compressed**
Compress the VCF file if not. VCF file that are compressed have a *.gz* extension

View the variants of a VCF file
```bash
bcftools view file.vcf.gz | less -S
```

View the variants of a VCF file **with no headers** 
```bash
bcftools view --no-header file.vcf.gz | less -S
```

View the variants of a VCF file and include only *PASS* variants in the filter
```bash
bcftools view --apply-filters 'PASS' file.vcf.gz | less -S
```

View the variants of a VCF and filter to include specific components of the INFO or FORMAT fields
```bash
bcftools view --include '<TAG> = <VALUE>' file.vcf.gz | less -S
```

View the variants of a VCF and filter to exclude specific components of the INFO or FORMAT fields
```bash
bcftools view --exclude '<TAG> = <VALUE>' file.vcf.gz | less -S
```

Compute Mendelian rules from a family trio, the sample names corresponding the the mother, father and proband need to be given in specific order.
```bash
bcftools +mendelian --trio M,F,P -m c file.vcf.gz
```
