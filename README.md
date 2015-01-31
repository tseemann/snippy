#Snippy

Rapid haploid SNP calling by Torsten Seemann

##Description

###Synopsis
Snippy finds SNPs between a haploid reference genome and your NGS sequence reads (handles reads >500bp long). It will find both substitutions and insertions/deletions (indels). It will use as many CPUs and RAM as you can give it on a single computer (tested to 64 cores). It is designed with speed in mind, and it may occasionally miss some complicated SNPs.

###Input
Snippy needs a reference genome in FASTA format (can be in multiple contigs) and Illumina read files in FASTQ (can be .gz compressed) format. The reads can be single-end, interleaved paired-end, or separate pair files.

###Output
Snippy provides the called SNPs in VCF, TSV, BED and GFF3 format. It can optionally produce a TXT file with stacked alignments for each SNP suited to manual examination. The index BAM and FASTA files are also kept by default.

###Etymology
The name Snippy is a combination of SNP (pronounced "snip") , snappy (meaning "quick") and Skippy the Bush Kangaroo (to represent its Australian origin)

##Usage

```
% snippy --threads 32 --outdir mysnps --ref Ecoli.gbk R1.fastq.gz R2.fastq.gz
Walltime used: 3 min, 42 sec
Results in folder: mysnps
```

```
% ls mysnps
snps.vcf snps.bed snps.gff snps.csv aln.bam ref.fa ...
```

```
% head -5 mysnps/snps.csv
CHROM   POS     TYPE    REF     ALT     EVIDENCE        FEATURES
chr      5958   snp     A       G       G:44            ECO_0001 dnaA replication protein DnaA
chr     35524   snp     G       T       T:73 G:1 C:1      
chr     45722   ins     ATT     ATTT    ATTT:43 ATT:1   ECO_0045  gyrA gyrase
chr    100541   del     CAAA    CAA     CAA:38 CAAA:1   ECO_0121       hypothetical protein
```

###License

Snippy is free software, released under the GPL (version 3).

###Requirements

* Perl >= 5.6
* BioPerl >= 1.6
* bwa mem >= 0.7.5 
* samtools >= 0.1.19+
* freebayes >= 0.9.9.2 
* vcflib (vcfstreamsort, vcfuniq, vcffirstheader)

