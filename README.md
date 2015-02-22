#Snippy

Rapid haploid SNP calling by Torsten Seemann

##Description

###Synopsis
Snippy finds SNPs between a haploid reference genome and your NGS sequence reads. It will find both substitutions (snps) and insertions/deletions (indels). It will use as many CPUsas you can give it on a single computer (tested to 64 cores). It is designed with speed in mind, and produces a clean set of output files.

###Etymology
The name Snippy is a combination of [SNP](http://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (pronounced "snip") , [snappy](http://www.thefreedictionary.com/snappy) (meaning "quick") and [Skippy the Bush Kangaroo](http://en.wikipedia.org/wiki/Skippy_the_Bush_Kangaroo) (to represent its Australian origin)

###Input Files
Snippy needs 
* a reference genome in FASTA or GENBANK format (can be in multiple contigs)
* sequence read files in FASTQ or FASTA format (can be .gz compressed) format

###Output Files

Extension | Description
----------|--------------
.tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants
.csv | A [comma-separated](http://en.wikipedia.org/wiki/Comma-separated_values) version of the .tab file
.html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file
.vcf | The variants in [VCF](http://en.wikipedia.org/wiki/Variant_Call_Format) format
.vcf.gz | Compressed .vcf file via [BGZIP](http://blastedbio.blogspot.com.au/2011/11/bgzf-blocked-bigger-better-gzip.html) 
.vcf.gz.tbi | Index for the .vcf.gz via [TABIX](http://bioinformatics.oxfordjournals.org/content/27/5/718.full)
.bed | The variants in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format
.gff | The variants in [GFF3](http://www.sequenceontology.org/gff3.shtml) format
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Note that multi-mapping and unmapped reads are not present.
.bam.bai | Index for the .bam file
.raw.vcf | The unfiltered variant calls from Freebayes
.log | A log file with the commands run and their outputs
.consensus.fa | A version of the reference genome with all variants instantiated


##Usage

```
% snippy --cpus 16 --outdir mysnps --ref EHEC.gbk --R1 FDA_R1.fastq.gz FDA_R2.fastq.gz
Walltime used: 3 min, 42 sec
Results folder: mysnps
Done.
```

```
% ls mysnps
reference/ snps.vcf snps.bed snps.gff snps.csv aln.bam  ...
```

```
% head -5 mysnps/snps.csv
CHROM   POS     TYPE    REF     ALT     EVIDENCE        FEATURES
chr      5958   snp     A       G       G:44 A:0        ECO_0001 dnaA replication protein DnaA
chr     35524   snp     G       T       T:73 G:1 C:1      
chr     45722   ins     ATT     ATTT    ATTT:43 ATT:1   ECO_0045  gyrA gyrase
chr    100541   del     CAAA    CAA     CAA:38 CAAA:1   ECO_0121       hypothetical protein
chr    241221   mnp     GA      CT      CT:39 CT:0      ECO_0121       hypothetical protein
chr    341819   complex GATC    AATA    GATC:28 AATA:0  
```

###License

Snippy is free software, released under the GPL (version 3).

###Requirements

* Perl >= 5.6
* BioPerl >= 1.6
* bwa mem >= 0.7.12 
* samtools >= 1.1
* freebayes >= 0.9.20 
* vcflib (vcfstreamsort, vcfuniq, vcffirstheader)

