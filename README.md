#Snippy
Rapid haploid variant calling and core SNP phylogeny

##Author
Torsten Seemann (@torstenseemann)

##Synopsis

Snippy finds SNPs between a haploid reference genome and your NGS sequence
reads.  It will find both substitutions (snps) and insertions/deletions
(indels).  It will use as many CPUs as you can give it on a single computer
(tested to 64 cores).  It is designed with speed in mind, and produces a
consistent set of output files in a single folder.  It can then take a set
of Snippy results using the same reference and generate a core SNP alignment
(and ultimately a phylogenomic tree).

##Quick Start
```
% snippy --cpus 16 --outdir mysnps --ref Listeria.gbk --R1 FDA_R1.fastq.gz --R1 FDA_R2.fastq.gz
<cut>
Walltime used: 3 min, 42 sec
Results folder: mysnps
Done.

% ls mysnps
snps.vcf snps.bed snps.gff snps.csv snps.tab snps.html 
snps.bam snps.txt reference/ ...

% head -5 mysnps/snps.tab
CHROM  POS     TYPE    REF   ALT    EVIDENCE        FTYPE STRAND NT_POS AA_POS LOCUS_TAG GENE PRODUCT EFFECT
chr      5958  snp     A     G      G:44 A:0        CDS   +      41/600 13/200 ECO_0001  dnaA replication protein DnaA missense_variant c.548A>C p.Lys183Thr
chr     35524  snp     G     T      T:73 G:1 C:1    tRNA  -   
chr     45722  ins     ATT   ATTT   ATTT:43 ATT:1   CDS   -                    ECO_0045  gyrA DNA gyrase
chr    100541  del     CAAA  CAA    CAA:38 CAAA:1   CDS   +                    ECO_0179      hypothetical protein
plas      619  complex GATC  AATA   GATC:28 AATA:0  
plas     3221  mnp     GA    CT     CT:39 CT:0      CDS   +                    ECO_p012  rep  hypothetical protein

% snippy-core --prefix core mysnps1 mysnps2 mysnps3 mysnps4 
Loaded 4 SNP tables.
Found 2814 core SNPs from 96615 SNPs.

% ls core.*
core.aln core.tab core.txt
```

#Installation

##Homebrew
Install [HomeBrew](http://brew.sh/) (Mac OS X) or [LinuxBrew](http://brew.sh/linuxbrew/) (Linux).

    brew tap homebrew/science
    brew tap chapmanb/cbl
    brew tap tseemann/homebrew-bioinformatics-linux
    brew install snippy
    snippy --help

##Source
This will install the latest version direct from Github. You'll need to add the ```bin``` directory to your PATH.

    cd $HOME
    git clone https://github.com/tseemann/snippy.git
    $HOME/snippy/bin/snippy --help

#Calling SNPs

##Input Requirements
* a reference genome in FASTA or GENBANK format (can be in multiple contigs)
* sequence read files in FASTQ or FASTA format (can be .gz compressed) format
* a folder to put the results in

##Output Files

Extension | Description
----------|--------------
.tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants
.csv | A [comma-separated](http://en.wikipedia.org/wiki/Comma-separated_values) version of the .tab file
.html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file
.vcf | The final annotated variants in [VCF](http://en.wikipedia.org/wiki/Variant_Call_Format) format
.vcf.gz | Compressed .vcf file via [BGZIP](http://blastedbio.blogspot.com.au/2011/11/bgzf-blocked-bigger-better-gzip.html) 
.vcf.gz.tbi | Index for the .vcf.gz via [TABIX](http://bioinformatics.oxfordjournals.org/content/27/5/718.full)
.bed | The variants in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format
.gff | The variants in [GFF3](http://www.sequenceontology.org/gff3.shtml) format
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Note that multi-mapping and unmapped reads are not present.
.bam.bai | Index for the .bam file
.raw.vcf | The unfiltered variant calls from Freebayes
.filt.vcf | The filtered variant calls from Freebayes
.log | A log file with the commands run and their outputs
.consensus.fa | A version of the reference genome with all variants instantiated
.aligned.fa | A version of the reference but with - for unaligned and N for depth < --minfrac (**does not have variants**)
.depth.gz | Output of ```samtools depth``` for the .bam file
.depth.gz.tbi | Index for the .depth.gz (_currently unused_)

##Columns in the TAB/CSV/HTML formats

Name | Description
-----|------------
CHROM | The sequence the variant was found in eg. the name after the ```>``` in the FASTA reference
POS | Position in the sequence, counting from 1
TYPE | The variant type: snp msp ins del complex
REF | The nucleotide(s) in the reference
ALT | The alternate nucleotide(s) supported by the reads
EVIDENCE | Frequency counts for REF and ALT

If you supply a Genbank file as the ```--reference``` rather than a FASTA file, Snippy will fill in these extra columns by using the genome annotation to tell you which feature was affected by the variant:

Name | Description
-----|------------
FTYPE | Class of feature affected: CDS tRNA rRNA ...
STRAND | Strand the feature was on: + - .
NT_POS | Nucleotide position of the variant withinthe feature / Length in nt
AA_POS | Residue position / Length in aa (only if FTYPE is CDS)
LOCUS_TAG | The ```/locus_tag``` of the feature (if it existed)
GENE | The ```/gene``` tag of the feature (if it existed)
PRODUCT | The ```/product``` tag of the feature (if it existed)
EFFECT | The ```snpEff``` annotated consequence of this variant

##Variant Types

Type | Name | Example
-----|------|-------------
snp  | Single Nucleotide Polymorphism |  A => T
mnp  | Multiple Nuclotide Polymorphism | GC => AT
ins  | Insertion | ATT => AGTT
del  | Deletion | ACGG => ACG
complex | Combination of snp/mnp | ATTC => GTTA

##The variant caller
The variant calling is done by [Freebayes](https://github.com/ekg/freebayes). However, Snippy uses a very simple model for reporting variants, relying on two main options:
* ```--mincov``` is the minimum number of reads covering the variant position.
* ```--minfrac``` is the minimum proportion of those reads which must differ from the reference.

By default Snippy uses ```--mincov 10 --minfrac 0.9``` which is reasonable for most cases, but for very high coverage data you may get mixed populations such as (REF:310 ALT:28). Snippy may use a more statistical approach in future versions like [Nesoni](https://github.com/Victorian-Bioinformatics-Consortium/nesoni) does.

#Correcting assembly errors

The _de novo_ assembly process attempts to reconstruct the reads into the original 
DNA sequences they were derived from. These reconstructed sequences are called 
_contigs_ or _scaffolds_. For various reasons, small errors can be introduced into
the assembled contigs which are not supported by the original reads used in the 
assembly process.

A common strategy is to align the reads back to the contigs to check for discrepancies.
These errors appear as variants (SNPs and indels). If we can _reverse_ these variants
than we can "correct" the contigs to match the evidence provided by the original reads.
Obviously this strategy can go wrong if one is not careful about _how_ the read alignment
is performed and which variants are accepted.

Snippy is able to help with this contig correction process. In fact, it produces a
`snps.consensus.fa` FASTA file which is the `ref.fa` input file provided but with the
discovered variants in `snps.vcf` applied! 

However, Snippy is not perfect and sometimes finds questionable variants. Typically
you would make a copy of `snps.vcf` (let's call it `corrections.vcf`) and remove those
lines corresponding to variants we don't trust. For example, when correcting Roche 454
and PacBio SMRT contigs, we primarily expect to find A/T homopolymer errors and hence
expect to see `ins` more than `snp` type variants. In this case you need to run the
correcting process manually using these steps:

```
% cd snippy-outdir
% cp snps.vcf corrections.vcf
% $EDITOR corrections.vcf
% bgzip -c corrections.vcf > corrections.vcf.gz
% tabix -p vcf corrections.vcf.gz
% vcf-consensus corrections.vcf.gz < ref.fa > corrected.fa
```


#Core SNP phylogeny

If you call SNPs for multiple isolates from the same reference, you can
produce an alignment of "core SNPs" which can be used to build a
high-resolution phylogeny (ignoring possible recombination).  A "core site"
is a genomic position that is present in _all_ the samples.  A core site can
have the same nucleotide in every sample ("monomorphic") or some samples can
be different ("polymorphic" or "variant").  If we ignore the complications
of "ins", "del" variant types, and just use variant sites, these are the "core SNP genome".

##Input Requirements
* a set of Snippy folders which used the same ``--ref`` sequence.

##Output Files

Extension | Description
----------|--------------
.aln | A core SNP alignment in the ```--aformat``` format (default FASTA)
.tab | Tab-separated columnar list of core SNP sites with alleles and annotations
.txt | Tab-separated columnar list of alignment/core-size statistics

#Information

##Etymology
The name Snippy is a combination of [SNP](http://en.wikipedia.org/wiki/Single-nucleotide_polymorphism) (pronounced "snip") , [snappy](http://www.thefreedictionary.com/snappy) (meaning "quick") and [Skippy the Bush Kangaroo](http://en.wikipedia.org/wiki/Skippy_the_Bush_Kangaroo) (to represent its Australian origin)

##License
Snippy is free software, released under the GPL (version 3).

##Issues
Please submit suggestions and bug reports here: https://github.com/tseemann/snippy/issues

##Requirements
* Perl >= 5.6
* BioPerl >= 1.6
* bwa mem >= 0.7.12 
* samtools >= 1.1
* GNU parallel > 2013xxxx
* freebayes >= 0.9.20 
* freebayes sripts (freebayes-parallel, fasta_generate_regions.py)
* vcflib (vcffilter, vcfstreamsort, vcfuniq, vcffirstheader)
* vcftools (vcf-consensus)
* snpEff >= 4.1

##Bundled binaries
For a modern Linux system (Ubuntu >= 12.04) and Mac OS X all the binaries, JARs and scripts are included. 

