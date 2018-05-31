[![Build Status](https://travis-ci.org/tseemann/snippy.svg?branch=master)](https://travis-ci.org/tseemann/snippy) [![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)

# Snippy
Rapid haploid variant calling and core SNP phylogeny

## Author
Torsten Seemann (@torstenseemann)

## Synopsis

Snippy finds SNPs between a haploid reference genome and your NGS sequence
reads.  It will find both substitutions (snps) and insertions/deletions
(indels).  It will use as many CPUs as you can give it on a single computer
(tested to 64 cores).  It is designed with speed in mind, and produces a
consistent set of output files in a single folder.  It can then take a set
of Snippy results using the same reference and generate a core SNP alignment
(and ultimately a phylogenomic tree).

## Quick Start
```
% snippy --cpus 16 --outdir mysnps --ref Listeria.gbk --R1 FDA_R1.fastq.gz --R2 FDA_R2.fastq.gz
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
core.aln core.tab core.nway.tab core.txt core.vcf
```

# Installation

## Conda
Install [Conda](https://conda.io/docs/) or [Miniconda](https://conda.io/miniconda.html):
```
conda install -c bioconda -c conda-forge snippy
```

## Homebrew
Install [HomeBrew](http://brew.sh/) (Mac OS X) or [LinuxBrew](http://linuxbrew.sh/) (Linux).
```
brew untap homebrew/science
brew tap brewsci/bio
brew install snippy # COMING SOON!
```

## Source
This will install the latest version direct from Github. 
You'll need to add Prokka's `bin` directory to your `$PATH`.
```
cd $HOME
git clone https://github.com/tseemann/snippy.git
$HOME/bin/snippy --help
```

# Check installation
Ensure you have the latest version:
```
snippy --version
```
Check that all dependencies are installed and working:
```
snippy --check
```

# Calling SNPs

## Input Requirements
* a reference genome in FASTA or GENBANK format (can be in multiple contigs)
* sequence read file(s) in FASTQ or FASTA format (can be .gz compressed) format
* a folder to put the results in

## Output Files

Extension | Description
----------|--------------
.tab | A simple [tab-separated](http://en.wikipedia.org/wiki/Tab-separated_values) summary of all the variants
.csv | A [comma-separated](http://en.wikipedia.org/wiki/Comma-separated_values) version of the .tab file
.html | A [HTML](http://en.wikipedia.org/wiki/HTML) version of the .tab file
.vcf | The final annotated variants in [VCF](http://en.wikipedia.org/wiki/Variant_Call_Format) format
.bed | The variants in [BED](http://genome.ucsc.edu/FAQ/FAQformat.html#format1) format
.gff | The variants in [GFF3](http://www.sequenceontology.org/gff3.shtml) format
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Note that multi-mapping and unmapped reads are not present.
.bam.bai | Index for the .bam file
.log | A log file with the commands run and their outputs
.aligned.fa | A version of the reference but with `-` at position with `depth=0` and `N` for `0 < depth < --mincov` (**does not have variants**)
.consensus.fa | A version of the reference genome with *all* variants instantiated
.consensus.subs.fa | A version of the reference genome with *only substitution* variants instantiated
.raw.vcf | The unfiltered variant calls from Freebayes
.filt.vcf | The filtered variant calls from Freebayes
.vcf.gz | Compressed .vcf file via [BGZIP](http://blastedbio.blogspot.com.au/2011/11/bgzf-blocked-bigger-better-gzip.html)
.vcf.gz.csi | Index for the .vcf.gz via `bcftools index`)

:warning: :x: Snippy 4.x does **NOT** produce the following files that Snippy 3.x did

Extension | Description
----------|--------------
.vcf.gz.tbi | Index for the .vcf.gz via [TABIX](http://bioinformatics.oxfordjournals.org/content/27/5/718.full)
.depth.gz | Output of `samtools depth -aa` for the `.bam` file
.depth.gz.tbi | Index for the `.depth.gz` file

## Columns in the TAB/CSV/HTML formats

Name | Description
-----|------------
CHROM | The sequence the variant was found in eg. the name after the ```>``` in the FASTA reference
POS | Position in the sequence, counting from 1
TYPE | The variant type: snp msp ins del complex
REF | The nucleotide(s) in the reference
ALT | The alternate nucleotide(s) supported by the reads
EVIDENCE | Frequency counts for REF and ALT

If you supply a Genbank file as the `--reference` rather than a FASTA
file, Snippy will fill in these extra columns by using the genome annotation
to tell you which feature was affected by the variant:

Name | Description
-----|------------
FTYPE | Class of feature affected: CDS tRNA rRNA ...
STRAND | Strand the feature was on: + - .
NT_POS | Nucleotide position of the variant withinthe feature / Length in nt
AA_POS | Residue position / Length in aa (only if FTYPE is CDS)
LOCUS_TAG | The `/locus_tag` of the feature (if it existed)
GENE | The `/gene` tag of the feature (if it existed)
PRODUCT | The `/product` tag of the feature (if it existed)
EFFECT | The `snpEff` annotated consequence of this variant (ANN tag in .vcf)

## Variant Types

Type | Name | Example
-----|------|-------------
snp  | Single Nucleotide Polymorphism |  A => T
mnp  | Multiple Nuclotide Polymorphism | GC => AT
ins  | Insertion | ATT => AGTT
del  | Deletion | ACGG => ACG
complex | Combination of snp/mnp | ATTC => GTTA

## The variant caller

The variant calling is done by
[Freebayes](https://github.com/ekg/freebayes).  

If you wish to force more traditional cutoffs you can use these two options:

* `--mincov` is the minimum number of reads covering the variant position.
* `--minfrac` is the minimum proportion of those reads which must differ from the reference.

:warning: Snippy versions prior to 4.x used `--mincov 10 --micfrac 0.9` but this was removed in
Snippy 4.x, which now relies primarily on the Freebayes statistical models to find homozygous
variants of high probability. This increases true-positives and helps capture variants
in extreme GC regions where Illumina coverage can drop below 10x.

## Looking at variants in detail with `snippy-vcf_report`

If you run Snippy with the `--report` option it will automatically run
`snippy-vcf_report` and generate a `snps.report.txt` which has a section
like this for each SNP in `snps.vcf`:
```
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
>LBB_contig000001:10332 snp A=>T DP=7 Q=66.3052 [7]

         10301     10311     10321     10331     10341     10351     10361
tcttctccgagaagggaatataatttaaaaaaattcttaaataattcccttccctcccgttataaaaattcttcgcttat
........................................T.......................................
,,,,,,  ,,,,,,,,,,,,,,,,,,,,,t,,,,,,,,,,t,,t,,,,,,,,,,,,,,,,g,,,,,,,g,,,,,,,,,t,
,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,, .......T..................A............A.......
.........................A........A.....T...........    .........C..............
.....A.....................C..C........CT.................TA.............
,a,,,,,a,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,t,t,,,g,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
,,,,,ga,,,,,,,c,,,,,,,t,,,,,,,,,,g,,,,,,t,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,
                            ............T.C..............G...............G......
                                                    ,,,,,,,g,,,,,,,,g,,,,,,,,,,,
                                                           g,,,,,,,,,,,,,,,,,,,,
```

If you wish to generate this report *after* you have run Snippy, you can
run it directly:
```
cd snippydir
snippy-vcf_report --cpus 8 --auto > snps.report.txt
```
If you want a HTML version for viewing in a web browser, use the `--html` option:
```
cd snippydir
snippy-vcf_report --html --cpus 16 --auto > snps.report.html
```
It works by running `samtools tview` for each variant, which can be very slow
if you have 1000s of variants. Using `--cpus` as high as possible is recommended.

# Core SNP phylogeny

If you call SNPs for multiple isolates from the same reference, you can
produce an alignment of "core SNPs" which can be used to build a
high-resolution phylogeny (ignoring possible recombination).  A "core site"
is a genomic position that is present in _all_ the samples.  A core site can
have the same nucleotide in every sample ("monomorphic") or some samples can
be different ("polymorphic" or "variant").  If we ignore the complications
of "ins", "del" variant types, and just use variant sites, these are the "core SNP genome".

## Input Requirements
* a set of Snippy folders which used the same `--ref` sequence.

## Output Files

Extension | Description
----------|--------------
.aln | A core SNP alignment in the `--aformat` format (default FASTA)
.full.aln | A whole genome SNP alignment (includes invariant sites)
.tab | Tab-separated columnar list of **core** SN sites with alleles and annotations
.nway.tab | Tab-separated columnar list of **all** SNP sites with alleles and annotations
.vcf | Multi-sample VCF file with genotype `GT` tags for all discovered alleles
.txt | Tab-separated columnar list of alignment/core-size statistics

## Options
* If you want to mask certain regions of the genome, you can provide a BED file
  with the `--mask` parameter. Any SNPs in those regions will be excluded. This
  is common for genomes like *M.tuberculosis* where pesky repetitive PE/PPE/PGRS
  gebes cause false positives.
* If you use the `snippy --cleanup` option the reference files will be deleted.
  This means `snippy-core` can not "auto-find" the reference. In this case you
  simply use `snippy-core --reference REF` to provide the reference in FASTA format.
* If you want to exclude the reference genome from the alignment, 
  use `snippy-core --noref`.

# Advanced usage

## Increasing speed when too many reads

Sometimes you will have far more sequencing depth that you need to call SNPs.
A common problem is a whole MiSeq flowcell for a single bacterial isolate,
where 25 million reads results in genome depth as high as 2000x. This makes
Snippy far slower than it needs to be, as most SNPs will be recovered with
50-100x depth. If you know you have 10 times as much data as you need,
Snippy can randomly sub-sample your FASTQ data:
```
# have 1000x depth, only need 100x so sample at 10%
snippy --subsample 0.1  ...
<snip>
Sub-sampling reads at rate 0.1
<snip>
```

## Only calling SNPs in particular regions

If you are looking for specific SNPs, say AMR releated ones in particular genes
in your reference genome, you can save much time by only calling variants there.
Just put the regions of interest into a BED file:
```
snippy --targets sites.bed ...
```

## Finding SNPs between contigs

Sometimes one of your samples is only available as contigs, without
corresponding FASTQ reads. You can still use these contigs with Snippy
to find variants against a reference. It does this by shredding the contigs
into 250 bp single-end reads at `2 &times; --mincov` uniform coverage.

To use this feature, instead of providing `--R1` and `--R2` you use the
`--ctgs` option with the contigs file:

```
% ls
ref.gbk mutant.fasta

% snippy --outdir mut1 --ref ref.gbk --ctgs mut1.fasta
Shredding mut1.fasta into pseudo-reads.
Identified 257 variants.

% snippy --outdir mut2 --ref ref.gbk --ctgs mut2.fasta
Shredding mut2.fasta into pseudo-reads.
Identified 413 variants.

% snippy-core mut1 mut2 
Found 129 core SNPs from 541 variant sites.

% ls
core.aln core.full.aln ...
```

This output folder is completely compatible with `snippy-core` so you can
mix FASTQ and contig based `snippy` output folders to produce alignments.

## Correcting assembly errors

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
and PacBio SMRT contigs, we primarily expect to find homopolymer errors and hence
expect to see `ins` more than `snp` type variants. 

In this case you need to run the correcting process manually using these steps:

```
% cd snippy-outdir
% cp snps.vcf corrections.vcf
% $EDITOR corrections.vcf
% bgzip -c corrections.vcf > corrections.vcf.gz
% tabix -p vcf corrections.vcf.gz
% vcf-consensus corrections.vcf.gz < ref.fa > corrected.fa
```

You may wish to _iterate_ this process by using `corrected.fa` as a new `--ref` for
a repeated run of Snippy. Sometimes correcting one error allows BWA to align things
it couldn't before, and new errors are uncovered.

Snippy may not be the best way to correct assemblies - you should consider
dedicated tools such as [PILON](http://www.broadinstitute.org/software/pilon/) 
or [iCorn2](http://icorn.sourceforge.net/), or adjust the 
Quiver parameters (for Pacbio data).

## Unmapped Reads

Sometimes you are interested in the reads which did *not* align to the reference genome.
These reads represent DNA that was novel to *your* sample which is potentially interesting.
A standard strategy is to *de novo* assemble the unmapped reads to discover these novel
DNA elements, which often comprise mobile genetic elements such as plasmids.

By default, Snippy does **not** keep the unmapped reads, not even in the BAM file.
If you wish to keep them, use the `--unmapped` option and the unaligned reads will
be saved to a compressed FASTQ file:

```
% snippy --outdir out --unmapped ....

% ls out/
snps.unmapped.fastq.gz ....
```

# Information

## Etymology

The name Snippy is a combination of
[SNP](http://en.wikipedia.org/wiki/Single-nucleotide_polymorphism)
(pronounced "snip") , [snappy](http://www.thefreedictionary.com/snappy)
(meaning "quick") and [Skippy the Bush
Kangaroo](http://en.wikipedia.org/wiki/Skippy_the_Bush_Kangaroo) (to
represent its Australian origin)

## License

Snippy is free software, released under the 
[GPL (version 2)](https://raw.githubusercontent.com/tseemann/snippy/master/LICENSE).

## Issues

Please submit suggestions and bug reports to the 
[Issue Tracker](https://github.com/tseemann/snippy/issues)

## Requirements

* Perl >= 5.12
* Perl Modules: Time::Piece (core with modern Perl), Bioperl >= 1.6
* bwa mem >= 0.7.12 
* samtools >= 1.7
* bcftools >= 1.7
* GNU parallel >= 2013xxxx
* freebayes >= 1.1 (freebayes, freebayes-parallel, fasta_generate_regions.py)
* vcflib >= 1.0 (vcfstreamsort, vcfuniq, vcffirstheader)
* snpEff >= 4.3
* samclip >= 0.2

## Bundled binaries

For Linux (compiled on Ubuntu 16.04 LTS) and macOS (compiled on High Sierra Brew) all
the binaries, JARs and scripts are included.

