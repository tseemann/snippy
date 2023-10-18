[![CI](https://github.com/tseemann/snippy/actions/workflows/main.yml/badge.svg)](https://github.com/tseemann/snippy/actions/workflows/main.yml)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html)
![Don't judge me](https://img.shields.io/badge/Language-Perl_5-steelblue.svg)

# Snippy
Rapid haploid variant calling and core genome alignment

## Author
[Torsten Seemann](https://twitter.com/torstenseemann)

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

% snippy-core --corefreq .95 --prefix core mysnps1 mysnps2 mysnps3 mysnps4 
Loaded 4 SNP tables.
Found 2814 core SNPs from 96615 SNPs.

% ls core.*
core.aln core.f.95.aln core.full.aln core.tab core.tab core.txt core.vcf
```

# Installation

## Conda
Install [Bioconda](https://bioconda.github.io/user/install.html) then:
```
conda install -c conda-forge -c bioconda -c defaults snippy
```

## Homebrew
Install [Homebrew](http://brew.sh/) (MacOS)
or [LinuxBrew](http://linuxbrew.sh/) (Linux) then:
```
brew install brewsci/bio/snippy
```

## Source
This will install the latest version direct from Github. 
You'll need to add Snippy's `bin` directory to your `$PATH`.
```
cd $HOME
git clone https://github.com/tseemann/snippy.git
$HOME/snippy/bin/snippy --help
```

# Check installation
Ensure you have the desired version:
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
.bam | The alignments in [BAM](http://en.wikipedia.org/wiki/SAMtools) format. Includes unmapped, multimapping reads. Excludes duplicates.
.bam.bai | Index for the .bam file
.log | A log file with the commands run and their outputs
.aligned.fa | A version of the reference but with `-` at position with `depth=0` and `N` for `0 < depth < --mincov` (**now have variants**)
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
The key parameters under user control are:

* `--mincov` - the minimum number of reads covering a site to be considered (default=10)
* `--minfrac` - the minimum proportion of those reads which must differ from the reference
* `--minqual` - the minimum VCF variant call "quality" (default=100)

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

## Options

* `--rgid` will set the Read Group (`RG`) ID (`ID`) and Sample (`SM`) in the BAM and VCF file.
If not supplied, it will will use the `--outdir` folder name for both `ID` and `SM`.

* `--mapqual` is the minimum mapping quality to accept in variant calling. BWA MEM using `60`
to mean a read is "uniquely mapped". 

* `--basequal` is minimum quality a nucleotide needs to be used in variant calling. We use
`13` which corresponds to error probability of ~5%. It is a traditional SAMtools value.

* `--maxsoft` is how many bases of an alignment to allow to be soft-clipped before discarding
the alignment. This is to encourage global over local alignment, and is passed to the
`samclip` tool.

* `--mincov` and `--minfrac` are used to apply hard thresholds to the variant calling
beyond the existing statistical measure.. The optimal values depend on your sequencing
depth and contamination rate. Values of 10 and 0.9 are commonly used.

* `--targets` takes a BED file and only calls variants in those regions. Not normally needed
unless you are only interested in variants in specific locii (eg. AMR genes) but are still
performing WGS rather than amplicon sequencing.

* `--contigs` allows you to call SNPs from contigs rather than reads. It shreds the contigs
into synthetic reads, as to put the calls on even footing with other read samples in a 
multi-sample analysis.


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

### Using `snippy-multi`

To simplify running a set of isolate sequences (reads or contigs)
against the same reference, you can use the `snippy-multi` script.
This script requires a *tab separated* input file as follows, and
can handle paired-end reads, single-end reads, and assembled contigs.
```
# input.tab = ID R1 [R2]
Isolate1	/path/to/R1.fq.gz	/path/to/R2.fq.gz
Isolate1b	/path/to/R1.fastq.gz	/path/to/R2.fastq.gz
Isolate1c	/path/to/R1.fa		/path/to/R2.fa
# single end reads supported too
Isolate2	/path/to/SE.fq.gz
Isolate2b	/path/to/iontorrent.fastq
# or already assembled contigs if you don't have reads
Isolate3	/path/to/contigs.fa
Isolate3b	/path/to/reference.fna.gz
```
Then one would run this to generate the output script.
The first parameter should be the `input.tab` file.
The remaining parameters should be any remaining
shared `snippy` parameters. The `ID` will be used for
each isolate's `--outdir`.
```
% snippy-multi input.tab --ref Reference.gbk --cpus 16 > runme.sh

% less runme.sh   # check the script makes sense

% sh ./runme.sh   # leave it running over lunch
```
It will also run `snippy-core` at the end to generate the
core genome SNP alignment files `core.*`.

## Output Files

Extension | Description
----------|--------------
.aln | A core SNP alignment in the `--aformat` format (default FASTA)
.full.aln | A whole genome SNP alignment (includes invariant sites)
.tab | Tab-separated columnar list of **core** SNP sites with alleles but NO annotations
.vcf | Multi-sample VCF file with genotype `GT` tags for all discovered alleles
.txt | Tab-separated columnar list of alignment/core-size statistics
.ref.fa | FASTA version/copy of the `--ref`
.self_mask.bed | BED file generated if `--mask auto` is used.

## Why is `core.full.aln` an alphabet soup?

The `core.full.aln` file is a FASTA formatted mutliple sequence alignment file.
It has one sequence for the reference, and one for each sample participating in
the core genome calculation.  Each sequence has the same length as the reference
sequence.

Character | Meaning
----------|-----------
`ATGC`    | Same as the reference
`atgc`    | Different from the reference
`-`       | Zero coverage in this sample **or** a deletion relative to the reference
`N`       | Low coverage in this sample (based on `--mincov`)
`X`       | Masked region of reference (from `--mask`)
`n`       | Heterozygous or poor quality genotype  (has `GT=0/1` or `QUAL < --minqual` in `snps.raw.vcf`)

You can remove all the "weird" characters and replace them with `N` using the included
`snippy-clean_full_aln`.  This is useful when you need to pass it to a tree-building
or recombination-removal tool:

```
% snippy-clean_full_aln core.full.aln > clean.full.aln
% run_gubbins.py -p gubbins clean.full.aln
% snp-sites -c gubbins.filtered_polymorphic_sites.fasta > clean.core.aln
% FastTree -gtr -nt clean.core.aln > clean.core.tree
```

## Options

* If you want to mask certain regions of the genome, you can provide a BED file
  with the `--mask` parameter. Any SNPs in those regions will be excluded. This
  is common for genomes like *M.tuberculosis* where pesky repetitive PE/PPE/PGRS
  genes cause false positives, or masking phage regions. A `--mask` bed file
  for *M.tb* is provided with Snippy in the `etc/Mtb_NC_000962.3_mask.bed`
  folder. It is derived from the XLSX file from https://gph.niid.go.jp/tgs-tb/
* If you use the `snippy --cleanup` option the reference files will be deleted.
  This means `snippy-core` can not "auto-find" the reference. In this case you
  simply use `snippy-core --reference REF` to provide the reference in FASTA format.

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

* perl >= 5.18
* bioperl >= 1.7
* bwa mem >= 0.7.12 
* minimap2 >= 2.0
* samtools >= 1.7
* bcftools >= 1.7
* bedtools >= 2.0
* GNU parallel >= 2013xxxx
* freebayes >= 1.1 (freebayes, freebayes-parallel, fasta_generate_regions.py)
* vcflib >= 1.0 (vcfstreamsort, vcfuniq, vcffirstheader)
* [vt](https://genome.sph.umich.edu/wiki/Vt) >= 0.5
* snpEff >= 4.3
* samclip >= 0.2
* seqtk >= 1.2
* snp-sites >= 2.0
* any2fasta >= 0.4
* wgsim >= 1.8 (for testing only - `wgsim` command)

## Bundled binaries

For Linux (compiled on Ubuntu 16.04 LTS) and macOS (compiled on High Sierra Brew) 
some of the binaries, JARs and scripts are included.

