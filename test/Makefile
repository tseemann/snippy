SHELL := /bin/bash
MAKEFLAGS += --no-builtin-rules
MAKEFLAGS += --no-builtin-variables

.SUFFIXES:
.DELETE_ON_ERROR:
.SECONDARY:
.ONESHELL:
.DEFAULT: all

CPUS = 1
SNIPPY = ../bin/snippy --cpus $(CPUS)
SNIPPYCORE = ../bin/snippy-core
R1 = reads_R1.fastq.gz
R2 = reads_R2.fastq.gz

GBK=example.gbk
REF=example.fna
BED=example.bed

all: core.txt
	wc -l ?/snps.tab

core.txt : A/snps.tab B/snps.tab C/snps.tab D/snps.tab
	$(SNIPPYCORE) --prefix core --ref $(REF) A B C D 2> core.log
	$(SNIPPYCORE) --prefix core-noref --ref $(REF) --noref A B C D 2> core-noref.log
	$(SNIPPYCORE) --prefix core-mask --ref $(REF) --mask $(BED) A B C D 2> core-mask.log

%/snps.tab : %/R1.fq %/R2.fq
	$(SNIPPY) --force --ref $(REF) --outdir $(@D) --R1 $(word 1,$^) --R2 $(word 2,$^) 2> $(@D)/snippy.log

%/R1.fq %/R2.fq : $(REF)
	mkdir -p $(@D)
	wgsim -S $(shell echo $(@D) | od -N 1 -i | head -n 1 | cut -c9-) \
		-h -r 0.02 -N 10000 -1 100 -2 100 -d 200 \
		$< $@ $(subst R2,R1,$@) > $(@D)/dwgsim.tab

%.fna : %.gbk
	seqret -auto -filter -osformat2 fasta < $< > $@

bigtest/snps.tab : genome.fa reads_R1.fq.gz reads_R2.fq.gz
	$(SNIPPY) --force --cpus 36 --outdir bigtest --ref $(word 1,$^) --R1 $(word 2,$^) --R2 $(word 3,$^)

genome.fa : 
	curl http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/cv/CVOW01.fasta.gz \
		| gunzip | sed -e 's/ENA|CVOW01000001|//' > $@

reads_R1.fq.gz :
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/006/SRR1157476/SRR1157476_1.fastq.gz \
		| gunzip | seqtk sample 0.2 - | pigz -9 > $@
	
reads_R2.fq.gz :
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/006/SRR1157476/SRR1157476_2.fastq.gz \
		| gzip -d -c | seqtk sample 0.2 - | pigz -9 > $@

clean :
	rm -fr ./bigtest core*.* ./{A..D}

	