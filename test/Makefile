
OPTS =
SNIPPY = ../bin/snippy $(OPTS)
R1 = reads_R1.fastq.gz
R2 = reads_R2.fastq.gz

.PHONY: test1 test2 test3 test4 data bigtest

quick: test3

all: test1 test2 test3 test4

test1:
	$(SNIPPY) --outdir test1 --force --ref example.fna --pe1 $(R1) --pe2 $(R2)

test2:
	$(SNIPPY) --outdir test2 --force --ref example.fna --R1 $(R1) --R2 $(R2)

test3:
	$(SNIPPY) --outdir test3 --force --ref example.gbk --R1 $(R1) --pe2 $(R2)

test4:
	$(SNIPPY) --outdir test4 --force --ref example.gbk --se $(R1)

bigtest:  genome.fa reads_R1.fq.gz reads_R2.fq.gz
	$(SNIPPY) --force --cpus 36 --outdir bigtest --ref $(word 1,$^) --R1 $(word 2,$^) --R2 $(word 3,$^)

genome.fa: 
	curl http://ftp.ebi.ac.uk/pub/databases/ena/wgs/public/cv/CVOW01.fasta.gz \
		| gunzip | sed -e 's/ENA|CVOW01000001|//' > $@

reads_R1.fq.gz:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/006/SRR1157476/SRR1157476_1.fastq.gz \
		| gunzip | seqtk sample 0.2 - | pigz -9 > $@
	
reads_R2.fq.gz:
	curl ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR115/006/SRR1157476/SRR1157476_2.fastq.gz \
		| gzip -d -c | seqtk sample 0.2 - | pigz -9 > $@
	

data:
	rm -f reads*
	ln -sf example.fna example.fa
	fq-simulate_illumina_reads.pl --ref example.fa --len 150 --cov 80 --indels --ambig --prefix reads
	rm -f reads.log example.fa
	pigz reads*.fastq
	
clean:
	rm -fr ./test?
	