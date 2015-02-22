
.PHONY: test1 test2 test3 test4 data

quick: test3

all: test1 test2 test3 test4

test1:
	../bin/snippy --outdir test1 --force --ref example.fna --pe1 reads_R1.fastq.gz --pe2 reads_R2.fastq.gz

test2:
	../bin/snippy --outdir test2 --force --ref example.fna --R1 reads_R1.fastq.gz --R2 reads_R2.fastq.gz

test3:
	../bin/snippy --outdir test3 --force --ref example.gbk --R1 reads_R1.fastq.gz --pe2  reads_R2.fastq.gz

test4:
	../bin/snippy --outdir test4 --force --ref example.gbk --se reads_R1.fastq.gz

data:
	rm -f reads*
	ln -sf example.fna example.fa
	fq-simulate_illumina_reads.pl --ref example.fa --len 150 --cov 80 --indels --ambig --prefix reads
	rm -f reads.log example.fa
	pigz reads*.fastq
	
clean:
	rm -fr ./test?
	