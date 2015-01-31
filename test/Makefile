
.PHONY: test1 test2 test3 test4 data

quick: test3

all: test1 test2 test3 test4

test1:
	../bin/snippy --outdir test1 --force --ref example.fna reads_R1.fastq.gz reads_R2.fastq.gz

test2:
	../bin/snippy --report --outdir test2 --force --ref example.fna reads_R1.fastq.gz reads_R2.fastq.gz

test3:
	../bin/snippy --outdir test3 --force --ref example.gbk reads_R1.fastq.gz reads_R2.fastq.gz

test4:
	../bin/snippy --report --outdir test4 --force --ref example.gbk reads_R1.fastq.gz

data:
	rm -f reads*
	ln -sf example.fna example.fa
	fq-simulate_illumina_reads.pl --ref example.fa --len 150 --cov 80 --indels --ambig --prefix reads
	rm -f reads.log example.fa
	pigz reads*.fastq
	
clean:
	rm -fr ./test?
	