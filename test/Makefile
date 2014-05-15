
all:
	rm -f reads*
	fq-simulate_illumina_reads.pl --ref example.fna --len 150 --cov 40 --indels --ambig --prefix reads
	rm -f reads.log
	gzip reads*.fastq
	
