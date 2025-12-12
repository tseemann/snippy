#!/usr/bin/env python3

#Written by Nick Loman (@pathogenomenick)

import sys
import csv
from operator import itemgetter
import argparse

def make_tag(snp):
	return "%s-%s-%s-%s-%s" % (snp['CHROM'], snp['POS'], snp['TYPE'], snp['REF'], snp['ALT'])

def read_snp_table(d):
	fh = open(d+'/snps.tab')
	snps = []
	snps = list(csv.DictReader(fh, dialect='excel-tab'))
	evidence = {}
	for snp in snps:
		tag = make_tag(snp)
		evidence[tag] = snp['EVIDENCE']
	return snps, evidence

def go(directories):
	evidence = {}	
	all_snps = []
	samples = []
	fields = csv.DictReader(open(directories[0]+'/snps.tab'), dialect='excel-tab').fieldnames
	for d in directories:
		snps, sample_evidence = read_snp_table(d)
		evidence[d] = sample_evidence
		all_snps.extend(snps)
		samples.append(d)

	fields.extend(samples)
	fields.append('number_samples_with_variant')

	all_snps.sort(key=itemgetter('POS'))

	seen = set()

	out = csv.DictWriter(sys.stdout, fieldnames=fields, dialect='excel-tab')
	out.writeheader()
	for snp in all_snps:
		tag = make_tag(snp)
		if tag in seen:
			continue
		seen.add(tag)

		evidence_count = 0
		for sample in samples:
			if tag in evidence[sample]:
				snp[sample] = evidence[sample][tag]
				evidence_count += 1

		snp['number_samples_with_variant'] = str(evidence_count)

		out.writerow(snp)


parser = argparse.ArgumentParser(description='Process a set of snippy SNP output files into a single report.')
parser.add_argument('directory', metavar='directory', nargs='+', help='list of snippy output directories')
args = parser.parse_args()

go(args.directory)

