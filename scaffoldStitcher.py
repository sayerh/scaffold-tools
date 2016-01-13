#!/usr/bin/env python
#  python scaffoldSticher.py input.fasta "upstream flanking text" "downstream flanking text"
# generates intermediate file for 

# loads fasta into memory as dictionary. 


import sys
from Bio import Seq # OO data structures supporting sequence manipulation
#from Bio import IUPAC # ACGT ACGU, YRW...etc.
from Bio import SeqIO # Input-Output, file adaptors
import re # Regular Expressions native library in python

#example fastaheader = ">IWGSC_CSS_2DL_scaff_4847783 dna:scaffold scaffold:IWGSC1.0+popseq:IWGSC_CSS_2DL_scaff_4847783:1:200:1"

#alphabet = IUPAC.unambiguous_dna

fasta = sys.argv[1]
fasta = open(fasta)
#insert a string of 100 n's between scaffolds
n=100
filler = "".join(["n"] * n)

records = SeqIO.parse(fasta,"fasta")

# Chromosome -> Scaffold name -> Sequence
chromosomes = dict() # dictionary of chromosome arms
# temporary structure
farecords = dict()


for record in records:
	chromosome = re.search( 'CSS_(...)_scaff', record.id ) # regular expression search in record.id
	if 'group' not in dir(chromosome):continue
	chromosome = chromosome.group(1)
	if chromosome not in chromosomes:
		chromosomes[chromosome] = dict()	
	chromosomes[chromosome][record.id] = record.seq

# could use fasta writer from SeqIO
#for chrom in chromosomes:
#	print ">%s_scaffolds" % chrom
#        print filler.join([str(seq) for seq in chromosomes[chrom].viewvalues()])
# JBrowse support for gtf is extended from gff may the force be with you.

# sort by Chromosome Number 1-7, Letter ABD, arm short long. Scaffolds by assigned number 9,10, etc. MT and PT come last. 

def chromSort(scaffa,scaffb):
	def sorthelper(a):
		return int(a.split("_")[-1])
	return sorthelper(scaffa) < sorthelper(scaffb)

sorted(chromosomes[chrom].viewitems(), key=lamda scaffold: int(scaffold[0].split("_")[-1]))

for chrom in chromosomes:
	start = 1
	for scaffold, seq in chromosomes[chrom].viewitems():
		start
		end = start + len(str(seq)) - 1		
		print chrom + "_scaffolds\t" + "scaffold_element" + "\t" + scaffold + "\t" + str(start) + "\t" + str(end) + "\t" + ".\t+\t.\tscaffold_id \"" + scaffold + "\""
		start = end + n + 1;
# blast scaffolds against new sequences to get scaffold names and locations or just record and print vcf
# to be continued...
