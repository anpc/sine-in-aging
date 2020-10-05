#!/usr/bin/env python3.6
import math
import tre
import sys
import gzip
import sines_io
import gene_lib
from gene_lib import log
from sines_io import fastq_zst_records
import math
from math import ceil

import os

# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord


### MAIN ###
[in_fname] = sys.argv[1:]
def write_fastq_to_stdout():
	log('About to transform ',in_fname)
	with gene_lib.open_any(in_fname, 'rt') as in_f1_handle:
		records = SeqIO.parse(in_f1_handle, format="fasta")
		
		for (rec) in records:
			string = str(rec.seq)
			res_seq = SeqRecord(Seq(string), id=rec.id, name=rec.name, description = rec.description, letter_annotations =
										{"phred_quality":[30 for i in range(len(string))]})
			Bio.SeqIO.write(res_seq, handle=sys.stdout, format='fastq')
	
write_fastq_to_stdout()
