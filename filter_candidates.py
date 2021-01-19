#!/usr/bin/env python3.8

import sys
import tre

import Bio
from Bio import SeqIO
#from Bio.Seq import Seq
#from Bio.Alphabet import IUPAC
#from Bio.SeqRecord import SeqRecord

import gene_lib
from gene_lib import log
from gene_lib import get_sine_forward

def filter_potential_sines(in_fname, sine_string, sine_header=67, maxerr=19, reverse_complement=False):
    """
    Finds candidate SINEs with a certain distance from a prefix length.
    To be used for preliminary screening (input for later steps).
    """
    with gene_lib.open_any(in_fname, 'rt') as in_file_handle:
        records = SeqIO.parse(in_file_handle, format="fastq")
        re = tre.compile(sine[:sine_header], tre.EXTENDED)
        fuzziness = tre.Fuzzyness(maxerr=maxerr)

        for rec in records:
            if reverse_complement:
                cur_seq = rec.seq.reverse_complement()
            else:
                cur_seq = rec.seq

            match = re.search(str(cur_seq), fuzziness)
            if match:
                # log(rec.seq)
                #sine_location = match.groups() #returns tuple of tuples (in this case: ((2,78), ) for example
                SeqIO.write(rec, sys.stdout, 'fastq')

# Writes to stdout, uncompresed
[sine_fname, header_len, max_error, reverse_complement, merged_input_fname] = sys.argv[1:]
if reverse_complement not in {"forward", "rc"}:
    raise ValueError('reverse_complement arg must be "forward" or "rc"')

log(f"About to screen {merged_input_fname} ({reverse_complement}) for {sine_fname} first {header_len} up to {max_error} err")
sine = gene_lib.get_sine_forward(sine_fname) #"B1.fasta"
filter_potential_sines(
    in_fname=merged_input_fname,
    sine_string=sine,
    sine_header=int(header_len),
    maxerr=int(max_error),
    reverse_complement=(reverse_complement == "rc"))
