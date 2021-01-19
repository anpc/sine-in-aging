#!/usr/bin/env python3.8
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

#current file generated with 15, 4
common_req = 20
init_err = 4
tot_err = 4
padding = 0
step = 10000

def merged_paired_ends(records1, records2):
    tot_good = 0
    tot_great = 0
    tot = 0
#    log('in merged_paired_ends',records1,records2)
    for (rec1,rec2) in zip(records1, records2):
        tot += 1
        str1 = str(rec1.seq)
        str2 = str(rec2.seq.reverse_complement())
#        log('-------------------------------------------\n matching ',str1,'\n',str2,'\n===================================================')
        end1 = str1[-common_req:]
        re = tre.compile(end1, tre.EXTENDED)
        # we expect small errors here
        res_seq = None
        match = re.search(str2, tre.Fuzzyness(maxerr = init_err))
        if match:
            tot_good += 1
            match_loc = match.groups()[0][0]
            to_search_len = match_loc + common_req
            fuzzyness = max(tot_err, ceil(0.1*to_search_len))
            re = tre.compile(str1[-to_search_len :], tre.EXTENDED)
            match_tot = re.search(str2, tre.Fuzzyness(maxerr = fuzzyness))
 #           log('step1: matched ',end1,' at',match_loc,' testing prefix ',str2[:to_search_len],'cost ',match.cost)
            if match_tot:
            #    if (tot_good % 100 == 0):
            #        log('fuzzyness = ', fuzzyness)
  #              log('step2: matched ',str1[-to_search_len:],' at',match_tot.groups()[0][0],' testing prefix ','cost ',match.cost)
                tot_great += 1
                # An arbitrary decision: take the common string from r2
                res_str = str1[ : -to_search_len] + str2
                # TODO: preserve qualities
                res_seq = SeqRecord(Seq(res_str), id=rec1.id, name=rec1.name, description = rec1.description, letter_annotations =
                                    {"phred_quality":[30 for i in range(len(res_str))]})
                if (tot_great % step == 0):
                    log('nicely matched ',str1,'\n',str2, 
                        to_search_len, match_tot.group(0), match.group(0), match_tot.cost, match.cost)
   #             log('result = ',str(res_seq.seq))
                yield res_seq
                continue

        res_str = str1 + ('N' * padding) + str2
        res_seq = SeqRecord(Seq(res_str), id=rec1.id, name=rec1.name,
                   description = rec1.description, letter_annotations = {"phred_quality":[30 for i in range(len(res_str))]}) 
        if  (tot % step == 0):
            log(tot, tot_good, tot_great)
           # log('matched ',str1,'\n',str2, len(str1), len(str2))
  #      log('result = ',str(res_seq.seq))
        yield res_seq


### MAIN ###
[in_fname1, in_fname2] = sys.argv[1:]
def write_merged_to_stdout():
    # We assume for now all records match. Generally, need to verify.
    log('About to merge ',in_fname1, in_fname2)
    with gene_lib.open_any(in_fname1, 'rt') as in_f1_handle:
        records1 = SeqIO.parse(in_f1_handle, format="fastq")
        with gene_lib.open_any(in_fname2, 'rt') as in_f2_handle:
            records2 = SeqIO.parse(in_f2_handle, format="fastq")
            merged = merged_paired_ends(records1, records2)
            Bio.SeqIO.write(merged, handle=sys.stdout, format='fastq')

write_merged_to_stdout()
