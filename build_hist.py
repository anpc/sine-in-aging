#!/usr/bin/env python3.8

import gzip
import io

import shutil
import tre
from itertools import product, repeat

from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import collections

import gene_lib

try:
    from itertools import izip as zip
except ImportError: # will be 3.x series
    pass

import bz2
import sys
from datetime import datetime
from tqdm import tqdm
import os
from gene_lib import *


#TODO: copied from organized_all.py - is there a better way? 
def barcode_parts(record, part_len):
    barcode = record.seq
    barcode_len = len(barcode)

    assert barcode_len % part_len == 0
    num_of_parts = int(barcode_len / part_len)

    for i in range(0,num_of_parts-1):
        rec_part = record[part_len * i : part_len * (i + 1)] #todo: insted rec_part write seq_part ?
        yield rec_part

def count_neighbors(main_dict, maxerr, rec):
    fuzziness = tre.Fuzzyness(maxerr=maxerr)
    key_size = 9
    str_barc = str(rec.seq)
    re = tre.compile(str_barc, tre.EXTENDED)
    barc_parts_list = barcode_parts(rec, key_size)
    
    neighbors = []
    for rec_part in barc_parts_list:
        str_barc_part = str(rec_part.seq)
        sec_dict = main_dict[str_barc_part]
        for key, val in sec_dict.items():
            if re.search(str(key), fuzziness):
                # val is a list of record id's having exactly the same barcode.
                # It suffices, because we just want tonot count the same barcode in different buckets. The barcode is
                # inserted with the same list of id's into every window
                if ((val[0] in neighbors) == False):
                    neighbors += val

    return len(neighbors)	

# This is copied and adapted from organized_all - may be a better way in terms of modularity
def write_hist_file():
    maxerr = 3
    
    file_ext = None
    for ext in ['.fastq', '.fastq.gz', '.fastq.bz2']:
        if candidate_barcodes.endswith(ext):
            file_ext = ext
            break
    assert file_ext != None, f"Unknown file extension in {candidate_barcodes}"
    file_base = out_dir + '/' + os.path.basename(candidate_barcodes[:-len(file_ext)])

    # We assume for now all records neighbors. Generally, need to verify.
    log('About to build hist for',candidate_barcodes)

    # histogram {number of neighbors: count of barcodes that had this # of neighbors}
    distribution_of_neighbors = collections.Counter()

    # code copied from organized_all.py
    print_step("Start gen_new_or_inherited for histogram: load dict")
    with open_compressed(global_dict, "rb") as handle_dict:
        global_dict_obj = pickle.load(handle_dict)

    with open_compressed(candidate_barcodes, "rt") as handle_read_initial_filtering:
        records = gene_records_parse(handle_read_initial_filtering)
        # track records inserted into the histogram        
        count = 0
        for rec in records:
            num_neighbors = count_neighbors(global_dict_obj, maxerr, rec)
            distribution_of_neighbors[num_neighbors] += 1 		
            count += 1		
            if count % 5000 == 0 :
                print_step(f"{candidate_barcodes} {count} completed")
                 

    # save histogram 
    with open_compressed( file_base + '_histogram', "wb") as handle_hist:
        pickle.dump(distribution_of_neighbors, handle_hist, protocol=pickle.HIGHEST_PROTOCOL)
    return

### MAIN ###
[candidate_barcodes, global_dict, out_dir] = sys.argv[1:]
write_hist_file()
