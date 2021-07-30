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
    """Returns non-overlapping windows."""
    barcode = record.seq
    barcode_len = len(barcode)

    assert barcode_len % part_len == 0
    num_of_parts = int(barcode_len / part_len)

    for i in range(0,num_of_parts-1):
        rec_part = record[part_len * i : part_len * (i + 1)] #todo: insted rec_part write seq_part ?
        yield rec_part

def count_neighbors(main_dict, maxerr, rec):
    """Counts how many neighbors `rec` has.
    
    Note that for parllelization this script may process a subset of `rec`
    but the full `main_dict` so this gives the correct total for this `rec`.
    """
    fuzziness = tre.Fuzzyness(maxerr=maxerr)
    key_size = 9
    str_barc = str(rec.seq)
    #TODO: we assume maxerr = 3 here. Appending the 3+1 Y's will ensure anchoring from the end - $ seems to be buggy
    re = tre.compile('^'+str_barc + '1234', tre.EXTENDED)
    barc_parts_list = barcode_parts(rec, key_size)
    
    neighbors = set()
    dist_hist = [0] * 4
    #print('considering id '+rec.id)
    for rec_part in barc_parts_list:
        str_barc_part = str(rec_part.seq)
        sec_dict = main_dict[str_barc_part]
        for key, val in sec_dict.items():
            m = re.search(str(key)+'1234', fuzziness)
            if m:
                # val is a set of record id's having exactly the same barcode.
                # It suffices, because we just want to not count the same barcode in different buckets. The barcode is
                # inserted with the same list of id's into every window
                for v in val:
                    if ((v in neighbors) == False):
                        dist_hist[m.cost] += 1
                     #   print('adding neighbor '+v)
                     #   print('with cost ',m.cost, (m.numdel, m.numins, m.numsub))
                    neighbors.add(v)

    if len(neighbors) > 1000:
        print(str_barc+' has ',dist_hist)
    return len(neighbors)

# This is copied and adapted from organized_all - may be a better way in terms of modularity
def write_hist_file():
    maxerr = 3
    
    file_ext = None
    for ext in ['.fastq', '.fastq.gz', '.fastq.bz2', '.fasta', '.fasta.gz', '.fasta.bz2']:
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
    with open(file_base + '_histogram.pickle', "wb") as handle_hist:
        pickle.dump(distribution_of_neighbors, handle_hist, protocol=pickle.HIGHEST_PROTOCOL)
        print("Saved histogram to", file_base + '_histogram.pickle')
    return

### MAIN ###
if __name__ == '__main__':
    [candidate_barcodes, global_dict, out_dir] = sys.argv[1:]
    write_hist_file()
