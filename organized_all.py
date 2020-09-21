import gzip
import io

import shutil
import tre
from itertools import product, repeat
from multiprocessing import Process

from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
#from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
<<<<<<< HEAD
from itertools import product, repeat
from multiprocessing import Process
import gene_lib
import Bio
import queue
=======

import gene_lib
>>>>>>> 43f3ca09184d4218470cf7f545d6841e1e2ec1ad


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

#log(open_any("/media/sf_gene/original/B1.fasta", "rt"))
#log(open_any("/media/sf_gene/original/wt-lung_R1_001.fastq.gz", "rt"))
#log(open_any("/media/sf_gene/unified/wt-lung_unified_001.fastq.bz2", "r"))
#log(open_any("/media/sf_gene/unified/wt-lung_unified_001.fastq.bz2.blabla", "rt"))


#sys.exit()

# creates gzip file from bz2 file:
def bz2_to_gzip(in_bz2, out_gzip):
    with bz2.BZ2File('/media/sf_gene/10k_data/unified_10k.fastq.bz2', 'rb') as \
            f_in, gzip.open('/media/sf_gene/10k_data/unified_10k.fastq.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

#creates gzip file from txt file:
def txt_to_gzip(in_txt, out_gzip):
    with open('/media/sf_gene/50Mill_consolidate.fastq', 'rb') as\
            f_in, gzip.open('/media/sf_gene/50Mill_consolidate.fastq.gz', 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)

############## part 1 ##############

def filter_potential_sines_and_locations_proc(q, re, fuzziness):
    while True:
        recs = q.get()
        # log(rec)

        if recs is None:
            break

        for rec in recs:
            match = re.search(str(rec.seq), fuzziness)
            if match:
                # log(rec.seq)
                sine_location = match.groups() #returns tuple of tuples (in this case: ((2,78), ) for example

                q.put((rec, sine_location))

def filter_potential_sines_and_locations_write(q, handle_write_sine, handle_write_loc):
    while not q.empty():
        (rec, sine_location) = q.get()

        gene_record_write(rec, handle_write_sine, 'fasta')
        handle_write_loc.write(",".join([str(i) for i in sine_location[0]]) + "\n")

    handle_write_sine.flush()
    handle_write_loc.flush()


# this function gets fastq unify file (unify of R1 and R2 fastq files) in addition it gets a sine file
# (contains one of the sines B1 or B2 or B4), sine_header=67 and maxerr=14.
# for every record in the unify file it checks if it contains the current sine (the one in the sine file)
# and create a new fastq file with all the records contains the sine,
# and another file contains the sines locations (tuple of (sine_start, sine_end) for each sine)
def filter_potential_sines_and_locations(in_file_unify, in_file_sine, out_file_with_sine, out_file_location, sine_header=67, maxerr=14):
    sine = gene_lib.get_sine_forward(in_file_sine)  #"B1.fasta"
    re = tre.compile(sine[:sine_header], tre.EXTENDED)
    fuzziness = tre.Fuzzyness(maxerr=maxerr)

    # Create slave processes
    procs = []
    for _ in range(multiprocessing.cpu_count() - 3):
        # Create a communication queue between this process and slave process
        q = GeneDQueue()
        
        # Create and start slave process
        p = Process(target=filter_potential_sines_and_locations_proc, args=(q, re, fuzziness))
        p.start()

        procs.append({
            'p': p,
            'q': q,
            'batch': [],
            'write_i': 0
        })

    with open_any(in_file_unify, "rt") as handle_read, \
         open_any(out_file_with_sine, "wt") as handle_write_sine,\
         open_any(out_file_location, "wt") as handle_write_loc:


        records = gene_records_parse(handle_read)
        rec_i = 0

        for rec in tqdm(records, miniters=100):
            # Simple round-robin between the slave processes
            proc = procs[rec_i % len(procs)]

            # Add a new record into a local batch array of slave process
            proc['batch'].append(rec)

            if len(proc['batch']) >= 20:
                # Get found potential sine from slave process queue
                #
                # Optimization: 
                # Don't check the slave queue every iteration, as the check slows down.
                # Moreover we won't get a potential sine for every record.
                if proc['write_i'] > 3:
                    filter_potential_sines_and_locations_write(proc['q'], handle_write_sine, handle_write_loc)
                    proc['write_i'] = 0
                else:
                    proc['write_i'] += 1

                # Put batch of new records into slave process queue
                proc['q'].put(proc['batch'])

                # Reset local batch of slave process
                proc['batch'] = []

            # Uncomment for testing a small amount of records
            # if rec_i == 100000:
            #     break

            rec_i += 1
        
        # Cleanup slave processes
        for proc in procs:
            # Get found potential sine from slave process queue, before last batch
            filter_potential_sines_and_locations_write(proc['q'], handle_write_sine, handle_write_loc)

            # Put last batch, if avaliable
            if len(proc['batch']):
                proc['q'].put(proc['batch'])
                proc['batch'] = []
            
            # Make slave proccess terminate
            proc['q'].put(None)

            # Wait for termination
            proc['p'].join()
            
            # Get found potential sine from slave process queue, very last time
            filter_potential_sines_and_locations_write(proc['q'], handle_write_sine, handle_write_loc)


# def filter_potential_sines_and_locations(in_file_unify, in_file_sine, out_file_with_sine, out_file_location, sine_header=67, maxerr=14):
#     sine = gene_lib.get_sine_forward(in_file_sine) #"B1.fasta"
#     re = tre.compile(sine[:sine_header], tre.EXTENDED)
#     fuzziness = tre.Fuzzyness(maxerr=maxerr)
#     with open_any(in_file_unify, "rt") as handle_read:
#         records = gene_records_parse(handle_read, "fastq")
#         with open_any(out_file_with_sine, "wt") as handle_write_sine,\
#             open_any(out_file_location, "wt") as handle_write_loc:
#             for rec in tqdm(records):
#                 match = re.search(str(rec.seq), fuzziness)
#                 if match:
#                     gene_record_write(rec, handle_write_sine, 'fastq')
#                     sine_location = match.groups() #returns tuple of tuples (in this case: ((2,78), ) for example
#                     handle_write_loc.write(",".join([str(i) for i in sine_location[0]]) + "\n")

# filter_potential_sines_and_locations('/media/sf_gene/10k_data/unified_10k.fastq.gz', 'B1.fasta',
#                                      '/media/sf_gene/10k_data/unif10k_withSine.fastq.gz',
#                                      '/media/sf_gene/10k_data/unif10k_sineLocation.fastq.gz')



# this function gets a len of barcode and two files:
# 1. fastq file contains all the records which it's sequence contains sine
# 2. txt file that every row in it contains tuple (x, y) such that x and y is the begin and the end of sine
# (respectively to the rows in file 1)
# It creates a new fastq file contains all the records from file 1 so that the sequence of each record
# represents the barcode of the sine in the corresponding row.
def filter_potential_sines_barcode(sine_barcode_len, in_file_sine, in_file_location, out_file_sine_prefix):
    with open_any(in_file_sine, "rt") as handle_read_sine,\
         open_any(in_file_location, "rt") as handle_read_location,\
         open_any(out_file_sine_prefix, "wt") as handle_write_barcode:
        records = gene_records_parse(handle_read_sine)
        for rec, location in zip(tqdm(records), handle_read_location):
            sine_location = location.split(",") #gets string and delimiter (',' in this case) returns list of strings
            sine_location = [int(i) for i in sine_location]
            if(sine_location[0] >= sine_barcode_len):
                new_rec = rec[sine_location[0] - sine_barcode_len: sine_location[0]]
                gene_record_write(new_rec, handle_write_barcode)

# filter_potential_sines_barcode(36, '/media/sf_gene/10k_data/unif10k_withSine.fastq.gz',
#                                '/media/sf_gene/10k_data/unif10k_sineLocation.fastq.gz',
#                                '/media/sf_gene/10k_data/unif10k_sineBarcode.fastq.gz')

############## END of part 1 ##############



############## part 2 ##############

#==== Initial filter====#
# filter_potential_new_sines_prefix
# this function do an initial filtering (with d=0) to the sines prefixes (barcodes..)
# in order to distinguish between new sines and inherited sines.
# it creates a new file containing the sine suspicious as new and another new file contains the
# inherited sines
def new_SINES_Initial_filter_rec(in_file_sine_prefix, out_file_potential_new, out_file_inherited):
    with open_any(in_file_sine_prefix, "rt") as handle_read_prefix,\
         open_any(out_file_potential_new, "wt") as handle_potential_new,\
         open_any(out_file_inherited, "wt") as handle_write_inherited:

        records = gene_records_parse(handle_read_prefix)
        new_dict = {}
        for rec in tqdm(records):
            prefix = str(rec.seq)
            if prefix in new_dict:
                new_dict[prefix].append(rec)
            else:
                new_dict[prefix] = [rec]

        for val in tqdm(new_dict.values()):
            if len(val) == 1:
                gene_record_write(val[0], handle_potential_new)
            else: #len(val) = 1
                for rec in val:
                    gene_record_write(rec, handle_write_inherited)

# new_SINES_Initial_filter_rec('/media/sf_gene/10k_data/unif10k_sineBarcode.fastq.gz',
#                              '/media/sf_gene/10k_data/unif10k_potentialNewSINE.fastq.gz',
#                              '/media/sf_gene/10k_data/unif10k_inheritedSINE.fastq.gz')


#==== list of barcode parts====#
# this function gets record and the length of record's part and
# returns list with all the parts of the record.
def barcode_parts(record, part_len):
    barcode = record.seq
    barcode_len = len(barcode)

    assert barcode_len % part_len == 0
    num_of_parts = int(barcode_len / part_len)

    for i in range(0,num_of_parts-1):
        rec_part = record[part_len * i : part_len * (i + 1)] #todo: insted rec_part write seq_part ?
        yield rec_part


#==== list of barcode parts====#
# this function gets record and the length of record's part and
# yield all the parts of length "part_len" exist in the record.
def barcode_wins(record, part_len):
    barcode = record.seq
    barcode_len = len(barcode)
    num_of_winds = barcode_len - part_len + 1
    for i in range(0,num_of_winds-1):
        rec_wind = record[i : part_len + i]
        yield rec_wind


#==== dictionary build====#
# This function build a main_dictionary that will be used in the second filtering step.
# The main_dictionary keys length is "main_key_len" = barcode_len / d+1
# (where d is the "maximum error" - above it, two barcodes will be considered as different).
# this key represents one possible part of a barcode (out of all the possible options to get DNA subsequence of main_key_len)
# The main_dictionary's value is another dictionary - it's keys (secondary key) are all the barcodes
# containing the main_key as a substring,
# and its value is the barcode id (related to it's record).
def build_dictionary(in_file_prefix, out_file_dict, sine_barcode_len = 36, maxerr = 3):#main_key_len=9):
    assert sine_barcode_len % (maxerr + 1) == 0
    main_key_len = int(sine_barcode_len / (maxerr + 1))
    main_dict = {}

    print_step("Start build_dictionary: product of 'ATCGN'")
    for tuple in tqdm(product('ATCGN', repeat=main_key_len)): #product returns iterator
        #log(tuple)
        main_key = "".join([str(x) for x in tuple]) #without str() ??
        #log(main_key)
        main_dict[main_key] = {}
    
    print_step("Start build_dictionary: fill with records")
    with open_any(in_file_prefix, "rt") as handle_read_prefix:
        records = gene_records_parse(handle_read_prefix)
        for rec in tqdm(records):
            # barcode_parts_list = barcode_parts(rec, main_key_len)
            # for rec_part in barcode_parts_list:
            str_barc = str(rec.seq)
            for rec_part in barcode_wins(rec, main_key_len):
                str_barc_part = str(rec_part.seq)
                sec_dict = main_dict[str_barc_part]
                sec_dict[str_barc] = rec.id
    
    print_step("Start build_dictionary: write to file")
    with open_any(out_file_dict, "wb") as handle_dict:
        pickle.dump(main_dict, handle_dict, protocol=pickle.HIGHEST_PROTOCOL)
        handle_dict.flush()
        
    print_step("Start build_dictionary: done")

#dict = build_dictionary('/media/sf_gene/10k_data/unif10k_sineBarcode.fastq.gz')




#==== dictionary build====#
#the same dictionary, only this one save in the "sec_dict" a list of all the id with the same barcode.
#in_file_prefix- the barcodes, out_file_dict- the dictionary-empty at first.
def build_dictionary_for_histogram(in_file_prefix, out_file_dict, sine_barcode_len = 36, maxerr = 3):#main_key_len=9):
	assert sine_barcode_len % (maxerr + 1) == 0
	main_key_len = int(sine_barcode_len / (maxerr + 1))
	main_dict = {}

	print_step("Start build_dictionary: product of 'ATCGN'")
	for tuple in tqdm(product('ATCGN', repeat=main_key_len)): #product returns iterator
		#log(tuple)
		main_key = "".join([str(x) for x in tuple]) #without str() ??
		#log(main_key)
		main_dict[main_key] = {}

	print_step("Start build_dictionary: fill with records")
	with open_any(in_file_prefix, "rt") as handle_read_prefix:
		records = gene_records_parse(handle_read_prefix)
		for rec in tqdm(records):
			# barcode_parts_list = barcode_parts(rec, main_key_len)
			# for rec_part in barcode_parts_list:
			str_barc = str(rec.seq)
			for rec_part in barcode_wins(rec, main_key_len):
				str_barc_part = str(rec_part.seq)
				sec_dict = main_dict[str_barc_part]
				if sec_dict.get(str_barc) is None:
					sec_dict[str_barc] = [rec.id]
				else:
					sec_dict[str_barc].append(rec.id)
					

	print_step("Start build_dictionary: write to file")
	with open_any(out_file_dict, "wb") as handle_dict:
		pickle.dump(main_dict, handle_dict, protocol=pickle.HIGHEST_PROTOCOL)
		handle_dict.flush()
		
	print_step("Start build_dictionary: done")

	#dict = build_dictionary('/media/sf_gene/10k_data/unif10k_sineBarcode.fastq.gz')




#====check barcode match to its inner dict barcodes====#
# this function gets the barcode, its id, a dictionary contains all the records we want to check match with,
# and a maximum error argument for the check
# it will returns True if the barcode matches at list one of the barcodes in the dictionary.
def is_match_barcodes(sec_dict, barcode_id, re, fuzziness):
    for key, val in sec_dict.items():
        if barcode_id != val:
            match = re.search(str(key), fuzziness)
            if match:
                return True
    return False


#====check barcode match to its inner dict barcodes====#
# this function gets the barcode, its id, a dictionary contains all the records we want to check match with,
# and a maximum error argument for the check
# it update the match list
def is_match_barcodes_hist(sec_dict, barcode_id, re, fuzziness, match, lenght):
	for key, val in sec_dict.items():
		if len(match) == lenght:
			return
		if re.search(str(key), fuzziness):
			if ((val[0] in match) == False):
				match.extend(val)
				
				#log(len(val),len(match))
				


#====check matching of two string by tre====#
# this function gets two strings and maximum error value and compare the strings using the tre functions
# If found a match returns True else return False
#     def is_match_tre(str1, str2, maxerr):
#         assert len(str1) == len(str2)
#         re = tre.compile(str1, tre.EXTENDED)
#         fuzziness = tre.Fuzzyness(maxerr=maxerr)
#         match = re.search(str2, fuzziness)
#         if match:
#             return True
#         return False

def new_SINES_filter_proc(q, main_dict, key_size, fuzziness):
    while True:
        recs = q.get()
        # log(rec)

        if recs is None:
            q.put(None)
            break

        for rec in recs:
            str_barc = str(rec.seq)
            re = tre.compile(str_barc, tre.EXTENDED)
            barc_parts_list = barcode_parts(rec, key_size)
            match = False
            
            for rec_part in barc_parts_list:
                if is_match_barcodes(main_dict[str(rec_part.seq)], rec.id, re, fuzziness):
                    match = True
                    break

            q.put((rec, match))
    
    log("Slave process exited")
	
	
# the same as the previous function,
# only here the match is a list of all the barcodes id that close to the barcode
<<<<<<< HEAD
def new_SINES_filter_proc_histogram(recs, main_dict, noDuplicate, key_size, fuzziness, distribution_of_neighbors, lenght):
	
	with open_any(noDuplicate, "wt") as handle_noDuplicate:
=======
def new_SINES_filter_proc_histogram(q, main_dict, key_size, fuzziness):
    while True:
        recs = q.get()
        # log(rec)

        if recs is None:
            q.put(None)
            break

        for rec in recs:
            str_barc = str(rec.seq)
            re = tre.compile(str_barc, tre.EXTENDED)
            barc_parts_list = barcode_parts(rec, key_size)
            match = []
            
            for rec_part in barc_parts_list:
                is_match_barcodes_hist(main_dict[str(rec_part.seq)], rec.id, re, fuzziness, match)
                    
                    

            q.put((rec, match))
    
    log("Slave process exited")	
>>>>>>> 43f3ca09184d4218470cf7f545d6841e1e2ec1ad

		count = 0
		for rec in recs:
			str_barc = str(rec.seq)
			re = tre.compile(str_barc, tre.EXTENDED)
			barc_parts_list = barcode_parts(rec, key_size)
			match = []
			
			for rec_part in barc_parts_list:
				is_match_barcodes_hist(main_dict[str(rec_part.seq)], rec.id, re, fuzziness, match, lenght)
					
			count = count + 1		
			if count % 100000 == 0 :
				print_step(count)
			
			if len(match) == 1:
				gene_record_write(rec, handle_noDuplicate)
				
			if(len(match)>= lenght):
				distribution_of_neighbors[lenght-1] = distribution_of_neighbors[lenght-1] + 1
			else:
				distribution_of_neighbors[len(match)-1] = distribution_of_neighbors[len(match)-1] + 1
				
		print("Slave process exited")	
		

def new_SINES_filter_write(q, handle_write_inherited, handle_write_new, wait_none=False):
    while not q.empty() or wait_none:
        obj = q.get()

        if obj is None:
            assert wait_none
            break

        (rec, match) = obj

        if match:
            gene_record_write(rec, handle_write_inherited)
        else:
            gene_record_write(rec, handle_write_new)

    # handle_write_inherited.flush()
    # handle_write_new.flush()


# 
def update_distribution(q, distribution_of_neighbors, wait_none: bool=False):
	while not q.empty() or wait_none:
		obj = q.get()

		if obj is None:
			assert wait_none
			break

		(rec, match) = obj
		if(len(match)>= 50):
			distribution_of_neighbors[49] = distribution_of_neighbors[49] + 1
		else:
			distribution_of_neighbors[len(match)-1] = distribution_of_neighbors[len(match)-1] + 1# Can create an exception, change the size
		


#====new sines filter====#
# this function do a second filtering (with d = maxerr) to the sines prefixes.
# in order to distinguish between new sines and inherited sines.
# it creates a 2 new files: one containing the new sines and the second containing the inherited sines.
def new_SINES_filter(in_file_initial_filtering, out_file_new_SINES, out_file_inherited_SINES,
                     main_dict, key_size=9, maxerr=3):

    fuzziness = tre.Fuzzyness(maxerr=maxerr)

    # Create slave processes
    procs = []
    for _ in range(multiprocessing.cpu_count() - 3):
        # Create a communication queue between this process and slave process
        q = GeneDQueue()
        
        # Create and start slave process
        p = Process(target=new_SINES_filter_proc, args=(q, main_dict, key_size, fuzziness))
        p.start()

        procs.append({
            'p': p,
            'q': q,
            'batch': [],
            'write_i': 0
        })

    with open_any(in_file_initial_filtering, "rt") as handle_read_initial_filtering,\
         open_any(out_file_new_SINES, "wt") as handle_write_new,\
         open_any(out_file_inherited_SINES, "wt") as handle_write_inherited:

        records = gene_records_parse(handle_read_initial_filtering)
        rec_i = 0
        for rec in tqdm(records):
            # Simple round-robin between the slave processes
            proc = procs[rec_i % len(procs)]
            # Add a new record into a local batch array of slave process
            proc['batch'].append(rec)

            if len(proc['batch']) >= 10:
                new_SINES_filter_write(proc['q'], handle_write_inherited, handle_write_new)

                # Put batch of new records into slave process queue
                proc['q'].put(proc['batch'])

                # Reset local batch of slave process
                proc['batch'] = []

            # Uncomment for testing a small amount of records
            # if rec_i == 100000:
            #     break

            rec_i += 1
        
        print_step("cleanup")
        
        # Cleanup slave processes
        for proc in procs:
            # Get found potential sine from slave process queue, before last batch
            new_SINES_filter_write(proc['q'], handle_write_inherited, handle_write_new)

        for proc in procs:
            # Put last batch, if avaliable
            if len(proc['batch']):
                proc['q'].put(proc['batch'])
                proc['batch'] = []
            
        for proc in procs:
            # Make slave proccess terminate
            proc['q'].put(None)
            
        for proc in procs:
            # Get found potential sine from slave process queue, very last time
            new_SINES_filter_write(proc['q'], handle_write_inherited, handle_write_new, wait_none=True)

        for proc in procs:
            # Wait for termination
            proc['p'].join()
            

#in_file_initial_filtering- the barcodes, main_dict- the dictionary, distribution_of_neighbors- list
def new_SINES_filter_for_histogram(in_file_initial_filtering, main_dict, noDuplicate, distribution_of_neighbors, lenght, key_size=9, maxerr=3):

	fuzziness = tre.Fuzzyness(maxerr=maxerr)

	# Create slave processes
	

	with open_any(in_file_initial_filtering, "rt") as handle_read_initial_filtering:

		records = gene_records_parse(handle_read_initial_filtering)
		#q = queue.Queue()
		new_SINES_filter_proc_histogram(records, main_dict, noDuplicate, key_size, fuzziness, distribution_of_neighbors, lenght)
		
		
		#update_distribution(q, distribution_of_neighbors)# למה?



# new_SINES_filter('/media/sf_gene/10k_data/unif10k_potentialNewSINE.fastq.gz',
#                  '/media/sf_gene/10k_data/unif10k_NewSINE.fastq.gz',
#                  '/media/sf_gene/10k_data/unif10k_inheritedSINE(2).fastq.gz', dict)



def SINES_new_or_inherited(in_file_dict,
                           in_file_initial_filtering, out_file_new_SINES, out_file_inherited_SINES):                     

    print_step("Start SINES_new_or_inherited: load dict")
    with open_any(in_file_dict, "rb") as handle_dict:
        dict = pickle.load(handle_dict)

    print_step("Start new_SINES_filter")
    new_SINES_filter(in_file_initial_filtering, out_file_new_SINES, out_file_inherited_SINES, dict)



# in_file_dict- the dictionary, in_file_initial_filtering - the barcodes, distribution_of_neighbors- list for counting the neighbors of barcods.
def SINES_histogram_of_neighbors(in_file_dict, in_file_initial_filtering,noDuplicate, distribution_of_neighbors, lenght):
	print_step("Start SINES_new_or_inherited for histogram: load dict")
	with open_any(in_file_dict, "rb") as handle_dict:
		dict = pickle.load(handle_dict)

	print_step("Start new_SINES_filter")
	new_SINES_filter_for_histogram(in_file_initial_filtering, dict, noDuplicate, distribution_of_neighbors,lenght)

#activate the last lines to create a graph
<<<<<<< HEAD
def print_histogram(distribution_of_neighbors, name, lenght, mode = 5):
	
	if mode == 5:
		str1 = name + "_distribution_of_neighbors"+str(lenght)+".txt"
		dist = open(str1, "wt")
		dist.write(','.join([str(elem) for elem in distribution_of_neighbors]))
	
		dist.close()
	percetage = [0.0]*lenght
	normal1 =  [0.0]*lenght
	normal2 = [0.0]*(lenght-1)
	for i in range(len(distribution_of_neighbors)):
			
		normal1[i] = distribution_of_neighbors[i]/(i+1)
		if i < (lenght-1):
			normal2[i] = distribution_of_neighbors[i]/(i+1)
					
	s = sum(distribution_of_neighbors)
	s1 = sum(normal1)
	s2 = sum(normal2)
		
	for i in range(lenght):
		percetage[i] = (distribution_of_neighbors[i]/s) * 100
		normal1[i] = (normal1[i]/s1) * 100
		if(i<(lenght -1)):
			normal2[i] = (normal2[i]/s2) * 100
	print("distribution_of_neighbors:", distribution_of_neighbors)
	print("percetage:", percetage)
	print("normal all: ", normal1)
	print("normal -1: ", normal2)
=======
def print_histogram(distribution_of_neighbors):
	log(distribution_of_neighbors)
	#indices = np.arange(len(distribution_of_neighbors))
	#plt.bar(indices, distribution_of_neighbors)
>>>>>>> 43f3ca09184d4218470cf7f545d6841e1e2ec1ad
	
	

	
	


#SINES_new_or_inherited('/media/sf_gene/10k_data/unif10k_sineBarcode.fastq.gz',
#                       '/media/sf_gene/10k_data/unif10k_potentialNewSINE.fastq.gz',
#                       '/media/sf_gene/10k_data/unif10k_NewSINE.fastq.gz',
#                       '/media/sf_gene/10k_data/unif10k_inheritedSINE(2).fastq.gz')

def run_part_1(in_file, B_file, out_dir):
	file_ext = None
	for ext in ['.fastq', '.fastq.gz', '.fastq.bz2']:
		if in_file.endswith(ext):
			file_ext = ext
			break

	assert file_ext != None, "Unknown file extension in %s" % (in_file)

	file_base = out_dir + '/' + os.path.basename(in_file[:-len(file_ext)])

	print_step("file_base = %s, file_ext = %s" % (file_base, file_ext))
	print_step()

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	print_step("Start filter_potential_sines_and_locations")
	filter_potential_sines_and_locations(in_file, B_file,
										 file_base + '_withSine' + file_ext,
										 file_base + '_sineLocation' + file_ext)

#lenght- size of histogram
def run_all(in_file, B_file, out_dir, mode = 3, lenght = 50):
	
		
	file_ext = None
	for ext in ['.fastq', '.fastq.gz', '.fastq.bz2']:
		if in_file.endswith(ext):
			file_ext = ext
			break

	assert file_ext != None, "Unknown file extension in %s" % (in_file)

	file_base = out_dir + '/' + os.path.basename(in_file[:-len(file_ext)])

	print_step("file_base = %s, file_ext = %s" % (file_base, file_ext))
	print_step()

	if not os.path.exists(out_dir):
		os.makedirs(out_dir)

	# part 0 - detect potential sines
	if (mode == 1):
		print_step("Start filter_potential_sines_and_locations")
		filter_potential_sines_and_locations(in_file, B_file,
											 file_base + '_withSine' + file_ext,
											 file_base + '_sineLocation' + file_ext)
		return

	# part 1 - detect barcodes of potential sines
	if (mode == 2): 
		print_step("Start filter_potential_sines_barcode")
		filter_potential_sines_barcode(36, file_base + '_withSine' + file_ext,
										   file_base + '_sineLocation' + file_ext,
										   file_base + '_sineBarcode' + file_ext)
		return 
		

	#create the dictionary
	if(mode == 4):
		print_step("Start build_dictionary")
		build_dictionary_for_histogram(file_base + '_sineBarcode' + file_ext,
								file_base + '_mainDict' + file_ext)

		return
	# part 2 - identify new sines
	
	if mode == 5:#after we have the dictionary
		distribution_of_neighbors = [0]*lenght
		print_step("Start SINES_new_or_inherited histogram")
		SINES_histogram_of_neighbors(file_base + '_mainDict' + file_ext,
									 file_base + '_sineBarcode' + file_ext,
									 file_base + '_NewSINE' + file_ext,
									 distribution_of_neighbors, lenght)
		print_histogram(distribution_of_neighbors, file_base, lenght)
		return
	

	#to print the distribution_of_neighbors data
	if mode == 6:
		
		file_handle = open(file_base + '_distribution_of_neighbors'+str(lenght)+'.txt', "rt")
		numbers = [0]*lenght
		for line in file_handle:
			line = line.rstrip('\n')
		
			numbers = list(line.split(","))
			numbers = [int(i) for i in numbers]
			break
		print_histogram(numbers, file_base, lenght, mode)
		return
		
	
	print_step("Start build_dictionary")
	build_dictionary(file_base + '_sineBarcode' + file_ext,
							file_base + '_mainDict' + file_ext)



	print_step("Start new_SINES_Initial_filter_rec")
	new_SINES_Initial_filter_rec(file_base + '_sineBarcode' + file_ext,
								 file_base + '_potentialNewSINE' + file_ext,
								 file_base + '_inheritedSINE' + file_ext)


	print_step("Start SINES_new_or_inherited")
	SINES_new_or_inherited(file_base + '_mainDict' + file_ext,
						   file_base + '_potentialNewSINE' + file_ext,
						   file_base + '_NewSINE' + file_ext,
						   file_base + '_inheritedSINE_2' + file_ext)

	print_step("DONE ALL!")
