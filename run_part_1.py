#!/usr/bin/env python3.6

from organized_all import run_all
from organized_all import run_part_1
import sys

# TODO decide what command line params we need
if len(sys.argv[1:]) == 5:
	[merged_input_fname, sine_name, output_dir, mode, lenght] = sys.argv[1:]
	run_all(merged_input_fname, sine_name, output_dir, int(mode), int(lenght))
else:
	[merged_input_fname, sine_name, output_dir, mode] = sys.argv[1:]
	run_all(merged_input_fname, sine_name, output_dir, int(mode))

# run_all("/media/sf_gene/unified/wt-lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B1.fasta",
#         "/media/sf_gene/B1-wt")

# run_all("/media/sf_gene/unified/wt-lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B2.fasta",
#         "/media/sf_gene/B2-wt")

# run_all("/media/sf_gene/unified/wt-lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B4.fasta",
#         "/media/sf_gene/B4-wt")

# run_all("/media/sf_gene2/unified_old/old_lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B1.fasta",
#         "/media/sf_gene2/B1-old")

# run_all("/media/sf_gene2/unified_old/old_lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B2.fasta",
#         "/media/sf_gene2/B2-old")

# run_all("/media/sf_gene2/unified_old/old_lung_unified_001.fastq.gz",
#         "/media/sf_gene/original/B4.fasta",
#         "/media/sf_gene2/B4-old")




#run_part_1(f"/media/sf_gene/original/wt-lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B1.fasta",
#        "/media/sf_gene/B1-wt_R1")

#run_part_1("/media/sf_gene/original/wt-lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B1.fasta",
#        "/media/sf_gene/B1-wt_R2")

#run_part_1("/media/sf_gene/original/wt-lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B2.fasta",
#        "/media/sf_gene/B2-wt_R1")

#run_part_1("/media/sf_gene/original/wt-lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B2.fasta",
#        "/media/sf_gene/B2-wt_R2")

#run_part_1("/media/sf_gene/original/wt-lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B4.fasta",
#        "/media/sf_gene/B4-wt_R1")

#run_part_1("/media/sf_gene/original/wt-lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B4.fasta",
#        "/media/sf_gene/B4-wt_R2")

#run_part_1("/media/sf_gene2/original_old/old_lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B1.fasta",
#        "/media/sf_gene2/B1-old_R1")

#run_part_1("/media/sf_gene2/original_old/old_lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B1.fasta",
#        "/media/sf_gene2/B1-old_R2")

#run_part_1("/media/sf_gene2/original_old/old_lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B2.fasta",
#        "/media/sf_gene2/B2-old_R1")

#run_part_1("/media/sf_gene2/original_old/old_lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B2.fasta",
#        "/media/sf_gene2/B2-old_R2")

#run_part_1("/media/sf_gene2/original_old/old_lung_R1_001.fastq.gz",
#        "/media/sf_gene/original/B4.fasta",
#        "/media/sf_gene2/B4-old_R1")

#run_part_1("/media/sf_gene2/original_old/old_lung_R2_001.fastq.gz",
#        "/media/sf_gene/original/B4.fasta",
#        "/media/sf_gene2/B4-old_R2")
