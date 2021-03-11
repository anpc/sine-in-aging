#!/usr/bin/env python3.8

import sys

from organized_all import build_dictionary_for_histogram

# TODO decide what command line params we need
[file_base] = sys.argv[1:]
build_dictionary_for_histogram(file_base + '_sineBarcode.fasta.gz',
                               file_base + '_mainDictHistogram.pickle.gz')
