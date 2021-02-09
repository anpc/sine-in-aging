#!/usr/bin/env python3.8

import sys

from organized_all import build_dictionary_for_histogram

# TODO decide what command line params we need
[file_base, in_file_ext, sine_name] = sys.argv[1:]
build_dictionary_for_histogram(file_base + '_sineBarcode' + in_file_ext,
                               file_base + '_mainDictHistogram.pickle.gz')
