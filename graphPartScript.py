#!/usr/bin/env python3.8

from organized_all_2 import SINES_graph_of_neighbors
import sys

if len(sys.argv[1:]) == 3:
	[dict_fname, barcode_name ,part] = sys.argv[1:]
	SINES_graph_of_neighbors(dict_fname, barcode_name, int(part))

