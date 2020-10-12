#!/usr/bin/env python3.6

from showResults import showResult
import sys

if len(sys.argv[1:]) == 2:
	[centers_fname, sineName_name] = sys.argv[1:]
	showResult(centers_fname, sineName_name)
