import gzip
import io
import gzip
import shutil
import tre
#import matplotlib.pyplot as plt
import networkx as nx
from networkx import Graph
from itertools import product, repeat
from multiprocessing import Process
import sys
from resource import getrusage, RUSAGE_SELF
import numpy as np
import nltk


from Bio import SeqIO
from Bio import Align
from Bio import pairwise2
from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC

import gene_lib

def showResult(file_centers,in_file_sine,sine_header=67, maxerr=19):
    sum = 0
    hist = {}
    sine = gene_lib.get_sine_forward(in_file_sine)  # "B1.fasta"
    re = tre.compile(sine[:sine_header], tre.EXTENDED)
    stringSine=sine
    print ('original sine',stringSine)
    fuzziness = tre.Fuzzyness(maxerr=maxerr)
    with open(file_centers, "r") as centerFile:
        for line in centerFile:
            currentLine = line.strip()
#            re2 = tre.compile(currentLine, tre.EXTENDED)
#            match = re2.search(stringSine, fuzziness)
            match = re.search(currentLine, fuzziness)
            sine_location=match.groups()
#            print (sine_location)
#            print ('current center', currentLine)
#            print ('match sine', str(sine[sine_location[0][0] :sine_location[0][1]]))
#            print ('current center',  nltk.edit_distance(sine[sine_location[0][0] :sine_location[0][1]],currentLine))
            hist[nltk.edit_distance(stringSine[sine_location[0][0] :sine_location[0][1]],currentLine)] = hist.get(nltk.edit_distance(stringSine[sine_location[0][0] :sine_location[0][1]],currentLine), 0) + 1
            sum = sum + nltk.edit_distance(stringSine[sine_location[0][0] :sine_location[0][1]],currentLine)
        print(sum/1000)
        print(sorted(hist.items()))

