import gzip
import bz2
from datetime import datetime
import sys, os
import itertools
from multiprocessing import connection, Lock
import pickle
import multiprocessing

# http://biopython.org/
import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

USE_BIO = False
BIO_TYPE = 'fastq'

def log(*args, **kwargs):
    """Like print() but to stderr.  Want to keep stdout clean for streaming."""
    print(*args, **kwargs, file=sys.stderr)

def print_step(line=''):
    """Print a line of text that begins with a time stamp."""
    log("[%s] %s" % (str(datetime.now())[:19], line))


# Get filename base and extenstion
#
# I.e. '/media/gene/my_genes.fastq.gz' => ('my_genes', '.fastq.gz')
def file_base_and_ext(filename):
    file_ext = None
    for ext in ['.fastq', '.fastq.gz', '.fastq.bz2']:
        if filename.endswith(ext):
            file_ext = ext
            break

    assert file_ext != None, "Unknown file extension in %s" % (filename)

    return (os.path.basename(filename[:-len(file_ext)]), file_ext)


# Open gz/bzip/fasta/fastq file with correct function according to filename ending
def open_compressed(filename, mode):
    if filename.endswith('.bz2'):
        # 't' flag is not supported in Python2 for bz2
        mode = mode.replace('t', '')
        return bz2.BZ2File(filename, mode)
    elif filename.endswith('.gz'):
        # mgzip: Allows multi-threading when using gzip
        # if sys.version_info >= (3, 0): #and 'r' in mode:
        #     import mgzip
        #     return mgzip.open(filename, mode, blocksize=10**5)
        # else:
            return gzip.open(filename, mode)
    elif filename.endswith('.fastq') or filename.endswith('.fasta'):
        return open(filename, mode)

    log("Warning: Unknown file extension in %s, treating as uncompressed" % (filename))
    return open(filename, mode)

# Replace Bio SeqRecord for performance
#
# A very simple record that is easy replacable in our code.
# Supports only needed fields from original SeqRecord
class GeneRecord:
    def __init__(self, id = None, seq = None, format = None):
        self.id = id
        self.seq = seq
        self.format = format

    def __len__(self):
        return len(self.seq)

    # Support rec[start:len] and rec[i]
    # Partly copied from SeqRec
    def __getitem__(self, index):
        if isinstance(index, int):
            return self.seq[index]
        elif isinstance(index, slice):
            if self.seq is None:
                raise ValueError("If the sequence is None, we cannot slice it.")

            return self.__class__(id = self.id,
                                  seq = self.seq[index],
                                  format = self.format)

        raise ValueError("Invalid index")


# Try to guess record format according to first line
#   - fasta: first line starts with '>'
#   - fastq: first line starts with '@'
def record_format_get(line_0):
    if line_0[0] == '>':
        return 'fasta'
    elif line_0[0] == '@':
        return 'fastq'
    
    assert False


# Replace Bio SeqIO.parse for performance
#
# file_handle - Handle after open/gzip.open/etc.
# format - 'fasta' or 'fastq', if None then will try to detect format
# 
# Performance testing for only reading the records from "fastq.gz" file, no operation on the records:
#   - Bio: 4,500 rec/sec
#	- Gene: 170,000 rec/sec
#	- Using mgzip (instead of gzip) doesn't change the numbers for this test
def gene_records_parse(file_handle, format = None):
    assert format == 'fasta' or format == 'fastq' or format == None

    # Fallback to Bio using the same API - for testing
    if USE_BIO:
        if format == None:
            format = BIO_TYPE
        for rec in SeqIO.parse(file_handle, format):
            yield rec

    line_i = 0
    rec = None

    for line in file_handle:
        line = line.rstrip('\n')
        
        if format == None and line_i == 0:
            format = record_format_get(line)

        # fasta format example:
        #   >EAS54_6_R1_2_1_413_324
        #   CCCTTCTTGTCTTCAGCGTTTCTCC
        if format == 'fasta':
            if line_i == 0:
                assert line[0] == '>'
                rec = GeneRecord()
                rec.id = line[1:].split()[0]
            elif line_i == 1:
                rec.seq = line
                line_i = -1

        # fastq format example:
        #   @EAS54_6_R1_2_1_413_324
        #   CCCTTCTTGTCTTCAGCGTTTCTCC
        #   +
        #   ;;3;;;;;;;;;;;;7;;;;;;;88
        elif format == 'fastq':
            if line_i == 0:
                assert line[0] == '@'
                rec = GeneRecord()
                rec.id = line[1:].split()[0]
            elif line_i == 1:
                rec.seq = line
            elif line_i == 2:
                assert line[0] == '+'
            elif line_i == 3:
                assert len(line) == len(rec.seq)
                line_i = -1
        
        # Done parsing a record
        if line_i == -1:
            rec.format = format
            yield rec
        
        line_i += 1


# Replace Bio SeqIO.write for performance
#
# rec - Record to write
# file_handle - Handle after open/gzip.open/etc.
# format - 'fasta' or 'fastq', if None then will try to detect format from record 
def gene_record_write(rec, file_handle, format = None):
    assert format == 'fasta' or format == 'fastq' or format == None

    # Fallback to Bio using the same API - for testing
    if USE_BIO:
        if format == None:
            format = BIO_TYPE
        SeqIO.write(rec, file_handle, format)

    if format == None and rec.format != None:
        format = rec.format

    assert format != None

    # fasta format example:
    #   >EAS54_6_R1_2_1_413_324
    #   CCCTTCTTGTCTTCAGCGTTTCTCC
    if format == 'fasta':
        file_handle.write('>' + rec.id + '\n')
        file_handle.write(str(rec.seq) + '\n')

    # fastq format example:
    #   @EAS54_6_R1_2_1_413_324
    #   CCCTTCTTGTCTTCAGCGTTTCTCC
    #   +
    #   #########################
    elif format == 'fastq':
        file_handle.write('@' + rec.id + '\n')
        file_handle.write(str(rec.seq) + '\n')
        file_handle.write('+\n')
        file_handle.write('#' * len(rec.seq) + '\n') # Dummy quality


# Double queue based on full-duplex pipes
# Master pipe: Master puts an object for the slave, and gets an object from the slave
# Slave pipe: Slave puts an object for the master, and gets an object from the master
class GeneDQueue(object):
    def __init__(self):
        self._master, self._slave = connection.Pipe(duplex=True)
        self._master_pid = os.getpid()
        self._my_pipe = None

        self.__set_resources()

        multiprocessing.util.register_after_fork(self, GeneDQueue.__set_resources)

    def __set_resources(self):
        _is_master = (self._master_pid == os.getpid())
        
        if _is_master:
            self._my_pipe = self._master
        else:
            self._my_pipe = self._slave

    def empty(self):
        return not self._my_pipe.poll()

    def get(self):
        res = self._my_pipe.recv_bytes()
        return pickle.loads(res)

    def put(self, obj):
        obj = pickle.dumps(obj)
        self._my_pipe.send_bytes(obj)

# this function gets sine fasta file and extracts from it the sine sequence only.
def get_sine_forward(sine_fname):
    """Only in direction given in file."""
    [sine_record] = SeqIO.parse(sine_fname, "fasta", alphabet=IUPAC.IUPACAmbiguousDNA())
    # TODO: If we return it as dumb string, why did we bother about the alphabet?
    # TODO: The reference SINEs do contain a couple ambiguous chars - N, Y. 
    return str(sine_record.seq)
