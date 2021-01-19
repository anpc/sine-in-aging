#!/usr/bin/env python3.8
import gzip
import sys
#import zstandard as zstd
import itertools
import time
import os

from gene_lib import log

def chunks(iterable, chunk_size):
    """
    Split input into fixed-sized chunk iterators, without buffering.

    >>> for chunk in chunks(['a', 'b', 'c', 'd', 'e', 'f'], 3):
    ...     print(list(chunk))
    ['a', 'b', 'c']
    ['d', 'e', 'f']
    >>> for chunk in chunks(['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h'], 3):
    ...     print(list(chunk))
    ['a', 'b', 'c']
    ['d', 'e', 'f']
    ['g', 'h']

    As we're gonna use it with chunk size >20GB (uncompressed), 
    we obviously don't want to buffer a whole chunk in memory!
    So each inner `chunk` we produce must be an iterator, 
    yielding values as we get them from input.

    As we don't know the length of the input until we consumed it all,
    the outer interace is also an iterator!
    It will not start a new chunk iterator until previous one is fully consumed.

    >>> def verbose_iterator(name, items):
    ...     for x in items:
    ...         print(f"{name} iterator yielding -> {repr(x)}")
    ...         yield x
    ...     print(f"{name} iterator finished.")
    
    >>> input = verbose_iterator("        input", ['a', 'b', 'c', 'd', 'e'])
    >>> chunks_iter = verbose_iterator("chunks", chunks(input, chunk_size=3))
    >>> for chunk_iter in chunks_iter:
    ...     for x in verbose_iterator("    inner", chunk_iter):
    ...         pass
    ...
            input iterator yielding -> 'a'
    chunks iterator yielding -> <itertools.chain object at ...>
        inner iterator yielding -> 'a'
            input iterator yielding -> 'b'
        inner iterator yielding -> 'b'
            input iterator yielding -> 'c'
        inner iterator yielding -> 'c'
        inner iterator finished.
            input iterator yielding -> 'd'
    chunks iterator yielding -> <itertools.chain object at ...>
        inner iterator yielding -> 'd'
            input iterator yielding -> 'e'
        inner iterator yielding -> 'e'
            input iterator finished.
        inner iterator finished.
    chunks iterator finished.
    """
    # Implementation stolen from https://stackoverflow.com/a/24527424/1413499
    # See there for explanation & several equivalent ways to write it.
    shared_input_iterator = iter(iterable)
    for first_item_of_chunk in shared_input_iterator:
        yield itertools.chain([first_item_of_chunk],
                              itertools.islice(shared_input_iterator, chunk_size - 1))

def split_recompress(in_fname, out_basename, power_of_ten, skip=[]):
    """
    skip - indexes of parts to skip (still need to decompress but not compress)
    """
    with gzip.open(in_fname, 'rb') as in_fastq:
        for i, chunk_iter in enumerate(chunks(in_fastq, chunk_size=10**power_of_ten)):
            out_fname = f'{out_basename}.part{i}e{power_of_ten}.fastq.gz'
            t0 = time.time()
            if i in skip:
                log('Skipping', out_fname)
                max(chunk_iter)  # consume iterator, discarding values
            else:
                log('Writing', out_fname)
                with gzip.open(out_fname + '.tmp', 'wb', compresslevel=2) as out_fastq_gz:
                #with open(out_fname, 'wb') as out_fastq_gz:
                    out_fastq_gz.writelines(chunk_iter)
                os.rename(out_fname + '.tmp', out_fname)
            t1 = time.time()
            log(t1 - t0, 'sec.')

#split_recompress('Old-lung/old_lung_R1_001.fastq.gz', 'Old-lung/old_lung_R1_001')

if __name__ == '__main__':
    [in_fname, out_basename, power_of_ten] = sys.argv[1:]
    split_recompress(in_fname, out_basename, int(power_of_ten))
