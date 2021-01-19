TODO document all stages

# Formats we use

[.fasta](https://en.wikipedia.org/wiki/FASTA_format) - a record header line, followed by AGCT characters.
But can contain other letters that represent reading uncertainty, e.g. Y = C, T or U.

[.fastq](https://en.wikipedia.org/wiki/FASTQ_format) - similar but also contains per-character *(Q)uality* info.

Our particular FASTQ inputs don't wrap lines, so contain exactly 4 lines per record.

.fastq.gz, .fastq.zst - FASTQ file, compressed with off-the-shelf compression tools (gzip, zstd).

We have a function `open_any()` that does de/compression automatically, depending on file extension.

# Install tools

On ubuntu linux, try running `./INSTALL/ALL.sh`.

Note it installs python libraries for a _single_ python version.  All our *.py scripts should use that version.

See INSTALL/README.md for details on install script.

# Dowload input from AWS S3 storage

Choose an organ whose samples you want to process:

```bash
aws s3 ls s3://endogene/mice_wgs/Aging/
aws s3 ls s3://endogene/mice_wgs/Aging/Old-liver/
mkdir Old-liver
aws s3 cp s3://endogene/mice_wgs/Aging/Old-liver/old_liver_R1_001.fastq.gz Old-liver/
```

Now repeat for `_R2_` and for both Old and Young samples — but if you're short on disk space, you'll want to next steps for each one so you can delete inputs...

## streaming

The command supports `-` as either source (meaning stdin) or destination (meaning stdout):
```bash
aws s3 cp s3://endogene/mice_wgs/Aging/Old-liver/old_liver_R1_001.fastq.gz - | gunzip --stdout | head --lines=20
```

# breaking into chunks

for example in p53-liver
```bash
./split_recompress.py p53-liver_R1_001.fastq.gz p53-liver_R1_001 8 && rm --verbose p53-liver_R1_001.fastq.gz
./split_recompress.py p52-liver_R2_001.fastq.gz p53-liver_R2_001 8 && rm --verbose p52-liver_R2_001.fastq.gz
...
```
(The `rm` is of course optional.)

The `8` parameter makes it cut input into chunks of 100 000 000 (1e8) lines - which are exactly 25 000 000 records - writing files named:
```
p53-liver_R1_001.part0e8.fastq.gz
p53-liver_R1_001.part1e8.fastq.gz
...
p53-liver_R1_001.part27e8.fastq.gz
...
```
(last chunk will be smaller.)

This allows following steps to be run from the middle and/or parallelized.

## streaming directly from AWS S3

This uses bash "Process substitution" syntax.  
See blog post [Pipes, process substitution and why should a biologist ever care
](http://manutamminen.info/posts/process_subst/).

```bash
./split_recompress.py <(aws s3 cp s3://endogene/mice_wgs/Aging/p53-liver/p53-liver_R2_001.fastq.gz -) p53-liver/p53-liver_R2_001 8
```
What this `<(...)` syntax does is open a pipe from `aws s3 cp ... -` to `split_recompress.py` process,
but not to stdin. It gets assigned some new file descriptor, say 63.
Bash then replaces the `<(...)` with a special file name like `/dev/fd/63`.

The neat thing is split_recompress.py isn't even aware anything special happened, it asks operating system to open the file name given to it, and gets back the pipe!

- One limitation is that `open_any()` function can't detect input compression from file extension (`/dev/fd/NN` has no extension).
  (But we can always inset `| gunzip --stdout` / `| gzip --stdout -2` process...)
  And split_recompress.py assumes input is always in .gz format.

## Tests
```bash
python3.6 -m doctest --option=ELLIPSIS split_recompress.py
```

# r1r2merge.py

To support streaming, this script no longer takes an output file name;
it now always writes to stdout, *without compression*.

To merge the whole input:

```bash
./r1r2merge.py p53-liver/p53-liver_R1_001.fastq.gz p53-liver/p53-liver_R2_001.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged.fastq.gz
```

But we can run it faster by parallelizing over chunks produced above:
You should replace 2nd argument of `seq` with number of chunks you have.
`--jobs=-2` means number of CPU cores you have minus 2.

```bash
seq 0 32 | parallel --jobs=-2 --bar --eta --joblog=p53-liver/merge.joblog '
  ./r1r2merge.py p53-liver/p53-liver_R1_001.part{}e8.fastq.gz p53-liver/p53-liver_R2_001.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged.part{}e8.fastq.gz &&
  rm --verbose p53-liver/p53-liver_R1_001.part{}e8.fastq.gz p53-liver/p53-liver_R2_001.part{}e8.fastq.gz
'
```

## Tests

```bash
./r1r2merge.py merge-test-R1.fastq merge-test-R2.fastq > merge-test-merged.fastq
diff --report-identical-files merge-test-merged.fastq merge-test-expected.fastq
```

# Filter (crude) candidates

Now we can find lines that are candidate (aka potential) SINEs.  
(Shouldn't do it before merging, as we'd miss cases where each read contains half of the SINE.)

The parameters we'd used in 2019 were lookig for first 67 chars of a SINE, allowing edit distance up to 14.

But it's better not to filter with exact params we want, but first take a **crude superset**!
Motivation: simply reading & decompressing the full ~60GB takes hours.  By allowing say edit distance up to 19, we're already reducing the input size by 2 orders of magnitude, and reading that takes less than a minute!  That means once we do this filtering, we can upload it to AWS S3 and later experiment with more precise thresholds.

```bash
seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B1forward.joblog './filter_candidates.py B1.fasta 67 19 forward p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B1forward-head67err19.part{}e8.fastq.gz'
seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B1rc.joblog './filter_candidates.py B1.fasta 67 19 rc p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B1rc-head67err19.part{}e8.fastq.gz'

seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B2forward.joblog './filter_candidates.py B2.fasta 67 19 forward p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B2forward-head67err19.part{}e8.fastq.gz'
seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B2rc.joblog './filter_candidates.py B2.fasta 67 19 rc p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B2rc-head67err19.part{}e8.fastq.gz'

seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B4forward.joblog './filter_candidates.py B4.fasta 67 19 forward p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B4forward-head67err19.part{}e8.fastq.gz'
seq 0 32 | parallel --jobs=-2 --eta --joblog=p53-liver/filter-B4rc.joblog './filter_candidates.py B4.fasta 67 19 rc p53-liver/p53-liver_merged.part{}e8.fastq.gz | gzip --stdout -2 > p53-liver/p53-liver_merged-candidates-B4rc-head67err19.part{}e8.fastq.gz'
```

# TODO document rest of process 

- Currently, run_part_1.py is the high level script. This is messy and should be refactored in the future
- To generate potential sines run mode = 1.
- To generate barcodes run mode = 3
- cat $(ls -v old_liver_merged*/old_liver_merged.part*e8_sineBarcode.fastq.gz) > old_liver_merged_sineBarcode.fastq.gz
