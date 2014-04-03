"""Understand the ENCODE3 pipeline comparison fastqs.

this was some supporting code to try and figure out
what was happening with our fastq files.

Unfortunately I had the bright idea of saving
this all into a repository long after I wrote it 
so i'm not quite sure how this works.
"""

from __future__ import print_function
from idfastq.fastq_summary import FastqSummary, read_pretty
import os

def main():
    summaries = load_summaries(['georgi.fastqs.txt', 'dcc.fastqs.txt'])
    distances = compute_distances(summaries)

    for fastq in distances:
        print(fastq.filename, ':')
        for best in distances[fastq][:3]:
            if best[0] < 1.0:
                print('  ', best[0], best[1].filename)
              
def load_summaries(filenames):
    summaries = []
    for filename in filenames:
        with open(filename, 'r') as instream:
            summaries.extend(read_pretty(instream))

    print('Total recorded fastqs:', len(summaries))
    return summaries

def make_md5_index(summaries):
    md5s = {}
    for fastq in summaries:
        md5s[fastq.md5sum] = fastq
    return md5s

def make_filename_index(summaries):
    filenames = {}
    for fastq in summaries:
        filenames[fastq.filename] = fastq
    return filenames
    
def compute_distances(summaries):
    fastqs = {}
    for fastq in summaries:
        distances = [(fastq.distance(x), x) for x in summaries if x != fastq]
        distances.sort()
        fastqs[fastq] = distances
    return fastqs        
    
def get_filename(x):
    _, name = os.path.split(x.filename)
    return name
    
if __name__ == "__main__":
    main()
