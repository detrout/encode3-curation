#!/usr/bin/env python3
"""Convert a gff/gtf file to an indexed HDF5 file.

It takes about 11 minutes to run on the ENCODE GENCODE v19 file,
but afterwords I can query for attributes by gene name in 10-15 ms.
"""
import argparse
import os
import collections
import logging

import pandas
import numpy
import time

logger = logging.getLogger('gff2table')

class AttributesParser:
    def __init__(self):
        self.index = 0
        self.terms = collections.OrderedDict()
        self.max_string = collections.OrderedDict()
        
    def __call__(self, cell):
        records = [ x for x in cell.split(';') if len(x) > 0 ]
        for term in records:
            term = term.strip()
            space = term.find(" ")
            name = term[:space]
            if name == 'ont':
                continue
            
            value = term[space+1:]
            if value == '"NULL"':
                value = None
            elif value[0] == '"':
                value = value[1:-1]
                prev = self.max_string.get(name, 0)
                self.max_string[name] = max(prev, len(value))
            else:
                value = int(value)
            column = self.terms.setdefault(name, {})
            if name == 'ont':
                column.setdefault(self.index, []).append(value)
            else:
                column[self.index] = value
        self.index += 1

        return len(records)
    
def parse_score(x):
    if x == '.':
        return numpy.nan
    else:
        return x

def parse_strand(x):
    if x == '+':
        return 1
    elif x == '-':
        return -1
    else:
        return numpy.nan

def parse_phase(x):
    if x == '.':
        return numpy.nan
    else:
        return int(x)

def convert_gff(inname, outname):
    logger.info("Converting %s to %s", inname, outname)
    attribute_parser = AttributesParser()
    tzero = time.monotonic()
    gtf = pandas.read_csv(
        inname, 
        sep='\t', 
        header=None,
        index_col=False,
        names=['chromosome', 'source', 'type', 'start', 'stop',
               'score', 'strand', 'frame',
               'attributes'],
        na_values='.',
        converters={
            'strand': parse_strand,
            'attributes': attribute_parser,
        },
    )
    tnow = time.monotonic()
    tprev = tnow
    logger.info("Parsed in {:.3} seconds".format(tnow-tzero))
    for name in ['chromosome', 'type']:
        attribute_parser.max_string[name] = max(gtf[name].map(len))

    # drop my synthetic column counting how many records are in the
    # variant column
    gtf.drop('attributes', inplace=True)
    columns = set(gtf.columns)
    gtf.dropna(axis=1, how='all', inplace=True)
    logger.info("Removed columns for no data: {}",
                ",".join(columns.difference(set(gtf.columns))))
    attributes = pandas.DataFrame(attribute_parser.terms)
    gff = gtf.merge(attributes, left_index=True, right_index=True)
    tnow = time.monotonic()
    logger.info("Merged table in {:.3} seconds".format(tnow-tprev))
    tprev = tnow
    
    store = pandas.HDFStore(outname, mode='w', complevel=9, complib='bzip2')
    name = 'v19_tRNAs_ERCC'
    store.append(name, gff,
                 min_itemsize = attribute_parser.max_string)
    tnow = time.monotonic()
    logger.info("Wrote table in {:.3} seconds".format(tnow-tprev))
    tprev = tnow
    store.create_table_index(name, optlevel=9, kind='full')
    tnow = time.monotonic()
    logger.info("Wrote index in {:.3} seconds".format(tnow-tprev))
    store.close()

def main(cmdline=None):
    parser = argparse.ArgumentParser()
    parser.add_arguments('-v', '--verbose', action='store_true', default=False)
    parser.add_arguments('-d', '--debug', action='store_true', default=False)
    parser.add_arguments('-o', '--output',
                         help='specify output name, defaults to name.h5')
    parser.add_arugments('filename', nargs=1,
                         help='specify input GFF/GTF filename')

    args = parser.parse_args(cmdline)

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARN)
        
    if args.output is None:
        name, ext = os.path.splitext(args.filename)
        args.output = name + '.h5'
        
    convert_gff(args.filename, args.outname)


if __name__ == '__main__':
    main()
