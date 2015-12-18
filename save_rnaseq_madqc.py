#!/usr/bin/env python3
import collections
import dateutil
import time
import pandas
import shelve

from htsworkflow.submission.encoded import ENCODED

def url_end(url):
    """Convert DCC object ID URL to just the object ID
    """
    url = [ x for x in url.split('/') if len(x) > 0 ]
    return url[-1]

def aliases_to_wold_id(aliases):
    """Convert Library Alias to Wold lab library id.
    
    And only do this for barbara-wold: aliases since I don't
    know what other groups IDs mean.
    """
    alias = ''
    if len(aliases) > 1:
        raise RuntimeError("You need to do more work. Several library aliases: %s", aliases)
    elif len(aliases) == 1 and aliases[0].startswith('barbara-wold:'):
        alias = aliases[0].replace('barbara-wold:', '')

    return alias

def starting_quantity(quantity, units):
    """Try to normalize starting quantity.
    
    Ignore the difference between "cell-equivalent" and "cells"
    Squish our 11 cell group into its intended 10 cell set.
    """
    if quantity is None:
        return None
    
    if units == 'cell-equivalent':
        units = 'cells'
    
    quantity = float(quantity)
    assert int(quantity) == quantity
    quantity = int(quantity)

    formatted = "{} {}".format(quantity, units)
    if formatted == '11 cells':
        formatted = '10 cells'
        
    return formatted


def make_experiment_df(experiments):
    """Create a dataframe from a blob of DCC experiment json.
    """
    rows = []
    for accession, detail in experiments.items():
        # For all the files in an experiment
        for f in detail['files']:
            # look at files that have quality control metrics attached
            if len(f['quality_metrics']) > 0:
                # sort qc metrics by date created
                sorted_qc = sorted(
                    [ qc for qc in f['quality_metrics'] ], 
                    key=lambda x: dateutil.parser.parse(x['date_created']))
                # choose the most recent run.
                most_recent_qc = sorted_qc[-1]

                if 'replicate' not in f:
                    print('Skipping:', f['@id'])
                    continue

                replicate = f['replicate']
                library = replicate['library']
                biosample = library['biosample']
                alias = aliases_to_wold_id(library['aliases'])
                bio_rep = replicate['biological_replicate_number']
                tech_rep = replicate['technical_replicate_number']

                # we only want to look at the MAD scores right now.
                # we're ignoring the STAR & RSEM scores and just
                # investigating Rafa's MAD QC output.
                # 
                # Also the MAD qc scores are symetric between replicates so we 
                # only need scores from one file. 
                # (the reason for requiring bio_rep and tech_rep == 1)
                if 'MadQualityMetric' in most_recent_qc['@type'] and bio_rep == 1 and tech_rep == 1:
                    record = collections.OrderedDict(
                        # long list of things to identify a particular experiment replicate
                        (('experiment', accession), 
                         ('lab', detail['lab']['title']),
                         ('rfa', detail['award']['rfa']),
                         ('description', detail['description']),
                         ('organism', url_end(biosample['organism'])),
                         ('biosample', biosample['biosample_term_name']),
                         ('biosample_lab', url_end(biosample['lab'])),
                         ('age', biosample['age'] + ' ' + biosample.get('age_units', '')),
                         ('starting', starting_quantity(library.get('nucleic_acid_starting_quantity'),
                                                        library.get('nucleic_acid_starting_quantity_units'))),
                         ('library_id', alias),
                         # the qc metrics
                         ('Pearson', most_recent_qc['Pearson correlation']),
                         ('Spearman', most_recent_qc['Spearman correlation']),
                         ('MAD', most_recent_qc['MAD of log ratios']),
                         ('SD', most_recent_qc['SD of log ratios']),
                        ))
                    rows.append(record)
    return pandas.DataFrame(rows, columns=record.keys())

def caching_encoded_experiment_loader(url, cache_name):
    """Cache DCC objects found with query in a python shelf

    Parameters:
      - query is the url from ENCODED
      - cache_name base name of the cache file

    Note: to refresh objects delete cache_name.db
    """
    server = ENCODED('www.encodeproject.org')
    server.load_netrc()
    
    shelf_name = cache_name
    shelf_db_name = shelf_name + '.db'

    experiments = shelve.DbfilenameShelf(shelf_name)
    query = server.get_json(url)
    tzero = time.monotonic()
    tprev = tzero
    progress = len(query['@graph']) // 10
    for i, record in enumerate(query['@graph']):
        accession = record['@id'][len('/experiments/'):-1]
        if accession not in experiments:
            experiments[accession] = server.get_json(record['@id'])
            
        if (i+1) % progress  == 0:
            tnow = time.monotonic()
            print("Reading {} of {} records in {} seconds".format(
                  (i+1),
                  len(query['@graph']), 
                  tnow - tprev))
            tprev = tnow
    print("Read {} records in {} seconds".format(len(query['@graph']), tnow-tzero))
    #experiments.close()

    #experiments = shelve.DbfilenameShelf(shelf_name)
    return experiments

def main():

    query_url = 'search/?type=experiment&assay_term_name=RNA-seq'
    cache_name = 'rnaseq-experiments.shelf'
    experiments = caching_encoded_experiment_loader(query_url, cache_name)
    df = make_experiment_df(experiments)
    df.to_csv('experiment-mad-qc.csv', index=False)

if __name__ == '__main__':
    main()
