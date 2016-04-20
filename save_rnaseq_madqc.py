#!/usr/bin/env python3
import collections
import dateutil
import time
import pandas
import shelve
import logging

from htsworkflow.submission.encoded import ENCODED

logger = logging.getLogger(__name__)

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


def isrsem(experiment_file):
    """Does this ENCODED file object refer to a RSEM file.
    """
    analysis = experiment_file.get('analysis_step_version')
    if analysis:
        software_versions = analysis.get('software_versions')
        if software_versions:
            for version in software_versions:
                software = version['software']
                if software['name'] == 'rsem':
                    return True

    return False


def get_replicate_tuple(replicate):
    """Return a tuple of the biological and technical replicate numbers
    """
    bio = replicate['biological_replicate_number']
    tech = replicate['technical_replicate_number']
    return (bio, tech)


RSEMInfo = collections.namedtuple(
    'FileInfo',
    ['date_created', 'file_id', 'library_id', 'experiment_id', 'spikes_used', 'href']
)

def find_rsem(files):
    """Find most recent RSEM urls for given a list of file objects.
    """
    best_reps = {}
    for file in files:
        if file['file_format'] == 'tsv' and file['output_type'] == 'gene quantifications':
            if not isrsem(file):
                continue
            replicate = file['replicate']
            rep_id = get_replicate_tuple(replicate)
            library = replicate['library']
            library_id = library['accession']
            experiment = replicate['experiment']
            experiment_id = experiment['accession']
            spikes = library['spikeins_used']

            spikes_used = [ url_end(x) for x in spikes ]

            file_info = RSEMInfo(
                file['date_created'], file['@id'], library_id, experiment_id, spikes_used, file['href']
            )
            best_reps.setdefault(rep_id, []).append(file_info)

    for rep_id in best_reps:
        reps = best_reps[rep_id]
        reps = sorted(reps)
        yield reps[-1]

def load_rsems(cache, experiment_keys, quantification='fpkm', limit=None):
    """Return r
    """
    score_col = rsem_quantification_to_column(quantification)

    keys = list(keys)
    total = len(keys)
    chunk = max(total // 10, 1)
    tzero = time.monotonic()
    tprev = tzero

    for i, experiment_id in enumerate(exerpiment_keys):
        experiment = cache[experiment_id]
        score = []
        for file in save_rnaseq_madqc.find_rsem(experiment['files']):
            url = 'https://www.encodeproject.org' + file.href
            score = pandas.read_csv(
                url, usecols=[0,score_col], sep='\t', index_col=0)
            score.columns = [file.library_id]
            scores.append(scores)

        if scores:
            yield (experiment_id, pandas.concat(scores, axis=1))

        if (i + 1) % chunk == 0:
            tnow = time.monotonic()
            logger.info("{} of {} in {:.2f} sec".format(
                i, total, tnow-tprev))
            tprev = tnow

        if limit and i > limit:
            return

def rsem_quantification_to_column(name):
    """Convert quantification name to RSEM column number
    """
    scores = {
        'length': 2,
        'effective_length': 3,
        'expected_count': 4,
        'TPM': 5,
        'FPKM': 6,
    }

    score_column = scores.get(column)
    if score_column is None:
        raise ValueError("Unrecognized column name")

    return score_column


def make_experiment_df(experiments):
    """Create experiment containing some descriptive and qc status.
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
                         ('spikeins_used', ','.join([url_end(x) for x in library['spikeins_used']])),
                         # the qc metrics
                         ('Pearson', most_recent_qc['Pearson correlation']),
                         ('Spearman', most_recent_qc['Spearman correlation']),
                         ('MAD', most_recent_qc['MAD of log ratios']),
                         ('SD', most_recent_qc['SD of log ratios']),
                        ))
                    rows.append(record)
    return pandas.DataFrame(rows, columns=record.keys())

def caching_encoded_experiment_loader(query_url, cache_name):
    """Cache DCC objects found with query in a python shelf

    Parameters:
      - query is the url from ENCODED
      - cache_name base name of the cache file

    Note: to refresh objects delete cache_name.db
    """
    shelf_name = cache_name
    #shelf_db_name = shelf_name + '.db'

    experiments = shelve.DbfilenameShelf(shelf_name)
    encoded_experiment_loader(query_url, experiments)
    return experiments

def encoded_experiment_loader(query_url, experiments=None):
    server = ENCODED('www.encodeproject.org')
    server.load_netrc()

    if experiments is None:
        experiments = {}

    query = server.get_json(query_url)
    tzero = time.monotonic()
    tnow = tzero
    tprev = tzero
    progress = len(query['@graph']) // 10
    for i, record in enumerate(query['@graph']):
        accession = record['@id'][len('/experiments/'):-1]
        if accession not in experiments:
            experiments[accession] = server.get_json(record['@id'])

        if progress != 0 and (i+1) % progress  == 0:
            tnow = time.monotonic()
            print("Reading {} of {} records in {} seconds".format(
                  (i+1),
                  len(query['@graph']),
                  tnow - tprev))
            tprev = tnow
    print("Read {} records in {} seconds".format(
        len(query['@graph']), tnow-tzero))

    return experiments


def main():
    query_url = 'search/?type=experiment&assay_term_name=RNA-seq'
    cache_name = 'rnaseq-experiments.shelf'
    experiments = caching_encoded_experiment_loader(query_url, cache_name)
    df = make_experiment_df(experiments)
    df.to_csv('experiment-mad-qc.csv', index=False)

if __name__ == '__main__':
    main()
