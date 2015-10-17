import os
import gzip
import sys
import time
import logging

from IPython import get_ipython

logger = logging.getLogger(__name__)


# maybe we should just install htsworkflow
def custom_paths():
    paths = [ os.path.expanduser('~/proj/htsworkflow'),
              ]
    for p in paths:
        if p not in sys.path:
            sys.path.append(p)
custom_paths()

from htsworkflow.submission.encoded import ENCODED
from htsworkflow.util.rdfhelp import get_model, dump_model, load_into_model
from htsworkflow.util.rdfjsonld import load_into_model as load_jsonld_into_model

if not 'DJANGO_SETTINGS_MODULE' in os.environ:
    os.environ['DJANGO_SETTINGS_MODULE'] = 'htsworkflow.settings.myrada'

ipython = get_ipython()

ipython.magic('load_ext rdfmagic')
ipython.magic('%addns htsw http://jumpgate.caltech.edu/wiki/LibraryOntology#')
ipython.magic('%addns library http://jumpgate.caltech.edu/library/')
ipython.magic('%addns flowcell http://jumpgate.caltech.edu/flowcell/')

ipython.magic('%addns biosample https://www.encodeproject.org/profiles/biosample.json#')
ipython.magic('%addns experiment https://www.encodeproject.org/profiles/experiment.json#')
ipython.magic('%addns library https://www.encodeproject.org/profiles/library.json#')
ipython.magic('%addns platform https://www.encodeproject.org/profiles/platform.json#')
ipython.magic('%addns replicate https://www.encodeproject.org/profiles/replicate.json#')
ipython.magic('%addns file https://www.encodeproject.org/profiles/file.json#')

ipython.magic('%addns experiments https://www.encodeproject.org/experiments/')
ipython.magic('%addns libraries https://www.encodeproject.org/libraries/')
ipython.magic('%addns files https://www.encodeproject.org/files/')


def build_geneid_to_gene_from_gtf(gencode):
    """Build a dictionary mapping from gene_id to gene_name.

    Arguments:
        gencode (str): compressed filename to read

    Returns:
        dictionary mapping gene_id to gene_name
    """
    logger.info("Loading %s", gencode)
    start = time.time()
    names = {}
    with gzip.GzipFile(gencode, 'r') as instream:
        for line in instream:
            line = line.decode('ascii')

            if line.startswith('#'):
                continue

            columns = line.split('\t')

            gene_id = None
            gene_name = None

            for item in columns[-1].split(';'):
                item = item.strip()
                if len(item) == 0:
                    continue
                item = item.split()
                if len(item) != 2:
                    print("Confused: {} {}".format(item, len(item)))
                name, value = item
                if name == 'gene_id':
                    gene_id = value[1:-1]
                elif name == 'gene_name':
                    gene_name = value[1:-1]

            if gene_id and gene_name:
                names[gene_id] = gene_name

    logger.info("loaded in %d seconds", time.time() - start)
    return names
