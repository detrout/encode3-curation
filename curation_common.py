import os
import sys

from IPython import get_ipython

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

ipython = get_ipython()

ipython.magic('load_ext rdfmagic')
ipython.magic('%addns htsw http://jumpgate.caltech.edu/wiki/LibraryOntology#')
ipython.magic('%addns library http://jumpgate.caltech.edu/library/')
ipython.magic('%addns flowcell http://jumpgate.caltech.edu/flowcell/')

ipython.magic('%addns biosample https://www.encodedcc.org/profiles/biosample.json#')
ipython.magic('%addns experiment https://www.encodedcc.org/profiles/experiment.json#')
ipython.magic('%addns library https://www.encodedcc.org/profiles/library.json#')
ipython.magic('%addns platform https://www.encodedcc.org/profiles/platform.json#')
ipython.magic('%addns replicate https://www.encodedcc.org/profiles/replicate.json#')

ipython.magic('%addns experiments https://www.encodedcc.org/experiments/')
