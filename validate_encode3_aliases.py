import gzip
import os
from io import BytesIO
import sys
import json
import jsonschema
import logging
import re

LOGGER = logging.getLogger(__name__)

import RDF
import pandas as pd
from IPython.display import display_html, clear_output

import requests
try:
    from pysam import Samfile
except ImportError as e:
    LOGGER.error("Pysam not available, bam reading wont work")

from htsworkflow.submission.ucsc import get_encodedcc_file_index
from htsworkflow.util.rdfjsonld import load_into_model
from htsworkflow.util.rdfhelp import get_model

from rdfmagic import load_source, make_temp_model, LibRdfResults, HTTPError


class CheckDCCWoldAlias:
    def __init__(self, server):
        self.model = get_model()
        self.server = server
        self.DCCIndex = LookupSubmittedFile()
        self.read_patterns = [
            b":(?P<fc>[A-Z0-9]*A[ACD]XX):(?P<lane>[\d]):.* (?P<end>[12]):[YN]",
            b":(?P<fc>[A-Z0-9]*A[ACD]XX):(?P<lane>[\d]):",
            b"_(?P<fc>[A-Z0-9]*AAXX)_(?P<lane>[\d])(_(?P<end>[\d]))?",
            b"_(?P<fc>FC[\d]+)_(?P<lane>[\d])",
                    ]
    def load_datasets(self, urls):
        """Load ENCODE3 datasets into our model.

        Once the core dataset classes are loaded
        also add in detailed information about the experiments.
        """
        self.load(urls)
        experiments = RDF.SPARQLQuery("""
            select distinct ?exp ?lib
    where {
    ?exp <http://submit.encodedcc.org/profiles/experiment.json#replicates> ?rep .
    ?rep <http://submit.encodedcc.org/profiles/replicate.json#library> ?lib .
    }""")
        self.load(( str(row['lib']) for row in experiments.execute(self.model)))

    def load(self, urls):
        for url in urls:
            LOGGER.info("Loading: %s", url)
            data = self.server.get_jsonld(url, Embed=False)
            load_into_model(self.model, data)

    def parse_read_id(self, read_id):
        """Parse different versiond of the read ID

        early flowcells were named FC<digits>
        modern flowcells are named <alphanum>AAXX

        some versions of the illumina pipeline didn't include the
        flowcell name, there were a few cases where I had added it in.
        Unfortunately I added it in in a different place from
        the official versions. So I took the easy way out and special
        cased those regexes.
        """
        for p in self.read_patterns:
            match = re.search(p, read_id)
            if match:
                return match.groupdict()

    def match_dcc_alias_and_our_files(self):
        """Check encode libraries that have aliases for matching woldlab library IDs.

        I know its a bit long, but basically this is generating a report.
        that shows the alias:library_id and either valid if all the metadata
        can be looked up, or reporting what metadata can be found if
        we don't have definitive results.
        """
        query = RDF.SPARQLQuery("""
    PREFIX experiment: <http://submit.encodedcc.org/profiles/experiment.json#>
    PREFIX replicate: <http://submit.encodedcc.org/profiles/replicate.json#>
    PREFIX library: <http://submit.encodedcc.org/profiles/library.json#>
    PREFIX file: <http://submit.encodedcc.org/profiles/file.json#>

    select distinct ?exp ?lib ?alias ?file ?format ?path ?submitted
    where {
    ?exp experiment:replicates ?rep .
    ?rep replicate:library ?lib .
    ?lib library:aliases ?alias .
    ?file file:replicate ?rep ;
            file:file_format ?format ;
            file:download_path ?path ;
            file:submitted_file_name ?submitted.
    }
    order by ?alias ?file
    """)

        prev_alias = None
        messages = []
        valid_library = False

        for row in query.execute(self.model):

            alias = str(row['alias'])
            _, library_id = alias.split(':')

            file_format = str(row['format'])
            path = str(row['path'])
            submitted = str(row['submitted'])
            header = None
            if alias != prev_alias:
                if prev_alias is not None:
                    display_file_library_ids(prev_alias, messages, valid_library)
                valid_library = False
                messages = []
                prev_alias = alias

            if not path:
                messages.append("No path: {}".format(str(row['file'])))
                continue

            href = 'http://encodedcc.sdsc.edu/warehouse/' + path
            read_id = self.read_flowcell_from_url(href, file_format, submitted)
            if not read_id:
                messages.append("Couldn't parse: {}".format(header))
                messages.append(self.DCCIndex.format_encode2_metadata(submitted))
                continue

            libraries = lookup_library_id_from_lims(read_id['fc'],
                                                    read_id['lane'])
            if not libraries:
                messages.append("Couldn't lookup: {}".format(
                    ','.join(read_id.items())))
                valid_library = False
                continue

            for l in libraries:
                if l == library_id:
                    valid_library = True
            messages.append('Library: {}'.format(','.join(libraries)))

        display_file_library_ids(alias, messages, valid_library)

    def read_flowcell_from_url(self, href, file_format, submitted=None):
        """lookup flowcell and lane using read id from a url
        """
        if not submitted:
            submitted = href

        if file_format == 'fastq':
            header = fastq_read_id(href)
        elif file_format == 'bam':
            header = bam_read_id(href)
        else:
            # unsported file type
            return

        if not header:
            return

        return self.parse_read_id(header)

    def load_fastq_headers(self, cache, file_href_to_wold_query):
        """given a sparql query result, open the fastq url and read the read id.

        Parameters
          - cache dictionary to hold the results

        The sparql query needs the attributes:
          - href
          - aliases
          - file_accession
        """
        for i, row in enumerate(file_href_to_wold_query):
            attributes = { k: str(v) for k, v in row.items() }
            href = attributes['href']
            if href not in cache:
                clear_output()
                print("loading: {} {}/{}".format(href, i, len(file_href_to_wold_query)))
                header = fastq_read_id(href)
                attributes['header'] = header
                cache[href] = attributes
                
        return cache

    
def display_file_library_ids(alias, results, valid_library):
    """Display formatted html for our report
    """
    if valid_library:
        display_html('<p>{} validated</p>'.format(alias), raw=True)
    else:
        display_html('<dl><dt>{}</dt><dd>{}</dd></dl>'.format(
            alias,
            '</br>'.join(results)), raw=True)


def lookup_library_id_from_lims(flowcell, lane):
    '''Given a flowcell & lane, go retrieve its library id.

    this involves downloading metadata from my LIMS and
    attempting to parse out the library id.

    It's returning a list because later libraries need to
    be identified by flowcell/lane/multiplex
    '''
    url = 'http://jumpgate.caltech.edu/flowcell/{flowcell}/{lane}'.format(
        flowcell=flowcell, lane=lane)
    m = make_temp_model()
    try:
        load_source(m, url)
    except HTTPException as e:
        print("Unable to open:", url)
        return []

    query = RDF.SPARQLQuery('''
prefix htsw: <http://jumpgate.caltech.edu/wiki/LibraryOntology#>

select ?flowcell ?lane_number ?library
where {
  ?flowcell a htsw:IlluminaFlowcell ;
            htsw:has_lane ?lane .
  ?lane htsw:library ?library ;
        htsw:lane_number ?lane_number .
}
''')
    libs = []
    for row in query.execute(m):
        library = str(row['library'])
        if library[-1] == '/':
            library = library[:-1]
        _, library = os.path.split(library)
        libs.append(library)

    return libs

class LookupSubmittedFile:
    """Lookup ENCODE3 submitted filenames in ENCODE2 metadata.

    The ENCODE3 imports linked to the UCSC ENCODE2 names.
    which isn't so useful for me to find anything.
    This class caches the files for different tracks, and
    can munge the ENCODE3 filename to match what my parser
    was returning.
    """
    def __init__(self):
        self.loaded = set()
        self.dcc_index = {}
        self.base_url = 'http://hgdownload-test.cse.ucsc.edu/goldenPath/'

    def update_cache(self, genome, track):
        if (genome, track) not in self.loaded:
            self.dcc_index.update(get_encodedcc_file_index(genome, track))
            self.loaded.add((genome, track))

    def lookup_filename(self, filename):
        parts = filename.split('/')
        if parts[0] in ('mm9','hg19'):
            genome = parts[0]
            track = parts[1]
            filename = parts[2]
            self.update_cache(genome, track)
            url = self.base_url + '{}/encodeDCC/{}/{}'.format(genome, track, filename)
            return self.dcc_index.get(url, None)

    def format_encode2_metadata(self, submitted):
        """Given a submitted filename return formatted metadata message.

        This restricts the report to the two most useful terms
        labExpId and md5sum (as when I looked I didn't see my filenames)
        """
        dcc2_metadata = self.lookup_filename(submitted)
        if dcc2_metadata:
            dcc2_attributes = []
            for term in ('labExpId', 'md5sum'):
                if term in dcc2_metadata:
                    dcc2_attributes.append('{}: {}'.format(term, dcc2_metadata[term]))
            return "Old metadata: {}".format(','.join(dcc2_attributes))
        else:
            return "Couldn't lookup file name"

def fastq_read_id(url):
    '''Read the first line (containing the read id) out of a remote fastq file'''
    data = requests.get(url, stream=True)

    block = BytesIO(next(data.iter_content(1024)))
    compressed = gzip.GzipFile(None, 'r', fileobj=block)
    header = compressed.readline().rstrip()
    return header

def bam_read_id(url):
    '''Read first read id out of a remote bam file.

    Note: requires a patched version of pysam
    '''
    stream = Samfile(url, 'rb')
    read = stream.next()
    return read.qname

