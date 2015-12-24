import pandas
import os
import collections

from htsworkflow.submission.encoded import ENCODED

def load_caltech1000x():
    # our lab doesn't have a nice TSV file posted 
    caltech_1000x = pandas.DataFrame([ (x[9:], y) for x,y in [
        ('gSpikein_ERCC-00130',  3612000,),
        ('gSpikein_ERCC-00004',   903000,),
        ('gSpikein_ERCC-00136',   225750,),
        ('gSpikein_ERCC-00108',   112875,),
        ('gSpikein_ERCC-00116',    56438,),
        ('gSpikein_ERCC-00092',    28219,),
        ('gSpikein_ERCC-00095',    14109,),
        ('gSpikein_ERCC-00131',    14109,),
        ('gSpikein_ERCC-00062',     7055,),
        ('gSpikein_ERCC-00019',     3527,),
        ('gSpikein_ERCC-00144',     3527,),
        ('gSpikein_ERCC-00170',     1764,),
        ('gSpikein_ERCC-00154',      882,),
        ('gSpikein_ERCC-00085',      882,),
        ('gSpikein_ERCC-00028',      441,),
        ('gSpikein_ERCC-00033',      220,),
        ('gSpikein_ERCC-00134',      220,),
        ('gSpikein_ERCC-00147',      110,),
        ('gSpikein_ERCC-00097',       55,),
        ('gSpikein_ERCC-00156',       55,),
        ('gSpikein_ERCC-00123',       28,),
        ('gSpikein_ERCC-00017',       14,),
        ('gSpikein_ERCC-00083',        3,),
        ('gSpikein_ERCC-00096',  1806000,),
        ('gSpikein_ERCC-00171',   451500,),
        ('gSpikein_ERCC-00009',   112875,),
        ('gSpikein_ERCC-00042',    56438,),
        ('gSpikein_ERCC-00060',    28219,),
        ('gSpikein_ERCC-00035',    14109,),
        ('gSpikein_ERCC-00025',     7055,),
        ('gSpikein_ERCC-00051',     7055,),
        ('gSpikein_ERCC-00053',     3527,),
        ('gSpikein_ERCC-00148',     1764,),
        ('gSpikein_ERCC-00126',     1764,),
        ('gSpikein_ERCC-00034',      882,),
        ('gSpikein_ERCC-00150',      441,),
        ('gSpikein_ERCC-00067',      441,),
        ('gSpikein_ERCC-00031',      220,),
        ('gSpikein_ERCC-00109',      110,),
        ('gSpikein_ERCC-00073',      110,),
        ('gSpikein_ERCC-00158',       55,),
        ('gSpikein_ERCC-00104',       28,),
        ('gSpikein_ERCC-00142',       28,),
        ('gSpikein_ERCC-00138',       14,),
        ('gSpikein_ERCC-00117',        7,),
        ('gSpikein_ERCC-00075',        2,),
        ('gSpikein_ERCC-00074',  1806000,),
        ('gSpikein_ERCC-00113',   451500,),
        ('gSpikein_ERCC-00145',   112875,),
        ('gSpikein_ERCC-00111',    56438,),
        ('gSpikein_ERCC-00076',    28219,),
        ('gSpikein_ERCC-00044',    14109,),
        ('gSpikein_ERCC-00162',     7055,),
        ('gSpikein_ERCC-00071',     7055,),
        ('gSpikein_ERCC-00084',     3527,),
        ('gSpikein_ERCC-00099',     1764,),
        ('gSpikein_ERCC-00054',     1764,),
        ('gSpikein_ERCC-00157',      882,),
        ('gSpikein_ERCC-00143',      441,),
        ('gSpikein_ERCC-00039',      441,),
        ('gSpikein_ERCC-00058',      220,),
        ('gSpikein_ERCC-00120',      110,),
        ('gSpikein_ERCC-00040',      110,),
        ('gSpikein_ERCC-00164',       55,),
        ('gSpikein_ERCC-00024',       28,),
        ('gSpikein_ERCC-00016',       28,),
        ('gSpikein_ERCC-00012',       14,),
        ('gSpikein_ERCC-00098',        7,),
        ('gSpikein_ERCC-00057',        2,),
        ('gSpikein_ERCC-00002',  1806000,),
        ('gSpikein_ERCC-00046',   451500,),
        ('gSpikein_ERCC-00003',   112875,),
        ('gSpikein_ERCC-00043',    56438,),
        ('gSpikein_ERCC-00022',    28219,),
        ('gSpikein_ERCC-00112',    14109,),
        ('gSpikein_ERCC-00165',     7055,),
        ('gSpikein_ERCC-00079',     7055,),
        ('gSpikein_ERCC-00078',     3527,),
        ('gSpikein_ERCC-00163',     1764,),
        ('gSpikein_ERCC-00059',     1764,),
        ('gSpikein_ERCC-00160',      882,),
        ('gSpikein_ERCC-00014',      441,),
        ('gSpikein_ERCC-00077',      441,),
        ('gSpikein_ERCC-00069',      220,),
        ('gSpikein_ERCC-00137',      110,),
        ('gSpikein_ERCC-00013',      110,),
        ('gSpikein_ERCC-00168',       55,),
        ('gSpikein_ERCC-00041',       28,),
        ('gSpikein_ERCC-00081',       28,),
        ('gSpikein_ERCC-00086',       14,),
        ('gSpikein_ERCC-00061',        7,),
        ('gSpikein_ERCC-00048',        2,),
    ]], columns=['spikes', 'Caltech 1000x'])
    caltech_1000x.set_index('spikes', inplace=True)
    return caltech_1000x

def load_caltech_single():
    caltech_cpc = load_caltech1000x() / 1000
    caltech_cpc.columns = ['Caltech Single']
    return caltech_cpc

# the Gravely table is messy and has header lines split over two lines.
# pools 12,13,14,15 have units 'actual concentration (nmol/ul)'
def load_graveley(name, href):
    pool = name[-2:]
    columns=['Control', 'Subpool', 'original Genbank', 'Source', 'Length', 'GC content', 'MW', '12', '13', '14', '15']
    graveley = pandas.read_csv(href, skiprows=11, header=None, sep='\t', names=columns, usecols=['Control', pool])
    graveley.columns = ['spike', name]
    graveley.set_index('spike', inplace=True)
    return graveley

def load_simple(name, href):
    """The gingeras group used a simple format with two columns and some # comments
    """
    table = pandas.read_csv(href, comment='#', sep='\t', header=None,
                            names=['spike', name])
    table.set_index('spike', inplace=True)
    return table
    
def load_spike_table():
    """Return reported concentrations for the ENCODE RNA-Seq spikes.
    """
    spikes = {
        'ENCSR133ALU':  load_caltech_single(),
        'ENCSR884LPM':  load_caltech1000x(),
        'ENCSR156CIL':  load_simple('Gingeras-Ambion Mix 1', 'https://www.encodeproject.org/documents/7ae25f97-ab0a-448b-b345-cfffd6368c90/@@download/attachment/AmbionMix1Concentrations.tsv'),
        'ENCSR470JZL':  load_simple('Gingeras-NIST Profile 13', 'https://www.encodeproject.org/documents/2c038562-a136-4df0-9dc8-fde86ab809a8/@@download/attachment/NistPool13Concentrations.tsv'),
        'ENCSR402QNO':  load_simple('Gingeras-NIST Profile 14', 'https://www.encodeproject.org/documents/f96a19a6-dc66-4023-993a-f4c28af3d8e2/@@download/attachment/NistPool14Concentrations.tsv'),
        'ENCSR449DXG':  load_graveley('Gravely-NIST Profile 15', 'https://www.encodeproject.org/documents/122fc48c-aec6-4113-9821-178a435a8646/@@download/attachment/NIST-ERCC-Pool15.tsv'),
    }
    spikes_order = ['ENCSR133ALU', 'ENCSR884LPM', 'ENCSR156CIL', 'ENCSR470JZL', 'ENCSR402QNO', 'ENCSR449DXG']

    allspikes = pandas.concat([spikes[s] for s in spikes_order], axis=1)
    return allspikes


    
