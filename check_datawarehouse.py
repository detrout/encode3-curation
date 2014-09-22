from __future__ import print_function
import requests
import gzip
import StringIO
import pysam
import ctypes

class BamHeader(ctypes.Structure):
    _fields_ = [('n_targets', ctypes.c_uint),
                ('ignore_sam_err', ctypes.c_int),
                ('l_text', ctypes.c_uint), 
                ('target_len', ctypes.POINTER(ctypes.c_uint32)), #uint32 *
                ('cigar_tab', ctypes.POINTER(ctypes.c_uint8)), #unit8 *
                ('target_name', ctypes.POINTER(ctypes.c_char_p)), #char **
                ('text', ctypes.c_char_p),
                ('sdict', ctypes.c_void_p)]

class Bam1Core(ctypes.Structure):
     _fields_ = [('tib', ctypes.c_int32),
                 ('pos', ctypes.c_int32),
                 ('qname_flags', ctypes.c_uint32),
                 ('cigar_flags', ctypes.c_uint32),
                 ('qseq', ctypes.c_int32),
                 ('mtid', ctypes.c_int32),
                 ('mpos', ctypes.c_int32),
                 ('isize', ctypes.c_int32)
                 ]
     
class Bam1(ctypes.Structure):
     _fields_ = [('core', Bam1Core),
                 ('l_data', ctypes.c_int32),
                 ('m_data', ctypes.c_int32),
                 ('data', ctypes.c_int8),
                 ('id', ctypes.c_uint64),
                 ]
     
htslib = ctypes.cdll.LoadLibrary('libhts.so.0')
htslib.sam_hdr_read.restype=ctypes.POINTER(BamHeader)
htslib.bam_read1.argtypes = (ctypes.c_int32, ctypes.POINTER(Bam1))
    
def fastq_header(url):
    data = requests.get(url, stream=True)
    
    block = StringIO.StringIO(data.iter_content(1024).next())
    compressed = gzip.GzipFile(None, 'r', fileobj=block)
    header = compressed.readline().rstrip()
    return header

def bam_header(url):
    fd = htslib.hts_open(url, 'rb')
    h = htslib.sam_hdr_read(fd)
    header = {}
    for line in StringIO.StringIO(h.contents.text):
        record = line.rstrip().split('\t')
        if record[0] == '@PG':
            return record[1:]

def bam_read(url):
    fd = htslib.hts_open(url, 'rb')
    read = Bam1()
    htslib.bam_read1(fd, ctypes.byref(read))
    
    
if __name__ == "__main__":
    bam='http://encodedcc.sdsc.edu/warehouse/2013/4/18/ENCFF001ICO.bam'
    bam_read(bam)
    #bam='../mcf7/accepted_hits.bam'
    #print(' '.join(bam_header(bam)))
    #fastq1 = 'http://encodedcc.sdsc.edu/warehouse/2013/4/18/ENCFF001ICW.fastq.gz'
    #print(fastq_header(fastq1))
