#!/usr/bin/env python3
import argparse
import pandas
import logging
import os

from mergequantifications import \
  (load_quantification, remove_never_expressed, gencode_only,
   expressed_in_at_least,
   add_name_from_gencode_gene_id,
  )

logger = logging.getLogger('human_single')

def main(cmdline=None):
    logging.basicConfig(level=logging.DEBUG)
    
    gm12878 = load_gm12878()
    cont_purk = load_quantification(
        os.path.expanduser('~/proj/single-cell/Hs_purkinje_single_FPKM.h5'),
        'cont_purk')
    asp_purk = load_quantification(
        os.path.expanduser('~/proj/single-cell/Hs_asp_purkinje_UMB5294_single_FPKM.h5'),
        'asp_purk')

    combined = pandas.concat([gm12878, cont_purk, asp_purk], axis=1)
    filtered = expressed_in_at_least(gencode_only(combined), level=3, atleast=1)

    gtf_cache = os.path.expanduser('~/proj/encode3-curation/human.vV19-tRNAs-ERCC.h5')
    filtered = add_name_from_gencode_gene_id(gtf_cache, filtered)
    filtered = remove_sno(filtered)
    
    logger.info('final shape: %s', filtered.shape)
    filtered.to_csv('all_human_gm12878_control_purk_asp_purk_spike-_trna-sno-.tsv', sep='\t')

def load_gm1278():
    gm = pandas.read_csv(os.path.expanduser('~diane/proj/Hs_GM12878.csv'))
    columns = [ x + '_gm12878' for x in gm.columns]
    gm.columns = columns
    return gm


def remove_sno(filtered):
    # note requires name
    assert 'NAME' in filtered.columns
    not_snos = [ not x.startswith('Snor') for x in filtered['NAME']]
    return filtered[not_snos]

if __name__ == '__main__':
    main()
