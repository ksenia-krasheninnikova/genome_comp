#!/hive/groups/recon/local/bin/python

import argparse
from collections import Counter
import ragout_synteny_blocks

def get_location(seq_id, chroms):
    for c in chroms:
        if seq_id == c.seq_id:
            return c.description
    raise Exception('No such entry id!')

def get_specie_region(seq_id, chroms):
    for c in chroms:
        if seq_id == c.seq_id:
            return c.description.split('.')
    raise Exception('No such entry id!')

def print_specie_bed(blocks, chroms, ref_specie):
    for b in blocks:
        entries = b.entries
        for e in entries:
            specie, region = get_specie_region(e.seq_id,chroms)
            if specie == ref_specie:
                print '\t'.join([region, str(e.start), str(e.end), str(b.id), '0', e.strand])

def print_bed(blocks, chroms):
    for b in blocks:
        entries = b.entries
        for e in entries:
            specie = get_location(e.seq_id,chroms)
            print '\t'.join([specie, str(e.start), str(e.end), str(b.id), '0', e.strand])

if  __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--specie')
    args = parser.parse_args()
    chroms = ragout_synteny_blocks.parse_chromosomes(args.file)
    blocks = ragout_synteny_blocks.parse_blocks(args.file)
    if args.specie:
        print_specie_bed(blocks, chroms, args.specie)
    else:
        print_bed(blocks, chroms)
