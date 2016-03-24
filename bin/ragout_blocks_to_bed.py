#!/hive/groups/recon/local/bin/python

import argparse

import synteny_blocks.model as model


def get_location(seq_id, chroms):
    for c in chroms:
        if seq_id == c.seq_id:
            return c.description
    raise Exception('No such entry id!')

'''
def get_specie_region(seq_id, chroms):
    for c in chroms:
        if seq_id == c.seq_id:
            return c.description.split('.')
    raise Exception('No such entry id!')
'''

#handle cases when there are '.'s inside of chromosome name
def get_specie_region(seq_id):
    a = seq_id.split('.')
    return a[0], '.'.join(a[1:])

def print_specie_bed(blocks, chroms, ref_specie):
    for b in blocks:
        entries = b.entries
        for e in entries:
            #specie, region = get_specie_region(e.seq_id,chroms)
            specie, region = get_specie_region(e.seq_id)
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
    chroms = model.parse_chromosomes(args.file)
    blocks = model.parse_blocks(args.file)
    if args.specie:
        print_specie_bed(blocks, chroms, args.specie)
    else:
        print_bed(blocks, chroms)
