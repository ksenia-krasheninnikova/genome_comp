#!/usr/bin/env python

import argparse
from synteny_blocks.model import MAF_Entry
from synteny_blocks.model import parse_bed

def intersect(maf_entries, bed_entries):
    intersected_bed_entries = []
    for m in maf_entries:
        for b in bed_entries:
            if m.genome == b.genome and m.chrom == b.chrom:
                #print m.genome, m.chrom
                if not (b.end < m.global_start or m.global_end < b.start) :
                    intersected_bed_entries.append(b)
    return intersected_bed_entries

def process(bed_file, maf_file):
    regions = parse_bed(bed_file)

    with open(maf_file) as maf:
        maf_entries = []
        for line in maf:
            line = line.strip()
            if not line or '#' in line:
                continue
            if line[0] == 'a':
                if maf_entries:
                    intersected_regions = intersect(maf_entries, regions)
                    if intersected_regions:
                        for e in intersected_regions:
                            print '##',
                            e.print_out()
                        print 'a'
                        for maf_e in maf_entries:
                            maf_e.print_out()
                        print
                    maf_entries = []
            else:
                whole_line = line
                line = line.split()
                line = line[1:]
                genome = line[0].split('.')[0]
                chrom = '.'.join(line[0].split('.')[1:])
                maf_entries.append(MAF_Entry(genome, chrom, int(line[1]), int(line[2]),\
                                             line[3], int(line[4]), line[4]))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('bed', help='bed file')
    parser.add_argument('in_maf', help='input maf file')

    args = parser.parse_args()
    process(args.bed, args.in_maf)
