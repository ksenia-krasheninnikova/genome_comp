#!/usr/bin/env python

import argparse
import utils
import model

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--circos_output', help='output the --species for circos plot')
    parser.add_argument('--species', nargs='+', help='species to check')
    parser.add_argument('--prefixes', nargs='+', help='prefixes for circos naming in species')
    parser.add_argument('--old_prefixes', nargs='+', help='prefixes to rename')
    args = parser.parse_args()
    chroms = model.parse_chromosomes(args.file)
    blocks, count_chrs = model.parse_blocks(args.file, True)
    blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
    if len(args.species) != 2:
        raise Exception("Can only draw circos plot for two species")
    if len(args.prefixes) != 2:
        raise Exception("Specify two species prefixes with --prefixes")
    utils.output_for_circos(blocks, args.species, args.prefixes, args.old_prefixes, args.circos_output)