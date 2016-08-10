#!/usr/bin/env python

import os
import argparse

import model
import utils
from blocks_to_paths_processor import BlocksToPathsProcessor
import rearrangements_type
import breakpoints_classifier

def print_out_genome_thread(entries):
    #with open(file_name,'w') as f:
        i = 0
        for c in entries:
            i += 1
            #f.write(str(i)+'\n')
            print i
            for e in c:
                 #f.write('seq_id: ' + str(e.seq_id) + ' block_id: ' + str(e.block_id) + ' strand: '\
                 #+ str(e.strand) + ' start: ' + str(e.start) + ' end: ' + str(e.end) + '\n')
                print 'seq_id: ' + str(e.seq_id) + ' block_id: ' + str(e.block_id) + ' strand: '\
                + str(e.strand) + ' start: ' + str(e.start) + ' end: ' + str(e.end)
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--report_transpositions', action='store_true', help='report transpositions in specie2 related to specie1')
    parser.add_argument('--report_translocations', action='store_true', help='report translocations in specie2 related to specie1')
    parser.add_argument('--report_reversals', action='store_true', help='report reversals in specie2 related to specie1')
    parser.add_argument('--report_duplications', action='store_true', help='search for duplications in each of --species')
    parser.add_argument('--count_breakpoints', action='store_true', help='print number of breakpoints in each of --species')
    parser.add_argument('--species', nargs='+', help='species to check')
    parser.add_argument('--classify_breakpoints', action='store_true', help='find out which species contain breakpoint')
    parser.add_argument('--ref_genome')
    parser.add_argument('--print_out_genomes', help='prints out genomes of --species in terms of blocks')
    parser.add_argument('--filter', help='filter blocks for regions mentioned in bed file')

    args = parser.parse_args()
    chroms = model.parse_chromosomes(args.file)
    blocks, count_chrs = model.parse_blocks(args.file, True)

    if args.filter:
        f_blocks = utils.filter_bed(blocks, args.filter)
    elif args.classify_breakpoints:
        if not args.ref_genome:
            raise Exception('reference genome is needed for breakpoints classification')
        breakpoints = breakpoints_classifier.run(blocks, args.ref_genome)
        for k in breakpoints.keys():
            k[0].print_out()
            k[1].print_out()
            print breakpoints[k]
    elif args.report_duplications:
            for sp in args.species:
                entries = utils.get_specie_entries(blocks, sp)
                entries = utils.thread_specie_genome(entries)
                for c in entries:
                    count_dup = 0
                    dup = rearrangements_type.check_duplications(c, blocks, sp)
                    for e in dup:
                        this_prev = e[0]
                        this_dup = e[1]
                        if not this_prev in map(lambda x:x[1], dup):
                            count_dup += 1
                        print 'duplication:',
                        this_dup.print_out()
                    if count_dup:
                        print 'overall duplications', count_dup

    elif args.report_translocations or args.report_transpositions or args.report_reversals\
             or args.report_duplications:
            if len(args.species) != 2:
                raise Exception("Can evaluate rearrangements only between two species")
            blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
            #get all the entries from specie1
            entries = utils.get_specie_entries(blocks, args.species[0])
            #sort entries by chromosomes for specie1
            specie1 = utils.thread_specie_genome(entries)
            #Tfor testing purposes
            #Ttest_path = '/Users/admin/projects/felidae_comp/analysis/synteny/solenodon/1000/'
            #Tprint_out_genome_thread(args.species[0],specie1,os.path.join(test_path,'tmp1'))
            entries = utils.get_specie_entries(blocks, args.species[1])
            #Tprint_out_genome_thread(args.species[1],utils.thread_specie_genome(entries),os.path.join(test_path,'tmp2'))
            #Texit()
            specie2 = utils.thread_specie_genome(entries)
            specie2_grouped = []
            #group entries in specie2 according to order of blocks on chromosomes
            #in specie1
            for sp in specie1:
                specie2_grouped.append([])
                for y in sp:
                    c = filter(lambda x: x.block_id == y.block_id, entries)
                    specie2_grouped[-1].append(c)
                    if not c:
                        print y.block_id
            specie2 = []
            cnt_empty = 0
            for e in specie2_grouped:
                p = BlocksToPathsProcessor.search_paths(e)
                if not p:
                    cnt_empty += 1
                specie2.append(p)
            print 'unresolved:', cnt_empty
            specie1,specie2 = utils.normalize(specie1, specie2)
            for c in specie2:
                if args.report_transpositions:
                    count_trp = 0
                    trp = rearrangements_type.check_transpositions(c)
                    for e in trp:
                        this_prev = e[0]
                        this_trp = e[1]
                        #count transposition only once if
                        #it occured in neighbouring blocks
                        if not this_prev in map(lambda x: x[1], trp):
                            count_trp += 1
                        print 'transposition:',
                        this_trp.print_out()
                    if count_trp:
                        print 'overall transpositions', count_trp
                if args.report_translocations:
                    trl = rearrangements_type.check_translocations(c)
                    for e in trl:
                        this_prev = e[0]
                        this_trl = e[1]
                        #if not this_prev in all_translocated_entries:
                        #    count_trl += 1
                        print 'whole chromosome:'
                        for x in c:
                            x.print_out()
                        print 'translocation:'
                        for x in e[1]:
                            x.print_out()
                    if trl:
                        print 'overall translocations:', len(trl)
                if args.report_reversals:
                    count_rev = 0
                    rev = rearrangements_type.check_reversals(c)
                    print 'whole chromosome'
                    for x in c:
                        x.print_out()
                    for e in rev:
                        this_prev = e[0]
                        this_rev = e[1]
                        #count reversal only once if
                        #it occured in neighbouring blocks
                        if not this_prev in map(lambda x: x[1], rev):
                            count_rev += 1
                        print 'reversal:',
                        this_rev.print_out()
                    if count_rev:
                        print 'overall reversals', count_rev
    elif args.print_out_genomes :
        blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
        for sp in args.species:
            entries = utils.get_specie_entries(blocks, sp)
            specie_genome = utils.thread_specie_genome(entries)
            print_out_genome_thread(sp)

    elif args.count_breakpoints:
        blocks = utils.filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
        for s in args.species:
            entries = utils.get_specie_entries(blocks, s)
            genome = utils.thread_specie_genome(entries)
            print sum(map(lambda x: len(x)-1, genome))

