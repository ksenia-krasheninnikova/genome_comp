#!/usr/bin/env python

import utils
import itertools

class Decision:
    true = True
    false = False
    noinfo = None

def get_map_block_id_to_entries(blocks):
    b2e = {}
    for b in blocks:
        b2e[b.id] = b.entries
    return b2e

def get_next_entry(e, whole_chrom):
    return whole_chrom[whole_chrom.index(e) + 1]


def get_map_entry_to_chromosome(threaded_genomes):
    e2c = {}
    for sp in threaded_genomes.keys():
        for ch in threaded_genomes[sp]:
            for e in ch:
                e2c[e] = ch
    return e2c

def get_set_entries(blocks):
    list_of_list_entries = map(lambda x: x.entries, blocks)
    set_of_entries = set(list(itertools.chain(*list_of_list_entries)))
    set_of_species = set(map(lambda x: x.get_specie(), set_of_entries))
    return set_of_species

def run(blocks, ref_specie) :
    species = get_set_entries(blocks)
    threaded_genomes = {}
    for sp in species:
            entries = utils.get_specie_entries(blocks, sp)
            threaded_genomes[sp] = utils.thread_specie_genome(entries)

    b2e_map = get_map_block_id_to_entries(blocks)
    e2c_map = get_map_entry_to_chromosome(threaded_genomes)

    ref_genome = threaded_genomes[ref_specie]
    breakpoints = {}
    #iterate over chromosomes in the reference genome
    for ref_chrom in ref_genome:
        #take reference chromosome
        for i in range(len(ref_chrom)-1):
            e_ref = ref_chrom[i]
            e_next = ref_chrom[i+1]
            breakpoints[(e_ref, e_next)] = []
            ref_block_id = e_ref.block_id
            #get entries of not reference species for blocks the same that ones on reference chromosome
            not_ref_entries = filter(lambda e_x: not e_x.get_specie() == ref_specie, b2e_map[ref_block_id])
            for sp in species:
                if sp == ref_specie:
                    continue
                not_ref_next_entries = set()
                is_breakpoint_this_specie = Decision.true
                for e_not_ref in filter(lambda e_x: e_x.get_specie() == sp, not_ref_entries):
                    e_not_ref_chrom = e2c_map[e_not_ref]
                    if len(e_not_ref_chrom) == 1 :
                        #can't judge about breakpoints here
                        #one block is the whole chromosome
                        is_breakpoint_this_specie = Decision.noinfo
                        break
                    if e_not_ref_chrom.index(e_not_ref) == len(e_not_ref_chrom) - 1:
                        #can't judge about breakpoints here
                        #looks like we got to the end of list of entries
                        is_breakpoint_this_specie = Decision.noinfo
                    else:
                        e_not_ref_next_entry = get_next_entry(e_not_ref, e_not_ref_chrom)
                        not_ref_next_entries.add(e_not_ref_next_entry)
                if is_breakpoint_this_specie == Decision.true :
                    for e_not_ref_next_entry in not_ref_next_entries:
                        if e_not_ref_next_entry.block_id == e_next.block_id:
                            is_breakpoint_this_specie = Decision.false
                breakpoints[(e_ref, e_next)].append((sp, is_breakpoint_this_specie))
    return breakpoints




