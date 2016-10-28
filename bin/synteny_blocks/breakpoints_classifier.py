#!/usr/bin/env python

import utils
import itertools
from collections import Counter, defaultdict

class Decision:
    true = True
    false = False
    noinfo = None

def get_set_entries(blocks):
    list_of_list_entries = map(lambda x: x.entries, blocks)
    set_of_entries = set(list(itertools.chain(*list_of_list_entries)))
    set_of_species = set(map(lambda x: x.get_specie(), set_of_entries))
    return set_of_species

def get_blocks_ids(blocks):
    return set(map(lambda x: x.id, blocks))

def create_indices(species, threaded_genomes):
    ind = defaultdict(list)
    for g in threaded_genomes.keys():
       chrs = threaded_genomes[g]
       for k in range(len(chrs)):
           for i in range(len(chrs[k])):
               prev_id = -1 if (i == 0) else chrs[k][i-1].block_id
               next_id = -2 if (i == len(chrs[k]) - 1) else chrs[k][i+1].block_id
               ind[(chrs[k][i].block_id, g)].append((prev_id,next_id)) 
    return ind

def run(blocks, print_table=False):
    species = sorted(list(get_set_entries(blocks)))
    threaded_genomes = {}
    for sp in species:
            entries = utils.get_specie_entries(blocks, sp)
            threaded_genomes[sp] = utils.thread_specie_genome(entries)
    index = create_indices(species, threaded_genomes)
    blocks_ids = get_blocks_ids(blocks)
    dupls_num = 0
    blocks_num = 0
    entries_num = 0
    if print_table:
        header = '\t'.join(['breakpoint block']+species)
        print header
    for b in blocks_ids:
        blocks_num += 1
        block_inds = filter(lambda x: x[0] == b, index)
        neighbours = []    
        species_status = {}
        if not print_table:
            for ind in block_inds:
                print index[ind]
        #just linearize two-dimensional data
        #[(prev1, next1), (prev2, next2)] -> [prev1, next1, prev2, next2]
        for ind in block_inds:
            entries_num += len(index[ind])
            if len(index[ind]) > 1:
                dupls_num += len(index[ind])
                index[ind] = [(-3,-3)]
                species_status[ind[1]] = 'DUP'
                #TODO solve duplications!!
            #beware of dupl!
            neighb = index[ind][0]    
            neighbours.append(neighb[0])
            neighbours.append(neighb[1])
        #sort by popularity in descending order
        #and leave only non-ending
        c = Counter(neighbours).most_common()
        c = filter(lambda x:x[0] != -2 and x[0] != -1 and x[0] != -3, c)
        #print c
        #if len is less or equal than two (most popular from left and from right), 
        #then breakpoint is likely to be caused by assembly incompleteness
        if len(c) > 2:
            c2 = Counter(dict(c)).most_common()
            if c[2][1] == c[1][1]:
                if not print_table:
                    print 'cant distinguish two most common!'
                    print
                for ind in block_inds:
                    if not ind[1] in species_status.keys():
                        #could not resolve breakpoint
                        species_status[ind[1]] = 'NA'
                continue
            first_common = c[0][0]
            second_common = c[1][0]
            nodef = set([-1,-2])
            allowable = set([-1,-2,first_common,second_common])
            br = False
            for ind in block_inds:
                #beware of dupl!
                prev,next = index[ind][0]
                if prev in nodef and next in nodef:
                    #the whole block is a full scaffold in the specie
                    species_status[ind[1]] = 'END'
                #    print 'possible breakpoint', ind[1], prev, '-', ind[0], '-', next
                    continue
                if prev in allowable and next in allowable:
                    species_status[ind[1]] = '-'
                    continue
                if not prev in allowable:
                    species_status[ind[1]] = 'BR'
                    if not print_table:
                        print 'breakpoint', ind[1], prev, '-', ind[0]
                        br = True
                if not next in allowable:
                    species_status[ind[1]] = 'BR'
                    if not print_table:
                        print 'breakpoint', ind[1], ind[0], '-', next
                        br = True
            #if br:
            #    print
        if not print_table:    
            print
        if print_table:
            l = str(b)
            if not species_status.keys():
                #in case breakpoint is caused by assembly incompleteness
                continue
            for e in species:
                if not e in species_status.keys():
                    l += '\t'+'not in block'
                else:
                    l+='\t'+species_status[e]
            print l
    if not print_table:
        print 'STAT Also:'
        print 'STAT number of blocks:', blocks_num
        print 'STAT number of entries:', entries_num
        print 'STAT number of dupls (among entries):', dupls_num
        print 'STAT rate of duplications:', float(dupls_num)/entries_num
        

    
