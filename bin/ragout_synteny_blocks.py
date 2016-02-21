#!/hive/groups/recon/local/bin/python

import argparse
import math
import itertools
from collections import Counter
import re
import os

SPLITTER = '-----------------------------------'

    
def find_graph_edges(entries):
    i = 0
    j = 1
    edges = []
    while i < len(entries) - 1:
        v = entries[i]
        v_next = entries[i+1]
        appended = False
        for e in v:
            for e_next in v_next:
                if e.seq_id == e_next.seq_id:
                    edges.append((e,e_next,i))
                    appended = True
        #in case we didn't found the next block on the same chromosome
        #add all the following blocks (on all the chromosomes)
        if not appended:
            for e_next in v_next:
                edges.append((e,e_next,i))
        i += 1
    ##in case we didn't filter out unsplitted chromosomes
    ##we add a loop
    if len(entries) == 1:
        edges.append((entries[0][0],entries[0][0],0))
    return edges

def dfs(v, v_prev, edges, level, max_level):
    new_paths = []
    e_level = []
    for e in edges:
        if e[2] == level:
            e_level.append(e) 
            if level == max_level:
                new_paths.append([(e[0],e[1])])
            else:
                for l in dfs(e[1], e[0], edges, level+1, max_level):
                    new_paths.append([(e[0],e[1])] + l)
    #there can be multiple paths in case of at some point the neighbouring blocks 
    #have multiple alternative chromosomes like this:
    #a b d
    #  c
    #here a,b,c,d - chromosomes, and alternative paths are a b d and a c d
    #also, b and c can belong to different locations on the same chromosome
    return new_paths

#entries - a chromosome consists in a list of blocks
#each block is also a list of possible alternative locations of this block
def search_paths(entries):
    '''
    #PRINT ENTRIES
    print 'entries:'
    for e in entries:
        for p in e:
            p.print_out()
        print '---'
    '''
    edges = find_graph_edges(entries)
    '''
    print 'edges:'
    for e in edges:
        print e[0].print_out(), e[1].print_out(), e[2]
        print
    '''
    threads = []
    for e in entries[0]:
        #-1 because max level number is len(entries)-1
        #-1 because number of edges is (number of entries)-1
        #in total -2
        max_level = len(entries)-2
        if len(entries) == 1:
            max_level = 0
        path = dfs(e, e, edges, 0, max_level)
        if len(path) > 1:
            print 'Alternative solutions!'
            print 'returning empty chromosome'
            return [] 
        thread = [path[0][0][0]]
        for p in path[0]:
            thread.append(p[1])
        threads.append(thread)
    if len(threads) > 1:
        print 'Alternative solutions!'
        print 'returning empty chromosome'
        return []
    return threads[0]

class Chromosome:
    def __init__(self,seq_id,size,description):
        self.seq_id = seq_id
        self.size = size
        self.description = description

    def get_specie(self):
        return self.description.split('.')[0] 

    def print_out(self):
        print self.seq_id, self.size, self.description

class Entry:
    def __init__(self,seq_id,strand,start,end,length):
        self.seq_id = seq_id
        self.strand = strand
        self.block_id = -1
        if start > end:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end
        self.length = length

    def set_block_id(self, block_id):
        self.block_id = block_id

    def get_specie(self):
        return self.seq_id.split('.')[0]

    def print_out(self):
        print 'seq_id:',self.seq_id, 'block_id:', self.block_id, 'strand:', self.strand, 'start:', self.start, 'end:', self.end, 'length:', self.length

class Block:
    def __init__(self, id,entries):
        self.id = id
        self.entries = entries
        for e in self.entries:
            e.set_block_id(id)

    def print_out(self):
        print self.id 
        for e in self.entries:
            e.print_out()

def parse_chromosomes(f):
    chroms = {}
    with open(f) as blocks_file:
        blocks_file.readline()
        for line in blocks_file:
            line = line.strip().split()
            if (len(line)!= 3) :
                return chroms
            chroms[int(line[0])] = Chromosome(int(line[0]), int(line[1]), line[2])
    return chroms

def parse_blocks(f, count_c=False):
    count_chrs = {}
    with open(f) as blocks_file:
        blocks_section = False
        entries = []
        blocks = []
        id = ''
        for line in blocks_file:
            line = line.strip()
            if 'Block #' in line:
                blocks_section = True
            if blocks_section:
                if 'Block #' in line:
                    if id:
                        blocks.append(Block(id, entries))
                    id = int(line[7:])
                    entries = []
                    continue
                if 'Seq_id' in line or SPLITTER in line: 
                    continue
                line = line.split()
                #seq_id = int(line[0])
                seq_id = line[0]
                entries.append(Entry(seq_id,line[1],int(line[2]),int(line[3]),int(line[4])))
                if count_c:
                    if  seq_id in count_chrs.keys():
                        count_chrs[seq_id] += 1
                    else:
                        count_chrs[seq_id] = 1
        blocks.append(Block(id, entries))
    if count_c:
        return blocks, count_chrs
    else:
        return blocks
                    
    
def filter_unsplitted_chromosomes(blocks, count_chrs, sps):
    upd_blocks = []
    for b in blocks:
        entries = b.entries
        upd_entries = []
        upd_species = set()
        for e in entries:
            #seq_id = e.seq_id 
            #specie = chroms[int(seq_id)].get_specie()
            specie = e.get_specie()
            if specie in sps:
                if count_chrs[e.seq_id] > 1:
                    upd_entries.append(e)
                    upd_species.add(specie)
        #also count duplications?
        if len(upd_entries) >= len(sps) and len(upd_species) == len(sps):
        #if len(upd_entries) == 1:
        #if upd_entries:
            upd_blocks.append(Block(b.id,upd_entries))
    return upd_blocks        
 
#traverses blocks and collects all the entries
#related to the specie
def get_specie_entries(blocks, specie):
    specie_entries = []
    for b in blocks:
        for e in b.entries:
            if e.get_specie() == specie:
                specie_entries.append(e)
    return specie_entries
            
def thread_specie_genome(specie_entries):
    genome = []
    chromosomes_names = set(map(lambda x: x.seq_id, specie_entries))
    chromosomes = []
    # group entries by chromosomes
    for c in chromosomes_names:
        c_entries = filter(lambda x: x.seq_id==c, specie_entries)
        chromosomes.append(c_entries)
    # sort entries in chromosomes by position
    sorted_chromosomes = []
    for c in chromosomes:
        sorted_c = sorted(c, key=lambda x: x.start)
        sorted_chromosomes.append(sorted_c)
    return sorted_chromosomes


def find_distance(e1,e2):
    if e1.seq_id != e2.seq_id:
        raise Exception('Can find nucleotide distance for blocks from different chromosomes!')
    #one is positive the other is negative
    return max(e2.start - e1.end, e1.start - e2.end)

#normalization means we revert all the negative-strand blocks of the chromosome in specie1
#and change the strand of the corresponding block in specie2
#this is needed in order to search for reversals only in specie2 related to specie1
def normalize(specie1, specie2):
    specie1_upd = []
    specie2_upd = []
    for i in range(len(specie1)):
        c1 = specie1[i]
        c2 = specie2[i]
        if not c1 or not c2:
            print 'skipping empty chromosome'
            #print c1
            #print c2
            continue
        for j in range(len(c1)):
            if c1[j].strand == '-':
                c1[j].strand = '+'
                if c2[j].strand == '-':
                    c2[j].strand = '+'
                elif c2[j].strand == '+':
                    c2[j].strand = '-' 
        specie1[i] = c1
        specie2[i] = c2
    return specie1,specie2

def get_order(c, seq_id):
    c = filter(lambda x: x.seq_id == seq_id, c)
    ends = map(lambda x: x.end, c)
    if sorted(ends) == ends:
        return '+'
    else: 
        return '-'


def check_order(list_a, sorted_list):
    transpositions = [] 
    for i in range(len(list_a)):
        if list_a[i] != sorted_list[i]:
            transpositions.append(sorted_list[i])
            list_a.remove(sorted_list[i])
            list_a.insert(i,sorted_list[i])
    return transpositions

'''
we need it for two reasons:
1. in order to count rearrangements we look at the previous entry and
check if it's rearranged. sometimes we can find multiple possible previous entries.
in this case we need to know which one is previous in this case
2. in order to be able to report the breakpoints
'''
def get_previous_entries(list_entries, c):
    rearrangement_ids = []
    for e in list_entries:
        rearrangement_ids.append(c.index(e))
    rearrangement_ids = map(lambda x: x-1, rearrangement_ids)
    rearrangement_prev = []
    for i in rearrangement_ids:
        if i == -1:
            rearrangement_prev.append(None)
        else:
            rearrangement_prev.append(c[i])
    return rearrangement_prev
 

'''    
transposition is a change of the genomic location on the same chromosome
idea: transposition is a change of order in a sorted list
collect to possible interpretations of transpositions:
in case sorted list in increasing or decreasing
(this corresponds to entries coming in different order along a chromosome)
choose the one interpretation that corresponds to less transpositions
returns list of (prev entry, rearranged entry) s. get_previous_entries
'''
def check_transpositions(c):
    transpositions = []
    seq_ids = map(lambda x: x.seq_id, c)
    seq_ids = set(seq_ids)
    #consider blocks separately for each chromosome
    for s in seq_ids:
        c_seq_id = filter(lambda x: x.seq_id == s, c)
        c_seq_id_sorted = sorted(c_seq_id, key = lambda x: x.start)
        if c_seq_id == c_seq_id_sorted or c_seq_id == c_seq_id_sorted[::-1]:
            continue
        transpositions_up = check_order(c_seq_id, c_seq_id_sorted)
        transpositions_down = check_order(c_seq_id, c_seq_id_sorted[::-1])
        if len(transpositions_up) <= len(transpositions_down):
            transpositions += transpositions_up
        else:
            transpositions += transpositions_down
    #find indices of entries that are previous to the transposed entries 
    tps_prev = get_previous_entries(transpositions, c)
    return zip(tps_prev, transpositions)

'''
translocation is a change of the genomic position to another chromosome
choose the one interpretation that corresponds to less translocations (by length)
generally the translocations can be counted without information about previous entries
because they are already grouped into lists of unbroken segments of blocks.
however we report translocations in the sema manner as reversals and trnaspositions
for uniformity and for being able to report the breakpoints if needed
'''
def check_translocations(c):
    #seq_ids = map(lambda x: x.seq_id, c)
    #seq_ids = set(seq_ids)
    c_seq_ids =[]
    for e in itertools.groupby(c, lambda x: x.seq_id):
        c_seq_ids.append(list(e[1]))
    #prev_seq_id = ''
    #for e in c:
    lengths = map(lambda y: sum(map(lambda x: math.fabs(int(x.end)-int(x.start)), y)), c_seq_ids)
    ls = zip(lengths, c_seq_ids)
    ls_sorted = sorted(ls, key=lambda x: x[0])
    translocations = map(lambda x: x[1], ls_sorted[:-1])
    trans_prev = []
    for tl in translocations:
        trans_prev.append(get_previous_entries(tl, c))
    return zip(trans_prev,translocations)

'''
reversal is the change of strand
here we use the normalization of two genomes. i.e. we brought the genome of specie1
to the form where all the entries have '+' strand, changing also the strand of 
the corresponding entries in specie2
every '-' is called reversal
but if the whole 'chromosome' is '-' than nothing is reversed
'''
def check_reversals(c):
    c_rev = filter(lambda x: x.strand == '-', c)
    reversals = []
    if len(c_rev) < len(c):
        reversals = c_rev
    if reversals:
        rev_prev = get_previous_entries(reversals,c)
        return zip(rev_prev,reversals)
    else:
        return []
#def check_duplications(blocks, specie):
#    l = get_specie_entries(blocks,specie)
#    c = Counter(map(lambda x: x.block_id, l))
#    dupl_block_ids = filter(lambda x: c[x] > 1, c.keys())
#    return filter(lambda x: x.block_id in dupl_block_ids, l)
'''
If block appears in genome several times then it's a duplication
returns species entries that belong to duplicated blocks
'''
def check_duplications(c, blocks, specie):
    l = get_specie_entries(blocks,specie)
    cnt = Counter(map(lambda x: x.block_id, l))
    dups = filter(lambda x: cnt[x.block_id] > 1, c)
    if dups:
        dup_prev = get_previous_entries(dups,c)
        return zip(dup_prev,dups)
    else:
        return []

def output_for_circos(blocks, species, prefixes, old_prefixes, output):
    old_prefix='|'.join(args.old_prefixes)
    with open(output,'w') as f:
        for b in blocks:
            e = b.entries
            e1 = filter(lambda x: x.seq_id.split('.')[0] == species[0], e)[0]
            e2 = filter(lambda x: x.seq_id.split('.')[0] == species[1], e)[0]
            id = [s.strip() for s in re.split(old_prefix, e1.seq_id.split('.')[1])][1]
            name1 = prefixes[0]+id
            id = [s.strip() for s in re.split(old_prefix, e2.seq_id.split('.')[1])][1]
            name2 = prefixes[1]+id
            if e1.strand == '+':
                start1 = e1.start
                end1 = e1.end
            else:
                end1 = e1.start
                start1 = e1.end
            if e2.strand == '+':
                start2 = e2.start
                end2 = e2.end
            else:
                end2 = e2.start
                start2 = e2.end
            f.write(name1 + ' ' + str(start1) + ' ' + str(end1) + ' ' + name2 + ' ' + str(start2) + ' ' + str(end2)+ '\n')

def print_out_genome_thread(species, entries, file_name):
    with open(file_name,'w') as f:
        i = 0
        for c in entries:
            i += 1
            f.write(str(i)+'\n')
            for e in c:
                f.write('seq_id: ' + str(e.seq_id) + ' block_id: ' + str(e.block_id) + ' strand: '\
                + str(e.strand) + ' start: ' + str(e.start) + ' end: ' + str(e.end) + '\n')
    
#def get_stat_prev(entries):
#    upd_entries = {}
#    for c in entries:
#        a = dict(zip(c,[None]+c[:-1]))
#        for e in a.keys():
#            if e in upd_entries.keys():
#                upd_entries[e] = upd_entries[e] + a[e]
#            else:
#                upd_entries[e] = a[e]
#        #upd_entries.update(dict(zip(c,[None]+c[:-1])))
#    return upd_entries

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--report_transpositions', action='store_true', help='report transpositions in specie2 related to specie1')
    parser.add_argument('--report_translocations', action='store_true', help='report translocations in specie2 related to specie1')
    parser.add_argument('--report_reversals', action='store_true', help='report reversals in specie2 related to specie1')
    parser.add_argument('--report_duplications', action='store_true', help='search for duplications in each of --species')
    parser.add_argument('--count_breakpoints', action='store_true', help='print number of breakpoints in each of --species')
    parser.add_argument('--circos_output', help='output the --species for circos plot')
    parser.add_argument('--species', nargs='+', help='species to check')
    parser.add_argument('--prefixes', nargs='+', help='prefixes for circos naming in species')
    parser.add_argument('--old_prefixes', nargs='+', help='prefixes to rename')
    args = parser.parse_args()
    chroms = parse_chromosomes(args.file)
    blocks, count_chrs = parse_blocks(args.file, True)
    if args.report_duplications:
        for sp in args.species:
            entries = get_specie_entries(blocks, sp)
            entries = thread_specie_genome(entries)
            for c in entries:
                count_dup = 0
                dup = check_duplications(c, blocks, sp)
                for e in dup:
                    this_prev = e[0]
                    this_dup = e[1]
                    if not this_prev in map(lambda x:x[1], dup):
                        count_dup += 1
                    print 'duplication:',
                    this_dup.print_out()
                if count_dup:
                    print 'overall duplications', count_dup
    blocks = filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
    if args.report_translocations or args.report_transpositions or args.report_reversals\
         or args.report_duplications:
        if len(args.species) != 2:
            raise Exception("Can evaluate rearrangements only between two species")
        #get all the entries from specie1
        entries = get_specie_entries(blocks, args.species[0])
        #sort entries by chromosomes for specie1
        specie1 = thread_specie_genome(entries)
        #Tfor testing purposes
        #Ttest_path = '/hive/groups/recon/projs/felidae_comp/bin'
        #Tprint_out_genome_thread(args.species[0],specie1,os.path.join(test_path,'tmp1'))
        entries = get_specie_entries(blocks, args.species[1])
        #Tprint_out_genome_thread(args.species[1],thread_specie_genome(entries),os.path.join(test_path,'tmp2'))
        specie2 = thread_specie_genome(entries)
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
        for e in specie2_grouped:
            specie2.append(search_paths(e))
        specie1,specie2 = normalize(specie1, specie2)
        for c in specie2:
            '''
            if args.report_duplications:
                count_dup = 0
                dup = check_duplications(c, blocks, args.species[1])
                for e in dup:
                    this_prev = e[0]
                    this_dup = e[1]
                    if not this_prev in map(lambda x:x[1], dup):
                        count_dup += 1
                    print 'duplication:',
                    this_dup.print_out()
                if count_dup:
                    print 'overall duplications', count_dup
                    '''
            if args.report_transpositions:
                count_trp = 0
                trp = check_transpositions(c)
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
                count_trl = 0
                trl = check_translocations(c)
                all_translocated_entries = map(lambda x: x[1], trl)
                all_translocated_entries = sum(all_translocated_entries, [])
                for e in trl:
                    this_prev = e[0]
                    this_trl = e[1]
                    if not this_prev in all_translocated_entries:
                        count_trl += 1
                    print 'translocation:'
                    for x in e[1]:
                        x.print_out()
                if count_trl: 
                    #count_trl == len(e[0]) 
                    #for explanation s. docs to check_translocations()
                    print 'overall translocations', count_trl
            if args.report_reversals:
                count_rev = 0
                rev = check_reversals(c)
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
    if args.count_breakpoints:
        for s in args.species:
            entries = get_specie_entries(blocks, s)
            genome = thread_specie_genome(entries)
            print sum(map(lambda x: len(x)-1, genome))
    if args.circos_output:
        if len(args.species) != 2:
            raise Exception("Can only draw circos plot for two species")
        if len(args.prefixes) != 2:
            raise Exception("Specify two species prefixes with --prefixes")
        output_for_circos(blocks, args.species, args.prefixes, args.old_prefixes, args.circos_output)

