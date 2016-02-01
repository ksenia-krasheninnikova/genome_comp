#!/hive/groups/recon/local/bin/python

import argparse
import math
import itertools

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
    #if len(entries) == 1:
    #    edges = [(entries[0],entries[0],0)]
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
        path = dfs(e, e, edges, 0, len(entries)-2)
        if len(path) > 1:
            print 'Alternative solutions!'
            return [] 
        thread = [path[0][0][0]]
        for p in path[0]:
            thread.append(p[1])
        threads.append(thread)
    if len(threads) > 1:
        print 'Alternative solutions!'
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
        for e in entries:
            #seq_id = e.seq_id 
            #specie = chroms[int(seq_id)].get_specie()
            specie = e.get_specie()
            if specie in sps:
                if count_chrs[e.seq_id] > 1:
                    upd_entries.append(e)
        #also count duplications?
        if len(upd_entries) >= len(sps):
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

'''
#use this to sort blocks by specie chromosomes
#def get_first_entry_for_specie_in_block(block, specie):
    for e in block.entries:
        if chroms[e.seq_id].description:
            return e.seq_id

#splits entries by chromosomes and sorts entries
#on each chromosome according to start pos
def sort_blocks_on_chromosomes(blocks, chroms, specie):
    sorted_blocks = sorted(blocks, key=lambda x: get_any_entry_for_specie(x, specie))
    acc_blocks = []
    chromosomes = []
    prev_id = -1
    for e in sorted_blocks:
        if get_any_entry_for_specie(e,specie) != prev_id:
            chromosomes.append(acc_blocks)
            acc_blocks = [e]
            prev_id = e.entries[0].seq_id
            continue
        acc_blocks.append(e)
    return chromosomes
'''

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
            print c1
            print c2
            continue
        print 'c1',c1
        print 'c2',c2
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
transposition is a change of the genomic location on the same chromosome
idea: transposition is a change of order in a sorted list
collect to possible interpretations of transpositions:
in case sorted list in increasing or decreasing
(this corresponds to entries coming in different order along a chromosome)
choose the one interpretation that corresponds to less transpositions
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
    return transpositions

'''
translocation is a change of the genomic position to another chromosome
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
    return map(lambda x: x[1], ls_sorted[:-1])
    '''
    for e in c_seq_ids:
        if prev_seq_id and prev_seq_id != e.seq_id:
            pos = c.index(e)
            rest_c = c[:pos]+c[pos+1:]
            rest_length = sum(map(lambda x: math.abs(int(x.end)-int(x.start)), rest_c))
            e_length = mas.abs(int(e.end) - int(e.start))
            if e_length > rest_length:
                translocations.append(e)
            else:
                translocations.append(rest_c)
        prev_seq_id = e.seq_id
    return translocations
    '''

'''
reversal is the change of strand
here we use the normalization of two genomes. i.e. we brought the genome of specie1
to the form where all the entries have '+' strand, changing also the strand of 
the corresponding entries in specie2
every '-' is called reversal
'''
def check_reversals(c):
    return filter(lambda x: x.strand == '-', c)

def identify_rearrangements_by_type(sp1, sp2) :
    #print 'homolog:', c2_homolog.seq_id, c2_homolog.block_id
    homologs = zip(sp1,sp2)
    for c1,c2 in homologs:
        if not c1 or not c2:
            print 'skipping empty chromosome'
            print c1
            print c2
            continue
        blocks = zip(c1,c2)
        #for x in c1:
        '''
        PRINT CHROMOSOMES 
        print 'c1'
        for e in c1:
            e.print_out(), 
        print 'c2'
        for e in c2: 
            e.print_out(),
        #    x.print_out()
        prev_e2 = ''
        reversal = False
        for e1,e2 in blocks:
            if prev_e2:
                if prev_e2.strand != e2.strand:
                    if not reversal:
                        reversal = True
                        #we consider each '-' as reversal because we suppose 
                        #that the genomes are normalized
                        if prev_e2.strand == '-':
                            print 'reversal', prev_e2.seq_id+':', prev_e2.end
                        elif  e2.strand == '-':
                            print 'reversal', e2.seq_id+':', e2.start    
                else:
                    reversal = False
                if prev_e2.seq_id != e2.seq_id:
                    print 'translocation start:', e2.seq_id+':',e2.start
                elif prev_e2.end > e2.end and get_order(c2,prev_e2.seq_id) == '+' or prev_e2.end < e2.end and get_order(c2,prev_e2.seq_id) == '-':
                    print 'transposition start:', e2.seq_id+':',e2.start
            prev_e2 = e2
        print '-----------' 
        '''
    # reversal: names and chromosomes and surroundings of blocks are the same, but signs of the blocks are different
    # translocation: different chromosomes
    # transposition: different order on the same chromosome
 
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    parser.add_argument('--report_breakpoints', action='store_true',help='report breakpoints among two --species')
    parser.add_argument('--species', nargs='+', help='species to check')
    args = parser.parse_args()
    chroms = parse_chromosomes(args.file)
    if args.report_breakpoints:
        blocks,count_chrs = parse_blocks(args.file, True)
        blocks = filter_unsplitted_chromosomes(blocks, count_chrs, args.species)
        #get all the entries from specie1
        entries = get_specie_entries(blocks, args.species[0])
        #sort entries by chromosomes for specie1
        specie1 = thread_specie_genome(entries)
        entries = get_specie_entries(blocks, args.species[1])
        #specie2 = thread_specie_genome(entries)
        specie2_grouped = []
        #group entries in specie2 according to order of blocks on chromosomes
        #in specie1
        for sp in specie1:
            specie2_grouped.append([])
            for y in sp:
                specie2_grouped[-1].append(filter(lambda x: x.block_id == y.block_id, entries))
        specie2 = []
        #this is for testing purpose
        #for i in range(len(specie1)):
        #    chr_sp1 = specie1[i]
        #    chr_sp2 = specie2_grouped[i]
        #    for j in range(len(chr_sp1)):
        #        block_sp1 = chr_sp1[j]
        #        block_sp1.print_out()
        #    print
        #    for j in range(len(chr_sp1)):
        #        blocks_sp2 = chr_sp2[j]
        #        for x in blocks_sp2:
        #            x.print_out()
        #    print '---'
        #exit()
        ##
        for e in specie2_grouped:
            specie2.append(search_paths(e))
        #search_paths(specie2_grouped[182])
        #specie2 = find_homologous_chromosomes(specie1, specie2)
        specie1,specie2 = normalize(specie1, specie2)
        #identify_rearrangements_by_type(specie1, specie2)
        for c in specie2:
            for e in check_transpositions(c):
                print 'transposition:',
                e.print_out()
            for e in check_translocations(c):
                print 'translocation:'
                for x in e:
                    x.print_out()
            for e in check_reversals(c):
                print 'reversal',
                e.print_out()
    else:
        blocks = parse_blocks(args.file)
        
