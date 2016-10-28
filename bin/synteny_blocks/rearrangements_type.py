from collections import Counter
import itertools
import math
import utils

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

def get_next_entries(list_entries, c):
    rearrangement_ids = []
    for e in list_entries:
        rearrangement_ids.append(c.index(e))
    rearrangement_ids = map(lambda x: x+1, rearrangement_ids)
    rearrangement_prev = []
    for i in rearrangement_ids:
        if i == len(c):
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
        print 'fc:'
        for x in c_seq_id:
            x.print_out() 
        print 'fc sorted according pt:'
        for x in c_seq_id_sorted:
            x.print_out()
        transpositions_up = check_order(c_seq_id, c_seq_id_sorted)
        transpositions_down = check_order(c_seq_id, c_seq_id_sorted[::-1])
        if len(transpositions_up) <= len(transpositions_down):
            transpositions += transpositions_up
        else:
            transpositions += transpositions_down
    #find indices of entries that are previous to the transposed entries
    #transpositions = sorted(transpositions, key = lambda x: x.start, reverse=True)
    tps_prev = get_previous_entries(transpositions, c)
    tps_next = get_next_entries(transpositions, c)
    '''
    for p,t,n in zip(tps_prev, transpositions, tps_next):
        if p:
            p.print_out()
        else:
            print None
        t.print_out()
        if n:
            n.print_out()
        else:
            print None
        print
    '''
    return zip(tps_prev, transpositions, tps_next)

'''
translocation is a change of the genomic position to another chromosome
choose the one interpretation that corresponds to less translocations (by length)
generally the translocations can be counted without information about previous entries
because they are already grouped into lists of unbroken segments of blocks.
however we report translocations in the same manner as reversals and transpositions
for uniformity and for being able to report the breakpoints if needed
'''
def check_translocations(c):
    c_seq_ids =[]
    for e in itertools.groupby(c, lambda x: x.seq_id):
        c_seq_ids.append(list(e[1]))
    lengths = map(lambda y: sum(map(lambda x: math.fabs(int(x.end)-int(x.start)), y)), c_seq_ids)
    ls = zip(lengths, c_seq_ids)
    ls_sorted = sorted(ls, key=lambda x: x[0])
    translocations = map(lambda x: x[1], ls_sorted[:-1])
    main_chrom = ls_sorted[-1][1][0].get_chrom()
    '''
    for e in ls_sorted:
        print e[0]
        for x in e[1]:
            x.print_out()
        print
    '''
    return main_chrom, translocations

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
    if len(c_rev) == len(c) :
        return []
    if c_rev:
        rev_prev = get_previous_entries(c_rev,c)
        return zip(rev_prev,c_rev)
    return []
'''
If block appears in genome several times then it's a duplication
returns species entries that belong to duplicated blocks
'''
def check_duplications(c, blocks, specie):
    l = utils.get_specie_entries(blocks,specie)
    cnt = Counter(map(lambda x: x.block_id, l))
    dups = filter(lambda x: cnt[x.block_id] > 1, c)
    if dups:
        dup_prev = get_previous_entries(dups,c)
        return zip(dup_prev,dups)
    else:
        return []
