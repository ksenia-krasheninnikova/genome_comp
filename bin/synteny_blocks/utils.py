#!/usr/bin/env python

import model
import re

def intersect(entry, bed_entries):
    intersected_bed_entries = []
    for b in bed_entries:
        if (entry.get_specie() == b.genome or b.genome == '') and entry.get_chrom() == b.chrom:
            if not (b.end < entry.start or entry.end < b.start) :
                intersected_bed_entries.append(b)
    return intersected_bed_entries

def filter_bed(blocks, bed):
    bed_entries = model.parse_bed(bed)
    for b in blocks:
        intersected_bed_entries = []
        for e in b.entries:
            intersected_bed_entries += intersect(e, bed_entries)
        if intersected_bed_entries:
            #bed_entries = map(lambda x: bed_entries.remove(x), intersected_bed_entries)
            for bed in intersected_bed_entries:
                print '#',
                bed.print_out()
            b.print_out()
            #if len(intersected_bed_entries) == 0:
            #    return


#traverses blocks and collects all the entries
#related to the specie
def get_specie_entries(blocks, specie):
    specie_entries = []
    for b in blocks:
        #i = 0
        for e in b.entries:
            if e.get_specie() == specie:
                #i += 1
                specie_entries.append(e)
        #if i > 1:
            #print i, specie
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
                #if count_chrs[e.seq_id] > 1:
                if count_chrs[e.seq_id] > 1:
                    upd_entries.append(e)
                    upd_species.add(specie)
        #also count duplications?
        #if so than only blocks when both chromosomes are split counted
        if len(upd_entries) >= len(sps) and len(upd_species) == len(sps):
        #if len(upd_entries) == 1:
        
        #if so then counted also those blocks that partly split but in some
        #species it can be the whole scaffold
        #if upd_entries:
            upd_blocks.append(model.Block(b.id, upd_entries))
    return upd_blocks

def filter_absent_species(blocks, sps):
    upd_blocks = []
    for b in blocks:
        entries = b.entries
        upd_species = []
        for e in entries:
            specie = e.get_specie()
            if specie in sps:
                upd_species.append(specie)
        if len(upd_species) == len(sps):
            upd_blocks.append(b)
    return upd_blocks

#must be fixed: in case of duplications in a genome
#there can be ambiguities in prev entries
def find_prev_block_in_specie(entry,specie):
    for c in specie:
        find = filter(lambda x: x.block_id == entry.block_id, c)
        if len(find) > 1:
            raise Exeption('duplicated entry!')
        if find:
            l = c.index(find[0])
            if l == 0:
                return None
            return c[l-1] 
    raise Exeption('No such block! ' + entry.block_id)

#must be fixed: in case of duplications in a genome
#there can be ambiguities in next entries
def find_next_block_in_specie(entry,specie):
    for c in specie:
        find = filter(lambda x: x.block_id == entry.block_id, c)
        if len(find) > 1:
            raise Exeption('duplicated entry!') 
        if find:
            l = c.index(find[0])
            if l == len(c)-1:
                return None
            return c[l+1] 
    raise Exeption('No such block! ' + entry.block_id)

def output_for_circos(blocks, species, prefixes, old_prefixes, output):
    old_prefix='|'.join(old_prefixes)
    with open(output,'w') as f:
        for b in blocks:
            e = b.entries
            e1 = filter(lambda x: x.seq_id.split('.')[0] == species[0], e)
            if e1:
                e1 = e1[0]
                id = [s.strip() for s in re.split(old_prefix, e1.seq_id.split('.')[1])][1]
                name1 = prefixes[0]+id
                if e1.strand == '+':
                    start1 = e1.start
                    end1 = e1.end
                else:
                    end1 = e1.start
                    start1 = e1.end
                f.write(name1 + ' ' + str(start1) + ' ' + str(end1)) 
            e2 = filter(lambda x: x.seq_id.split('.')[0] == species[1], e)
            if e2:
                e2 = e2[0]
                id = [s.strip() for s in re.split(old_prefix, e2.seq_id.split('.')[1])][1]
                name2 = prefixes[1]+id
                if e2.strand == '+':
                    start2 = e2.start
                    end2 = e2.end
                else:
                    end2 = e2.start
                    start2 = e2.end
                if e1:
                    f.write(' ')
                f.write(name2 + ' ' + str(start2) + ' ' + str(end2)+ '\n')
            elif e1:
                f.write('\n')
