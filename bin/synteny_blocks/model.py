#!/usr/bin/env python

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

    def equals(self, e):
        return self.seq_id == e.seq_id and self.strand == e.strand\
                and self.start == e.start and self.end == e.end\
                and self.length == e.length

    def equals_to_list(self, e_list):
        for x in e_list:
            if not self.equals(x):
                return False
        return True

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

    def get_species(self):
        return set(map(lambda x: x.get_specie(), self.entries))

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

SPLITTER = '-----------------------------------'

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

