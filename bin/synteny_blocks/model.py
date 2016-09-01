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

    def get_chrom(self):
        return self.seq_id.split('.')[1]

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

class BED_Entry:
    def __init__(self, genome, chrom, start, end):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.end = end

    def to_string(self):
        return  '_'.join(map(str,[self.genome+'.'+self.chrom, self.start, self.end]))

    def print_out(self):
        print ' '.join(map(str,[self.genome+'.'+self.chrom, self.start, self.end]))

def parse_bed(bed_file, margin=0):
    beds = []
    with open(bed_file) as bed:
        for line in bed:
            line = line.strip().split()
            delim = line[0].find('.')
            genome = line[0][:delim]
            chrom = line[0][delim+1:]
            #genome,chrom = line[0].split('.')
            beds.append(BED_Entry(genome,chrom,int(line[1])-margin,int(line[2])+margin))
            #beds.append(BED_Entry('',line[0],int(line[1]),int(line[2])))
    return beds

class MAF_Entry:
    def __init__(self, genome, chrom, start, length, strand, all_length, seq):
        self.genome = genome
        self.chrom = chrom
        self.start = start
        self.length = length
        self.strand = strand
        self.all_length = all_length
        self.seq = seq
        self.__init_global_coords()

    def __init_global_coords(self):
        if self.strand == '+':
            self.global_start = self.start
            #self.global_end = self.start + self.length + 1
            self.global_end = self.start + self.length
        elif self.strand == '-':
            #self.global_end = self.all_length - self.start - 1
            self.global_end = self.all_length - self.start 
            self.global_start = self.all_length - self.start - self.length

    def print_out(self):
        #print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.start,\
        #                        self.length, self.strand, self.all_length, self.seq]))
        print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.global_start,\
                                self.global_end, self.strand, self.all_length, self.seq]))

    #this prints canonical maf entry
    def print_out_local_coords(self):
        print ' '.join(map(str,['s', self.genome + '.' + self.chrom, self.start,\
                                self.length, self.strand, self.all_length, self.seq]))

def check_maf_for_no_overlaps(maf_entries):
    #i don't keep any storage for processed entries
    #in order not to load the whole maf into memory
    for e1 in maf_entries:
        for e2 in maf_entries:
            if e1 == e2 or e1.genome != e2.genome or e1.chrom != e2.chrom:
                continue
            if e1.start < e2.start and e1.start + e1.length > e2.start:
                print 'e2 overlaps e1:'
                e1.print_out_local_coords()
                e2.print_out_local_coords()
            if e1.global_start < e2.global_start and e1.global_end > e2.global_start: 
                print 'e2 overlaps e1 global:'
                e1.print_out()
                e2.print_out()
            if e2.start < e1.start and e2.start + e2.length > e1.start:
                print 'e1 overlaps e2:'
                e1.print_out_local_coords()
                e2.print_out_local_coords()
            if e2.global_start < e1.global_start and e2.global_end > e1.global_start: 
                print 'e1 overlaps e2 global:'
                e1.print_out()
                e2.print_out()


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
                if seq_id.find('random') > -1 or seq_id.find('hap') > -1:
                    continue
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

