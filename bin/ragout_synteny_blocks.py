#!/hive/groups/recon/local/bin/python

import argparse


SPLITTER = '--------------------------------------------------------------------------------'

class Chromosome:
    def __init__(self,seq_id,size,description):
        self.seq_id = seq_id
        self.size = size
        self.description = description

    def print_out(self):
        print self.seq_id, self.size, self.description

class Entry:
    def __init__(self,seq_id,strand,start,end,length):
        self.seq_id = seq_id
        self.strand = strand
        if start > end:
            self.start = end
            self.end = start
        else:
            self.start = start
            self.end = end
        self.length = length

    def print_out(self):
        print self.seq_id, self.strand, self.start, self.end, self.length

class Block:
    def __init__(self, id,entries):
        self.id = id
        self.entries = entries

    def print_out(self):
        print self.id, self.entries

def parse_chromosomes(f):
    chroms = []
    with open(f) as blocks_file:
        blocks_file.readline()
        for line in blocks_file:
            line = line.strip().split()
            if (len(line)!= 3) :
                return chroms
            chroms.append(Chromosome(int(line[0]), int(line[1]), line[2]))
    return chroms

def parse_blocks(f):
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
                entries.append(Entry(int(line[0]),line[1],int(line[2]),int(line[3]),int(line[4])))
    return blocks
                    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='blocks_coords.txt')
    args = parser.parse_args()
    chroms = parse_chromosomes(args.file)
    for e in chroms[:10]:
        e.print_out()
    blocks = parse_blocks(args.file)
    for e in blocks[:10]:
        e.print_out()
