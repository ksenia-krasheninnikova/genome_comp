#!/usr/bin/env python

import os
import sys
import matplotlib
#matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from model import Entry

def parse_genome(genome_path):
    chrs = []
    lengths = []
    names = set()
    with open(genome_path) as f:
        for line in f:
            line = line.strip().split()
            if len(line) == 1:
                chrs.append([])
                lengths.append(0)
            else:
                e= Entry(line[1],line[5],int(line[7]),int(line[9]),int(line[9])-int(line[7]))
                e.set_block_id(int(line[3]))
                chrs[-1].append(e)
                lengths[-1] += e.length
                names.add(e.seq_id.split('.')[1])
    return chrs,lengths

#genomes is the list of chrs (chrs is a list of chromosomes of a genome)
def karyoplot(genome, lengths, metadata={}, part=1):
    fig = plt.figure(1,figsize=(8, 6))
    ax = fig.add_subplot(111, aspect='equal')
    color = 'blue'
    breath = 0.05
#    block_length_param = 0.0000001
    block_length_param = 5e-08
    x_dist = 0.05
    y_dist = 0.05
    x_start = x_dist
    for c,l in zip(genome,lengths):
        y_start = y_dist
        for e in c:
            #y_size = float(e.length) / l * chrom_length_param 
            y_size = float(e.length) * block_length_param 
            r = Rectangle((x_start, y_start), breath, y_size, facecolor = color, alpha=None)
            y_start += y_size
            print r
            ax.add_patch(r)
        center_x = x_start + breath/2.0
        radius = breath/2.0
        theta1 = 0.0
        theta2 = 180.0
        print center_x, y_start, y_dist
        w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
        w2 = Wedge((center_x, y_dist), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
        ax.add_patch(w1)
        ax.add_patch(w2)
        x_start += breath + x_dist
        break

    fig.savefig('test.png', dpi=90)

'''
        #Plot semicircles at the beginning and end of the chromosomes
        center_x = x_start + (x_end-x_start)/2.0
        radius = (x_end-x_start)/2.0
        theta1 = 0.0
        theta2 = 180.0
        w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
        w2 = Wedge((center_x, y_end), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
        ax.add_patch(w1)
        ax.add_patch(w2)
        ax.plot([x_start, x_start], [y_start, y_end], ls='-', color='black')
        ax.plot([x_end, x_end], [y_start, y_end], ls='-', color='black')

        #Plot metadata
        if chromosome in metadata:
            for md in metadata[chromosome]:
                ax.plot([x_end + (DIM*0.015)], [y_start + (y_end-y_start) * (md/chromosome_length)], '.', color='black')

        ax.text(center_x, y_end - (DIM * 0.07), chromosome)
'''



if __name__ == '__main__':
    #import urllib
    #fn = 'karyotype_hg19.txt'
    #url = 'http://pastebin.com/raw.php?i=6nBX6sdE'
    #if not os.path.exists(fn):
    #    print 'Downloading %s to local file: %s' % (url, fn)
    #    with open(fn, 'w') as k_file:
    #        f = urllib.urlopen(url)
    #        k_file.write(f.read())

    print 'loading genomes..'
    g,l = parse_genome(sys.argv[1])
    print 'plotting...'
    karyoplot(g,l, {}, part=1)



