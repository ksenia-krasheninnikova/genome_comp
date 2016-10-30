#!/usr/bin/env python

import os
import sys
import matplotlib
import argparse
#matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from model import Entry
from color_cache import cache

def parse_genome(genome_path):
    chrs = []
    #lengths = []
    names = set()
    with open(genome_path) as f:
        for line in f:
            line = line.strip().split()
            if len(line) == 1:
                chrs.append([])
                #lengths.append(0)
            else:
                e= Entry(line[1],line[5],int(line[7]),int(line[9]),int(line[9])-int(line[7]))
                e.set_block_id(int(line[3]))
                chrs[-1].append(e)
                #lengths[-1] += e.length
                names.add(e.seq_id.split('.')[1])
    return chrs

def check_if_homolog_to_reference_chrom(c, reference_block_id2colors_sizes):
    homologous_blocks_len = 0
    c_hom = []
    for e in c:
        if e.block_id in reference_block_id2colors_sizes.keys():
            homologous_blocks_len += e.length
            c_hom.append((e,True))
        else:
            c_hom.append((e,False))
    return float(homologous_blocks_len) / sum(map(lambda x: x.length, c)) >= 0.5, c_hom

#genomes is the list of chrs (chrs is a list of chromosomes of a genome)
def karyoplot(ref_genome, genomes):
    fig = plt.figure(1,figsize=(8, 6))
    ax = fig.add_subplot(111, aspect='equal')
    default_color = 'black'
    breath = 0.05
#    block_length_param = 0.0000001
    #block_length_param = 7.5e-09
    block_length_param = 5e-08
    x_dist = 0.05
    y_dist = 0.05
    x_start = x_dist
    for c in ref_genome:
        reference_block_id2colors_sizes = {}
        cache_names = cache.keys()
        color_i = 0
        y_start = y_dist
        for e in c:
            if color_i == len(cache_names):
                i=0
                #raise Exception('No more colors!')
            y_size = float(e.length) * block_length_param 
            r = Rectangle((x_start, y_start), breath, y_size, facecolor = cache[cache_names[color_i]])
            #remember the colors in reference
            reference_block_id2colors_sizes[e.block_id] = (cache[cache_names[color_i]],y_size)
            y_start += y_size
            #print r
            ax.add_patch(r)
            color_i += 1
        center_x = x_start + breath/2.0
        radius = breath/2.0
        theta1 = 0.0
        theta2 = 180.0
        w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
        w2 = Wedge((center_x, y_dist), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
        ax.add_patch(w1)
        ax.add_patch(w2)
        
        for g in genomes:
            x_start += breath + x_dist
            y_start = y_dist
            for c in g:
                if_hom,c_hom = check_if_homolog_to_reference_chrom(c, reference_block_id2colors_sizes)
                if if_hom:
                    #w2 = Wedge((center_x, y_start), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
                    #ax.add_patch(w2)
                    for e,status in c_hom:
                        if status: #True if homologous
                            r_color, r_y_size = reference_block_id2colors_sizes[e.block_id]
                            r = Rectangle((x_start, y_start), breath, r_y_size, facecolor = r_color)
                            y_start += r_y_size
                            ax.add_patch(r)
                        else:
                            print 'else'
                            y_size = float(e.length) * block_length_param 
                            r = Rectangle((x_start, y_start), breath, y_size, facecolor = default_color)
                            y_start += r_y_size
                            ax.add_patch(r)
                    #radius = breath/2.0
                    #center_x = x_start + breath/2.0
                    #w1 = Wedge((center_x, y_start), radius, theta1, theta2, width=0.00001, facecolor='white', edgecolor='black')
                    #y_start += 2*radius
                    #ax.add_patch(w1)
                
                    y_start += radius
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
    parser = argparse.ArgumentParser()
    parser.add_argument('reference')
    parser.add_argument('others', nargs='+')
    args = parser.parse_args()
    print 'loading ref genome..'
    rg = parse_genome(args.reference)
    print 'loading other genomes..'
    gs = []
    for g in args.others:
       x =  parse_genome(g)
       gs.append(x)
    print 'plotting...'
    karyoplot(rg,gs)



