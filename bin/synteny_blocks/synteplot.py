#!/usr/bin/env python

import os
import sys
import matplotlib
import argparse
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle
from matplotlib.collections import PatchCollection
import matplotlib.pyplot as plt
from model import Entry
from color_cache import cache

length_threshold = 1000000

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

def check_if_homolog_to_reference_chrom(c, ref_blocks):
    homologous_blocks_len = 0
    c_hom = []
    for e in c:
        if e.block_id in ref_blocks and e.length >= length_threshold:
            homologous_blocks_len += e.length
            c_hom.append((e,True))
        else:
            c_hom.append((e,False))
    return float(homologous_blocks_len) / sum(map(lambda x: x.length, c)) >= 0.5, c_hom

class ReferenceBlockDrawing:
    def __init__(self,e, color, y_size, y_start):
        self.e = e
        self.color = color
        self.y_size = y_size
        self.y_start = y_start

#genomes is the list of chrs (chrs is a list of chromosomes of a genome)
def karyoplot(ref_genome, genomes,names, folder):
    default_color = 'black'
    breath = 0.05
    block_length_param = 5e-08
    x_dist = 0.05
    y_dist = 0.05
    plt.hold(False) 
    for c in ref_genome:
        chr_name = c[0].seq_id.split('.')[1]
        print 'processing',chr_name
        fig = plt.figure(1,figsize=(5, 7))
        ax = fig.add_subplot(111, aspect='equal')
        ax.set_ylim([0,1.5])
        x_start = x_dist
        block_id2drawing = {}
        cache_names = cache.keys()
        color_i = 0
        y_start = y_dist
        for e in c:
            if e.length < length_threshold:
                continue
            if color_i == len(cache_names):
                #i = 0
                raise Exception('No more colors!')
            y_size = float(e.length) * block_length_param 
            r = Rectangle((x_start, y_start), breath, y_size, facecolor = cache[cache_names[color_i]])
            print r
            #remember the colors in reference
            block_id2drawing[e.block_id] = ReferenceBlockDrawing(e, cache[cache_names[color_i]], y_size, y_start)
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
        
        
        for g,name in zip(genomes,names):
            x_start += breath + x_dist
            y_start = y_dist
            rectangles = []
            for c in g:
                if_hom,c_hom = check_if_homolog_to_reference_chrom(c, block_id2drawing.keys())
                if if_hom:
                    for e,status in c_hom:
                        if status: #True if homologous
                            ref_drawing = block_id2drawing[e.block_id]
                            r = Rectangle((x_start, y_start), breath, ref_drawing.y_size, facecolor = ref_drawing.color)
                            rectangles.append(r)
                            y_start += ref_drawing.y_size
                            ax.add_patch(r)
                        else:
                            if e.length < length_threshold:
                                continue
                            #print 'else'
                            y_size = float(e.length) * block_length_param 
                            r = Rectangle((x_start, y_start), breath, y_size, facecolor = default_color)
                            rectangles.append(r)
                            y_start += y_size
                            ax.add_patch(r)
                
                    y_start += radius
            name = ' '.join(name.split('_'))
            ax.text(x_start+breath/2.0, y_start, name, va='bottom', rotation='vertical')

        x_start += breath + x_dist                     
        ax.set_axis_off()
        ax.set_title(chr_name)
        fig.savefig(os.path.join(folder,chr_name+'.png'), dpi=90)
        fig.clf()
        plt.cla()
        plt.close()
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
    parser.add_argument('path',help='folder for the output')
    parser.add_argument('reference')
    parser.add_argument('others', nargs='+')
    args = parser.parse_args()
    print 'loading ref genome..'
    rg = parse_genome(args.reference)
    print 'loading other genomes..'
    gs = []
    names = []
    for g in args.others:
       names.append(g)
       x =  parse_genome(g)
       gs.append(x)
    print 'plotting...'
    karyoplot(rg,gs,names,args.path)



