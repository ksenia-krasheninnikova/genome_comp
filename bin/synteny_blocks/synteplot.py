#!/usr/bin/env python

import os
import sys
import matplotlib
import argparse
matplotlib.use('Agg')
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.patches import Circle, Wedge, Polygon, Rectangle, FancyArrowPatch
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
    #return float(homologous_blocks_len) / sum(map(lambda x: x.length, c)) >= 0.5, c_hom
    return homologous_blocks_len > 0, c_hom

class ReferenceBlockDrawing:
    def __init__(self,e, color, y_size, y_start):
        self.e = e
        self.color = color
        self.y_size = y_size
        self.y_start = y_start

def draw_main(ax, c, chr_name, x_start,  breath = 0.05, block_length_param = 5e-08, x_dist = 0.05, y_dist = 0.05, link_height = 0.05) :
        if chr_name == 'chrC1':
            ax.set_ylim([0,3])
        else:
            ax.set_ylim([0,2])
        block_id2drawing = {}
        cache_names = cache.keys()
        color_i = 0
        y_start = y_dist
        #link means belonging to the same scaffold
        if_draw_link = False
        for e in c:
            if e.length < length_threshold:
                #if_draw_link = False
                continue
            if color_i == len(cache_names):
                #i = 0
                raise Exception('No more colors!')
            if if_draw_link:
                link_loc  = x_start + breath/2.0
                link = FancyArrowPatch((link_loc,y_start),(link_loc,y_start+link_height),arrowstyle='-')
                ax.add_patch(link)
            y_start+=link_height
            if_draw_link = True
            y_size = float(e.length) * block_length_param 
            r = Rectangle((x_start, y_start), breath, y_size, facecolor = cache[cache_names[color_i]])
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
        w2 = Wedge((center_x, y_dist+link_height), radius, theta2, theta1, width=0.00001, facecolor='white', edgecolor='black')
        ax.add_patch(w1)
        ax.add_patch(w2)
        last_y = y_start + 0.3
        return ax, block_id2drawing

def draw_homologous(ax, block_id2drawing, g, x_start, name,unhomologous_len, homologous_len, unhomologous_num, homologous_num,\
                    breath = 0.05, x_dist = 0.05, y_dist = 0.05, link_height = 0.05): 
            x_start += breath + x_dist
            y_start = y_dist
            rectangles = []
            radius = breath/2.0
            for c in g:
                if_hom,c_hom = check_if_homolog_to_reference_chrom(c, block_id2drawing.keys())
                if if_hom:
                    #for e,status in c_hom:
                    for i in range(len(c_hom)):
                        e = c_hom[i][0]
                        homologous_len += e.end - e.start
                        homologous_num += 1
                        status = c_hom[i][1]
                        if status: #True if homologous
                            ref_drawing = block_id2drawing[e.block_id]
                            #r = Rectangle((x_start, y_start), breath, ref_drawing.y_size, facecolor = ref_drawing.color)
                            r = Rectangle((x_start, ref_drawing.y_start), breath, ref_drawing.y_size, facecolor = ref_drawing.color)
                            cur_y = ref_drawing.y_start + ref_drawing.y_size 
                            if cur_y > y_start:
                                y_start = cur_y
                            #if i > 0 and c_hom[i-1][1]:
                            prev_hom = filter(lambda x: x[1] == True, c_hom[:i])
                            if i > 0 and prev_hom:
                                prev_hom = prev_hom[-1]
                                link_loc  = x_start + breath/2.0
                                #prev_block_drawing = block_id2drawing[c_hom[i-1][0].block_id] 
                                prev_block_drawing = block_id2drawing[prev_hom[0].block_id] 
                                if prev_block_drawing.y_start < ref_drawing.y_start:
                                    #link down
                                    link = FancyArrowPatch((link_loc,ref_drawing.y_start),(link_loc,ref_drawing.y_start-link_height),arrowstyle='-')
                                    ax.add_patch(link)
                                else: #link up
                                    link = FancyArrowPatch((link_loc,prev_block_drawing.y_start),(link_loc,prev_block_drawing.y_start-link_height),arrowstyle='-')
                                    ax.add_patch(link)
                            rectangles.append(r)
                            ax.add_patch(r)
                    else:
                        e = c_hom[i][0]
                        unhomologous_len = e.end - e.start
                        unhomologous_num += 1
            name = ' '.join(name.split('_'))
            ax.text(x_start+breath/2.0, y_start+2*radius, name, va='bottom', rotation='vertical')
            return ax, x_start, unhomologous_len, homologous_len, unhomologous_num, homologous_num


#genomes is the list of chrs (chrs is a list of chromosomes of a genome)
def karyoplot(ref_genome, genomes,names, folder):
    default_color = 'black'
    breath = 0.05
    block_length_param = 5e-08
    x_dist = 0.05
    y_dist = 0.05
    plt.hold(False) 
    for c in ref_genome:
        fig = plt.figure(1,figsize=(5, 7))
        ax = fig.add_subplot(111, aspect='equal')
        x_start = x_dist
        chr_name = c[0].seq_id.split('.')[1]
        print 'processing',chr_name
        ax, block_id2drawing = draw_main(ax, c, chr_name, x_start) 
        unhomologous_len = 0 
        homologous_len = 0 
        unhomologous_num = 0
        homologous_num = 0
        for g,name in zip(genomes,names):
            ax, x_start, unhomologous_len, homologous_len, unhomologous_num, homologous_num = \
                draw_homologous(ax, block_id2drawing, g, x_start, name, unhomologous_len, homologous_len, unhomologous_num, homologous_num)
        print 'length of homologous', homologous_len
        print 'number of homologous', homologous_num
        print 'length of unhomologous', unhomologous_len, float(unhomologous_len)/(homologous_len+0.01), 'of homologous'
        print 'number of unhomologous', unhomologous_num
        x_start += breath + x_dist                     
        ax.set_axis_off()
        ax.set_title(chr_name)
        fig.savefig(os.path.join(folder,chr_name+'.png'), dpi=90)
        fig.clf()
        plt.cla()
        plt.close()


if __name__ == '__main__':
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



