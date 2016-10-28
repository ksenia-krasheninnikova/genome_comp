#!/usr/bin/env python

"""
Environment for handling synteny blocks relations
"""
class BlocksToPathsProcessor:

    @staticmethod
    def __find_graph_edges(entries):
        i = 0
        j = 1
        edges = []
        while i < len(entries) - 1:
            v = entries[i]
            v_next = entries[i+1]
            for e in v:
                appended = False
                for e_next in v_next:
                    if e.seq_id == e_next.seq_id:
                        edges.append((e,e_next,i))
                        #print i
                        #e.print_out()
                        #e_next.print_out()
                        appended = True
                #in case we didn't found the next block on the same chromosome
                #add all the following blocks (on all the chromosomes)
                if not appended:
                    for e_next in v_next:
                        edges.append((e,e_next,i))
            i += 1
        ##in case we didn't filter out unsplitted chromosomes
        ##we add a loop
        if len(entries) == 1:
            edges.append((entries[0][0],entries[0][0],0))
        return edges

    @staticmethod
    def __dfs(v, edges, level, max_level):
        new_paths = []
        for e in edges:
            #no check that the next vertex equals the end of the previous edge
            if e[2] == level:
                if level == max_level:
                    new_paths.append([(e[0],e[1])])
                else:
                    for l in BlocksToPathsProcessor.__dfs(e[1], edges, level+1, max_level):
                        '''
                        print 'l:'
                        for x in l:
                            x[0].print_out()
                            x[1].print_out()
                            print '-----'
                        '''
                        new_paths.append([(e[0],e[1])] + l)
        #there can be multiple paths in case of at some point the neighbouring blocks
        #have multiple alternative chromosomes like this:
        #a b d
        #  c
        #here a,b,c,d - chromosomes, and alternative paths are a b d and a c d
        #also, b and c can belong to different locations on the same chromosome
        '''
        for x in new_paths:
            for y in x:
                if y[0].block_id==993 or y[1].block_id==993:
                    print 'new paths:'
                    for a in new_paths:
                        for b in a:
                            b[0].print_out()
                            b[1].print_out()
                            print
                        print

                    return new_paths
        '''
        return new_paths

    '''
    in case a block contains one entry in specie1 and two entries in specie2
    it will result in two alternative paths (actually a bubbles) that we will try to merge
    by sorting according the positions of these two entries along the chromosome in specie2
    '''
    @staticmethod
    def try_merge(path) :
        entries_paths = []
        for e in path:
            entries_paths.append([e[0][1]])
            entries_paths[-1] += map(lambda x: x[1], e[1:])
        zipped_path = zip(*entries_paths)
        merged_path = []
        #for e in zipped_path:
        #    print e
        for e in zipped_path :
            #in case all are equal
            if e[0].equals_to_list(e[1:]):
                merged_path.append(e[0])
            else:
                #if all are on the same chromosome
                if len(set(map(lambda x: x.seq_id, e))) == 1:
                    a = sorted(e, key=lambda x: x.start)
                    merged_path += a
                else:
                    print 'cant merge!'
                    e[0].print_out()
                    for z in e[1:]:
                        z.print_out()
                    print '----'
                    return path
        return [merged_path]


    #entries - a chromosome consists in a list of blocks
    #each block is also a list of possible alternative locations of this block
    @staticmethod
    def search_paths(entries):
        '''
        #PRINT ENTRIES
        print 'entries:'
        for e in entries:
            for p in e:
                p.print_out()
            print '---'
        '''
        edges = BlocksToPathsProcessor.__find_graph_edges(entries)
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
            max_level = len(entries)-2
            if len(entries) == 1:
                max_level = 0
            path = BlocksToPathsProcessor.__dfs(e, edges, 0, max_level)
            if len(path) > 1:
                path = BlocksToPathsProcessor.try_merge(path)
                if len(path) == 1:
                    return path[0]
                else:
                    print 'Alternative solutions!'
                    '''
                    print len(path)
                    i = 0
                    for p in path:
                        i += 1
                        print i
                        for x in p:
                            x[0].print_out()
                            x[1].print_out()
                            print '----'
                    '''
                    print 'returning empty chromosome'
                    return []
            thread = [path[0][0][0]]
            for p in path[0]:
                thread.append(p[1])
            threads.append(thread)
        if len(threads) > 1:
            print 'Alternative solutions!'
            print 'returning empty chromosome'
            return []
        return threads[0]
