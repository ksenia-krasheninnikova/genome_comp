#!/usr/bin/env python

import argparse

def find_all(string, pattern):
    return [i for i in range(len(string)) if string.startswith(pattern, i)]

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='breakpoints table')
    parser.add_argument('--filter_num', help='filter rows containg number of breakpoints')
    args = parser.parse_args()

    with open(args.file) as f:
        header = f.readline()
        print header.strip()
        species = header.strip().split()[1:]
        for line in f:
            line = line.strip()
            if args.filter_num:
                n = int(args.filter_num)
                oc = find_all(line, 'BR')
                if len(oc) >= n:
                    print line
                    

