#!/usr/bin/env python

import argparse


def evaluate_weighted_coverage(header, vals):
    cum = 0.0
    prev_seq_len = 0.0
    for cov, seq_len in zip(header[::-1], vals[::-1]):
        cum += cov*(seq_len - prev_seq_len)
        prev_seq_len = seq_len
    return cum / vals[0]

def report_coverage_as_duplicated_rate(vals):
    print 'aligned unique:', str((vals[0] - vals[1])/vals[0]*100)+'%'
    print 'aligned with coverage at least 2X:', str(vals[1]/vals[0]*100)+'%'

def handle_cov(cov_file):
    with open(cov_file) as f:
        header = f.readline()
        header = map(lambda x: x.strip(), header.split(','))
        header = header[1:]
        header = map(lambda x: x.replace('sitesMapping',''), header)
        header = map(lambda x: x.replace('Times',''), header)
        header = map(int, header)
        result_map = {}
        for line in f:
            line = map(lambda x: x.strip(), line.split(','))
            genome = line[0]
            vals = map(float, line[1:])
            eval = evaluate_weighted_coverage(header, vals)
            print 'weighted coverage of', genome+':', eval
            report_coverage_as_duplicated_rate(vals)
            result_map[genome] = eval
        return result_map


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('file', help='coverage file from halStats')

    args = parser.parse_args()
    handle_cov(args.file)