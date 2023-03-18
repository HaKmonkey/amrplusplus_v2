#!/usr/bin/env python3

from argparse import ArgumentParser
from os import listdir, getcwd

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--indir', required=True)
    parser.add_argument('--manifest', required=True)
    return parser.parse_args()

def direction(read):
    if '_R1.fastq.gz' in read: return 0
    else: return 1

def get_read_pairs(indir):
    read_pairs = {}
    for f in listdir(indir):
        sample_id = f.split('_R1.fastq.gz')[0].split('_R2.fastq.gz')[0]
        read_pairs.setdefault(sample_id, ['',''])
        read_pairs[sample_id][direction(f)] = f'{getcwd()}/{indir}/{f}'
    return read_pairs

def write_manifest(manifest, read_pairs):
    with open(manifest, 'w') as o:
        o.write('sample-id,absolute-filepath,direction\n')
        [o.write(f'{k},{read_pairs[k][0]},forward\n{k},{read_pairs[k][1]},reverse\n') for k in read_pairs]

def main():
    args = parse_args()
    read_pairs = get_read_pairs(args.indir)
    write_manifest(args.manifest, read_pairs)

if __name__ == '__main__':
    main()