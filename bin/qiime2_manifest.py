#!/usr/bin/env python3

from argparse import ArgumentParser
from os import listdir, path

def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--indir', required=True)
    parser.add_argument('--manifest', required=True)
    return parser.parse_args()

def direction(read):
    if '_R1_001.fastq.gz' in read: return 0
    else: return 1

def get_read_pairs(indir):
    read_pairs = {}
    for f in listdir(indir):
        sample_id = f.split('_R1_001.fastq.gz')[0].split('_R2_001.fastq.gz')[0]
        read_pairs.setdefault(sample_id, ['',''])
        read_pairs[sample_id][direction(f)] = path.abspath(f)
    return read_pairs

def write_manifest(manifest, read_pairs):
    with open(manifest, 'w') as o:
        o.write('sample-id\tforward-absolute-filepath\treverse-absolute-filepath\n')
        [o.write(f'{k}\t{read_pairs[k][0]}\t{read_pairs[k][1]}\n') for k in read_pairs]

def main():
    args = parse_args()
    read_pairs = get_read_pairs(args.indir)
    write_manifest(args.manifest, read_pairs)

if __name__ == '__main__':
    main()