#!/usr/bin/python3
from argparse import ArgumentParser
from glob import glob
import os
def main(cmdline=None):
    parser = ArgumentParser()
    parser.add_argument('-g', '--genome-dir', required=True)
    parser.add_argument('-l', '--library-id')
    parser.add_argument('-o', '--output-dir', default='.')
    args = parser.parse_args(cmdline)

    make_comment_file(args.library_id, args.genome_dir, args.output_dir)

def make_comment_file(lib_id, genome_dir, output_dir):
    pathname = os.path.join(output_dir, 'COfile.txt')
    with open(pathname, 'wt') as outstream:
        if lib_id is not None:
            outstream.write('@CO\tLIBID:{}\n'.format(lib_id))
        for comment_file in glob(os.path.join(genome_dir, '*_bamCommentLines.txt')):
            with open(comment_file, 'rt') as instream:
                for line in instream:
                    outstream.write(line)


if __name__ == '__main__':
    main()
