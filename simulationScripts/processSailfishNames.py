from __future__ import print_function
import argparse
import os
import sys
import re

def main(args):
    ifile = args.input
    if (not os.path.isfile(ifile)):
        print('Could not find input file {}'.format(ifile), file=sys.stderr)
        sys.exit(1)

    lines = open(ifile).readlines()
    headers = [l for l in lines if l.startswith('#')]
    data = [l for l in lines if not l.startswith('#')]

    with open(args.output, 'w') as ofile:
        for l in headers:
            ofile.write(l)
        for l in data:
            toks = l.rstrip('\n').split('\t')
            name = re.split('\s+', toks[0])[0]
            newline = '\t'.join([name] + toks[1:])
            ofile.write("{}\n".format(newline))

    print("successfully modified file", file=sys.stderr)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process crazy transcript names')
    parser.add_argument('--input', type=str, required=True, help='the quant.sf file')
    parser.add_argument('--output', type=str, required=True, help='where the processed file should be placed')
    args = parser.parse_args()
    main(args)
