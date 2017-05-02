import argparse
import os
import sys

parser = argparse.ArgumentParser(prog=os.path.basename(sys.argv[1]), description='Process EDM file with pythia.', prefix_chars="+")

def rootfile(filename):
    if not filename.endswith(".root"):
        raise ValueError("{} should end with .root".format(filename))
    return filename

parser.add_argument("infile", type=rootfile)
parser.add_argument("firstevent", type=int)
parser.add_argument("lastevent", type=int)
parser.add_argument("outfile", type=rootfile)

args = parser.parse_args(sys.argv[2:])

infile = args.infile
firstevent = args.firstevent
lastevent = args.lastevent
outfile = args.outfile

if not os.path.exists(infile):
    raise ValueError("{} should exist".format(infile))
if os.path.exists(outfile):
    raise ValueError("{} should not exist".format(outfile))
