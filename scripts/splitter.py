#!/usr/bin/env python

import os
from os.path import exists
import sys
import argparse
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input multi-fasta to be split')
parser.add_argument('-d', '--out_dir', type=str, required=True,
                    help='Output directory')                    

if(exists(in_seq) == False):
  sys.exit('File not found)
if(exists(out_dir) == False):
  os.mkdir(out_dir)
file_list=[]

# split fasta file
with open(args.in_seq, 'r') as handle:
    for record in SeqIO.parse(handle, "fasta"):
        file_name = (args.out_dir+"/"+record.name.split(sep="#")[0]+".fasta")
        file_list.append(record.name.split(sep="#")[0]+".fasta")
        SeqIO.write(record, file_name, "fasta-2line")
# write file list
with open((args.out_dir+"/split_file_list.txt"), 'w') as fp:
    for item in file_list:
        # write each item on a new line
        fp.write("%s\n" % item)

