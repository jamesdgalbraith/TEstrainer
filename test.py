#!/usr/bin/env python

import os
import string
import statistics
import argparse
import re
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--in_seq', type=str, required=True,
#                     help='Input MSA to be trimmed')
# args = parser.parse_args()

# set names
seq_name='Echis_carinatus_rnd-1_family-0.fasta'
in_seq_path="data/mafft/"+seq_name
out_seq_path='data/TEtrim/'+seq_name

# read sequence
align = AlignIO.read(in_seq_path, "fasta")

# get name of first sequence (the consensus)
final_id = align[0].id
og_con = SeqRecord(seq= align[0].seq.ungap("-"), id=align[0].id, name=align[0].id)
align = align[1:len(align)]

# remove single base pair insertions
def single_trim(aln_in):
  # make empty alignment
  good=aln_in[:,0:0]
  for x in range(aln_in.get_alignment_length()):
    # extract columns with more than 1 base pair
    if len(aln_in) - aln_in[:, x].count("-") > 1:
      good=good+aln_in[:,x:x+1]
  return(good)

bp_trimmed=single_trim(align)
SeqIO.write(bp_trimmed, ('data/TEtrim/bp_trimmed_'+seq_name),"fasta-2line")

a = (bp_trimmed[:, x].count("a") + bp_trimmed[:, x].count("A"))/len(bp_trimmed)
t = (bp_trimmed[:, x].count("t") + bp_trimmed[:, x].count("T"))/len(bp_trimmed)
c = (bp_trimmed[:, x].count("c") + bp_trimmed[:, x].count("C"))/len(bp_trimmed)
g = (bp_trimmed[:, x].count("g") + bp_trimmed[:, x].count("G"))/len(bp_trimmed)
gap = bp_trimmed[:, x].count("-")
statistics.mean()
