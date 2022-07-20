import argparse
parser = argparse.ArgumentParser(description='Trim MSA and create consensus sequence')
parser.add_argument('-i', '--in_seq', help='in_sequence')
args = parser.parse_args()

import string
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
# from plotnine import ggplot, aes, geom_line
# import pandas as pd
# import numpy as np

# set names
seq_name=args.in_seq
# seq_name="data/mafft/Echis_carinatus_rnd-1_family-117.fasta"
in_seq_path="data/mafft/"+seq_name
out_seq_path='data/TEtrim/'+seq_name
# read sequence
align = AlignIO.read(in_seq_path, "fasta")
# get name of first sequence (the consensus)
final_id = align[0].id
# remove consensus
align = align[1:len(align)]
# make containers
coverage = []
acceptable = []
iden = []
con_string = ""
# loop to make consensus
for x in range(align.get_alignment_length()):
  # ensure coverage of > 2 sequences
  position_x = len(align) - align[:, x].count("-")
  coverage.append(position_x)
  if position_x > 2:
    acceptable.append(x)
    # count ratio of each base
    a = (align[:, x].count("a") + align[:, x].count("A"))/position_x
    t = (align[:, x].count("t") + align[:, x].count("T"))/position_x
    c = (align[:, x].count("c") + align[:, x].count("G"))/position_x
    g = (align[:, x].count("g") + align[:, x].count("G"))/position_x
    iden.append(max(a, t, c, g))
    # make consensus using most common base
    if(max(a, t, c, g) >= 0.5):
      pos_count = {'A': a, 'T': t, 'C': c, 'G': g}
      con_string+=max(pos_count, key=pos_count.get)
    else:
      con_string+="N"
# convert consensus to sequence

coverage

con_seq=Seq(con_string)
# name and describe sequence
con_seq_record=SeqRecord(con_seq, id=final_id, description="")
# write sequence to file
SeqIO.write(con_seq_record, out_seq_path, "fasta")

