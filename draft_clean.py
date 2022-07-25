#!/usr/bin/env python

import os
import string
import statistics
import argparse
import re
from math import ceil
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
from plotnine import ggplot, aes, geom_line, scale_x_continuous

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed')
args = parser.parse_args()

# set names
seq_name=re.sub('.*/', '', 'data/TEtrim_mafft/Echis_carinatus_rnd-1_family-2.fasta')
in_seq_path="data/TEtrim_mafft/"+seq_name
out_path='data/TEtrim/'
window_size=3

# read sequence
align = AlignIO.read(in_seq_path, "fasta")

# remove single base pair insertions
def single_trim(aln_in):
  # make empty alignment
  good=aln_in[:,0:0]
  for x in range(aln_in.get_alignment_length()):
    # extract columns with more than 1 base pair
    if len(aln_in) - aln_in[:, x].count("-") > 1:
      good=good+aln_in[:,x:x+1]
  return(good)

# find blocks
def find_block(aln_in):
  
  gaps=[]
  for x in range(aln_in.get_alignment_length()):
    gaps.append(aln_in[:, x].count("-")/len(aln_in))
  
  # step along in windows of 5, calculating average number of gaps
  gap_analysis=[]
  for x in range(ceil(len(gaps)/window_size)):
    gap_analysis.append(statistics.mean(gaps[(x*window_size):(((x+1)*window_size))]))
  
  # find starting good block
  for x in range(len(gap_analysis)-1):
    if (gap_analysis[x] < 0.5) & (gap_analysis[x+1] < 0.5):
      aln_st=x*window_size
      break
  
  # find ending good block
  for x in reversed(range(len(gap_analysis)-1)):
    if (gap_analysis[x] < 0.5) & (gap_analysis[x-1] < 0.5):
      aln_end=(x+1)*window_size
      break
  
  return(aln_in[:,aln_st:aln_end])

inner_trimmed=end_trimmed[:,x:x]
for x if range(len(end_trimmed)):
  if end_trimmed[:, x].count("-")/len(end_trimmed) < 0.5:
    inner_trimmed.extend(end_trimmed[:,x])
  
x=0


end_trimmed[:,0]
end_trimmed[0:1].id
# make consensus
def con_maker(aln_in):
  con_seq=str()
  for x in range(aln_in.get_alignment_length()):
    a = aln_in[:, x].count("a") + aln_in[:, x].count("A")
    t = aln_in[:, x].count("t") + aln_in[:, x].count("T")
    c = aln_in[:, x].count("c") + aln_in[:, x].count("C")
    g = aln_in[:, x].count("g") + aln_in[:, x].count("G")
    gap = aln_in[:, x].count("-")
    if (len(aln_in) - gap) > 2:
      pos_count = {'A': a, 'T': t, 'C': c, 'G': g}
      if len([k for k, v in pos_count.items() if v == max(pos_count.values())])>1:
        con_seq+='n'
      else:
        con_seq+=max(pos_count, key=pos_count.get)
  con_seq=Seq(con_seq)
  return(con_seq)

final_id="Echis_carinatus_rnd-1_family-2#DNA/PIF-Harbinger"
bp_trimmed=single_trim(align)
end_trimmed=find_block(bp_trimmed)
consensus=SeqRecord(line_con(end_trimmed),id=final_id,description=final_id)
SeqIO.write(end_trimmed, ('data/TEtrim_con_2/ends_trimmed_'+seq_name),"fasta")
SeqIO.write(consensus, ('data/TEtrim_con_2/consensus_'+seq_name),"fasta")


