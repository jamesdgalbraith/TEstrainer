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

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_seq', type=str, required=True,
                    help='Input MSA to be trimmed')
args = parser.parse_args()

# set names
seq_name=re.sub('.*/', '', args.in_seq)
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
SeqIO.write(bp_trimmed, ('data/TEtrim/'+seq_name),"fasta-2line")

# create unaligned sequences
with open(('data/TEtrim_unaln/temp_'+seq_name), "w") as o:
  for record in SeqIO.parse(('data/TEtrim/'+seq_name), "fasta-2line"):
    record.seq = record.seq.ungap("-")
    SeqIO.write(record, o, "fasta-2line")

# make consensus, remove gaps, convert x to n
consensus_seq=AlignInfo.SummaryInfo(align).gap_consensus(threshold=0.25).ungap("-").replace('X', 'n')
consensus_seq=SeqRecord(consensus_seq,id=final_id,description=final_id)
SeqIO.write(consensus_seq, ('data/TEtrim_con/'+seq_name),"fasta")

### run blast ###
os.system('blastn -query data/TEtrim_con/'+seq_name+' -subject data/TEtrim_unaln/temp_'+seq_name+' -outfmt "6 qseqid sseqid qcovs" -task dc-megablast | uniq > data/TEtrim_blast/'+seq_name+'.tsv')

# read blast
df = pd.read_table(('data/TEtrim_blast/'+seq_name+'.tsv'), names=['qseqid', 'sseqid', 'qcovs'])

# determine acceptable passed on coverage >50%
acceptable = list(df.query("qcovs>50")['sseqid'])

# create unaligned fasta of acceptables
with open(('data/TEtrim_unaln/'+seq_name), "w") as o:
  for record in SeqIO.parse(('data/TEtrim_unaln/temp_'+seq_name), "fasta-2line"):
    if record.id in acceptable:
      SeqIO.write(record, o, "fasta-2line")

# remove temportary unaligned file
os.remove(('data/TEtrim_unaln/temp_'+seq_name)) 

### run mafft ###
os.system('mafft --localpair data/TEtrim_unaln/'+seq_name+' > data/TEtrim_mafft/'+seq_name)

# read new alignment, make consensus
align_2=AlignIO.read(('data/TEtrim_mafft/'+seq_name), "fasta")

# make consensus from allignment, (optimised for lines)
def line_con(aln_in):
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

con_2=SeqRecord(line_con(align_2),id=final_id,description=final_id)
SeqIO.write(con_2, ('data/TEtrim_con_2/'+seq_name),"fasta")
