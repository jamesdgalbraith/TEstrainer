import string
import statistics
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
# from plotnine import ggplot, aes, geom_line
import pandas as pd
# import numpy as np

# set names
seq_name="Echis_carinatus_rnd-1_family-117.fasta"
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
  print(x)
  position_x = len(align) - align[:, x].count("-")
  coverage.append(position_x)
  a = (align[:, x].count("a") + align[:, x].count("A"))/len(align)
  t = (align[:, x].count("t") + align[:, x].count("T"))/len(align)
  c = (align[:, x].count("c") + align[:, x].count("G"))/len(align)
  g = (align[:, x].count("g") + align[:, x].count("G"))/len(align)
  gap = (align[:, x].count("-"))/len(align)
  iden.append(max(a, t, c, g))

windows=(align.get_alignment_length()-6)
sweet=[]
x=0
for x in range(windows):
  print((statistics.mean(coverage[x:(x+6)]) > 4) & (statistics.mean(iden[x:(x+6)]) > 0.5))
  if((statistics.mean(coverage[x:(x+6)]) > 4) & (statistics.mean(iden[x:(x+6)]) > 0.5)):
    sweet.append(x)

sweet

calculator = DistanceCalculator('identity', )
dismat = calculator.get_distance(align)
new_align = []
for x in range(20):
  if statistics.mean(dismat[x]) < 0.2:
    new_align.append(align[x,:])


