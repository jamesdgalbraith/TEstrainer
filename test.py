import string
import statistics
import os
from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import SeqIO, AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import pairwise2
# from plotnine import ggplot, aes, geom_line
import pandas as pd
import numpy as np

os.getenv('PATH')

# set names
seq_name="Echis_carinatus_rnd-1_family-117.fasta"
in_seq_path="data/mafft/"+seq_name
out_seq_path='data/TEtrim/'+seq_name
# read sequence
align = AlignIO.read(in_seq_path, "fasta")

# get name of first sequence (the consensus)
final_id = align[0].id
og_con = SeqRecord(seq= align[0].seq.ungap("-"), id=align[0].id, name=align[0].id)
align = align[1:len(align)]

# remove single base pair insertions
def single_trim(align):
  # make empty alignment
  good=align[:,0:0]
  for x in range(align.get_alignment_length()):
    # extract columns with more than 1 base pair
    if (len(align) - align[:, x].count("-"))/len(align) > 1/len(align):
      good=good+align[:,x:x+1]
  return(good)
bp_trimmed=single_trim(align)
SeqIO.write(bp_trimmed, ('data/TEtrim/'+seq_name),"fasta-2line")

# create unaligned sequences
with open(('data/TEtrim_unaln/temp_'+seq_name), "w") as o:
  for record in SeqIO.parse(('data/TEtrim/'+seq_name), "fasta-2line"):
    record.seq = record.seq.ungap("-")
    SeqIO.write(record, o, "fasta-2line")

# make consensus, remove gaps, convert x to n
consensus_seq=AlignInfo.SummaryInfo(align).gap_consensus(threshold=0.25)
consensus_seq=consensus_seq.ungap("-")
consensus_seq=consensus_seq.replace('X', 'n')
consensus_seq=SeqRecord(consensus_seq,id=final_id,description=final_id)
SeqIO.write(consensus_seq, ('data/TEtrim_con/'+seq_name),"fasta")

### run blast ###
# blastn -query data/TEtrim_con/Echis_carinatus_rnd-1_family-117.fasta -subject data/TEtrim_unaln/temp_Echis_carinatus_rnd-1_family-117.fasta -outfmt "6 qseqid sseqid qcovs" -task dc-megablast | uniq > data/TEtrim_blast/Echis_carinatus_rnd-1_family-117.tsv

os.system('blastn')

# read blast
df = pd.read_table("data/TEtrim_blast/Echis_carinatus_rnd-1_family-117.tsv", names=['qseqid', 'sseqid', 'qcovs'])

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






# calculate distances
calculator=DistanceCalculator('identity')
distmat=calculator.get_distance(bp_trimmed)




# make df of dismat
for x in range(len(dismat)):
  if x==0:
    disdf=pd.DataFrame(dismat[x])
    disdf.index=dismat.names
  else:
    holder=pd.DataFrame(dismat[x])
    holder.index=dismat.names
    disdf=pd.concat([disdf,holder], axis=1)
disdf.columns = list(range(len(align)))




for disLis in dismat:
  print(dismat.names)
  statistics.mean(disLis)
      
dismatDF = pd.DataFrame(distances)

# remove consensus

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

window_width=7

def window_counter(x):
  window=align[:, x:x+window_width]
  for i in range(window_width):
    a = (window[:, i].count("a") + window[:, i].count("A"))/len(window)
    t = (window[:, i].count("t") + window[:, i].count("T"))/len(window)
    c = (window[:, i].count("c") + window[:, i].count("G"))/len(window)
    g = (window[:, i].count("g") + window[:, i].count("G"))/len(window)
    gap = (window[:, i].count("-"))/len(window)
    if i==0:
      window_stats=np.array([a,t,c,g,gap])
    else:
      window_stats=np.vstack((window_stats, np.array([a,t,c,g,gap])))
  bp_average=[statistics.mean(window_stats[:,0]),
  statistics.mean(window_stats[:,1]),
  statistics.mean(window_stats[:,2]),
  statistics.mean(window_stats[:,3])]
  coverage_average=statistics.mean(window_stats[:,4])
  return coverage_average

for x in range(int(align.get_alignment_length()/window_width)-1):
  window_counter(x)

for x in range(align.get_alignment_length()):
  print(x)
  if x==0:
    pos=[]
  if (align[:, x].count("-"))/len(align)>=0.15:
    pos.append(x)

columns_as_strings=
for y in pos:
  if y == pos[0]:
    columns_as_strings=Seq(Seq(align[:, y]))
  else:
    columns_as_strings.append(Seq(align[:, y]))
align[:, (0):1]
for y in pos:
  outputMSA = tempMSA[:,:] + align[:, (y-1):y]
  tempMSA = outputMSA

SeqRecord(Seq==outputMSA,record=[rec.id for rec in align])



id(align)

Align.MultipleSeqAlignment()

[Seq(align[:, y]),Seq(align[:, y])]

MultipleSeqAlignment(columns_as_strings)




windows=(align.get_alignment_length()-6)

for x in range(windows):
  print((statistics.mean(coverage[x:(x+6)]) > 4) & (statistics.mean(iden[x:(x+6)]) > 0.5))
  if((statistics.mean(coverage[x:(x+6)]) > 4) & (statistics.mean(iden[x:(x+6)]) > 0.5)):
    sweet.append(x)

sweet






positions=pos
inputMSA=align
msa_array = np.array([list(rec) for rec in inputMSA], dtype=str)

align_contig_ids = [rec.id for rec in align]
new_numpy_array = msa_array[:,positions]

seq_record_list = []
counter = 0
for x in new_numpy_array:
    dna_str = "".join(x)
    record_temp = SeqRecord(Seq(dna_str),id=align_contig_ids[counter],description="")
    seq_record_list.append(record_temp)
    counter = counter + 1

filtered_msa = MultipleSeqAlignment(seq_record_list)
print(filtered_msa)
