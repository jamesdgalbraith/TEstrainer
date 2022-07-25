#!/usr/bin/env python

import os
import sys
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
parser.add_argument('-n', '--iteration', type=str, required=True,
                    help='iteration number')
parser.add_argument('-t', '--threads', type=str,
                    help='Threads to use', default=1)
parser.add_argument('-f', '--flank', type=int,
                    help='Length of flanks used', default=1500)
parser.add_argument('-d', '--debug', type=str,
                    help='Print debug messages', default='FALSE')
parser.add_argument('-w', '--window', type=int,
                    help='Window size for blocks', default=5)

args = parser.parse_args()
# args = argparse.Namespace(in_seq="data/run_1/mafft/Echis_carinatus_rnd-1_family-117.fasta", iteration="1", threads="1", flanks="1500", debug='FALSE', window_size=5)

# function for removing single base pair insertions
def single_trim(aln_in):
  # make empty alignment
  good=aln_in[:,0:0]
  for x in range(aln_in.get_alignment_length()):
    # extract columns with more than 1 base pair
    if len(aln_in) - aln_in[:, x].count("-") > 1:
      good=good+aln_in[:,x:x+1]
  return(good)

# function for making first consensus from alignment
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

# set names
seq_name=re.sub('.*/', '', args.in_seq)
in_seq_path='data/run_'+args.iteration+'/mafft/'+seq_name
out_seq_path='data/run_'+args.iteration+'/TEtrim/'+seq_name


# read sequence
if args.debug == 'TRUE':
  print('Reading og alignment')
align = AlignIO.read(in_seq_path, "fasta")

# get name of first sequence (the consensus)
final_id = align[0].id
print(final_id)
og_con = SeqRecord(seq= align[0].seq.ungap("-"), id=align[0].id, name=align[0].id)
SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_con/og_'+seq_name),"fasta")
align = align[1:len(align)]

# cancel if less than 3 sequences, write to file for troubleshooting
if args.debug == 'TRUE':
  print('Reading determining number of sequences')
if(len(align)<3):
  SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit((seq_name+" contains less than 3 sequences"))
print(seq_name+" contains at least 3 sequences")

# single bp trim
bp_trimmed=single_trim(align)
SeqIO.write(bp_trimmed, ('data/run_'+args.iteration+'/TEtrim_bp/trimmed_'+seq_name),"fasta")

# create unaligned sequences
if args.debug == 'TRUE':
  print('Creating unaligned sequences')
with open(('data/run_'+args.iteration+'/TEtrim_unaln/temp_'+seq_name), "w") as o:
  for record in SeqIO.parse(('data/run_'+args.iteration+'/TEtrim_bp/trimmed_'+seq_name), "fasta"):
    record.seq = record.seq.ungap("-")
    SeqIO.write(record, o, "fasta-2line")

### run blast ###
if args.debug == 'TRUE':
  print('Reading determining number of sequences')
os.system('blastn -query data/run_'+args.iteration+'/TEtrim_con/og_'+seq_name+' -subject data/run_'+args.iteration+'/TEtrim_unaln/temp_'+seq_name+' -outfmt "6 qseqid sseqid qcovs" -task blastn | uniq > data/run_'+args.iteration+'/TEtrim_blast/'+seq_name+'.tsv')

# read initial blast, determine acceptable passed on coverage >50% mean of coverage
if args.debug == 'TRUE':
  print('Reading initial blast')
df = pd.read_table(('data/run_'+args.iteration+'/TEtrim_blast/'+seq_name+'.tsv'), names=['qseqid', 'sseqid', 'qcovs'])
mean_covs=(statistics.mean(df.qcovs)/2)
acceptable=list(df.query("(qcovs>@mean_covs) or (qcovs>50)")['sseqid'])

# Check at least 3
if len(acceptable) < 3:
  SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit((seq_name+" contains less than 3 acceptable sequences"))

# create unaligned fasta of acceptables
if args.debug == 'TRUE':
  print('Creating list of unaligned acceptable sequences')
with open(('data/run_'+args.iteration+'/TEtrim_unaln/unaln_'+seq_name), "w") as o:
  for record in SeqIO.parse(('data/run_'+args.iteration+'/TEtrim_unaln/temp_'+seq_name), "fasta-2line"):
    if record.id in acceptable:
      SeqIO.write(record, o, "fasta-2line")

### run mafft ###
if args.debug == 'TRUE':
  print('Running initial mafft')
os.system('mafft --quiet --thread '+args.threads+' --localpair data/run_'+args.iteration+'/TEtrim_unaln/unaln_'+seq_name+' > data/run_'+args.iteration+'/TEtrim_mafft/mafft_'+seq_name)

# read new alignment, make consensus
if args.debug == 'TRUE':
  print('Reading new alignment and making consensus')
align_2=AlignIO.read(('data/run_'+args.iteration+'/TEtrim_mafft/mafft_'+seq_name), "fasta")
con_2=SeqRecord(con_maker(align_2),id=final_id,description=final_id)
SeqIO.write(con_2, ('data/run_'+args.iteration+'/TEtrim_con/cleaned_'+seq_name),"fasta")

# Run final blast
if args.debug == 'TRUE':
  print('Final blast')
os.system('blastn -query data/run_'+args.iteration+'/TEtrim_con/og_'+seq_name+' -subject data/run_'+args.iteration+'/TEtrim_con/cleaned_'+seq_name+' -outfmt "6 length qlen slen qcovs" -task dc-megablast -out data/run_'+args.iteration+'/TEtrim_blast/check_'+seq_name+'.tsv')

if os.path.getsize('data/run_'+args.iteration+'/TEtrim_blast/check_'+seq_name+'.tsv') == 0:
  SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit(('New consensus of '+seq_name+" has no homology to original sequence"))

# read in first line of blast
df2 = pd.read_table(('data/run_'+args.iteration+'/TEtrim_blast/check_'+seq_name+'.tsv'), names=['length', 'qlen', 'slen', 'qcovs']).head(1)

# determine values from blast table
length=df2.values[0,0]
qlen=df2.values[0,1]
slen=df2.values[0,2]
qcovs=df2.values[0,3]

if qcovs <= 80:
  SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit(('New consensus of '+seq_name+" covers < 80% than original sequence"))
elif (slen <= (0.9 * qlen)):
  SeqIO.write(og_con, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit(('New consensus of '+seq_name+" is shorter than original sequence"))
elif ((slen - qlen) <= (0.5 * args.flank)):
  SeqIO.write(con_2, ('data/run_'+args.iteration+'/TEtrim_complete/'+seq_name),"fasta")
  sys.exit(('New consensus of '+seq_name+" is good, no further work needed"))
else:
  SeqIO.write(con_2, ('data/run_'+args.iteration+'/TEtrim_further/'+seq_name),"fasta")
  sys.exit(('New consensus of '+seq_name+" is good, needs further extension."))
