#!/usr/bin/env python

import argparse
import sys
from os.path import exists

parser = argparse.ArgumentParser()
parser.add_argument('-d', '--directory', type=str, required=True,
                    help='Directory housing TEstrainer run')
parser.add_argument('-r', '--iteration', type=str, required=True,
                    help='Interation number of TEstrainer curation')
parser.add_argument('-s', '--seq_name', type=str, required=True,
                    help='Sequence being prepared for alignment')
parser.add_argument('-g', '--genome', type=str, required=True,
                    help='Path to genome sequence')
parser.add_argument('-f', '--flank', type=int, default=1500,
                    help='Length of flank to be extended')
parser.add_argument('-n', '--no_seq', type=int, default=20,
                    help='Number of sequences to use for alignment')
parser.add_argument('-D', '--debug', action='store_true',
                    help='Set for full messaging')
args = parser.parse_args()

def file_check(file_name, debug):
  if(exists(file_name) == False):
    if(debug == False):
      exit()
    else:
      sys.exit((file_name+" not found"))

file_check(args.directory, args.debug)
file_check(args.genome, args.debug)
file_check((args.directory+'/run_'+args.iteration+'/initial_seq/'+args.seq_name), args.debug)
file_check((args.directory+'/run_'+args.iteration+'/initial_blast/'+args.seq_name+".out"), args.debug)

import string
from os import system
import statistics
import re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import pyranges as pr

def size_check(var_name, size):
  if(len(var_name) < size):
      exit()

# read in starting seq
start_seq = SeqIO.read((args.directory+"/run_"+args.iteration+"/raw/"+args.seq_name), "fasta")

# read in blast table and filter
blast_df = pd.read_table((args.directory+'/run_'+args.iteration+'/initial_blast/'+args.seq_name+'.out'), names=['qseqid', 'seqnames', 'pident', 'length', 'qstart', 'qend', 'qlen', 'start', 'end', 'slen', 'evalue', 'bitscore', 'qcovs'])
blast_df = blast_df.query('pident >= 70 & qcovs >= 50').copy()
size_check(blast_df, 3)

# determine coordinates properly direction
rev = blast_df['start'] > blast_df['end']
fwd = blast_df['start'] < blast_df['end']
blast_df.loc[rev, ['start', 'end']] = (blast_df.loc[rev, ['end', 'start']].values)
blast_df.loc[fwd, ['start']] = blast_df.loc[fwd, ['start']] - 1
blast_df.loc[rev, ['end']] = blast_df.loc[rev, ['end']] +1
 
# select 50 best, add flanks and convert to ranges
blast_df = blast_df.sort_values(by = 'bitscore', ascending=False)
blast_df = blast_df.iloc[:50]
blast_df['start'] = blast_df['start'] - args.flank
blast_df['end'] = blast_df['end'] + args.flank
blast_df.loc[blast_df['start'] < 1, 'start'] = 0
blast_df.loc[blast_df['end'] > blast_df['slen'], 'end'] = blast_df['slen']
blast_df = blast_df.sort_values(by = ['seqnames', 'start'])
blast_gr = pr.from_dict({"Chromosome": blast_df.seqnames, "Start": blast_df.start, "End": blast_df.end, "Bitscore": blast_df.bitscore, "slen" : blast_df.slen})

# reduce/merge ranges
best_hits_df = blast_gr.cluster(strand=False).df.groupby(['Cluster']).agg({'Chromosome':'first', 'Start':'min', 'End':'max', 'Bitscore':'max'})[['Chromosome','Start','End','Bitscore']].reset_index()
best_hits_py = pr.from_dict({"Chromosome": best_hits_df.Chromosome, "Start": best_hits_df.Start, "End": best_hits_df.End, "Bitscore": best_hits_df.Bitscore})


# get sequence
best_hits_seq = pr.get_fasta(best_hits_py, path=args.genome)
# convert df to strings
best_hits_df = best_hits_df.astype('str')
best_hits_df['seq'] = best_hits_seq
best_hits_df['title'] = best_hits_df.Chromosome+":"+best_hits_df.Start+"-"+best_hits_df.End

# write sequences to file
with open((args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_1'), "w") as o:
  for x in range(len(best_hits_seq)):
    out_seq_name = (best_hits_df.title[x]+"#"+best_hits_df.Bitscore[x])
    next_seq = SeqRecord(seq=Seq(best_hits_df.seq[x]),id=out_seq_name,name=out_seq_name, description="")
    SeqIO.write(next_seq, o, "fasta-2line")

# blast against og sequence,read in and filter
system("blastn -task dc-megablast -query "+args.directory+'/run_0/og/'+args.seq_name+" -subject "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1 -outfmt \"6 sseqid pident qcovs\" -out "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1.out")
check_df = pd.read_table((args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_1.out"), names=['seqnames', 'pident', 'qcovs'])
check_df = check_df.query('pident >= 70 & qcovs >= 50')
size_check(check_df, 3)

# select one instance of each sequence
check_df = check_df.groupby(['seqnames']).head(n=1).reset_index()
check_df[['seqnames', 'bitscore']] = check_df['seqnames'].str.split('#', n=1, expand = True)
size_check(check_df, 3)

# sort by bitscore and select no_seq (default 20) highest bitscores
check_df = check_df.astype({'bitscore': 'float'})
check_df = check_df.sort_values(by = 'bitscore', ascending=False)

# select accurate and write ready for self blast
correct_df = best_hits_df.loc[best_hits_df.title.isin(check_df.seqnames)].reset_index()
with open((args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_2'), "w") as o:
  for x in range(len(correct_df)):
    next_seq = SeqRecord(seq=Seq(correct_df.seq[x]),id=correct_df.title[x],name=correct_df.title[x], description="")
    SeqIO.write(next_seq, o, "fasta-2line")

# self blast
system("blastn -task dc-megablast -query "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_2 -subject "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_2 -outfmt \"6 qseqid sseqid length pident qstart qend bitscore\" -out "+args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_2.out")
self_blast_df = pd.read_table((args.directory+"/run_"+args.iteration+"/self_search/"+args.seq_name+"_check_2.out"), names=['Chromosome', 'sseqid', 'length', 'pident', 'Start', 'End', 'bitscore'])
size_check(self_blast_df, 3)

# fix coordinates
self_blast_df.Start = self_blast_df.Start - 1

# filter self hits, small hits
self_blast_df = self_blast_df[self_blast_df.Chromosome != self_blast_df.sseqid].copy()
self_blast_df = self_blast_df[self_blast_df.length >= (0.5*len(start_seq))].copy()
size_check(self_blast_df, 3)

# calculate and filter quantiles
self_blast_df['q1'] = self_blast_df.groupby('Chromosome').Start.transform(lambda x: x.quantile(0.1))
self_blast_df['q9'] = self_blast_df.groupby('Chromosome').End.transform(lambda x: x.quantile(0.9))
self_blast_df = self_blast_df[(self_blast_df.Start.astype(float) >= self_blast_df.q1) & (self_blast_df.End.astype(float) <= self_blast_df.q9)].copy()

# create and reduce/merge ranges
self_blast_gr = pr.from_dict({"Chromosome": self_blast_df.Chromosome, "Start": self_blast_df.Start, "End": self_blast_df.End})
self_hits_trimmed_df = self_blast_gr.cluster(strand=False).df.groupby(['Cluster']).agg({'Chromosome':'first', 'Start':'min', 'End':'max'})[['Chromosome','Start','End']].reset_index()
self_hits_trimmed_py = pr.from_dict({"Chromosome": self_hits_trimmed_df.Chromosome, "Start": self_hits_trimmed_df.Start, "End": self_hits_trimmed_df.End})
size_check(self_hits_trimmed_py, 3)
 
# get trimmed sequence
self_hits_trimmed_seq = pr.get_fasta(self_hits_trimmed_py, path=(args.directory+'/run_'+args.iteration+'/self_search/'+args.seq_name+'_check_2'))
self_hits_trimmed_df = self_hits_trimmed_df.astype('str')
self_hits_trimmed_df['seq'] = self_hits_trimmed_seq
self_hits_trimmed_df['name'] = (self_hits_trimmed_df.Chromosome+":"+self_hits_trimmed_df.Start+"-"+self_hits_trimmed_df.End)

# write trimmed fasta to file
with open((args.directory+'/run_'+args.iteration+'/to_align/'+args.seq_name), "w") as o:
  SeqIO.write(start_seq, o, "fasta-2line")
  for x in range(len(self_hits_trimmed_df)):
    next_seq = SeqRecord(seq=Seq(self_hits_trimmed_df.seq[x]),id=self_hits_trimmed_df.name[x],name=self_hits_trimmed_df.name[x], description="")
    SeqIO.write(next_seq, o, "fasta-2line")
