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
import numpy as np
import pyranges

directory = 'TS_Mixia_osmundae-families.fa_0846'
iteration = "2"
seq_name = 'Mixia_osmundae-families.fa'
seq_no = "1"
flank = 1500
genome_path = "seq/Mixia_osmundae_v1.0.fasta"

# read in blast table and filter
blast_df = pd.read_table((directory+'/run_'+iteration+'/initial_blast/'+seq_name+'_seq_'+seq_no+'.fasta.out'), names=['qseqid', 'seqnames', 'pident', 'length', 'qstart', 'qend', 'qlen', 'start', 'end', 'slen', 'evalue', 'bitscore', 'qcovs'])
blast_df = blast_df.query('pident >= 70 & qcovs >= 50')
# determine coordinates properly direction
blast_df = blast_df.assign(strand="+")
m = blast_df['start'] > blast_df['end']
blast_df.loc[m, ['start', 'end']] = (blast_df.loc[m, ['end', 'start']].values)
blast_df.loc[m, ['strand']] = "-"
# select 50 best, add flanks and convert to ranges
blast_df = blast_df.sort_values(by = 'bitscore', ascending=False)
blast_df = blast_df.iloc[:50]
blast_df.start = blast_df.start - flank
blast_df.end = blast_df.end + flank
blast_df.loc[blast_df['start'] < 1, 'start'] = 1
blast_df.loc[blast_df['end'] > blast_df['slen'], 'end'] = blast_df['slen']
blast_df = blast_df.sort_values(by = ['seqnames', 'start'])
pd.DataFrame({'seqnames' : blast_df.seqnames, "start" : blast_df.start, "end" : blast_df.end, "strand" : blast_df.strand, "slen" : blast_df.slen, "bitscore" : blast_df.bitscore})
gr = pyranges.from_dict({"Chromosome": blast_df.seqnames, "Start": blast_df.start, "End": blast_df.end, "Strand": blast_df.strand, "Bitscore": blast_df.bitscore, "slen" : blast_df.slen})
# reduce/merge ranges
best_hits_df = gr.cluster(strand='TRUE').df.groupby(['Cluster']).agg({'Chromosome':'first', 'Start':'min', 'End':'max', 'Strand':'first', 'Bitscore':'max'})[['Chromosome','Start','End','Strand','Bitscore']].reset_index()
best_hits_py = pyranges.from_dict({"Chromosome": best_hits_df.Chromosome, "Start": best_hits_df.Start, "End": best_hits_df.End, "Strand": best_hits_df.Strand, "Bitscore": best_hits_df.Bitscore})
# get sequence
best_hits_seq = pyranges.get_fasta(best_hits_py, path=genome_path)
# convert df to strings
best_hits_df = best_hits_df.astype('str')

# write sequences to file
with open((directory+'/run_'+iteration+'/self_search/'+seq_name+'_seq_'+seq_no+'_check1.fasta'), "w") as o:
  for x in range(len(best_hits_seq)):
    seqname = (best_hits_df.Chromosome[x]+":"+best_hits_df.Start[x]+"-"+best_hits_df.End[x]+"("+best_hits_df.Strand[x]+")_"+best_hits_df.Bitscore[x])
    next_seq = SeqRecord(seq=Seq(best_hits_seq[x]),id=seqname,name=seqname, description="")
    SeqIO.write(next_seq, o, "fasta-2line")

"blastn -task dc-megablast -query "+directory+'/run_0/'+seq_name+" -subject "+directory+"/run_"+iteration+"/self_search/"seq_name"_check1.fasta -outfmt \"6 qseqid sseqid qcovs\""
