#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from re import sub
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--in_file', type=str, required=True,
                    help='Input file')
parser.add_argument('-o', '--out_dir', type=str,
                    help='Output directory', default = '')
args = parser.parse_args()

# check input is embl
if args.in_file[-5:] != '.embl':
  print('Input file must be in embl format with .embl extension')
  exit()

# make out folder
if args.out_dir == '':
  args.out_dir=args.in_file[0:-5]+'_fastas'
if os.path.exists(args.out_dir) == False:
  os.makedirs(args.out_dir)

# set initial variables
ID=''
new_seq=''
subclass=''
family=''
nm=''
os_embl=''
oc=''
seq_df=pd.DataFrame(data={'ID': [ID], 'Subclass': [subclass], 'Family': [family], 'Name': [nm], 'Species': [os_embl], 'Taxonomy': [oc]})

# remove compiled fasta if already existant
if os.path.exists(sub('.embl', '.fasta', args.in_file)):
  os.remove(sub('.embl', '.fasta', args.in_file)) 

# if file contains data,
with open(args.in_file, 'r') as embl:
  with open(sub('.embl', '.fasta', args.in_file), "a") as compiled:
    for lines in embl:
      if lines.split()[0] == "ID":
        # get repeat name
        # need to remove _;
        ID = sub(';', '', sub('_;', '', lines.split()[1]))
      elif lines.split()[0] == 'NM':
        nm=lines.split()[1]
      elif lines.split()[0] == "OS":
        os_embl+='_'.join(lines.split()[1:-1]) # species
      elif lines.split()[0] == "OC":
        oc+=''.join(lines.split()[1:-1])
      elif lines.split()[0] == 'CC':
        if ' Type: ' in lines:
          subclass = lines.split()[2]
        if ' SubType: ' in lines:
          family_line=lines.split()
          if len(family_line) == 3:
            family=family_line[2]
          else:
            family = ''
      elif lines[0:2] == '  ':
        new_seq+=''.join(lines.split()[0:-1])
      elif lines[0:2] == '//':
        if family == '':
          if nm == '':
            seq_id=ID+'#'+subclass
          else:
            seq_id=nm+'#'+subclass
        else:
          if nm == '':
            seq_id=ID+"#"+subclass+'/'+family
          else:
            seq_id=nm+"#"+subclass+'/'+family
        # create and write file
        out_seq=SeqRecord(seq=Seq(new_seq),id=seq_id,description=seq_id, name=seq_id)
        o=args.out_dir+'/'+ID+'.fasta'
        SeqIO.write(sequences=out_seq, handle=o, format="fasta") # write to individual file
        SeqIO.write(out_seq, compiled, "fasta") # write to compiled file
        seq_df=pd.concat([seq_df, pd.DataFrame(data={'ID': [ID], 'Subclass': [subclass], 'Family': [family], 'Name': [nm], 'Species': [os_embl], 'Taxonomy': [oc]})])
        # reset repeat sequence and classification
        ID=''
        new_seq=''
        subclass=''
        family=''
        nm=''
        oc=''
        os_embl=''


    
    
# write df to file
seq_df.to_csv(sub('.embl', '.tsv', args.in_file), sep='\t')
