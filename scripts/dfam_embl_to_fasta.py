#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from re import sub
import pandas as pd
# import argparse
import os

# parser = argparse.ArgumentParser()
# parser.add_argument('-i', '--in_file', type=str, required=True,
#                     help='Input file')
# args = parser.parse_args()

in_file = "seq/Dfam/Dfam_curatedonly.embl"
out_dir = "seq/Dfam/curatedonly_fastas/"

if os.path.exists(out_dir) == False:
  os.makedirs(out_dir)

# set initial variables
ID=''
new_seq=''
subclass=''
family=''
nm=''
os=''
oc=''

# if file contains data,
with open(in_file, 'r') as embl:
  for lines in embl:
    if lines.split()[0] == "ID":
      # get repeat name
      # need to remove _;
      ID = sub(';', '', sub('_;', '', lines.split()[1]))
    elif lines.split()[0] == 'NM':
      nm=lines.split()[1]
    elif lines.split()[0] == "OS":
      os+='_'.join(lines.split()[1:-1]) # species
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
      o=out_dir+'/'+ID+'.fasta'
      SeqIO.write(sequences=out_seq, handle=o, format="fasta")
      # reset repeat sequence and classification
      ID=''
      new_seq=''
      subclass=''
      family=''
      nm=''
      oc=''
      os=''

pd.DataFrame(data={'ID': [ID], 'Subclass': [subclass], 'Family': [family], 'Name': [nm], 'Species': [os], 'Taxonomy': [oc]})
