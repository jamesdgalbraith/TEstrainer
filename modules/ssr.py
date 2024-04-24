import os
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pyranges as pr
import re

def initial_trf(work_dir, seq_name):
    raw_dir=work_dir+"raw/"
    trf_dir=work_dir+"trf/"
    complete_dir=work_dir+"complete/"
    
    # Run trf
    trf_cmd = "trf "+raw_dir+seq_name+".fasta 2 7 7 80 10 50 500 -d -h -ngs > "+trf_dir+seq_name+".trf"
    os.system(trf_cmd)

    # Read in starting seq
    start_seq = SeqIO.read(raw_dir+seq_name+".fasta", "fasta")
    te_len=len(start_seq.seq)

    # Check if sequence is largely tandem repeat
    # Parse trf
    trf_df = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])
    with open(trf_dir+seq_name+".trf", 'r') as trf:
        for line in trf:
            if line.startswith('@'):
                seqnames = re.sub('@', '', line.split()[0])
            elif float(line.split()[3]) > 5:
                trf_df = pd.concat([trf_df, pd.DataFrame({'Chromosome':[seqnames], 'Start':[int(line.split()[0])], 'End':[int(line.split()[1])]})])

    # Merge trf pr
    trf_pr=pr.PyRanges(df=trf_df).merge()

    # Calculate total length of satellite sequences
    # TEs >90% consider satellites
    if trf_pr.length/te_len > 0.9:
        with open(complete_dir+seq_name+".fasta", "w") as o:
            SeqIO.write(start_seq, o, "fasta")
        return()
    # TEs >50% consider TEs needing trimming
    elif (trf_pr.length/te_len > 0.5) & (len(start_seq.seq) > 500):
        # make pr object of whole repeat
        te_pr=pr.from_dict({'Chromosome': [seqnames], 'Start': [0], 'End': [te_len]})
        # get inverse or trf, i.e. non tandem repeat section
        trf_inverse_pr=te_pr.subtract(trf_pr)
        # convert to df
        trf_inverse_df=trf_inverse_pr.as_df()
        # calculate width
        trf_inverse_df['Width']=trf_inverse_df.End-trf_inverse_df.Start+1
        # select longest non TR section and convert to ranges
        non_trf_longest_df=trf_inverse_df[trf_inverse_df.Width==max(trf_inverse_df.Width)].head(n=1)
        non_trf_longest_pr=pr.PyRanges(df=non_trf_longest_df)
        # get sequence of non TR section, write to file
        non_trf_seq=SeqRecord(seq=Seq(str(start_seq.seq)[int(non_trf_longest_df['Start']):int(non_trf_longest_df['End'])]), id=start_seq.id, name=start_seq.name)
        with open(complete_dir+seq_name+".fasta", "w") as o:
            SeqIO.write(non_trf_seq, o, "fasta-2line")
        return()
    return(seq_name)

def trf():
    return()

def mreps():
    return()

def sassr():
    return()

def ssr_check():
    return()