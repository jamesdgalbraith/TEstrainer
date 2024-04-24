import os
from Bio import SeqIO
import pandas as pd
import pyranges as pr
import re
import shutil

def prep(iteration: int, work_dir: str):
    beat_dir = work_dir+"/beat/run_"+iteration+"/"
    subdir = [ "raw/", "trf/", "initial_blast/", "self_search/", "to_align/", "mafft/", "consensus/", "complete/" ]
    for d in subdir:
        os.makedirs(beat_dir+d)
    if iteration > 1:
        # copy previous rounds consensus needing more work into directory 
        shutil.copytree(work_dir+"/beat/run_"+str(iteration)-1+"/consensus/", work_dir+"/beat/run_"+iteration+"/raw/")

def trf(good_list, work_dir: str, seq_name: str, complete_dir: str):
    raw_dir=work_dir+"/raw/"
    trf_dir=work_dir+"/trf/"

    # Run trf
    trf_cmd = "trf "+raw_dir+seq_name+".fasta 2 7 7 80 10 50 500 -d -h -ngs > "+trf_dir+seq_name+".trf"
    os.system(trf_cmd)
    if os.path.getsize(trf_dir+seq_name+".trf") == 0:
        good_list.append(seq_name)
    else:
        # Read in starting seq
        start_seq = SeqIO.read(raw_dir+seq_name+".fasta", "fasta")
        te_len=len(start_seq.seq)

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
            # get sequence of non TR section, write to file
            non_trf_seq=SeqRecord(seq=Seq(str(start_seq.seq)[int(non_trf_longest_df['Start']):int(non_trf_longest_df['End'])]), id=start_seq.id, name=start_seq.name)
            with open(complete_dir+seq_name+".fasta", "w") as o:
                SeqIO.write(non_trf_seq, o, "fasta-2line")
        else:
            good_list.append(seq_name)

def compare(previous_seq_path: str, current_seq_path: str, flank: int):
    blast_cmd = os.system("blastn -task dc-megablast -query "+previous_seq_path+" -subject "+current_seq_path+" -evalue 1e-5 -outfmt \"6 qseqid sseqid length pident qstart qend qlen sstart send slen\" -out "+current_seq_path+".out")
    os.system(blast_cmd)
    # Read in blast results
    self_blast = pd.read_csv(current_seq_path+".out", sep="\t", header=None, names=['qseqid',"sseqid","length","pident","qstart","qend","qlen","sstart","send","slen"])
    self_blast_filtered = self_blast.loc[(self_blast.length >=500) & (self_blast.pident >= 90) & (self_blast.sstart > self_blast.send)]
    # Check how far increased
    flank_5_change = min(self_blast_filtered["qstart"])
    flank_3_change = min(self_blast_filtered["qlen"] - self_blast_filtered["qend"])

    # Alter flanks if if necessary
    if flank_5_change<=flank/2:
        flank_5=0
    else:
        flank_5=flank
    if flank_3_change<=flank/2:
        flank_3=0
    else:
        flank_3=flank
    return flank_5, flank_3

def blast_extend(seq_name, flank_5, flank_3, genome_db, work_dir, complete_dir):
    if flank_5 == 0 and flank_3 == 0:
        shutil.copyfile(work_dir+"raw/"+seq_name+".fasta", complete_dir+seq_name+".fasta")
    else:
        blast_cmd = os.system("blastn -task dc-megablast -query "+work_dir+"raw/"+seq_name+".fasta -db "+genome_db+" -evalue 1e-5 -outfmt \"6 qseqid sseqid pident length qstart qend qlen sstart send slen evalue bitscore qcovs\" -out "+beat_dir+"initial_blast/"+seq_name+".out")
        os.system(blast_cmd)
        # Filter table
        # Rearrange Blast table
        # Extend flanks
        # Get sequences
        # Write consensus and sequences to file

def mafft(work_dir, sequence_list, iteration, threads):
    beat_dir = work_dir+"/beat/run_"+iteration
    os.makedirs(beat_dir+"/mafft")
    for seq_name in sequence_list:
        mafft_cmd = "mafft --thread "+threads+" --quiet --localpair --adjustdirectionaccurately "+beat_dir+"/to_align/"+seq_name+" > "+beat_dir+"/mafft/"+seq_name
        os.system(mafft_cmd)
    return()

def trim():
    return()

def further_check():
    more_needed = True
    if more_needed is True:
        return(True)
    else:
        return(False)
    
def compile(work_dir):
    # function to compile completed sequences
    os.makedirs(work_dir+"/complete_beat")
    seq_list = []
    return(seq_list)