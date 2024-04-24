from multiprocessing import Process, Manager
import os
from Bio import SeqIO
import argparse
import re
import pandas as pd

def splitter(sequence, split_dir):
    file_list=[]
    if(os.path.exists(split_dir) == False):
        os.makedirs(split_dir)
    with open(sequence, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            file_name = (split_dir+"/"+record.name.split(sep="#")[0]+".fasta")
            file_list.append(record.name.split(sep="#")[0])
            SeqIO.write(record, file_name, "fasta")
    return(file_list)

def initial_trf(good_list, work_dir, seq_name):
    raw_dir=work_dir+"/split/"
    trf_dir=work_dir+"/trf/"
    trf_out_dir=work_dir+"/trf_parsed/"

    # Run trf
    trf_cmd = "trf "+raw_dir+seq_name+".fasta 2 7 7 80 10 50 500 -d -h -ngs > "+trf_dir+seq_name+".trf"
    os.system(trf_cmd)

    if os.path.getsize(trf_dir+seq_name+".trf") > 0:
        trf_df = pd.DataFrame(columns=['Chromosome', 'Start', 'End'])
        with open(trf_dir+seq_name+".trf", 'r') as trf:
            for line in trf:
                if line.startswith('@'):
                    seqnames = re.sub('@', '', line.split()[0])
                elif float(line.split()[3]) > 5:
                    trf_df = pd.concat([trf_df, pd.DataFrame({'Chromosome':[seqnames], 'Start':[int(line.split()[0])], 'End':[int(line.split()[1])]})])
        trf_df.to_csv(trf_out_dir+seq_name+".tsv", sep='\t')
    else:
        good_list.append(seq_name)

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', type=str, default="",
                        help='Directory housing TEstrainer run')
    parser.add_argument('-l', '--library', type=str, required=True,
                        help='Sequence being prepared for alignment')
    parser.add_argument('-t', '--threads', type=int, default=1,)
    args = parser.parse_args()

    split_list = splitter(args.library, args.directory+"/split")
    
    if os.path.isdir(args.directory+"/trf_parsed") is False:
        os.makedirs(args.directory+"/trf_parsed")

    with Manager() as manager:
        good_list = manager.list()
        processes = []

        for i in split_list:
            p = Process(target=initial_trf, args=(good_list,args.directory,i))  # Passing the list
            p.start()
            processes.append(p)

        for p in processes:
            p.join()

        good_list = list(good_list) 
    print(good_list)