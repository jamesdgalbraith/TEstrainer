#!/bin/python3

# make necessary directories
def path_setup(out_dir):
    if(exists(out_dir+'/split') == False):
        mkdir(out_dir+'/split')

# split fasta file
def splitter(in_seq, out_dir):
    file_list=[]
    # split fasta
    with open(in_seq, 'r') as handle:
        for record in SeqIO.parse(handle, "fasta"):
            file_name = (out_dir+"/split/"+record.name.split(sep="#")[0]+".fasta")
            file_list.append(record.name.split(sep="#")[0]+".fasta")
            SeqIO.write(record, file_name, "fasta-2line")
    # write file list
    with open((out_dir+"/split/"+sub('.*/', '', in_seq)+"_split.txt"), 'w') as fp:
        for item in file_list:
            # write each item on a new line
            fp.write("%s\n" % item)
    # append full path to file list
    append_str=out_dir+'/split/'
    full_file_list = [append_str + sub for sub in file_list]
    return(full_file_list)

# run rpstblastn on single fasta file with custom blast headers
def rpstblastn(blast_headers, seq_path):
    print(seq_path)
    rps_cmd='rpstblastn -query '+seq_path+' -db /Users/jgalbrai/Databases/Cdd/Cdd -out '+seq_path+'.rps.out -outfmt \"6 '+blast_headers+'\" -evalue 0.01 -num_threads 1'
    system(rps_cmd)

if __name__ == "__main__":

    from argparse import ArgumentParser
    from os import mkdir, system
    from os.path import exists
    from re import sub
    from Bio import SeqIO
    from multiprocessing import Pool
    from functools import partial
    import sys

    parser = ArgumentParser()
    parser.add_argument('-i', '--in_seq', type=str, required=True,
                        help='FASTA file containing in sequences')
    parser.add_argument('-o', '--out_path', type=str, required=True,
                        help='Path to file to write and compile rps out into')
    parser.add_argument('-b', '--blast_headers', type=str, default='qseqid qstart qend qlen slen length evalue bitscore stitle',
                        help='BLAST outfmt 6 variables to use')
    parser.add_argument('-t', '--num_threads', type=int, default=4,
                        help='Number of cores to use (default 4)')
    args = parser.parse_args()

    if(exists(args.in_seq) == False):
        sys.exit('In sequence file not found')

    path_setup(args.out_dir)

    file_list=splitter(args.in_seq, args.out_dir)
    
    with Pool(processes=args.num_threads) as pool:
        rpstblastn_func = partial(rpstblastn, args.blast_headers)
        pool.map(rpstblastn_func, file_list)

    # compile rpsblast output with header
    with open(args.out_path, 'w') as tsv:
        tsv.write(sub(' ', '\t', args.blast_headers)+'\n')
        for file in file_list:
            with open(file+'.out', 'r') as rps:
                for line in rps:
                    tsv.write(line)