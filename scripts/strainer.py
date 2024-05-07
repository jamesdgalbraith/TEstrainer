#!/bin/python3

# make necessary directories
def path_setup(out_dir):
    if(exists(out_dir+'/split') == False):
        Path(out_dir+'/split').mkdir(parents=True, exist_ok=True)

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
def rpstblastn(blast_headers, cdd_database, seq_path):
    from os import system
    rps_cmd='rpstblastn -query '+seq_path+' -db '+cdd_database+' -out '+seq_path+'.rps.out -outfmt \"6 '+blast_headers+'\" -evalue 0.01 -num_threads 1'
    system(rps_cmd)

def library_strainer(reference_path, additional_path, rps_out, in_seq_path, out_dir):
    import pandas as pd
    # Read in database of acceptable hits
    acceptable_df = pd.read_table(reference_path)
    additional_df = pd.read_table(additional_path)
    acceptable_domains_set = set(list(acceptable_df['ref']) + list(additional_df['ref']))
    # Read in rps blast results, filter small hits, split stitle
    rps_out_df = pd.read_table(rps_out)
    # Filter outs hits <50% of either database or query length
    rps_out_df = rps_out_df.query('length/slen >= 0.5 | length >= qlen/6').copy()
    rps_out_df = rps_out_df.reset_index(drop=True)
    rps_out_df[['ref', 'abbrev', 'full']] = rps_out_df['stitle'].str.split(', ', n=2, expand = True)
    
    # Split into hits in acceptable and hits not in acceptable
    contains_acceptable_df = rps_out_df[rps_out_df['ref'].isin(acceptable_domains_set)].copy()
    contains_not_acceptable_df = rps_out_df[~rps_out_df['ref'].isin(acceptable_domains_set)].copy()
    # Split out only acceptable and only not in acceptable
    only_acceptable = contains_acceptable_df[~contains_acceptable_df['qseqid'].isin(contains_not_acceptable_df['qseqid'] )].copy()
    only_not_acceptable = contains_not_acceptable_df[~contains_not_acceptable_df['qseqid'].isin(contains_acceptable_df['qseqid'] )].copy()
    # Split out and combine chimeras (contain hits in both acceptable and not in acceptable)
    chimeric_acceptable = contains_acceptable_df[contains_acceptable_df['qseqid'].isin(contains_not_acceptable_df['qseqid'] )].copy()
    chimeric_acceptable['state'] = 'acceptable'
    chimeric_not_acceptable = contains_not_acceptable_df[contains_not_acceptable_df['qseqid'].isin(contains_acceptable_df['qseqid'] )].copy()
    chimeric_not_acceptable['state'] = 'not_acceptable'
    chimeric = pd.concat([chimeric_acceptable.reset_index(drop = True), chimeric_not_acceptable.reset_index(drop = True)]).reset_index(drop = True)
    # Make sets
    only_acceptable_set = set(only_acceptable['qseqid'])
    only_not_acceptable_set = set(only_not_acceptable['qseqid'])
    chimeric_set = set(chimeric['qseqid'])
    # sort and write chimeric rps data to file
    chimeric = chimeric.sort_values(by=['qseqid', 'qstart'])
    if rps_out.count("/") > 0:
        chimeric_path=out_dir+"/chimeric_"+rps_out.split("/")[-1]
    else:
        chimeric_path=out_dir+"chimeric_"+rps_out
    chimeric.to_csv(chimeric_path, sep="\t", index=False)

    # Split library into clean and dirty
    out_seq_path=out_dir+"/"+str(in_seq_path.split("/")[-1])
    with open(in_seq_path, 'r') as handle:
        with open(out_seq_path+'.clean', 'w') as clean, \
            open(out_seq_path+'.dirty', 'w') as dirty, \
            open(out_seq_path+'.chimeric', 'w') as chimeric:
            for record in SeqIO.parse(handle, "fasta"):
                if record.id in only_not_acceptable_set:
                    SeqIO.write(record, dirty, "fasta")
                elif record.id in chimeric_set:
                    SeqIO.write(record, clean, "fasta")
                    SeqIO.write(record, chimeric, "fasta")
                else:
                    SeqIO.write(record, clean, "fasta")
    return(only_not_acceptable_set)

def strain_gff(in_gff, out_dir, unacceptable_list):
    
    # Modify names of dirty list
    dirty_list=[]
    for entry in unacceptable_list:
        dirty_list.append(sub("#.*", "", entry).lower())
    
    in_gff_name = in_gff.split('/')[-1]
    out_gff = out_dir+'/'+ '.'.join(in_gff_name.split('.')[0:-1]) + ".strained." + in_gff_name.split('.')[-1]

    # Write new copy of gff
    with open(out_gff, 'w') as out_path:
        with open(in_gff, 'r') as record:
            for line in record:
                family = sub('ID=', '', line.split('\t')[8].split(';')[2]).lower()
                if family not in dirty_list:
                    out_path.write(line)
    
if __name__ == "__main__":

    from argparse import ArgumentParser
    from pathlib import Path
    from os.path import exists
    from os import remove, rmdir
    from re import sub
    from Bio import SeqIO
    from multiprocessing import Pool
    from functools import partial
    import sys
    import tqdm

    parser = ArgumentParser()
    parser.add_argument('-i', '--in_seq', type=str, required=True,
                        help='FASTA file containing in sequences')
    parser.add_argument('-o', '--out_dir', type=str, required=True,
                        help='Path to file to write and compile rps out into')
    parser.add_argument('-b', '--blast_headers', type=str, default='qseqid qstart qend qlen slen length evalue bitscore stitle',
                        help='BLAST outfmt 6 variables to use')
    parser.add_argument('-r', '--reference', type=str, default='data/acceptable_domains_2.tsv',
                        help='Reference tsv of TE typical protein domains')
    parser.add_argument('-a', '--additional_reference', type=str, default='data/additional_domains_2.tsv',
                        help='Additional custom reference tsv of TE typical protein domains')
    parser.add_argument('-g', '--in_gff', type=str,
                        help='Path to gff to strain')
    parser.add_argument('-s', '--strain', action='store_true',
                        help='Set to run strainer')
    parser.add_argument('-t', '--num_threads', type=int, default=4,
                        help='Number of cores to use (default 4)')
    parser.add_argument('-d', '--database', type=str, default = '/ceph/software/databases/cdd/Cdd',
                        help="Full path to NCBI CDD database (default setup for Club Ashworth)")
    args = parser.parse_args()

    if(exists(args.in_seq) == False):
        sys.exit('In sequence file not found')

    # if(args.strain is True):
    #     if(args.in_gff == '' or exists(args.in_gff) == False):
    #         sys.exit('If running strainer GFF must be provided')

    print('Path check')
    path_setup(args.out_dir)

    file_list=splitter(args.in_seq, args.out_dir)
    
    print('Performing RPSTBLAST')
    with Pool(processes=args.num_threads) as pool:
        rpstblastn_func = partial(rpstblastn, args.blast_headers, args.database)
        max_ = len(file_list)
        with tqdm.tqdm(total=max_) as pbar:
            for _ in pool.imap_unordered(rpstblastn_func, file_list):
                pbar.update()

    # compile rpsblast output with header
    print('Compiling RPSTBLAST output')
    with open(args.in_seq+'.rps.out', 'w') as tsv:
        tsv.write(sub(' ', '\t', args.blast_headers)+'\n')
        for file in file_list:
            with open(file+'.rps.out', 'r') as rps:
                for line in rps:
                    tsv.write(line)
            remove(file)
            remove(file+'.rps.out')
    
    # strain library
    if(args.strain is True):
        print('Straining library')
        only_not_acceptable = library_strainer(args.reference, args.additional_reference, args.in_seq+'.rps.out', args.in_seq, args.out_dir)
        
        if args.in_gff is None:
            print('Path to GFF not be provided, skipping step.')
        elif exists(args.in_gff) is False:
            print('GFF not found, skipping step.')
        else:
            print('Straining gff')
            strain_gff(args.in_gff, args.out_dir, only_not_acceptable)

    # Remove split file and folder
    remove(args.out_dir+"/split/"+sub('.*/', '', args.in_seq)+"_split.txt")
    rmdir(args.out_dir+"/split/")

    