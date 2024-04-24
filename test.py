from modules import beat
import os, argparse
from Bio import SeqIO
import time
import re

def variable_check_setup(work_dir, iterations, library_path, genome_path, flank, no_seq, coverage, pident, threads):
    library_name = ".".join(re.sub(".*/", "", library_path).split(".")[0:-1])

    if(work_dir == ""):
        work_dir = "output/TS_"+library_name+"_"+str(int(time.time()))
    if(os.path.exists(work_dir) == False):
        os.makedirs(work_dir)    
    return(work_dir, iterations, library_path, library_name, genome_path, flank, no_seq, coverage, pident, threads)

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

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--directory', type=str, default="",
                        help='Directory housing TEstrainer run')
    parser.add_argument('-r', '--iterations', type=int, default=10,
                        help='Interation number of TEstrainer curation')
    parser.add_argument('-l', '--library', type=str, required=True,
                        help='Sequence being prepared for alignment')
    parser.add_argument('-g', '--genome_path', type=str, required=True,
                        help='Path to genome sequence')
    parser.add_argument('-t', '--threads', type=int, default=1,
                        help='Number of sequences to use for alignment')
    parser.add_argument('-f', '--flank', type=int, default=1000,
                        help='Length of flank to be extended')
    parser.add_argument('-c', '--coverage', type=float, default=0.5,
                        help='Length of flank to be extended')
    parser.add_argument('-p', '--pident', type=int, default=70,
                        help='Length of flank to be extended')
    parser.add_argument('-n', '--no_seq', type=int, default=20,
                        help='Number of sequences to use for alignment')
    args = parser.parse_args()

    # Check necessary parsed variables and make normal variables
    work_dir, iterations, library_path, library_name, genome_path, flank, no_seq, coverage, pident, threads = variable_check_setup(args.directory, args.iterations, args.library, args.genome_path, args.flank, args.no_seq, args.coverage, args.pident, args.threads)

    # Set genome name
    genome_db = genome_path.split("/")[-1]
    
    # Make blast database
    os.makedirs(work_dir+"/beat/")
    print("Making Database")
    makedb_cmd = "makeblastdb -in "+genome_path+ " -dbtype nucl -out "+work_dir+"/beat/"+genome_db
    os.system(makedb_cmd)

    #  Make directory for completed beat sequences
    complete_dir = work_dir+"/beat/complete/"
    os.makedirs(complete_dir)

    # Split starting library
    sequence_list = splitter(library_path, work_dir+"/beat/run_1/raw/")

    for i in iterations:
        print("Starting iteration "+str(i))
        beat.prep(i, work_dir)
