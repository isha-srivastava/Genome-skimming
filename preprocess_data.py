import Bio
import sys
from Bio import SeqIO

def preprocess_fasta(inputfile):
    infile = SeqIO.parse(open(inputfile), 'fasta')
    bp_1 = 0 # base pairs
    seq = ""

    for fasta in infile:
        name, sequence = fasta.id, str(fasta.seq)
        bp_1 += len(sequence)
        seq += sequence

    print("Base pairs of simulated genome: ", bp_1)
    return seq

def write_fasta(outfile, ind):
    f = open(outfile, "w+")
    f.write('>seq1\n') # Only 1 sequence
    for i in range(0, len(ind)):
        #if i != 0 and i % 70 == 0: # limit set at 70
        #    f.write("\n")
        f.write(ind[i])
    #f.write("\n")
    f.close()
    
if __name__ == "__main__":
    
    if len(sys.argv) != 3: 
        print("Missing parameters")
        sys.exit(2)
    
    ind_10000 = preprocess_fasta(sys.argv[1])
    write_fasta(sys.argv[2], ind_10000)
    
    print("File written\n")