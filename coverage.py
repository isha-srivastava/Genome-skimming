from Bio import SeqIO
import sys


def read_ref():
    
    ''' Function to read reference genome '''
    
    ref_file = 'Priacma_serrata.fa'
    fasta_seq = SeqIO.parse(open(ref_file),'fasta')
    base_pairs = 0
    ref_genome = ""
    for fasta in fasta_seq:
        name, sequence = fasta.id, str(fasta.seq)
        base_pairs += len(sequence)
        ref_genome += sequence

    # print(base_pairs, "base pairs")
    return ref_genome


def read_ind(file):
    
    ''' Function to read input fasta file for each individual '''
    
    f = open(file, 'r')
    temp = ""
    for line in f.readlines():
        temp += line
        
    genes = ['A', 'C', 'G', 'T'] 
    ind = ""
    for i in range(len(temp)):
        if temp[i] in genes:
            ind += temp[i]
            
    # print(len(ind))
    return ind


def mut_index(ind, ref_genome):
    
    ''' Function to keep list of indices where mutations have been introduced '''
    
    mut = []
    count = 0
    for i in range(len(ref_genome)):
        if ref_genome[i] != ind[i]:
            mut.append(i)
            count += 1
    
    # print(count)
    # print(len(mut))
    return mut


def read_aln(file):
    
    ''' Function to read alignment file at given coverage '''
    
    i = 0
    f = open(file, 'r')
    indices = ""
    for line in f.readlines():
        i += 1
        if (i > 4 and (i + 1) % 3 == 0):
            # print(i, line)
            indices += (line[::-1][0:10])[::-1]
    
    return indices


def positions(idx):
    
    ''' Function to return processed indices from the alignment file '''
    
    pos_ = []
    for i in range(len(idx)):
        # 1st char = tab
        if idx[i][0] == '\t':
            temp = idx[i][1:]
            idx[i] = temp
        # 2nd char = tab
        elif idx[i][1] == '\t':
            temp = idx[i][2:]
            idx[i] = temp
        # 3rd char = tab
        elif idx[i][2] == '\t':
            temp = idx[i][3:]
            idx[i] = temp
        # 4th char = tab
        elif idx[i][3] == '\t':
            temp = idx[i][4:]
            idx[i] = temp
        # 5th char = tab
        elif idx[i][4] == '\t':
            temp = idx[i][5:]
            idx[i] = temp
        # 6th char = tab
        elif idx[i][5] == '\t':
            temp = idx[i][6:]
            idx[i] = temp
    
        if(idx[i] != ''): # find the relative start pos and the strand (+/-) 
            rel_index = int(idx[i][0 : len(idx[i]) - 2])
            strand = idx[i][-1]
    
        # read_len = 70
        # total_len = 9520930
        if strand == "+":
            # print("Index: ", rel_index, rel_index + 70)
            temp = []
            temp.append(rel_index)
            temp.append(rel_index + 70)
            pos_.append(temp)
        elif strand == "-":
            # print("Index: ", (9520930 - (rel_index + 70)), (9520930 - rel_index))
            temp = []
            temp.append(9520930 - (rel_index + 70))
            temp.append(9520930 - rel_index)
            pos_.append(temp)
            
    return pos_


def write_pos(pos_, file):
    
    ''' Function to write positions to a file for later usage '''
    
    f = open(file, 'w+')
    for i in range(len(pos_)):
        f.write("%d\t%d\n" % (pos_[i][0], pos_[i][1]))
        
    print("File written")

    
def compute_cov(mut, pos_):
    
    ''' Function to compute the number of sites covered by genome skim data '''
    
    covered = 0
    cov = []
    for i in range(len(pos_)):
        beg = pos_[i][0]
        end = pos_[i][1]
        for j in range(len(mut)):
            current = mut[j]
            if(current >= beg and current < end):
                # print(current, " is in range\n")
                mut[j] = -1 # set value unusable for next iteration
                covered += 1
                cov.append(j)
    
    return covered, cov


if __name__ == "__main__":
    
    if len(sys.argv) != 4: 
        print("Missing parameters")
        sys.exit(2)
    
    ind = read_ind(sys.argv[2])
    ref_genome = read_ref()
    # print("Input files read . . ")
    mut = mut_index(ind, ref_genome)
    idx = read_aln(sys.argv[3])
    idx_split = idx.strip().split('\n')
    # print("Alignment file read . . ")
    pos = positions(idx_split)
    # print("Alignment file processed . . ")
    covered, cov = compute_cov(mut, pos)
    print("Coverage: ", sys.argv[1])
    print("Total polymorphic sites: ", len(mut))
    print("Sites covered: ", covered)
    print("Ratio: ", covered / len(mut))
  


    