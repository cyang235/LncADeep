import sys
import numpy as np

np.seterr(all='ignore')

def sigmoid(x):
    return 1.0 / (1.0 + np.exp(-x))

def d_sigmoid(x):
    return sigmoid(x) * (1 - sigmoid(x))

def GetFasta(inputfile):
    '''Get sequence from input, seqid=the first part without space'''
    try:
        f = open(inputfile, 'r')    # input file
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    tmpseq = ''
    seqlist = []
    seqID = []

    for line in f.readlines():
        line = line.strip()
        if not len(line):
            continue
        elif line[0] == '>':
            seqID.append(line.split()[0][1:])
            if tmpseq != '':
                seqlist.append(tmpseq)
            tmpseq = ''
        else:
            tmpseq += line.upper()
    seqlist.append(tmpseq)      ## append the last sequence
    f.close()

    return [seqID, seqlist]


def RevComp(seq):
    '''reverse complement of sequence'''
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'T', 'c': 'G', 'g': 'C', 't': 'A' }
    reverse_complement = "".join(complement.get(base, base) for base in reversed(seq))
    return reverse_complement


def Codon2AA2(codon):
    '''convert codon to aa'''
    if codon == "TTT" or codon == "TTC":
        return 'F'
    elif codon == 'TTA' or codon == 'TTG' or codon == 'CTT' or codon == 'CTA' or codon == 'CTC' or codon == 'CTG':
        return 'L'
    elif codon == 'ATT' or codon == 'ATC' or codon == 'ATA':
        return 'I'
    elif codon == 'ATG':
        return 'M'
    elif codon == 'GTA' or codon == 'GTC' or codon == 'GTG' or codon == 'GTT':
        return 'V'
    elif codon == 'GAT' or codon == 'GAC':
        return 'D'
    elif codon == 'GAA' or codon == 'GAG':
        return 'E'
    elif codon == 'TCA' or codon == 'TCC' or codon == 'TCG' or codon == 'TCT':
        return 'S'
    elif codon == 'CCA' or codon == 'CCC' or codon == 'CCG' or codon == 'CCT':
        return 'P'
    elif codon == 'ACA' or codon == 'ACG' or codon == 'ACT' or codon == 'ACC':
        return 'T'
    elif codon == 'GCA' or codon == 'GCC' or codon == 'GCG' or codon == 'GCT':
        return 'A'
    elif codon == 'TAT' or codon == 'TAC':
        return 'Y'
    elif codon == 'CAT' or codon == 'CAC':
        return 'H'
    elif codon == 'CAA' or codon == 'CAG':
        return 'Q'
    elif codon == 'AAT' or codon == 'AAC':
        return 'N'
    elif codon == 'AAA' or codon == 'AAG':
        return 'K'
    elif codon == 'TGT' or codon == 'TGC':
        return 'C'
    elif codon == 'TGG':
        return 'W'
    elif codon == 'CGA' or codon == 'CGC' or codon == 'CGG' or codon == 'CGT':
        return 'R'
    elif codon == 'AGT' or codon == 'AGC':
        return 'S'
    elif codon == 'AGA' or codon == 'AGG':
        return 'R'
    elif codon == 'GGA' or codon == 'GGC' or codon == 'GGG' or codon == 'GGT':
        return 'G'
    # stop codon
    elif codon == 'TAA' or codon == 'TAG' or codon == 'TGA':
        return 'J'
    else:
        return 'Z'     ## IUPAC Ambiguity Codes 


def SixMer2AA(seq):
    '''Convert 6mer to 2 AA'''

    return Codon2AA2( seq[0:3] ) + Codon2AA2( seq[3:6] )

