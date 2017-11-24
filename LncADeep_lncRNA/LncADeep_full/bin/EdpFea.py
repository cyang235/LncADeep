import os,sys
import numpy as np

from utils import GetFasta
from utils import Codon2AA2 

from HmmFea import RunHMM
from HmmFea import ReadHmm
from HmmFea import GenerateTrans

from Hexamer import ReadLogScore 
from Hexamer import HexamerScore2

## AA list 
_AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

## 2mer list
_DNA = ['A', 'C', 'G', 'T']
_Kmer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        _Kmer_list.append(dna1+dna2)

## 3mer list
_3mer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            _3mer_list.append(dna1+dna2+dna3)

## IUPAC code
_IUPAC = {'A':'A','C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'M':'AC', 'K':'GT', 'S':'CG', 'W':'AT', 'H':'ACT', \
          'B':'CGT', 'V':'ACG', 'D':'AGT', 'N':'ACGT'}

def IUPAC_2mer(seq):
    '''Return a list of all possible 2mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            kmer_list.append(dna1+dna2)
    return kmer_list

def IUPAC_3mer(seq):
    '''Return a list of all possible 3mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            for dna3 in _IUPAC[seq[2]]:
                if Codon2AA2(dna1+dna2+dna3) != "J":
                    kmer_list.append(dna1+dna2+dna3)
    return kmer_list

def SixMer2AA(seq):
    '''Convert 6mer to 2 AA'''

    return Codon2AA2( seq[0:3] ) + Codon2AA2( seq[3:6] )



def GetFeature(seq, seqid, HMMdict1, logscore_dict, ORFFea=True, EdpFea=True, 
               FicFea=True, ORFkmerFea=True, UTRkmerFea=True, HmmFea=True, HexFea=True):

    # HMMER index
    if HmmFea:
        if HMMdict1.has_key(seqid):
            hmmfeature = HMMdict1[seqid][1] + "\t" + HMMdict1[seqid][2] + "\t" + HMMdict1[seqid][3]
        else:
            hmmfeature = "0\t0\t0"
    else:
        hmmfeature = "" 

    # search ORF
    seq = seq.strip()
    ORF, UTR5, UTR3 = GetORF_UTR(seq)
    transcript_len = len(seq)

    # ORF features
    if ORFFea:
        ORF_fea = str( len(ORF) ) + "\t" + str( len(ORF) * 1.0 / transcript_len )
    else:
        ORF_fea = ""

    # The EDP of ORF
    if EdpFea :
        if len(ORF) < 6:
            EDP_fea = GetEDP_noORF()
        else:
            EDP_fea = GetEDP(ORF, transcript_len)
    else:
        EDP_fea = ""

    if ORFkmerFea :
        Kmer_EDP_fea = GetKmerEDP(ORF)
    else:
        Kmer_EDP_fea = ""

    # UTR feature
    if UTRkmerFea :
        UTR5_fea = str(len(UTR5) * 1.0 / len(seq)) + "\t" + str(GetGC_Content(UTR5))
        UTR3_fea = str(len(UTR3) * 1.0 / len(seq)) + "\t" + str(GetGC_Content(UTR3))
    else:
        UTR5_fea = ""
        UTR3_fea = ""

    # Hexamer feature
    ave_hexamer_score = HexamerScore2(ORF, logscore_dict)
    if HexFea :
        hex_score = str(ave_hexamer_score)
    else:
        hex_score = ""

    # Fickett feature
    if FicFea :
        A_pos_fea = GetBasePositionScore(seq, 'A')
        C_pos_fea = GetBasePositionScore(seq, 'C')
        G_pos_fea = GetBasePositionScore(seq, 'G')
        T_pos_fea = GetBasePositionScore(seq, 'T')
        base_ratio = GetBaseRatio(seq)

        fickett_fea = "\t".join([str(A_pos_fea), str(C_pos_fea), str(G_pos_fea), str(T_pos_fea), str(base_ratio[0]), str(base_ratio[1]), str(base_ratio[2]), str(base_ratio[3])])
    else:
        fickett_fea = ""

    feature = "\t".join([EDP_fea, ORF_fea, hmmfeature, Kmer_EDP_fea, UTR5_fea, UTR3_fea, hex_score, fickett_fea])

    return feature


def GetBasePositionScore(seq, base):
    '''compute base(A, C, G, T) position score based
       on Fickett's methods
    '''
    base_num = [0, 0, 0]
    tmp = len(seq) / 3

    for i in range(0, tmp):
        for j in range(0, 3):
            if seq[j + 3*i] == base:
                base_num[j] += 1

    base_pos_score = max(base_num) * 1.0 / (min(base_num) + 1)

    return base_pos_score


def GetBaseRatio(seq):
    '''calculate the A, C, G, T percentage of ORF'''
    A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9 
    for i in range(0, len(seq)):
        if seq[i] == 'A':
            A += 1
        elif seq[i] == 'C':
            C += 1
        elif seq[i] == 'G':
            G += 1
        elif seq[i] == 'T':
            T += 1

    A = A * 1.0 / (len(seq) + 4e-9)
    C = C * 1.0 / (len(seq) + 4e-9)
    G = G * 1.0 / (len(seq) + 4e-9)
    T = T * 1.0 / (len(seq) + 4e-9)

    return [A, C, G, T]


def GetGC_Content(seq):
    '''calculate GC content of sequence'''
    A, C, G, T = 1e-9, 1e-9, 1e-9, 1e-9 
    for i in range(0, len(seq)):
        if seq[i] == 'A':
            A += 1
        elif seq[i] == 'C':
            C += 1
        elif seq[i] == 'G':
            G += 1
        elif seq[i] == 'T':
            T += 1

    GC = (G + C) / (A + C + G + T)

    return GC 


def GetORF_UTR(seq):
    '''Get ORF and UTR from sequence'''
    STP = {}
    STP[0] = []
    STP[1] = []
    STP[2] = []
    STP[0].append(0)
    STP[1].append(1)
    STP[2].append(2)
    AAnum = len(seq) / 3

    for i in range(0, 3):
        for j in range(0, AAnum):
            tmp = seq[(i+3*j):(i+3*(j+1))]
            if tmp == 'TAG' or tmp == 'TAA' or tmp == 'TGA':
                STP[i].append(i+3*j)

    ORF = {}

    for i in range(0,3):
        if len(STP[i]) < 2:
            continue
        for j in range(1, len(STP[i])):
            tmpN = (STP[i][j] - STP[i][j-1])/3
            for k in range(0, tmpN):
                tmpS = seq[ (STP[i][j-1] + 3*k):(STP[i][j-1] + 3*(k+1)) ]
                if tmpS == 'ATG':
                    ORF[3*k + STP[i][j-1]] = STP[i][j]
                    break

    # longest ORF
    ORFseq = []
    ORFlen = []
    ORFstart = []
    ORFend = []
    for (k,v) in ORF.items():
        ORFseq.append(seq[k:(v-3)])
        ORFlen.append(v - k)
        ORFstart.append(k)
        ORFend.append(v)

    # ORF doesn't exist
    if ORF:
        ORF_l = ORFseq[np.argmax(ORFlen)]
        UTR5 = ''
        UTR3 = ''
        if len(seq[ 0:ORFstart[np.argmax(ORFlen)] ]) > 0:
            UTR5 = seq[ 0:ORFstart[np.argmax(ORFlen)] ]
        if len(seq[ORFend[np.argmax(ORFlen)]:]) > 0:
            UTR3 = seq[ORFend[np.argmax(ORFlen)]:]
        return [ORF_l, UTR5, UTR3]
    else:
        return ['', '', '']



def GetEDP_noORF():
    '''return the default values'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9 

    sum_codon = 1e-9 * 20 

    H = 0.0
    for (k,v) in Codon.items():
        Codon[k] /= sum_codon
        Codon[k] = -Codon[k] * np.log2(Codon[k])
        H += Codon[k]

    EDP = {}
    for (k,v) in Codon.items():
        EDP[k] = Codon[k] / H

    outline = ''
    for i in range(20):
        outline += str(0) + "\t"

    return outline 


def GetEDP(seq, transcript_len):
    '''get features including: ORF length, ORF ratio, ORF EDP of codon'''

    Codon = {}
    for aa in _AA_list:
        Codon[aa] = 1e-9 

    sum_codon = 1e-9 * 20 

    if(len(seq) > 3):
        num = len(seq) / 3
        for i in range(0,num) :
            if Codon2AA2( seq[i*3:(i+1)*3] ) == "J":
                continue
            ## IUPAC codon
            elif Codon2AA2( seq[i*3:(i+1)*3] ) == "Z":
                tmp_kmer_list = IUPAC_3mer(seq[i*3:(i+1)*3])
                for tmp_kmer in tmp_kmer_list:
                    Codon[ Codon2AA2(tmp_kmer) ] += 1.0 / len(tmp_kmer_list)
                sum_codon += 1.0
            else:
                Codon[ Codon2AA2( seq[i*3:(i+1)*3] ) ] += 1.0 
                sum_codon += 1.0

        H = 0.0
        for (k,v) in Codon.items():
            Codon[k] /= sum_codon
            Codon[k] = -Codon[k] * np.log2(Codon[k])
            H += Codon[k]

        EDP = {}
        for (k,v) in Codon.items():
            EDP[k] = Codon[k] / H

        outline = ''
        for (k,v) in EDP.items():
            outline += str(v) + "\t"
        
        return outline 


def GetKmerEDP(seq):

    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16 

    if(len(seq) > 3):
        for i in range(0,len(seq)-1) :
            ## IUPAC kmer
            if seq[ i:(i+2) ] not in _Kmer_list:
                tmp_kmer_list = IUPAC_2mer(seq[i:(i+2)])
                for tmp_kmer in tmp_kmer_list:
                    Kmer[ tmp_kmer ] += 1.0 / len(tmp_kmer_list)
            else:
                Kmer[ seq[ i:(i+2)] ] += 1.0 
            sum_Kmer += 1.0

        H = 0.0
        for (k,v) in Kmer.items():
            Kmer[k] /= sum_Kmer
            Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
            H += Kmer[k]

        EDP = {}
        for (k,v) in Kmer.items():
            EDP[k] = Kmer[k] / H

        outline = ''
        for (k,v) in EDP.items():
            outline += str(v) + "\t"

        return outline

    else:
        return GetKmerEDP_Default()


def GetKmerEDP_Default():
    '''kmer entropy density'''
    Kmer = {}
    for aa in _Kmer_list:
        Kmer[aa] = 1e-9

    sum_Kmer = 1e-9 * 16 

    H = 0.0
    for (k,v) in Kmer.items():
        Kmer[k] /= sum_Kmer
        Kmer[k] = -Kmer[k] * np.log2(Kmer[k])
        H += Kmer[k]

    EDP = {}
    for (k,v) in Kmer.items():
        EDP[k] = Kmer[k] / H

    outline = ''
    for (k,v) in EDP.items():
        outline += str(v) + "\t"

    return outline
