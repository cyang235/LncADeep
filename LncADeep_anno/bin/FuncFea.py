#!/usr/bin/env python
import os
import sys
import numpy as np
import EdpFea

from utils import GetFasta
from Hexamer import ReadLogScore

from string import maketrans



__author__ = "Cheng Yang"
__copyright__ = "Copyright (c) 2016 Cheng Yang"
__license__ = "The MIT License (MIT)"
__version__ = "1.0"
__maintainer__ = "Cheng Yang"
__email__ = "ycheng@gatech.edu"

# AA list
_AA_list = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',
            'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y']

# 7 AA groups list
_AA_group = {'A': 'a', 'G': 'a', 'V': 'a', 'I': 'b', 'L': 'b', 'F': 'b', 'P': 'b',
             'Y': 'c', 'M': 'c', 'T': 'c', 'S': 'c', 'H': 'd', 'N': 'd', 'Q': 'd', 'W': 'd',
             'R': 'e', 'K': 'e', 'D': 'f', 'E': 'f', 'C': 'g'}

# 3mer of AA
_AA = ['a', 'b', 'c', 'd', 'e', 'f', 'g']
_AA_3mer = []
for aa1 in _AA:
    for aa2 in _AA:
        for aa3 in _AA:
            _AA_3mer.append(aa1 + aa2 + aa3)


# 2mer list
_DNA = ['A', 'C', 'G', 'T']
_Kmer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        _Kmer_list.append(dna1 + dna2)

# 3mer list
_3mer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            _3mer_list.append(dna1 + dna2 + dna3)


# 4mer list
_4mer_list = []
for dna1 in _DNA:
    for dna2 in _DNA:
        for dna3 in _DNA:
            for dna4 in _DNA:
                _4mer_list.append(dna1 + dna2 + dna3 + dna4)


# IUPAC code
_IUPAC = {'A': 'A', 'C': 'C', 'G': 'G', 'T': 'T', 'R': 'AG', 
          'Y': 'CT', 'M': 'AC', 'K': 'GT', 'S': 'CG', 'W': 'AT', 'H': 'ACT',
          'B': 'CGT', 'V': 'ACG', 'D': 'AGT', 'N': 'ACGT'}


def IUPAC_2mer(seq):
    '''Return a list of all possible 2mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            kmer_list.append(dna1 + dna2)
    return kmer_list


def IUPAC_3mer(seq):
    '''Return a list of all possible 3mers of the sequence'''

    kmer_list = []
    for dna1 in _IUPAC[seq[0]]:
        for dna2 in _IUPAC[seq[1]]:
            for dna3 in _IUPAC[seq[2]]:
                if Codon2AA2(dna1 + dna2 + dna3) != "J":
                    kmer_list.append(dna1 + dna2 + dna3)
    return kmer_list


def ReadPair(pair_file):
    '''Read lncRNA_protein pairs'''
    try:
        f = open(pair_file, 'r')
    except (IOError, ValueError) as e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    pair_list = []
    for line in f.readlines():
        line = line.strip()
        pair_list.append(line)

    f.close()

    return pair_list


def ReadLncRNA(lncRNA_fa):
    '''Read lncRNAs fasta'''

    SeqID, SeqList = GetFasta(lncRNA_fa)

    # lncRNA sequence
    lncRNA_dict = {}

    intap = "U"
    outap = "T"
    transtap = maketrans(intap, outap)

    for seqid, seq in zip(SeqID, SeqList):
        lncRNA_dict[seqid] = seq.translate(transtap)

    return lncRNA_dict


def ReadProtein(pro_fa):
    '''Read proteins fasta'''

    SeqID, SeqList = GetFasta(pro_fa)

    # lncRNA sequence
    protein_dict = {}

    for seqid, seq in zip(SeqID, SeqList):
        protein_dict[seqid] = seq

    return protein_dict


def GetRNAfea(NN_seq, logscore_dict=None):
    '''feature for rna'''
    nn_kmer_list = kmer_list(NN_seq, 4)
    nn_kmer_freq = NN_kmer_freq(nn_kmer_list)
    nn_edp_fea = GetEDP(nn_kmer_freq)

    rna_lncfea = EdpFea.GetFeature(NN_seq, "tmp", logscore_dict)

    return [nn_edp_fea, rna_lncfea]


def GetPROfea(AA_seq):
    '''feature for protein'''
    aa_kmer_list = kmer_list(AA_seq, 3)
    aa_kmer_freq = AA_kmer_freq(aa_kmer_list)
    aa_edp_fea = GetEDP(aa_kmer_freq)

    return aa_edp_fea


def GetBothFeature(NN_seq, AA_seq, logscore_dict=None):
    '''get features from NN and AA, for NN, consider two parts'''

    feature_1 = GetFeature(NN_seq, AA_seq)
    feature_2 = EdpFea.GetFeature(NN_seq, "tmp", logscore_dict)

    return feature_1 + "\t" + feature_2


def GetFeature(NN_seq, AA_seq):
    '''kmer features for NN and AA'''

    nn_kmer_list = kmer_list(NN_seq, 4)
    nn_kmer_freq = NN_kmer_freq(nn_kmer_list)
    nn_edp_fea = GetEDP(nn_kmer_freq)

    aa_kmer_list = kmer_list(AA_seq, 3)
    aa_kmer_freq = AA_kmer_freq(aa_kmer_list)
    aa_edp_fea = GetEDP(aa_kmer_freq)

    return nn_edp_fea + "\t" + aa_edp_fea


def GetFeatureName(NN_seq, AA_seq):
    '''kmer features names for NN and AA'''

    nn_kmer_list = kmer_list(NN_seq, 4)
    nn_kmer_freq = NN_kmer_freq(nn_kmer_list)
    nn_edp_fea_name = ",".join(nn_kmer_freq.keys())

    aa_kmer_list = kmer_list(AA_seq, 3)
    aa_kmer_freq = AA_kmer_freq(aa_kmer_list)
    aa_edp_fea_name = ",".join(aa_kmer_freq.keys())

    return "class," + nn_edp_fea_name + "," + aa_edp_fea_name


def Structure_fea(fea_file):
    '''read sturcture features'''
    try:
        f = open(fea_file, 'r')
    except (IOError, ValueError) as e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    fea_dict = {}

    for line in f.readlines():
        line = line.strip()
        if len(line.split()) < 5:
            continue
        rna = line.split()[0]
        fea = "\t".join(line.split()[1:])
        fea_dict[rna] = fea

    f.close()

    return fea_dict



def kmer_list(seq, k):
    '''return the kmer list of seq'''

    kmer_list = []

    for i in range(len(seq) - k + 1):
        kmer_list.append(seq[i:(i + k)])

    return kmer_list


def NN_kmer_freq(kmer_list):
    '''return the kmer frequency of nucleotide'''

    NN_4mer_freq = {}

    for nn in _4mer_list:
        NN_4mer_freq[nn] = 1e-12

    for k in NN_4mer_freq.keys():
        NN_4mer_freq[k] += kmer_list.count(k) * 1.0 / len(kmer_list)

    return NN_4mer_freq


def AA2Group(kmer):
    '''return the group of AA kmer of protein sequence'''
    new_kmer = ''
    for aa in kmer:
        if _AA_group.has_key(aa):
            new_kmer += _AA_group[aa]
        else:
            new_kmer += aa

    return new_kmer


def AA_kmer_freq(kmer_list):
    '''return the kmer frequency of AA of protein sequence'''

    for i in range(len(kmer_list)):
        kmer_list[i] = AA2Group(kmer_list[i])

    AA_kmer_freq = {}

    for aa in _AA_3mer:
        AA_kmer_freq[aa] = 1e-12

    for k in AA_kmer_freq.keys():
        AA_kmer_freq[k] += kmer_list.count(k) * 1.0 / len(kmer_list)

    return AA_kmer_freq


def GetEDP(kmer_freq):
    '''return the EDP of kmer_freq'''

    H = 0.0
    for k, v in kmer_freq.items():
        kmer_freq[k] = -v * np.log2(v)
        H += kmer_freq[k]

    EDP = {}
    for k, v in kmer_freq.items():
        EDP[k] = v / H
        if EDP[k] < 1e-7:
            EDP[k] = 0

    # print EDP
    outline = ''
    for (k, v) in EDP.items():
        outline += str(v) + "\t"

    return outline.strip()

