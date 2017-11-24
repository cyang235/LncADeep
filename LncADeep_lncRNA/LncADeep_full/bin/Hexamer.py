#!/usr/bin/env python

import os,sys
import numpy as np

from utils import GetFasta
from utils import Codon2AA2

## AA list 
_AA_list = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

## Di-codon list
_Di_Codon_list = []
for aa1 in _AA_list:
    for aa2 in _AA_list:
        _Di_Codon_list.append(aa1+aa2)

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

## 6mer list
_6mer_list = []
for mer1 in _3mer_list:
    for mer2 in _3mer_list:
        _6mer_list.append(mer1+mer2)


## IUPAC code
_IUPAC = {'A':'A','C':'C', 'G':'G', 'T':'T', 'R':'AG', 'Y':'CT', 'M':'AC', 'K':'GT', 'S':'CG', 'W':'AT', 'H':'ACT', \
          'B':'CGT', 'V':'ACG', 'D':'AGT', 'N':'ACGT'}

def IUPAC_2mer(seq):
    '''Return a list of all possible 2mers of the sequence
       Or randomly choose one from the list?
    '''

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



def HexamerGenerate(seq, logscore_dict):
    '''Generate hexamer'''

    hexamer_list = []
    hexamer_score_list = []

    if(len(seq) > 3):
        num = len(seq) / 3
        for i in range(0,num-1) :
            tmp = seq[ i*3:(i+2)*3] 
            hexamer_list.append(tmp)
            if logscore_dict.has_key(tmp):
                hexamer_score_list.append( logscore_dict[tmp] )
            else:
                hexamer_score_list.append(0)

    #return [hexamer_list, hexamer_score_list]
    return hexamer_score_list


def ReadLogScore(logscore_file):
    '''return a dict of logscore'''

    logscore_dict = {}

    with open(logscore_file, "rU") as fl:
        for line in fl.readlines():
            line = line.strip()
            logscore_dict[ line.split()[0] ] = float(line.split()[1])

        return logscore_dict



def HexamerScore2(seq, logscore_dict):
    '''compute the hexamer score given a sequence'''

    total = 0.0
    log_score = 0.0

    for k in HexamerGenerate(seq, logscore_dict):
        log_score += k

        total += 1

    try:
        return log_score/total
    except:
        return -1

