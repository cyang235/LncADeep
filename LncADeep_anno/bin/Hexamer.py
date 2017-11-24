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


def HexamerFrequency(fasta):
    '''count the hexamer usage as features'''

    HexamerCount = {}
    for k in _6mer_list:
        HexamerCount[k] = 0.0

    SeqID, SeqList = GetFasta(fasta) 

    totalcount = 0.0
    for seq in SeqList:
        ORF = GetORF(seq)
        if(len(ORF) > 3):
            num = len(ORF) / 3
            for i in range(0,num-1) :
                totalcount += 1.0
                tmp = ORF[ i*3:(i+2)*3] 
                if HexamerCount.has_key(tmp):
                    HexamerCount[tmp] += 1.0

    for k,v in HexamerCount.items():
        HexamerCount[k] = v / totalcount

    return HexamerCount 


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

    return hexamer_score_list


def ReadLogScore(logscore_file):
    '''return a dict of logscore'''

    logscore_dict = {}

    with open(logscore_file, "rU") as fl:
        for line in fl.readlines():
            line = line.strip()
            logscore_dict[ line.split()[0] ] = float(line.split()[1])

        fl.close()

        return logscore_dict


def MSSL(array, length):
    '''compute the maximum subarray with at least #length numbers'''

    if len(array) < length:
        print "The array is too short!"
        return -1

    SS = 0
    for i in range(length):
        SS += array[i]

    # the MM best value
    MM = SS

    # the global best value
    best = MM

    # the index of current, start, end
    cur = 0
    start = 0
    end = length 

    for k in range(length, len(array)):

        SS = SS + array[k] - array[k - length]

        if MM + array[k] > SS:
            MM = MM + array[k]
        else:
            MM = SS
            cur = k - length + 1

        if MM > best:
            best = MM
            start = cur
            end = k + 1

    return [start, end, best]


def MLC(seq, logscore_dict):
    '''find the most likely coding regions'''

    hexamer_score_list = [ HexamerGenerate(seq[i:], logscore_dict) for i in range(3)]

    start_index = [0,0,0]
    end_index = [0,0,0]
    best = [0,0,0]

    for i in range(3):
        start_index[i], end_index[i], best[i] = MSSL(hexamer_score_list[i], 3)

    phase = np.argmax(best)
    start = start_index[ phase ]
    end = end_index[ phase ]

    mlc_start = start*3 + phase
    mlc_end = (end+1)*3 + phase

    mlc_seq = seq[ (start*3 + phase):((end+1)*3 + phase) ]
    average_hexamer_score = best[phase] / (end - start)

    return [mlc_seq, average_hexamer_score, mlc_start, mlc_end]


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


def GetORF(seq):

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

    ORFseq = []
    ORFlen = []
    ORFstart = []
    ORFend = []
    for (k,v) in ORF.items():
        ORFseq.append(seq[k:(v+3)])
        ORFlen.append(v - k)
        ORFstart.append(k)
        ORFend.append(v)

    if ORF:
        ORF_l = ORFseq[np.argmax(ORFlen)]
        return ORF_l
    else:
        return ''


def main():
    args = sys.argv

    if len(args) < 2:
        print "Usage: this.py nfasta cfasta outfile"
    else:
        nfasta = args[1]
        cfasta = args[2]
        outfile = args[3]
        Hexamer_train(nfasta, cfasta, outfile)


if __name__ == "__main__":
    main()
