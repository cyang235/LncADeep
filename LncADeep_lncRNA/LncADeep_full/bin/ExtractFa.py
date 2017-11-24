#!/usr/bin/env python

__author__ = "Cheng Yang"
__copyright__ = "Copyright (c) 2016 Cheng Yang" 
__license__ = "The MIT License (MIT)" 
__version__ = "1.0"
__maintainer__ = "Cheng Yang"
__email__ = "ycheng@gatech.edu"

import sys,os
from utils import GetFasta

def Extract(infa, results, outfa):
    '''Extract lncRNA sequence'''

    SeqID, SeqList = GetFasta(infa)

    try:
        fr = open(results, 'rU')    # results file
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    idlist = set()
    for line in fr.readlines()[1:]:
        line = line.strip()
        if line.split("\t")[2] == "Noncoding":
            idlist.add( line.split("\t")[0] )
    fr.close()

    try:
        fo = open(outfa, 'w')    # output file
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e)
        sys.exit(1)

    for seqid, seq in zip(SeqID, SeqList):
        if seqid in idlist:
            fo.write(">" + seqid + "\n")
            fo.write(seq + "\n")
    fo.close()
