#!/usr/bin/env python

import json
import sys,os
import numpy as np
import random
import subprocess


def Read_mrmr(mrmr_file, fea_num=100):
    '''read mrmr_file and return the index of keeping features'''

    fm = open(mrmr_file, "rU")

    keep_fea = [0]  ## the label
    mrmr = fm.readlines()

    start = fea_num + 9
    end = start + fea_num

    for line in mrmr[start:end]:
        line = line.strip()
        keep_fea.append(int(line.split()[1]))

    fm.close()
    return keep_fea

