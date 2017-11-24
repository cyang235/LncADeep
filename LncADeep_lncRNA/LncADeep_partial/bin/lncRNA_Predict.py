import json
import sys,os
import numpy as np

from itertools import repeat
import multiprocessing as mul

import network
from utils import *

import EdpFea
from HmmFea import RunHMM
from HmmFea import ReadHmm
from HmmFea import GenerateTrans
from Hexamer import ReadLogScore 

from ExtractFa import Extract

def Normalize(data, mean, stdvar):
    '''Normalize input data'''
    data = ( data - mean ) / stdvar
    return data

def predict(filename=None, output_prefix=None, para_list_file=None, para_prefix=None, log_hexamer=None, 
        species="human", thread=1, HMMthread=8):

    ## make dir for prediction results
    tmpdir = os.path.join(os.path.abspath('.'), output_prefix + "_LncADeep_lncRNA_results")
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)
    print tmpdir

    exedir = os.path.dirname(os.path.abspath(__file__)) 

    ## translate fasta to protein and HMMER
    tmp_protein_file = os.path.join(tmpdir, output_prefix + "_protein.fasta")
    GenerateTrans(filename, tmp_protein_file)
    tmp_hmmer_prefix = os.path.join(tmpdir, output_prefix + "_hmmer")
    tmp_hmmer_out = os.path.join(tmpdir, output_prefix + "_hmmer.out2")
    pfam = os.path.join(exedir, "../../src/Pfam-A.hmm")
    RunHMM(tmp_protein_file, tmp_hmmer_prefix, pfam, thread=HMMthread)

    HMMdict1 = ReadHmm(tmp_hmmer_out)

    ## write results
    output = os.path.join(tmpdir, output_prefix + "_LncADeep.results") 
    try:
        fout = open(output, "w")
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e) 
        sys.exit(1)

    outlncRNA = os.path.join(tmpdir, output_prefix + "_LncADeep.predicted.lncRNA.fa") 

    ## default hexamer score
    if log_hexamer is None:
        log_hexamer = os.path.join(exedir, "../parameters/" + species + "/" + species + ".hexamer.logscore")
    logscore_dict = ReadLogScore(log_hexamer)

    ## a list of all parameter files
    ## default para_list_file
    if para_list_file is None:
        para_list_file = os.path.join(exedir, "../parameters/" + species + "/release.par.list")
    net_file_list, norm_file_list = ParaList(para_list_file)

    ## "para_prefix" refers to the dir which stores parameter files
    ## default para_prefix
    if para_prefix is None:
        para_prefix = os.path.join(exedir, "../parameters/" + species + "/")
    ## print para_prefix

    ## a list of all parameters
    mean_list = []
    stdvar_list = []
    net_list = []

    for net_file, norm_file in zip(net_file_list, norm_file_list):
        mean_list.append(load_para(para_prefix + norm_file)[0])
        stdvar_list.append(load_para(para_prefix + norm_file)[1])
        net_list.append(network.LoadNN(para_prefix + net_file))

    SeqID, SeqList = GetFasta(filename)

    fout.write("Transcript_ID\tMajorityVoteNum\tIndex\t21_ModelScores\n")

    para_tuple_list = zip( SeqID, SeqList, repeat(HMMdict1), repeat(logscore_dict), repeat(mean_list), repeat(stdvar_list), repeat(net_list) )

    ## single thread
    if thread == 1:
        for tmp_para_tuple in para_tuple_list:
            tmp_result = MajorityVoting(tmp_para_tuple)
            fout.write("\t".join(tmp_result) + "\n")

    ## multiple
    elif thread > 1:
        pool = mul.Pool(processes = thread)
        results = pool.map(MajorityVoting, para_tuple_list)
        pool.close()
        pool.join()

        for result in results:
            fout.write("\t".join(result) + "\n")

    fout.close()

    os.remove(tmp_protein_file)
    os.remove(tmp_hmmer_out)

    Extract(filename, output, outlncRNA)


def MajorityVoting(para_tuple):
    '''majority vote'''

    seqid = para_tuple[0]
    seq = para_tuple[1]
    HMMdict1 = para_tuple[2]
    logscore_dict = para_tuple[3]
    mean_list = para_tuple[4]
    stdvar_list = para_tuple[5]
    net_list = para_tuple[6]
    
    feature = EdpFea.GetFeature(seq, seqid, HMMdict1, logscore_dict) 
    if feature is not None:
        data_raw = np.fromstring(feature,sep=' ')

        pred_list = []
        label_vote = 0 
        for i in range(len(net_list)):
            data = Normalize(data_raw, mean_list[i], stdvar_list[i])
            data = data.reshape(data.shape[0], 1) 
            pred = net_list[i].propup(data)
            pred_list.append(str(pred[0][0]))
            if pred >= 0.5:
                label_vote += 1

        if label_vote >= 11:
            label_pred = "Coding"
        else:
            label_pred = "Noncoding"

        return [seqid, str(label_vote), label_pred, "\t".join(pred_list)]


def ParaList(para_list):
    '''load various ratio parameters'''
    net_list = []
    norm_list = []
    with open(para_list, "rU") as fp:
        for line in fp.readlines():
            line = line.strip()
            if line.find("net") != -1:
                net_list.append(line)
            if line.find("norm") != -1:
                norm_list.append(line)

    return [net_list, norm_list]


def load_para(filename):
    '''load the normalization parameters'''
    try:
        f = open(filename, 'r')
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e) 
        sys.exit(1)

    data = json.load(f)
    f.close()
    mean = data["mean"]
    stdvar = data["stdvar"]

    return [mean, stdvar]

