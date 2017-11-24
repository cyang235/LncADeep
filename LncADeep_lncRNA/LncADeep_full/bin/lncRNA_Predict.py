import json
import sys,os
import numpy as np

from itertools import repeat
import multiprocessing as mul

import EdpFea
import network

from utils import *
from HmmFea import RunHMM
from HmmFea import ReadHmm
from HmmFea import GenerateTrans
from Hexamer import ReadLogScore 

from ExtractFa import Extract

def Normalize(data, mean, stdvar):
    '''Normalize input data'''
    data = ( data - mean ) / stdvar
    return data

def predict(filename=None, net_para_file=None, norm_para_file=None, output_prefix=None, 
            log_hexamer=None, species="human", thread=1, HMMthread=8):

    ## make dir for prediction files 
    tmpdir = os.path.join(os.path.abspath('.'), output_prefix + "_LncADeep_lncRNA_results")
    exedir = os.path.dirname(os.path.abspath(__file__)) 

    print tmpdir
    if not os.path.exists(tmpdir):
        os.mkdir(tmpdir)

    ## translate fasta to protein and HMMER
    tmp_protein_file = os.path.join(tmpdir, output_prefix + "_protein.fasta")
    GenerateTrans(filename, tmp_protein_file)
    tmp_hmmer_prefix = os.path.join(tmpdir, output_prefix + "_hmmer")
    tmp_hmmer_out = os.path.join(tmpdir, output_prefix + "_hmmer.out2")
    pfam = os.path.join(exedir, "../../src/Pfam-A.hmm")
    RunHMM(tmp_protein_file, tmp_hmmer_prefix, pfam, thread=HMMthread)

    HMMdict1 = ReadHmm(tmp_hmmer_out)

    if species == "human":
        if log_hexamer is None:
            log_hexamer = os.path.join(exedir, "../parameters/human.hexamer.logscore")
        logscore_dict = ReadLogScore(log_hexamer)

        if net_para_file is None:
            net_para_file = os.path.join(exedir, "../parameters/human.net.para")
        if norm_para_file is None:
            norm_para_file = os.path.join(exedir, "../parameters/human.norm.para")
    elif species == "mouse":
        if log_hexamer is None:
            log_hexamer = os.path.join(exedir, "../parameters/mouse.hexamer.logscore")
        logscore_dict = ReadLogScore(log_hexamer)

        if net_para_file is None:
            net_para_file = os.path.join(exedir, "../parameters/mouse.net.para")
        if norm_para_file is None:
            norm_para_file = os.path.join(exedir, "../parameters/mouse.norm.para")

    ## write results
    output = os.path.join(tmpdir, output_prefix + "_LncADeep.results") 
    try:
        fout = open(output, "w")
    except (IOError,ValueError) as e:
        print >>sys.stderr, str(e) 
        sys.exit(1)

    outlncRNA = os.path.join(tmpdir, output_prefix + "_LncADeep.predicted.lncRNA.fa") 

    ## Load normalize parameters
    [mean, stdvar] = load_para(norm_para_file)

    ## Load NN
    net = network.LoadNN(net_para_file)

    SeqID, SeqList = GetFasta(filename)

    fout.write("Transcript_ID\tScore\tIndex\n")

    para_tuple_list = zip( SeqID, SeqList, repeat(HMMdict1), repeat(logscore_dict), repeat(mean), repeat(stdvar), repeat(net) )

    ## single thread
    if thread == 1:
        for tmp_para_tuple in para_tuple_list:
            tmp_result = Scoring(tmp_para_tuple)
            fout.write("\t".join(tmp_result) + "\n")

    ## multiple
    elif thread > 1:
        pool = mul.Pool(processes = thread)
        results = pool.map(Scoring, para_tuple_list)
        pool.close()
        pool.join()

        for result in results:
            fout.write("\t".join(result) + "\n")

    fout.close()

    os.remove(tmp_protein_file)
    os.remove(tmp_hmmer_out)

    Extract(filename, output, outlncRNA)



def Scoring(para_tuple):
    '''score'''

    seqid = para_tuple[0]
    seq = para_tuple[1]
    HMMdict1 = para_tuple[2]
    logscore_dict = para_tuple[3]
    mean = para_tuple[4]
    stdvar = para_tuple[5]
    net = para_tuple[6]
    
    feature = EdpFea.GetFeature(seq, seqid, HMMdict1, logscore_dict) 
    if feature is not None:
        data = np.fromstring(feature,sep=' ')
        data = Normalize(data, mean, stdvar)
        data = data.reshape(data.shape[0], 1) 
        pred = net.propup(data)

        if pred >= 0.5:
            label_pred = "Coding"
        else:
            label_pred = "Noncoding"

        return [seqid, str(pred[0][0]), label_pred]


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

