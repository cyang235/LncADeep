#!/usr/bin/env python

__author__ = "Cheng Yang"
__copyright__ = "Copyright (c) 2016 Cheng Yang"
__license__ = "The MIT License (MIT)"
__version__ = "1.0"
__maintainer__ = "Cheng Yang"
__email__ = "ycheng@gatech.edu"


import random
import sys
import os
import numpy as np
import json
import subprocess

from StringIO import StringIO

from Hexamer import ReadLogScore
from utils import GetFasta

from FuncFea import GetBothFeature
from FuncFea import ReadLncRNA
from FuncFea import ReadPair
from FuncFea import ReadProtein
from FuncFea import Structure_fea
from FuncFea import GetPROfea
from FuncFea import GetRNAfea

from SelectFeatures import Read_mrmr

import pandas as pd
from pandas import DataFrame
from pandas import Series

from PPI import SampleFreq
from PPI import ReactomeSampleFreq
from PPI import RunMCL
from PPI import ReadBackFreq
from PPI import ReadUniprotAnnoMeta
from PPI import ReadReactomeBackFreq
from PPI import ReadReactomeAnnoMeta

seed = 1234567

np.random.seed(seed)


def GenAllTestFea(pair_file=None, rna_fea_1=None, rna_fea_2=None,
               pro_fea=None, protein_s_dict=None, rna_s_dict=None):
    """generate all features for all pairs"""

    lncRNA_prot_pair = ReadPair(pair_file)

    pair_features = ""
    keep_pair = []
    for t_pair in lncRNA_prot_pair:
        t_lncRNA = t_pair.split()[0]
        t_protein = t_pair.split()[1]

        if not (rna_fea_1.has_key(t_lncRNA) and rna_fea_2.has_key(t_lncRNA)
                and pro_fea.has_key(t_protein)):
            continue

        if not (rna_s_dict.has_key(t_lncRNA) and protein_s_dict.has_key(t_protein)):
            continue

        tmpfeature1 = rna_fea_1[t_lncRNA] + "\t" + \
            pro_fea[t_protein] + "\t" + rna_fea_2[t_lncRNA]
        tmpfeature2 = rna_s_dict[t_lncRNA] + "\t" + protein_s_dict[t_protein]
        tmpfeature = tmpfeature1 + "\t" + tmpfeature2

        pair_features += tmpfeature + "\n"
        keep_pair.append(t_pair)

    return [keep_pair, pair_features]


def GeneratePairs(rna_file, pro_file, out_prefix):
    '''generate all rna_protein pairs'''

    rnaID = GetFasta(rna_file)[0]
    proID = GetFasta(pro_file)[0]

    # generate can
    pair_file_list = []
    rna_id_list = []

    for rnaid in rnaID:
        tmp_pair = out_prefix + "." + rnaid.split("|")[0]
        with open(tmp_pair, "w") as fo:
            for proid in proID:
                fo.write(rnaid + " " + proid + "\n")
            pair_file_list.append(tmp_pair)
            rna_id_list.append(rnaid)
        fo.close()

    return [rna_id_list, pair_file_list]


def SeparatePairs(test_pair, out_prefix):
    '''generate all rna_protein pairs'''

    ft = open(test_pair, "rU")

    # generate can
    pair_file_list = []
    rna_id_list = []
    rnaid_proid = {}

    for line in ft.readlines():
        line = line.strip()
        rnaid = line.split()[0]
        proid = line.split()[1]

        if rnaid_proid.has_key(rnaid):
            rnaid_proid[rnaid].add(proid)
        else:
            rnaid_proid[rnaid] = set()
            rnaid_proid[rnaid].add(proid)
    ft.close()

    for k, v in rnaid_proid.items():
        tmp_pair = os.path.join(out_prefix, k)
        with open(tmp_pair, "w") as fo:
            for proid in v:
                fo.write(k + " " + proid + "\n")
            pair_file_list.append(tmp_pair)
            rna_id_list.append(k)
        fo.close()

    return [rna_id_list, pair_file_list]


def SavePara(para, para_file):
    '''save parameters'''
    try:
        f = open(para_file, "w")
    except (IOError, ValueError) as e:
        print >> sys.stderr, e
        sys.exit(1)

    json.dump(para, f)
    f.close()


def LoadPara(para_file):
    '''load parameters'''
    try:
        f = open(para_file, "rU")
    except (IOError, ValueError) as e:
        print >> sys.stderr, e
        sys.exit(1)

    data = json.load(f)
    all_mean = data['mean']
    all_std = data['stdvar']
    f.close()

    return [all_mean, all_std]


def GetParameters(all_mean, all_std):
    '''Return the normalization parameters'''
    data = {'mean': [i.tolist() for i in all_mean],
            'stdvar': [j.tolist() for j in all_std]
            }

    return data


def RNA_StructureScore(rna_file, out_prefix=None):
    '''to compute structure features of rnas and proteins'''

    rna_out = os.path.join(out_prefix, "rna_score")

    #####################
    # lncRNA structure
    rna_file_part = os.path.join(out_prefix, "tmp.rna.file.")
    rna_file_list = os.path.join(out_prefix, "tmp.filelist")

    rnaID, rnaSeq = GetFasta(rna_file)
    i = 0
    for rnaid, rnaseq in zip(rnaID, rnaSeq):
        f_tmp = open(rna_file_part + str(i), "w")
        i += 1
        f_tmp.write(">" + rnaid + "\n")
        f_tmp.write(rnaseq + "\n")
        f_tmp.close()

    file_list_cmd = "ls " + out_prefix + " |grep tmp.rna.file > " + rna_file_list
    #print file_list_cmd
    subprocess.call(file_list_cmd, shell=True)

    exedir = os.path.dirname(os.path.abspath(__file__))
    RNAScore2 = os.path.join(exedir, "../src/RNAScore2")

    with open(rna_file_list, "rU") as fr:
        for tmp in fr.readlines():
            tmp = tmp.strip()
            tmpfile = os.path.join(out_prefix, tmp)
            tmpout = os.path.join(out_prefix, tmp + ".r_score")
            rna_cmd = RNAScore2 + " -i " + tmpfile + " -o " + tmpout + " -l 250 -r"
            #print rna_cmd
            subprocess.call(rna_cmd, shell=True)

            combine_cmd = "cat " + tmpout + " >> " + rna_out
            #print combine_cmd
            subprocess.call(combine_cmd, shell=True)
            os.remove(tmpfile)
            os.remove(tmpout)
        fr.close()
    #####################

    os.remove(rna_file_list)

    return rna_out


def Protein_StructureScore(protein_file, out_prefix=None):
    '''to compute structure features of rnas and proteins'''

    protein_out = os.path.join(out_prefix, "protein_score")

    #####################
    # protein structure
    protein_file_part = os.path.join(out_prefix, "tmp.protein.file.")
    protein_file_list = os.path.join(out_prefix, "tmp.filelist")

    exedir = os.path.dirname(os.path.abspath(__file__))
    stride_dat = os.path.join(exedir, "../src/stride.dat")
    stride_cmd = "cp " + stride_dat + " " + os.path.abspath('.')
    tmp_stride_dat = os.path.join(os.path.abspath('.'), "stride.dat")
    subprocess.call(stride_cmd, shell=True)

    proID, proSeq = GetFasta(protein_file)
    i = 0
    for proid, proseq in zip(proID, proSeq):
        f_tmp = open(protein_file_part + str(i), "w")
        i += 1
        f_tmp.write(">" + proid + "\n")
        f_tmp.write(proseq + "\n")
        f_tmp.close()

    file_list_cmd = "ls " + out_prefix + \
        " |grep tmp.protein.file > " + protein_file_list
    #print file_list_cmd
    subprocess.call(file_list_cmd, shell=True)

    exedir = os.path.dirname(os.path.abspath(__file__))
    RNAScore2 = os.path.join(exedir, "../src/RNAScore2")

    with open(protein_file_list, "rU") as fp:
        for tmp in fp.readlines():
            tmp = tmp.strip()
            tmpfile = os.path.join(out_prefix, tmp)
            tmpout = os.path.join(out_prefix, tmp + ".pro_score")
            protein_cmd = RNAScore2 + " -i " + tmpfile + " -o " + tmpout + " -p"
            #print protein_cmd
            subprocess.call(protein_cmd, shell=True)

            combine_cmd = "cat " + tmpout + " >> " + protein_out
            #print combine_cmd
            subprocess.call(combine_cmd, shell=True)
            os.remove(tmpfile)
            os.remove(tmpout)

        fp.close()
    #####################

    os.remove(protein_file_list)
    os.remove(tmp_stride_dat)

    return protein_out


def GenAnnoEDPfeature(rna_file, pro_file, logscore_dict=None):
    '''generate rna, protein features'''

    rna_fea_1 = {}
    rna_fea_2 = {}
    pro_fea = {}

    rna_ID, rna_Seq = GetFasta(rna_file)

    for rna_id, rna_seq in zip(rna_ID, rna_Seq):
        nn_edp_fea, rna_lncfea = GetRNAfea(rna_seq, logscore_dict)
        rna_fea_1[rna_id] = nn_edp_fea
        rna_fea_2[rna_id] = rna_lncfea

    exedir = os.path.dirname(os.path.abspath(__file__))
    pro_fea_file = os.path.join(exedir, "../src/Swiss-Uniprot.human.protein.seq_fea")
    fp1 = open(pro_fea_file, "rU")
    for line in fp1.readlines():
        line = line.strip()
        pro_fea[line.split()[0]] = "\t".join(line.split()[1:])
    fp1.close()

    return [rna_fea_1, rna_fea_2, pro_fea]


def GenEDPfeature(rna_file, pro_file, logscore_dict=None):
    '''generate rna, protein features'''

    rna_fea_1 = {}
    rna_fea_2 = {}
    pro_fea = {}

    rna_ID, rna_Seq = GetFasta(rna_file)
    pro_ID, pro_Seq = GetFasta(pro_file)

    for rna_id, rna_seq in zip(rna_ID, rna_Seq):
        nn_edp_fea, rna_lncfea = GetRNAfea(rna_seq, logscore_dict)
        rna_fea_1[rna_id] = nn_edp_fea
        rna_fea_2[rna_id] = rna_lncfea

    for pro_id, pro_seq in zip(pro_ID, pro_Seq):
        aa_edp_fea = GetPROfea(pro_seq)
        pro_fea[pro_id] = aa_edp_fea

    return [rna_fea_1, rna_fea_2, pro_fea]


def Annotate(test_pair=None, rna_file=None, out_prefix="tmp", para_num=33, mrmr_num=110):
    '''Annotate lncRNAs'''

    # testing/prediction working directory
    Test_workdir = os.path.join(
        os.path.abspath('.'), out_prefix + "_LncADeep_anno_results")
    if not os.path.exists(Test_workdir):
        os.mkdir(Test_workdir)

    # executive file directory
    exedir = os.path.dirname(os.path.abspath(__file__))

    #####################################
    # default annotation files 
    pro_file = os.path.join(exedir, "../src/Swiss-Uniprot.human.protein.fa") 
    protein_struct = os.path.join(exedir, "../src/Swiss-Uniprot.human.protein.struct_fea") 

    backfreq = os.path.join(exedir, "../src/Uniprot_KEGG_pathway_background.dat") 
    ko_freq = ReadBackFreq(backfreq=backfreq)

    react_backfreq = os.path.join(exedir, "../src/UniProt2Reactome_Background.dat") 
    react_freq = ReadReactomeBackFreq(backfreq=react_backfreq)

    uniprot_anno_file = os.path.join(exedir, "../src/Uniprot_KEGG.ANNO.meta")
    protein_anno_dict = ReadUniprotAnnoMeta(uniprot_anno_file=uniprot_anno_file)

    react_uniprot_anno_file = os.path.join(exedir, "../src/Human_UniProt2Reactome_filtered.txt")
    react_protein_anno_dict = ReadReactomeAnnoMeta(uniprot_anno_file=react_uniprot_anno_file)

    # default files for functional modules
    hippie_pairs = os.path.join(exedir, "../src/hippie_current.PPI.pairs") 
    #####################################

    #####################################
    # normalization parameter file
    para_dir = os.path.join(exedir, "../parameters") 
    para_prefix = ["6204_release_" + str(i+1) for i in range(para_num)]

    # default network parameters
    test_pair = os.path.join(Test_workdir, out_prefix + "_pairs")
    rna_id_list, pair_file_list = GeneratePairs(rna_file, pro_file, test_pair)
    #####################################

    ###################################
    # generate structure score
    rna_struct = RNA_StructureScore(rna_file, out_prefix=Test_workdir)
    ###################################

    #####################################
    # load files
    # Hexamer
    log_hexamer = os.path.join(exedir, "../src/Gencode.Refseq.logscore.logscore")
    logscore_dict = ReadLogScore(log_hexamer)

    # structure feature
    rna_s_dict = Structure_fea(rna_struct)
    protein_s_dict = Structure_fea(protein_struct)
    os.remove(rna_struct)

    # kmer edp features
    rna_fea_1, rna_fea_2, pro_fea = GenAnnoEDPfeature(rna_file, pro_file, logscore_dict)
    #####################################

    #####################################
    # generate testing features for each lncRNA
    # and output for each lncRNA
    for rna_id, pair_file in zip(rna_id_list, pair_file_list):
        # multiple tests for each lncRNA
        rna_protein_vote = {} 

        keep_pair, pair_features = GenAllTestFea(pair_file=pair_file, 
                        rna_fea_1=rna_fea_1, rna_fea_2=rna_fea_2, pro_fea=pro_fea, 
                        protein_s_dict=protein_s_dict, rna_s_dict=rna_s_dict)

        for tmp_para_prefix in para_prefix:
            norm_name = os.path.join(para_dir, tmp_para_prefix + "_norm.json")
            mrmr_file = os.path.join(para_dir, tmp_para_prefix + "_mrmr.out")

            # read mrmr_file
            keep_fea = Read_mrmr(mrmr_file, fea_num=mrmr_num)

            prediction_results = RunDNN(test_fea=pair_features, norm_name=norm_name, para_dir=para_dir, 
                    para_prefix=tmp_para_prefix, keep_fea=keep_fea)

            for pair,pred in zip(keep_pair, prediction_results):
                tmp_pro = pair.split()[1].split("|")[1]
                if pred >= 0.5:
                    if rna_protein_vote.has_key(tmp_pro):
                        rna_protein_vote[tmp_pro] += 1
                    else:
                        rna_protein_vote[tmp_pro] = 1
                else:
                    if not rna_protein_vote.has_key(tmp_pro):
                        rna_protein_vote[tmp_pro] = 0

        # output results
        outfile = os.path.join(Test_workdir, out_prefix +
                               "_" + rna_id.split("|")[0] + ".interaction.Results")
        outfile_protein = os.path.join(
            Test_workdir, out_prefix + "_" + rna_id.split("|")[0] + ".interacting.proteins")

        fo = open(outfile, "w")
        fop = open(outfile_protein, "w")
        fo.write("Pair\tIndex\tScore\n")

        for k, v in rna_protein_vote.items():
            if v > para_num * 0.5:
                pred_out = "Interacting"
                fop.write(k + "\n")
            else:
                pred_out = "Non-interacting"

            fo.write(rna_id + "\t" + k + ":\t" + pred_out + "\t" + str(v * 1.0/para_num) + "\n")

        fo.close()
        fop.close()

        ###############################
        # Annotation #1: KEGG and Reactome enrichment
        ###############################
        # the list of significant KEGG pathways
        outfile_KEGG_list = os.path.join(
            Test_workdir, out_prefix + "_" + rna_id.split("|")[0])
        SampleFreq(outfile_protein, ko_freq=ko_freq,
                   protein_anno_dict=protein_anno_dict, outprefix=outfile_KEGG_list)

        outfile_React_list = os.path.join(
            Test_workdir, out_prefix + "_" + rna_id.split("|")[0])
        ReactomeSampleFreq(outfile_protein, react_freq=react_freq,
                   protein_anno_dict=react_protein_anno_dict, outprefix=outfile_React_list)

        ###############################
        # Annotation #2: function module
        ###############################
        # the figure of functional modules
        outfile_module_list = os.path.join(
            Test_workdir, out_prefix + "_" + rna_id.split("|")[0])
        RunMCL(outfile_protein, hippie_pairs=hippie_pairs,
               outprefix=outfile_module_list)

        os.remove(pair_file)


def RunDNN(test_fea=None, norm_name=None, para_dir=None, para_prefix=None, keep_fea=None, mrmr_num=110):
    '''get the prediction results for features'''

    from keras.models import Sequential, load_model
    from keras.models import model_from_json
    from keras.layers import Dense, Dropout, Activation, Merge
    from keras.optimizers import SGD
    from keras.callbacks import ModelCheckpoint

    # generate the list of features to be used
    keep_fea = sorted([i-1 for i in keep_fea[1:]]) + [j+646 for j in range(80)]

    test_df = pd.read_table(StringIO(test_fea), header=None)
    test_df = test_df[keep_fea]
    column_names = [ i for i in range(mrmr_num + 80)]     # rename the columns
    test_df.columns = column_names

    ##################################
    # normalize datasets
    all_mean, all_std = LoadPara(norm_name)
    test_df = (test_df - all_mean[1:]) / all_std[1:]
    ##################################

    ##################################
    # process data
    test_Kmer = test_df.ix[:, 0:(mrmr_num - 1)]
    test_RNA_ss = test_df.ix[:, mrmr_num:(mrmr_num + 9)]
    test_RNA_hydro = test_df.ix[:, (mrmr_num + 10):(mrmr_num+19)]
    test_RNA_van = test_df.ix[:, (mrmr_num + 20):(mrmr_num + 29)]
    test_PRO_ss = test_df.ix[:, (mrmr_num + 30):(mrmr_num + 39)]
    test_PRO_hydro1 = test_df.ix[:, (mrmr_num + 40):(mrmr_num + 49)]
    test_PRO_hydro2 = test_df.ix[:, (mrmr_num + 50):(mrmr_num + 59)]
    test_PRO_van1 = test_df.ix[:, (mrmr_num + 60):(mrmr_num + 69)]
    test_PRO_van2 = test_df.ix[:, (mrmr_num + 70):(mrmr_num + 79)]
    ##############################################

    ##############################################
    # aggregate data
    test_Kmer = np.array(test_Kmer)
    test_ss = np.concatenate((test_RNA_ss, test_PRO_ss), axis=1)
    test_hydro = np.concatenate(
        (test_RNA_hydro, test_PRO_hydro1, test_PRO_hydro2), axis=1)
    test_van = np.concatenate(
        (test_RNA_van, test_PRO_van1, test_PRO_van2), axis=1)
    ##############################################

    ##############################################
    # load models

    ##########################
    # kmer branch
    kmer_json_name = os.path.join(
        para_dir, para_prefix + "_kmer_model.json")
    kmer_para_name = os.path.join(para_dir, para_prefix + "_kmer_model.h5")

    kmer_json_file = open(kmer_json_name, 'r')
    kmer_model_json = kmer_json_file.read()
    kmer_json_file.close()
    kmer_loaded_model = model_from_json(kmer_model_json)
    kmer_loaded_model.load_weights(kmer_para_name)

    # evaluate loaded model on test data
    kmer_loaded_model.compile(
        loss='binary_crossentropy', optimizer='rmsprop', metrics=['accuracy'])

    # output prediction for next cycle
    kmer_pred_test = kmer_loaded_model.predict([test_Kmer])

    test_Kmer_new = np.concatenate((test_Kmer, kmer_pred_test), axis=1)
    ##########################

    ##########################
    # structure branch
    str_json_name = os.path.join(
        para_dir, para_prefix + "_struct_model.json")
    str_para_name = os.path.join(
        para_dir, para_prefix + "_struct_model.h5")

    str_json_file = open(str_json_name, 'r')
    str_model_json = str_json_file.read()
    str_json_file.close()
    str_loaded_model = model_from_json(str_model_json)
    str_loaded_model.load_weights(str_para_name)

    # evaluate loaded model on test data
    str_loaded_model.compile(loss='binary_crossentropy',
                             optimizer='rmsprop', metrics=['accuracy'])

    # output prediction for next cycle
    structure_pred_test = str_loaded_model.predict(
        [test_ss, test_hydro, test_van])
    test_Struct_new = np.concatenate(
        (test_ss, test_hydro, test_van, structure_pred_test), axis=1)
    ##########################

    ##########################
    # merged model

    # No. 1 model
    json_name = os.path.join(para_dir, para_prefix + "_model.json.1")
    para_name = os.path.join(para_dir, para_prefix + "_model.h5.1")

    json_file = open(json_name, 'r')
    model_json = json_file.read()
    json_file.close()
    loaded_model = model_from_json(model_json)
    loaded_model.load_weights(para_name)

    # evaluate loaded model on test data
    loaded_model.compile(loss='binary_crossentropy',
                         optimizer='rmsprop', metrics=['accuracy'])

    # output prediction for next cycle
    prediction_test = loaded_model.predict(
        [test_Kmer_new, test_Struct_new])

    # No. 2 model
    for ii in range(5):
        test_Kmer_new = np.concatenate(
            (test_Kmer_new, prediction_test), axis=1)
        test_Struct_new = np.concatenate(
            (test_Struct_new, prediction_test), axis=1)

        json_name = os.path.join(
            para_dir, para_prefix + "_model.json." + str(ii + 2))
        para_name = os.path.join(
            para_dir, para_prefix + "_model.h5." + str(ii + 2))

        json_file = open(json_name, 'r')
        model_json = json_file.read()
        json_file.close()
        loaded_model = model_from_json(model_json)
        loaded_model.load_weights(para_name)

        # evaluate loaded model on test data
        loaded_model.compile(loss='binary_crossentropy',
                             optimizer='rmsprop', metrics=['accuracy'])

        # output prediction for next cycle
        prediction_test = loaded_model.predict(
            [test_Kmer_new, test_Struct_new])
    ##########################

    return prediction_test


def Predict(test_pair=None, rna_file=None, pro_file=None, para_num=33, mrmr_num=110,
            para_dir=None, para_prefix=None, out_prefix="tmp"):
    '''Predict lncRNA-protein interactions'''

    from keras.models import Sequential, load_model
    from keras.models import model_from_json
    from keras.layers import Dense, Dropout, Activation, Merge
    from keras.optimizers import SGD
    from keras.callbacks import ModelCheckpoint

    # testing/prediction working directory
    Test_workdir = os.path.join(
        os.path.abspath('.'), out_prefix + "_LncADeep_anno_results")
    if not os.path.exists(Test_workdir):
        os.mkdir(Test_workdir)

    # executive file directory
    exedir = os.path.dirname(os.path.abspath(__file__))

    ###################################
    # Hexamer
    log_hexamer = os.path.join(exedir, "../src/Gencode.Refseq.logscore.logscore")
    logscore_dict = ReadLogScore(log_hexamer)

    # generate structure score
    rna_struct = RNA_StructureScore(rna_file, out_prefix=Test_workdir)
    ###################################

    ###################################
    # default parameters
    if pro_file is None:
        pro_file = os.path.join(exedir, "../src/Swiss-Uniprot.human.protein.fa") 
        protein_struct = os.path.join(exedir, "../src/Swiss-Uniprot.human.protein.struct_fea") 
        # kmer edp features
        rna_fea_1, rna_fea_2, pro_fea = GenAnnoEDPfeature(rna_file, pro_file, logscore_dict)
    else:
        protein_struct = Protein_StructureScore(pro_file, out_prefix=Test_workdir)
        # kmer edp features
        rna_fea_1, rna_fea_2, pro_fea = GenEDPfeature(rna_file, pro_file, logscore_dict)

    if test_pair is None:
        test_pair = os.path.join(Test_workdir, out_prefix + "_pairs")
        rna_id_list, pair_file_list = GeneratePairs(
            rna_file, pro_file, test_pair)
    else:
        rna_id_list, pair_file_list = SeparatePairs(test_pair, Test_workdir)
    ###################################

    #####################################
    # normalization parameter file
    para_dir = os.path.join(exedir, "../parameters") 
    para_prefix = ["6204_release_" + str(i+1) for i in range(para_num)]
    #####################################

    #####################################
    # load files
    # structure feature
    rna_s_dict = Structure_fea(rna_struct)
    protein_s_dict = Structure_fea(protein_struct)
    os.remove(rna_struct)
    if protein_struct.find("Swiss-Uniprot.human.protein.struct_fea") == -1:
        os.remove(protein_struct)
    #####################################

    #####################################
    # generate testing features for each lncRNA
    # and output for each lncRNA
    for rna_id, pair_file in zip(rna_id_list, pair_file_list):
        # multiple tests for each lncRNA
        rna_protein_vote = {} 

        keep_pair, pair_features = GenAllTestFea(pair_file=pair_file, 
                        rna_fea_1=rna_fea_1, rna_fea_2=rna_fea_2, pro_fea=pro_fea, 
                        protein_s_dict=protein_s_dict, rna_s_dict=rna_s_dict)

        for tmp_para_prefix in para_prefix:
            norm_name = os.path.join(para_dir, tmp_para_prefix + "_norm.json")
            mrmr_file = os.path.join(para_dir, tmp_para_prefix + "_mrmr.out")

            # read mrmr_file
            keep_fea = Read_mrmr(mrmr_file, fea_num=mrmr_num)

            prediction_results = RunDNN(test_fea=pair_features, norm_name=norm_name, para_dir=para_dir, 
                    para_prefix=tmp_para_prefix, keep_fea=keep_fea)

            for pair,pred in zip(keep_pair, prediction_results):
                tmp_pro = pair.split()[1].split("|")[1]
                if pred >= 0.5:
                    if rna_protein_vote.has_key(tmp_pro):
                        rna_protein_vote[tmp_pro] += 1
                    else:
                        rna_protein_vote[tmp_pro] = 1
                else:
                    if not rna_protein_vote.has_key(tmp_pro):
                        rna_protein_vote[tmp_pro] = 0

        # output results
        outfile = os.path.join(Test_workdir, out_prefix +
                               "_" + rna_id.split("|")[0] + ".interaction.Results")
        outfile_protein = os.path.join(
            Test_workdir, out_prefix + "_" + rna_id.split("|")[0] + ".interacting.proteins")

        fo = open(outfile, "w")
        fop = open(outfile_protein, "w")
        fo.write("Pair\tIndex\tScore\n")

        for k, v in rna_protein_vote.items():
            if v > (para_num / 2):
                pred_out = "Interacting"
                fop.write(k + "\n")
            else:
                pred_out = "Non-interacting"

            fo.write(rna_id + "\t" + k + ":\t" + pred_out + "\t" + str(v * 1.0 / para_num) + "\n")

        fo.close()
        fop.close()

        os.remove(pair_file)

