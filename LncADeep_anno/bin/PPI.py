#!/usr/bin/env python

import sys
import os
import subprocess

__author__ = "Cheng Yang"
__copyright__ = "Copyright (c) 2016 Cheng Yang"
__license__ = "The MIT License (MIT)"
__version__ = "1.0"
__maintainer__ = "Cheng Yang"
__email__ = "ycheng@gatech.edu"


exedir = os.path.dirname(os.path.abspath(__file__))

###############################################################
# for KEGG and Reactome pathway enrichment analysis


def ReadBackFreq(backfreq=None):
    '''read background frequency of KEGG'''

    if backfreq is None:
        backfreq = os.path.join(
            exedir, "../src/Uniprot_KEGG_pathway_background.dat")
    fb = open(backfreq)

    ko_path_num = {}
    for line in fb.readlines():
        line = line.strip()
        ko_path_num[line.split()[0]] = line.split()[1]

    fb.close()

    return ko_path_num


def ReadReactomeBackFreq(backfreq=None):
    '''read background frequency of Reactome'''

    if backfreq is None:
        backfreq = os.path.join(
            exedir, "../src/UniProt2Reactome_Background.dat")
    fb = open(backfreq)

    path_num = {}
    for line in fb.readlines():
        line = line.strip()
        path_num[line.split()[0]] = line.split()[1]

    fb.close()

    return path_num


def ReadUniprotAnnoMeta(uniprot_anno_file=None):
    '''read Uniprot_KEGG.ANNO.meta'''

    if uniprot_anno_file is None:
        uniprot_anno_file = os.path.join(
            exedir, "../src/Uniprot_KEGG.ANNO.meta")
    fu = open(uniprot_anno_file, "rU")

    protein_anno_dict = {}

    for line in fu.readlines():
        line = line.strip()
        pro_id = line.split("\t")[0]
        anno = set(line.split("\t")[1].split(","))
        protein_anno_dict[pro_id] = anno

    fu.close()

    return protein_anno_dict


def ReadReactomeAnnoMeta(uniprot_anno_file=None):
    '''read Human_UniProt2Reactome_filtered.txt'''

    if uniprot_anno_file is None:
        uniprot_anno_file = os.path.join(
            exedir, "../src/Human_UniProt2Reactome_filtered.txt")
    fu = open(uniprot_anno_file, "rU")

    protein_anno_dict = {}

    for line in fu.readlines():
        line = line.strip()
        pro_id = line.split("\t")[0]
        if protein_anno_dict.has_key(pro_id):
            protein_anno_dict[pro_id].add("\t".join(line.split("\t")[1:]))
        else:
            protein_anno_dict[pro_id] = set()
            protein_anno_dict[pro_id].add("\t".join(line.split("\t")[1:]))

    fu.close()

    return protein_anno_dict


def SampleFreq(inter_file, ko_freq=None, protein_anno_dict=None, outprefix=None):
    '''return the sample frequency for KEGG pathway analysis'''

    fi = open(inter_file, "rU")

    universe_num = len(protein_anno_dict.keys())

    sample_num = 0
    sample_anno = {}
    for line in fi.readlines():
        line = line.strip()
        if line in protein_anno_dict.keys():
            sample_num += 1
            for tmp_anno in protein_anno_dict[line]:
                if sample_anno.has_key(tmp_anno):
                    sample_anno[tmp_anno] += 1
                else:
                    sample_anno[tmp_anno] = 1
    fi.close()

    out_sample_freq = outprefix + ".ko.freq"
    fo = open(out_sample_freq, "w")

    for k, v in sample_anno.items():
        ko_id = "PATH:ko" + k.split()[0]
        ko_anno = " ".join(k.split()[1:-1])
        fo.write("\t".join([ko_id, str(v), str(ko_freq[ko_id]), str(
            sample_num), str(universe_num), ko_anno]) + "\n")

    fo.close()

    out_kegg_pathway = outprefix + ".KEGG.pathway"
    KEGG_rscript = os.path.join(exedir, "../src/KEGG.enrichment.R")
    KEGG_cmd = "Rscript " + KEGG_rscript + " " + \
        out_sample_freq + " " + out_kegg_pathway
    print KEGG_cmd
    subprocess.call(KEGG_cmd, shell=True)
    os.remove(out_sample_freq)


def ReactomeSampleFreq(inter_file, react_freq=None, protein_anno_dict=None, outprefix=None):
    '''return the sample frequency for Reactome pathway analysis'''

    fi = open(inter_file, "rU")

    universe_num = len(protein_anno_dict.keys())

    sample_num = 0
    sample_anno = {}
    for line in fi.readlines():
        line = line.strip()
        if line in protein_anno_dict.keys():
            sample_num += 1
            for tmp_anno in protein_anno_dict[line]:
                if sample_anno.has_key(tmp_anno):
                    sample_anno[tmp_anno] += 1
                else:
                    sample_anno[tmp_anno] = 1
    fi.close()

    out_sample_freq = outprefix + ".reactome.freq"
    fo = open(out_sample_freq, "w")

    for k, v in sample_anno.items():
        react_id = k.split("\t")[0]
        react_anno = k.split("\t")[2]
        fo.write("\t".join([react_id, str(v), str(react_freq[react_id]), str(
            sample_num), str(universe_num), react_anno]) + "\n")

    fo.close()

    out_react_pathway = outprefix + ".Reactome.pathway"
    Reactome_rscript = os.path.join(exedir, "../src/Reactome.enrichment.R")
    Reactome_cmd = "Rscript " + Reactome_rscript + " " + \
        out_sample_freq + " " + out_react_pathway
    print Reactome_cmd
    subprocess.call(Reactome_cmd, shell=True)
    os.remove(out_sample_freq)

###############################################################


###############################################################
# for functional module analysis
def RunMCL(inter_proteins, hippie_pairs=None, outprefix=None):
    '''run MCL to find protein modules'''

    fh = open(hippie_pairs, "rU")
    fi = open(inter_proteins, "rU")

    predicted_inter_proteins = set()
    predicted_inter_pairs = []

    for line in fi.readlines():
        line = line.strip()
        predicted_inter_proteins.add(line)
    fi.close()

    ppi_out = outprefix + ".PPI"
    f_pairs = open(ppi_out, "w")

    for line in fh.readlines():
        line = line.strip()
        id1 = line.split()[0]
        id2 = line.split()[1]

        if id1 in predicted_inter_proteins and id2 in predicted_inter_proteins:
            f_pairs.write(line + "\n")
    fh.close()
    f_pairs.close()

    mcl_exe = os.path.join(exedir, "../src/mcl")
    mcl_out = outprefix + ".mcl"
    mcl_cmd = mcl_exe + " " + ppi_out + " --abc -o " + mcl_out
    print mcl_cmd
    subprocess.call(mcl_cmd, shell=True)

    graph_file_list = iGraphPairs(mcl_out, ppi_out, outprefix)

    for idx, graph_file in enumerate(graph_file_list):
        tmp_out = outprefix + ".function.module." + str(idx) + ".pdf"
        igraph_rscript = os.path.join(exedir, "../src/module.igraph.R")
        tmp_cmd = "Rscript " + igraph_rscript + " " + graph_file + " " + tmp_out
        subprocess.call(tmp_cmd, shell=True)


def iGraphPairs(mcl_out, ppi_pairs, outprefix):
    """generate protein modules for igraph"""
    fm = open(mcl_out, "rU")
    fp = open(ppi_pairs, "rU")

    ppi_pairs_list = fp.readlines()
    fp.close()
    output_list = []

    for idx, line in enumerate(fm.readlines()):
        # only keep the top 10 modules
        if idx > 9:
            break

        modules = set(line.split())
        output = outprefix + ".igraph.ppi." + str(idx)
        output_list.append(output)
        fo = open(output, "w")

        for i in range(len(ppi_pairs_list)):
            p_line = ppi_pairs_list[i].strip()
            print p_line
            id1 = p_line.split()[0]
            id2 = p_line.split()[1]

            if (id1 in modules) and (id2 in modules) and (id1 != id2):
                fo.write(p_line + "\n")
        fo.close()

    fm.close()

    return output_list

###############################################################
