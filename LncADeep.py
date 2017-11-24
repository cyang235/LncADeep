#!/usr/bin/env python

import argparse
import sys,os
import LncADeep_lncRNA.LncADeep_partial.bin.lncRNA_Predict as lncRNA_Predict_partial
import LncADeep_lncRNA.LncADeep_full.bin.lncRNA_Predict as lncRNA_Predict_full

from LncADeep_anno.bin.lncRNA_func_model_release import Predict as lncRNA_func_Predict
from LncADeep_anno.bin.lncRNA_func_model_release import Annotate as lncRNA_func_Annotate


def main():
    parser = argparse.ArgumentParser(usage='%(prog)s [options]', description='An ab initio lncRNA identification and functional annotation tool based on deep learning')

    parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')

    parser.add_argument('-MODE', '--MODE', action='store', dest='mode', choices=['lncRNA', 'anno'], type=str, default="lncRNA",
                        help='(Required) The mode used for lncRNA identification or functional annotation. If "lncRNA" is chosen, LncADeep will identify lncRNAs. If "anno" is chosen, LncADeep will predict lncRNA-protein interactions and annotate lncRNA functions. Default is "lncRNA"')

    parser.add_argument('-o', '--out', action='store', dest='out_prefix', 
                        help='(Required) The output prefix of results')


    ## lncRNA identification options
    parser.add_argument('-f', '--fasta', action='store', dest='fasta_file', 
                        help='(Required for lncRNA identification) Sequence file in FASTA format to be predicted')

    parser.add_argument('-m', '--model', action='store', dest='model', choices=['full', 'partial'], type=str, default="partial",
                        help='(Optional for lncRNA identification) The model used for lncRNA identification, default is "partial"')

    parser.add_argument('-s', '--species', action='store', dest='species', choices=['human', 'mouse'], type=str, default="human",
                        help='(Optional for lncRNA identification) The species used for lncRNA identification, default is "human"')

    parser.add_argument('-th', '--thread', action='store', dest='thread', default=1, type=int,
                        help='(Optional for lncRNA identification) Use multi-thread for predicting, default is 1')

    parser.add_argument('-HMM', '--HMMthread', action='store', dest='HMMthread', default=8, type=int,
                        help='(Optional for lncRNA identification) The thread number of using HMMER, default is 8')


    ## functional annotation options
    parser.add_argument('-l', '--lncRNA', action='store', dest='rna_file', help='(Required for functional annotation) lncRNA sequence file in FASTA format')

    parser.add_argument('-p', '--protein', action='store', dest='protein_file', help='(Optional for functional annotation) protein sequence file in FASTA format')

    parser.add_argument('-a', '--annotation', action='store', dest='annotation', choices=[1, 0], default=1, type=int,
        help='(Optional for functional annotation) To annotate lncRNA functions. If "1" is selected, LncADeep will annotate the functions for lncRNAs, otherwise LncADeep will only give the interacting proteins for lncRNAs. The default is "1".')

    parser.add_argument('-r', '--pair', action='store', dest='pair_file', default=None, 
            help = "(Optional for functional annotation) The lncRNA-protein pairs to be predicted. If this option is selected, LncADeep will only output interacting proteins.")


    args = parser.parse_args()

    if args.mode == "lncRNA":

        if args.model == "partial":
            if args.fasta_file and args.out_prefix:
                lncRNA_Predict_partial.predict(filename=args.fasta_file, output_prefix=args.out_prefix, species=args.species, thread=args.thread, HMMthread=args.HMMthread)
            else:
                parser.parse_args(['-h'])
        elif args.model == "full":
            if args.fasta_file and args.out_prefix:
                lncRNA_Predict_full.predict(filename=args.fasta_file, output_prefix=args.out_prefix, species=args.species, thread=args.thread, HMMthread=args.HMMthread)
            else:
                parser.parse_args(['-h'])

    elif args.mode == "anno":

        if args.pair_file is not None:
            args.annotation = 0

        if args.annotation == 1:
            if args.rna_file and args.out_prefix:
                lncRNA_func_Annotate(test_pair=args.pair_file, rna_file=args.rna_file, out_prefix=args.out_prefix)
            else:
                parser.parse_args(['-h'])
        else: 
            if args.rna_file and args.out_prefix:
                lncRNA_func_Predict(test_pair=args.pair_file, rna_file=args.rna_file, pro_file=args.protein_file, 
                                    out_prefix=args.out_prefix)
            else:
                parser.parse_args(['-h'])

    else:
        parser.parse_args(['-h'])



if __name__ == '__main__':
    main()
