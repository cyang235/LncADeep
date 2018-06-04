LncADeep: An ab initio lncRNA identification and functional annotation tool based on deep learning
===================================================================


* [Introduction](#introduction)
* [Requirements](#requirements)
    * [LncADeep.py](#lncadeeppy)
* [Version](#version)
* [Installation](#installation)
    * [Prerequisites](#prerequisites)
    * [Install LncADeep from zipped file](#install-lncadeep-from-zipped-file)
    * [Install LncADeep using git](#install-lncadeep-using-git)
* [Usage](#usage)
    * [LncADeep.py](#lncadeeppy-1)
* [Examples](#examples)
    * [For lncRNA identification.](#for-lncrna-identification)
    * [For lncRNA functional annotaion.](#for-lncrna-functional-annotaion)
* [Citation](#citation)
* [Contact](#contact)




Introduction
------------

Long noncoding RNAs (lncRNAs) play important biological roles and have been implicated in human diseases. To characterize lncRNAs, identifying and annotating lncRNAs is necessary. Here, we propose a novel lncRNA identification and functional annotation tool named **LncADeep**. First, **LncADeep** identifies lncRNAs by integrating sequence intrinsic and homology features based on deep belief networks. Second, **LncADeep** predicts lncRNA-protein interactions using sequence and structure features based on deep neural networks. Third, since accurate lncRNA-protein interactions can help to infer the functions of lncRNAs, **LncADeep** conducts KEGG and Reactome pathway enrichment analysis and functional module detection with the predicted interacting proteins of lncRNAs. Case studies show that LncADeep's annotations for lncRNAs comply with their known functions. As a tool for lncRNA identification and functional annotation based on deep learning, **LncADeep** has outperformed state-of-the-art tools on predicting lncRNAs and lncRNA-protein interactions, and can automatically provide informative functional annotations for lncRNAs.

**LncADeep** is freely available for non-commercial use at [http://cqb.pku.edu.cn/ZhuLab/lncadeep](http://cqb.pku.edu.cn/ZhuLab/lncadeep) or [https://github.com/cyang235/LncADeep](https://github.com/cyang235/LncADeep).


Requirements
-------------

### LncADeep.py

For **lncRNA identification**.

+ [Python version >= 2.7.6](https://www.python.org/)
+ [numpy version >= 1.8.2](http://www.numpy.org/)
+ [HMMER package 3.1b2](http://hmmer.org/download.html)


For predicting **lncRNA-protein interactions** and annotating **lncRNA functions**.

+ [Python, version >= 2.7.6](https://www.python.org/)
+ [numpy, version >= 1.8.2](http://www.numpy.org/)
+ [pandas, version >= 0.18.0](http://pandas.pydata.org/)
+ [theano, version >= 0.8.2](https://github.com/Theano/Theano)
+ [Keras, version 1.2.2](https://github.com/fchollet/keras/releases)
+ [h5py, version >= 2.5.0](http://www.h5py.org/)
+ [R, version >= 3.3.2](https://www.r-project.org/)
+ [iGraph R package](http://igraph.org/r/)
+ [MCL package](https://micans.org/mcl/)


Version
-------
+ LncADeep 1.0 (Tested on Linux_64, including CentOS 6.5 and Ubuntu 16.04)


Installation
------------

### Prerequisites

Please install **numpy, theano, pandas, Keras, h5py**, and **iGraph** according to their manuals. The following are examples for installing these prerequisites.

**numpy, theano, pandas**, and **h5py** are python packages, which can be installed with ``pip``, for example:

    # we use python v2.7.13
    # our machine is implemented with the following versions

    pip install numpy       # numpy v1.13.1
    pip install Theano      # Theano v0.9.0
    pip install pandas      # pandas v0.20.3
    pip install h5py        # h5py v2.7.0


**iGraph** is an R package, which can be installed with:

    # we use R v3.3.2
    # Download and install the package

    install.packages("igraph")


**Keras** is a Python Deep Learning library.
For **Keras**, please be noted that we use [Keras v1.2.2](https://github.com/fchollet/keras/releases), and we use **theano** as its backend. Please edit the file ``~/.keras/keras.json`` and change the backend. For example,

    # Download Keras v1.2.2
    wget https://github.com/fchollet/keras/archive/1.2.2.tar.gz

    # unpack the zipped file
    tar xzvf 1.2.2.tar.gz

    # install Keras v1.2.2
    cd keras-1.2.2
    python setup.py install

    # edit ~/.keras/keras.json and change the backend
    {
        "image_dim_ordering": "th",
        "epsilon": 1e-07,
        "floatx": "float32",
        "backend": "theano"
    }


HMMER and MCL package have been **included** in LncADeep package. You **don't need** to install them yourself.
After installing the above prerequisites, you can now install **LncADeep**.

### Install LncADeep from zipped file 

download the zipped file 

    wget http://cqb.pku.edu.cn/ZhuLab/LncADeep/LncADeep_v1.0.tgz 
    
unpack the zipped file 

    tar xzvf LncADeep_v1.0.tgz
    
change directory to LncADeep

    cd LncADeep_v1.0

configure and add directory to the PATH, and you are done!

    chmod +x configure
    ./configure

    source $HOME/.bash_profile
    

### Install LncADeep using git

clone LncADeep package

    git clone https://github.com/cyang235/LncADeep.git
    
change directory to LncADeep

    cd LncADeep

download Pfam 29.0 database

    wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam29.0/Pfam-A.hmm.gz
    gzip -d Pfam-A.hmm.gz
    mv Pfam-A.hmm ./LncADeep_lncRNA/src/

    # the Pfam-A.hmm need be put in directory /path to LncADeep/LncADeep_lncRNA/src/
    
configure and add directory to the PATH, and you are done!

    chmod +x configure
    ./configure

    source $HOME/.bash_profile
    

Usage
-----

### LncADeep.py

An ab initio lncRNA identification and functional annotation tool based on deep learning

	usage: LncADeep.py [options]

	An ab initio lncRNA identification and functional annotation tool based on
	deep learning

	optional arguments:
	  -h, --help            show this help message and exit
	  -v, --version         show program's version number and exit
	  -MODE {lncRNA,anno}, --MODE {lncRNA,anno}
							(Required) The mode used for lncRNA identification or
							functional annotation. If "lncRNA" is chosen, LncADeep
							will identify lncRNAs. If "anno" is chosen, LncADeep
							will predict lncRNA-protein interactions and annotate
							lncRNA functions. Default is "lncRNA"
	  -o OUT_PREFIX, --out OUT_PREFIX
							(Required) The output prefix of results
	  -f FASTA_FILE, --fasta FASTA_FILE
							(Required for lncRNA identification) Sequence file in
							FASTA format to be predicted
	  -m {full,partial}, --model {full,partial}
							(Optional for lncRNA identification) The model used
							for lncRNA identification, default is "partial"
	  -s {human,mouse}, --species {human,mouse}
							(Optional for lncRNA identification) The species used
							for lncRNA identification, default is "human"
	  -th THREAD, --thread THREAD
							(Optional for lncRNA identification) Use multi-thread
							for predicting, default is 1
	  -HMM HMMTHREAD, --HMMthread HMMTHREAD
							(Optional for lncRNA identification) The thread number
							of using HMMER, default is 8
	  -l RNA_FILE, --lncRNA RNA_FILE
							(Required for functional annotation) lncRNA sequence
							file in FASTA format
	  -p PROTEIN_FILE, --protein PROTEIN_FILE
							(Optional for functional annotation) protein sequence
							file in FASTA format
	  -a {1,0}, --annotation {1,0}
							(Optional for functional annotation) To annotate
							lncRNA functions. If "1" is selected, LncADeep will
							annotate the functions for lncRNAs, otherwise LncADeep
							will only give the interacting proteins for lncRNAs.
							The default is "1".
	  -r PAIR_FILE, --pair PAIR_FILE
							(Optional for functional annotation) The lncRNA-
							protein pairs to be predicted. If this option is
							selected, LncADeep will only output interacting
							proteins.



Examples
--------

The files for example are stored at directory `/path to LncADeep/data`

### For lncRNA identification.

1. Identify lncRNAs using **model** for transcripts including **full- and partial-length**

        python LncADeep.py -MODE lncRNA -f ./data/LncADeep_lncRNA/lncRNA_mRNA_test.fa -o test
    
    The output files will be generated at directory `test_LncADeep_lncRNA_results`. 

2. To use **multi-processes** (e.g., 4) for lncRNA identification

        python LncADeep.py -MODE lncRNA -f ./data/LncADeep_lncRNA/lncRNA_mRNA_test.fa -o test \
                                -th 4

3. To use the **model** for **full-length** transcripts, please use the following command

        python LncADeep.py -MODE lncRNA -f ./data/LncADeep_lncRNA/lncRNA_mRNA_test.fa -o test \
                                -m full

4. LncADeep has been trained on the datasets of two species, including "human" and "mouse", the default model is "human". To use the **model** trained on **mouse full- and partial-length** transcripts, please use the following command. 

        python LncADeep.py -MODE lncRNA -f ./data/LncADeep_lncRNA/lncRNA_mRNA_test.fa -o test \
                                -s mouse

5. To use the **model** trained on **mouse full-length** transcripts, please use the following command. 

        python LncADeep.py -MODE lncRNA -f ./data/LncADeep_lncRNA/lncRNA_mRNA_test.fa -o test \
                                -m full -s mouse

6. LncADeep accepts nucleotide **FASTA sequence** as _input_, e.g.:

        >RNA_id_1
        GGAAACGGCCGTGGGCATTTTGGTGTATTTTTATTCAACTTTGAAAGACATATTTTATTTTTACACATTTTATTTTATACAGTA
        TAGACATACATATGCATACACGCCTCCTCTCATGACATTAAACTTTTGCACAACTTCACAATTGTAAATGATCACAGAAAAATG
        CCTCAAAATGAATGTATCATATCCTAGCCCCACCACTTAACCTCTCTGTGCCTCAGTTTTCTCCTCTGTAAAACGGGGATAATA
        ATAGTATCTACTTTATAAGTTGCTTGTAAGGGTTCAATGTGATTATGGTGTGAATGTGGGAAGCGCTCAGAAAGTATCATTTTC
        ATTATTATTAGAACTATTATTCCTTAATTGCAAACATTTAAATTCTAATTTTAT
        
        >RNA_id_2
        CATCTCTTTCCTTCTCAGGAAATTTTATACATTGTCAATTATTCCTTCTCTCTAACTTCAACCTCGCCTTCTTTGCTGAGTCTG
        ACCCATCAACAGTTAAACATGATCAAGTCTTCCGATTTAAAAGTCCCTCTTTCTTGACACAGCTCATTTATAGCCAAACTTCTT
        TCTGAAGAGTAGTCTACATTCATTTTCTTTTTCTCCCTCACTTCTGATAATATTGAACCAACTCCATTTTAGTTTCTGTCCCTA
        TCATTCCTCTAAATTGATTAAGGTCTCCAGAATATTCCTCTGTATTTACGGGCATTATTCACTGCTCTTCTTATTTGACTACTC
        AGCAAGCATTTAACTTTTGATCAGTTTTTCCTTAAAATACTTTACTTGGCTTCCTTGACATCATGGTTTTTGTTCAGATCTCTG
        TGGTTATTTCTGTCTCCTTTGCTGCCTTCTCCTCTTGGTCCTTG
        
        # LncADeep will predict whether `RNA_id_1` and `RNA_id_2` are lncRNAs. 
        # More example input can be found at directory `/path to LncADeep/data/LncADeep_lncRNA`     
     

### For lncRNA functional annotaion.

1. To predict lncRNA-protein **interactions** and **annotate the functions** of lncRNAs. 

        python LncADeep.py -MODE anno -l ./data/LncADeep_anno/ENST00000424518.5.fa -o test
    
        # Here, LncADeep will predict the interactions between given lncRNAs and 20,121 reviewed proteins
        # and then annotate the functions of lncRNAs with their predicted interacting proteins.
        # The output files will be generated at directory `test_LncADeep_anno_results`

2. To predict lncRNA-protein **interactions**.

        python LncADeep.py -MODE anno -l ./data/LncADeep_anno/ENST00000424518.5.fa -o test -a 0

        # Here, LncADeep will predict the interactions between given lncRNAs and 20,121 reviewed proteins

3. To predict lncRNA-protein **interactions** for **given pairs**.

        python LncADeep.py -MODE anno -l ./data/LncADeep_anno/ENST00000424518.5.fa -o test \
                    -r ./data/LncADeep_anno/pair.dat -p ./data/LncADeep_anno/protein.fa

        # Here, LncADeep will predict the interactions between given lncRNAs and proteins for given pairs
        # Users are required to provide the lncRNA and protein sequences in FASTA format 
        # and lncRNA-protein pairs in text format, see below

4. lncRNA-protein **pairs** in text format as input, e.g.:

        ENST00000424518.5|ENSG00000228630.5|OTTHUMG00000152934.1|OTTHUMT00000328662.1|HOTAIR-001|HOTAIR|2421| sp|P27361|MK03_HUMAN
        ENST00000424518.5|ENSG00000228630.5|OTTHUMG00000152934.1|OTTHUMT00000328662.1|HOTAIR-001|HOTAIR|2421| sp|P53779|MK10_HUMAN
        ENST00000424518.5|ENSG00000228630.5|OTTHUMG00000152934.1|OTTHUMT00000328662.1|HOTAIR-001|HOTAIR|2421| sp|Q15049|MLC1_HUMAN
        ENST00000424518.5|ENSG00000228630.5|OTTHUMG00000152934.1|OTTHUMT00000328662.1|HOTAIR-001|HOTAIR|2421| sp|Q9UHC1|MLH3_HUMAN
        ENST00000424518.5|ENSG00000228630.5|OTTHUMG00000152934.1|OTTHUMT00000328662.1|HOTAIR-001|HOTAIR|2421| sp|P0DMT0|MLN_HUMAN

        # LncADeep will predict the interactions for the above lncRNA-protein pairs. 
        # Users are also required to provide the lncRNA and protein FASTA sequence files.
        # More example input can be found at directory `/path to LncADeep/data/LncADeep_anno`





Citation
--------
[Cheng Yang, Longshu Yang, Man Zhou, Haoling Xie, Chengjiu Zhang, May D Wang, Huaiqiu Zhu; LncADeep: An ab initio lncRNA identification and functional annotation tool based on deep learning, Bioinformatics, 2018, bty428](https://doi.org/10.1093/bioinformatics/bty428)


Contact
-------
Please direct your questions to: Dr. Huaiqiu Zhu, [hqzhu@pku.edu.cn](hqzhu@pku.edu.cn)
 
