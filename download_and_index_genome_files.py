#! /usr/bin/env python3

## version 1, 29 Nov 2019 by Xin Liu

import argparse
import sys
import os
import re
import time
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc

def get_args():
    global fi
    global prop
    global misc
    global properties_file
    global genome_name
    fi=fileutils()

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script downloading genome files')
    parser.add_argument('-p', '--properties_file', type=str, help='''Please provide the properties file, 
                               which including workdir''', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                               which is provided by genome_list.txt''', required=True)
    args=parser.parse_args()
    # check args
    fi.check_exist(args.properties_file)   
    properties_file=args.properties_file
    prop=properties(properties_file)
    misc=misc()
    misc.check_genome_avl(prop.get_attrib("available_genomes"), args.genome_name)
    # define variables        
    genome_name = args.genome_name
    print ("properties_file:",properties_file)
    print ("genome_name:", genome_name)

def dl_index_genome_files():
    ## ref_fasta_download 
    misc.download(workdir, "fasta", genome_name, "dna")
    ## gff3 download
    misc.download(workdir, "gff3", genome_name, "")
    time.sleep(3)
    REF_FASTA = "{}/{}.fasta".format(workdir, genome_name)
    command("samtools faidx {}".format(REF_FASTA)).run_comm(0)
    command("bwa index {}".format(REF_FASTA)).run_comm(0)
    ref_fasta_key=re.findall(".*/(.*?).fa",REF_FASTA)[0]
    command("bowtie2-build {} {}".format(REF_FASTA,ref_fasta_key)).run_comm(0)
    ref_dict_fpath=re.findall("(.*?).fa",REF_FASTA)[0]+".dict"
    if not os.path.isfile(ref_dict_fpath):
        command("{} CreateSequenceDictionary REFERENCE={} OUTPUT={}"
                 .format(picard,REF_FASTA,re.findall("(.*?).fa",REF_FASTA)[0]+".dict")).run_comm(0)

def initiate():
    print ("initiating...")
    global workdir
    subdir="ref"
    workdir=prop.workdir+"/"+subdir
    fi.create_processing_dir(workdir)

def execute():
    print ("executing...")
    fi.change_dir(workdir)
    dl_index_genome_files()

if __name__ == '__main__':
    get_args()
    print ("\n","Properties attributes:")
    print (prop.__dict__)
    getVar = lambda searchList, ind: [searchList[i] for i in ind]   
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    print (os.path.realpath(__file__)+" DONE")
