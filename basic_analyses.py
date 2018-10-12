#! /usr/bin/python

import argparse
import re
import sys
import time
import os
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict
from Bio import SeqIO

fi=fileutils()
misc=misc()

def check_fastq_postfix(fastq_f):
    if not re.search(".fq",fastq_f) and not re.search(".fastq",fastq_f):
        misc.my_exit("ERROR: fastq files postfix. Please use the common fastq file postfix '.fq' or '.fastq' for the fastq files.")       
    else:
        if re.search(".fq",fastq_f):
            print re.findall("^.*/(.*?).(fq)$",fastq_f)[0]
            return re.findall("^.*/(.*?).(fq)$",fastq_f)[0]
        else:
            print re.findall("^.*/(.*?).(fastq)$",fastq_f)[0]
            return re.findall("^.*/(.*?).(fastq)$",fastq_f)[0]

def get_args():    
    global properties_file
    global g_name_str
    global prop
    global min_homo
    global in_fastq1
    global in_fastq2
    global fastq1_key
    global fastq2_key
    global fastq1_postfix
    global fastq2_postfix
    global qc_sw
    global if_dedupQ
    global if_dedupM
    global prefix

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script provides basic analyses for the next generation sequences, including fastq file quality control, 
                                                    assembly by using spades, reference_mapping, and SNP calling using samtools 
                                                    for the genome available in the genome_list''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-pre', '--prefix', type=str,help='Please provide the prefix for the output file.',
                        required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name which is provided by genome_list.txt''',
                        required=True)
    parser.add_argument('-qc_sw', '--qc_software', type=str, help='''Please provide the quality control software, 
                        otherwise the default trim_galore will be used''', required=False)
    parser.add_argument('-deQ', '--dedupQ', action='store_true',help='if remove all the exact read duplications', default=False)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-deM', '--dedupM', action='store_true',help='if true, remove duplications after mapping using samtools', default=False)

    # check args
    args = parser.parse_args()
    fi.check_exist(args.properties_file)
    properties_file=args.properties_file
    prop=properties(properties_file)
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()):
        misc.my_exit("{} is not available, please try another genome".format(args.genome_name))     
    fi.check_exist(args.fastq1)
    (fastq1_key,fastq1_postfix)=check_fastq_postfix(args.fastq1)
    if args.fastq2 is not None:
        fi.check_exist(args.fastq2)
        (fastq2_key,fastq2_postfix)=check_fastq_postfix(args.fastq2)

    # define variables     
    properties_file=args.properties_file    
    prop=properties(properties_file)
    prefix=args.prefix
    g_name_str=args.genome_name
    in_fastq1=args.fastq1
    in_fastq2=args.fastq2
    qc_sw=args.qc_software
    if qc_sw is None:
        qc_sw="trim_galore"
    if_dedupQ=args.dedupQ   
    if_dedupM=args.dedupM

    print "properties_file:",properties_file
    print "prefix:",prefix
    print "gname:",g_name_str
    print "fastq1:",in_fastq1
    print "fastq2:",in_fastq2
    print "qc_software:",qc_sw
    print "assembly_software:spades"
    print "dedupQ:",if_dedupQ
    print "dedupM:",if_dedupM

def run_analysis():
    quality_control=prop.get_attrib("quality_control")
    assembly=prop.get_attrib("assembly")
    reference_mapping=prop.get_attrib("reference_mapping")
    snp_samtools=prop.get_attrib("snp_samtools")
    ### quality_control
    subdir="quality"
    fastq_dir=prop.workdir+"/quality/out"
    if if_dedupQ:
        if_dedupQ_str=" -de"
        if_dedupQ_postfix=".dedup"
    else:
        if_dedupQ_str=""
        if_dedupQ_postfix=""
    if in_fastq2 is not None:
        command("{} -p {} -fq1 {} -fq2 {} -pre {}{}".format(quality_control,properties_file,in_fastq1,in_fastq2,prefix,if_dedupQ_str)).run_comm(0)
    else:
        command("{} -p {} -fq1 {} -pre {}{}".format(quality_control,properties_file,in_fastq1,prefix,if_dedupQ_str)).run_comm(0)        
    ### assembly
    if in_fastq2 is not None:
        fastq1="{}/{}_{}_val_1{}.fq".format(fastq_dir,prefix,fastq1_key,if_dedupQ_postfix)
        fastq2="{}/{}_{}_val_2{}.fq".format(fastq_dir,prefix,fastq2_key,if_dedupQ_postfix)
        command("{} -p {} -fq1 {} -fq2 {} -pre {}".format(assembly,properties_file,fastq1,fastq2,prefix)).run_comm(0)
    else:
        fastq1="{}/{}_{}_trimmed{}.fq".format(fastq_dir,prefix,fastq1_key,if_dedupQ_postfix)
        command("{} -p {} -fq1 {} -pre {}".format(assembly,properties_file,fastq1,prefix)).run_comm(0)
    ### reference_mapping
    if if_dedupM:
        if_dedupM_str=" -de"
        if_dedupM_postfix=".dedup"
    else:
        if_dedupM_str=""
        if_dedupM_postfix=""
    if in_fastq2 is not None:
        command("{} -p {} -r {} -fq1 {} -fq2 {} -pre {}{}".format(reference_mapping,properties_file,fasta_fpath,
                fastq1,fastq2,prefix,if_dedupM_str)).run_comm(0)
    else:
        command("{} -p {} -r {} -fq1 {} -pre {}{}".format(reference_mapping,properties_file,fasta_fpath,fastq1,prefix,if_dedupM_str)).run_comm(0)        
    ### snp_samtools
    bam_dir=prop.workdir+"/reference_mapping/out"
    command("{} -p {} -r {} -bam {}/{}{}.bam -pre {}".format(snp_samtools,properties_file,fasta_fpath,bam_dir,
            prefix,if_dedupM_postfix,prefix)).run_comm(0)
        
def execute():
    print "executing..."
    global fasta_fpath
    fasta_fpath=os.getcwd()+"/"+misc.download("fasta",g_name_str,"dna")
    run_analysis()


if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #execute the main part of the program
    execute()
