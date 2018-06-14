#! /analysis/xin/parasite/bin/python/bin/python 

import argparse
import sys
from utilities import command
from utilities import fileutils
from utilities import properties
import os
import re

def get_args():    
    global properties_file
    global fastq1
    global fastq2
    global qc_sw
    global rm_dup_sw
    global prefix
    global fi
    global if_dedup
  
    fi=fileutils()
    default_qc_sw="trim_galore"
    rm_dup_sw="clumpify"
 
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file.', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the second fastq file.', required=False)
    parser.add_argument('-qc_sw', '--qc_software', type=str, help='''Please provide the quality control software, 
        otherwise the default trim_galore will be used''', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
        required=True)  
    parser.add_argument('-de', '--dedup', action='store_true',help='if remove all the exact read duplications', default=False)
    # check args
    args = parser.parse_args()
    fi.check_files_exist([args.properties_file, args.fastq1])
    if args.fastq2 is not None:
	fi.check_exist(args.fastq2)     

    # define variables
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    qc_sw=args.qc_software    
    if qc_sw is None:
        qc_sw = default_qc_sw
    prefix=args.prefix
    if_dedup=args.dedup
    
    print "properties_file:",properties_file
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "qc_software:",qc_sw
    print "prefix:",prefix
    print "dedup:",if_dedup
   

def run_trim_galore(fastqfiles):            
    global fastq1_name
    global fastq1_name_base
    global fastq2_name
    global fastq2_name_base
    global fqout1
    global fqout2
    fastq1_name=os.path.basename(fastq1)
    fastq1_name_base=fastq1_name.rstrip(".fastq")
    if fastq2 is None:
        fqout1=os.path.join(workdir,fastq1_name_base+"_trimmed.fq")
    else:
        fqout1=os.path.join(workdir,fastq1_name_base+"_val_1.fq")
        fastq2_name=os.path.basename(fastq2)
        fastq2_name_base=fastq2_name.rstrip(".fastq")
        fqout2=os.path.join(workdir,fastq2_name_base+"_val_2.fq")

    #check whether qc_software path was defined in the property file
    qc_sw_path=prop.get_attrib(qc_sw) 
    
    #pair fastq files    
    if len(fastqfiles)==2:
        command(qc_sw_path+" --paired -q 20 "+fastqfiles[0]+" "+fastqfiles[1]).run_comm(0)
    #single fastq files    
    elif len(fastqfiles)==1:
        command(qc_sw_path+" -q 20 "+fastqfiles[0]).run_comm(0)

    
def deduplication():
    global fq_dedup_out1
    global fq_dedup_out2
    rm_dup_sw_path=prop.get_attrib(rm_dup_sw)
    fq_dedup_out1=fqout1.replace(".fq",".dedup.fq")
    command("{} in={} out={} dedupe=t".format(rm_dup_sw_path,fqout1,fq_dedup_out1)).run_comm(0)
    
    if fastq2 is not None:
        fq_dedup_out2=fqout2.replace(".fq",".dedup.fq")    
        command("{} in={} out={} dedupe=t".format(rm_dup_sw_path,fqout2,fq_dedup_out2)).run_comm(0)

def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global type
        
    subdir="quality"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_to_destdir(fastq1,indir)
    type='single'    
    if fastq2 is not None:
        type='paired'
        fi.copy_file_to_destdir(fastq2,indir) 
        
def execute():
    print "executing..."
    fi.change_dir(workdir)
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_trim_galore(fastqfiles)
    if if_dedup:
        print("removing exact duplicates...")  
        deduplication()


def post_process():
    print "post_processing..."
    if fastq2 is not None:
        report2=os.path.join(workdir,fastq2_name+"_trimming_report.txt")
        fi.copy_file_add_prefix(fqout2,outdir,prefix+"_")  
        fi.copy_file_add_prefix(report2,outdir,prefix+"_") 
        if if_dedup :
           fi.copy_file_add_prefix(fq_dedup_out2,outdir,prefix+"_")

    report1=os.path.join(workdir,fastq1_name+"_trimming_report.txt")
    fi.copy_file_add_prefix(fqout1,outdir,prefix+"_")
    fi.copy_file_add_prefix(report1,outdir,prefix+"_")
    if if_dedup :
	fi.copy_file_add_prefix(fq_dedup_out1,outdir,prefix+"_")


if __name__ == '__main__':
    get_args()
    global prop
    prop=properties(properties_file)
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print os.path.realpath(__file__)+" DONE"
     
