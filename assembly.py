#! /analysis/xin/parasite/bin/python/bin/python

import argparse
import sys
import os
from utilities import command
from utilities import fileutils
from utilities import properties


def get_args():

    global properties_file
    global fastq1
    global fastq2
    global assem_sw
    global prefix
    global fi
    default_assem_sw='spades'
    fi=fileutils()
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-s', '--assembly_software', type=str, help='''Please provide the assembly software: 
			spades or it will be used by default''', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
        required=True)  

    # check args
    args = parser.parse_args()
    fi.check_files_exist([args.properties_file,args.fastq1]) 
    if args.fastq2 is not None:
        fi.check_exist(args.fastq2)

    # define variables
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    assem_sw=args.assembly_software    
    if assem_sw is None:
        assem_sw=default_assem_sw
    else:
        if assem_sw not in ['spades']:
            print "ERROR: assembly software can ONLY be spades!"
            parser.print_help()
            sys.exit(1)

    prefix=args.prefix    
    
    print "properties_file:",properties_file
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "assembly_software:",assem_sw
    print "prefix:",prefix


def run_assem_sw(fastqfiles):            
    #check whether qc_software path was defined in the property file
    assem_sw_path=prop.get_attrib( assem_sw) 
    #define command     
    if len(fastqfiles)==2:
        command(assem_sw_path+" -1 {} -2 {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],fastqfiles[1],workdir)).run_comm(0)
    elif len(fastqfiles)==1:
        command(assem_sw_path+" -s {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],workdir)).run_comm(0)


def initiate():
    print "initiating..."
    global indir
    global outdir
    global workdir

    subdir="assembly"
    workdir=prop.workdir+"/"+subdir
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_to_destdir(fastq1,indir)  
    if fastq2 is not None:
        fi.copy_file_to_destdir(fastq2,indir)
    

def execute():
    print "executing..."
    fi.change_dir(workdir)
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_assem_sw(fastqfiles)


def post_process():
    print "post_processing..."
    fi.copy_file_add_prefix(workdir+"/scaffolds.fasta",outdir,prefix+"_")
    fi.copy_file_add_prefix(workdir+"/spades.log",outdir,prefix+"_")

    
    
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
