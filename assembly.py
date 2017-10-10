#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

import argparse
import sys
from utilities import command
from utilities import fileutils
from utilities import properties


def get_args():

    global properties_file
    global fastq1
    global fastq2
    global assem_sw
    global prefix
    default_assem_sw='spades'
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-s', '--assembly_software', type=str, help='''Please provide the assembly software: spades or velvet, 
        otherwise spades will be used''', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
        required=True)  
    
    args = parser.parse_args()
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    assem_sw=args.assembly_software    
    if assem_sw is None:
        assem_sw=default_assem_sw
    else:
        if assem_sw not in ['spades', 'velvet']:
            print "Please choose assembly software ONLY from spades and velvet!"
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
    try:
        assem_sw_path=getattr(prop, assem_sw) 
    except AttributeError:
        print "ERROR: {} path was not defined in {}. Please define it and re-run.".format(assem_sw,properties_file)        
        sys.exit(1)        
    #define command     
    if len(fastqfiles)==2:
        comm=assem_sw_path+" -1 {} -2 {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],fastqfiles[1],workdir)   
    elif len(fastqfiles)==1:
        comm=assem_sw_path+" -s {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],workdir)
    else:
        print "ERROR: wrong number of fastq files " + fastqfiles.size()
    #run command  
    returncode, stdout, stderr=command(comm).run(3600)
    if returncode!=0:
        print "ERROR: {} failed.\n stderr: {}\n stdout: {}".format(comm,stderr,stdout)    
        sys.exit(1) 
    print returncode, stdout, stderr
    print prop.workdir
    

def initiate():
    print "initiating..."
    global indir
    global outdir
    global fi
    global out_prefix
    global run_type     
    global workdir
    
    workdir=prop.workdir
    fi=fileutils()
    indir=workdir+"/assembly/in/"
    outdir=workdir+"/assembly/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_into_dest(fastq1,indir)  
    run_type='single'
    if fastq2 is not None:
        run_type='paired'
        fi.copy_file_into_dest(fastq2,indir)
    out_prefix="{}_{}_{}_".format(prefix.upper(),assem_sw.lower(),run_type)
    

def execute():
    print "executing..."
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_assem_sw(fastqfiles)


def post_process():
    print "post_processing..."
    fi.copy_file_add_prefix(workdir+"/scaffolds.fasta",outdir,out_prefix)
    fi.copy_file_add_prefix(workdir+"/spades.log",outdir,out_prefix)

    
    
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
   
     
