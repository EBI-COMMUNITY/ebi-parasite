#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

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
    global prefix
    default_qc_sw="trim_galore"
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file.', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the second fastq file.', required=False)
    parser.add_argument('-qc_sw', '--qc_software', type=str, help='''Please provide the quality control software, 
        otherwise the default trim_galore will be used''', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
        required=True)  
    
    args = parser.parse_args()
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    qc_sw=args.qc_software    
    if qc_sw is None:
        qc_sw = default_qc_sw
    prefix=args.prefix
   
    
    print "properties_file:",properties_file
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "qc_software:",qc_sw
    print "prefix:",prefix
   

def run_trim_galore(fastqfiles):            
    #check whether qc_software path was defined in the property file
    try:
        qc_sw_path=getattr(prop, qc_sw) 
    except AttributeError:
        print "ERROR: {} path was not defined in {}. Please define it and re-run.".format(qc_sw,properties_file)        
        sys.exit(1)        
    #pair fastq files    
    if len(fastqfiles)==2:
        comm=qc_sw_path+" --paired -q 20 "+fastqfiles[0]+" "+fastqfiles[1]
    #single fastq files    
    elif len(fastqfiles)==1:
        comm=qc_sw_path+" -q 20 "+fastqfiles[0]

    comm_obj=command(comm)
    returncode, stdout, stderr=comm_obj.run(3600)
    print returncode, stdout, stderr
    print prop.workdir
   
    
def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global fi
    global type
    global out_prefix
        
    fi=fileutils()
    workdir=prop.workdir
    indir=workdir+"/quality/in/"
    outdir=workdir+"/quality/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_into_dest(fastq1,indir)
    type='single'    
    if fastq2 is not None:
        type='paired'
        fi.copy_file_into_dest(fastq2,indir) 
    out_prefix="{}_{}_{}_".format(prefix.upper(),qc_sw.lower(),type)   

        
def execute():
    print "executing..."
    os.chdir(workdir)
    fastqfiles=[]
    fastaq_dir=workdir+"/quality/in/"
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_trim_galore(fastqfiles)


def post_process():
    print "post_processing..."
    fastq1_name=os.path.basename(fastq1) 
    
    if fastq2 is None: 
        fqout1=os.path.join(workdir,fastq1_name+"_trimmed.fq")       
    else:        
        fqout1=os.path.join(workdir,fastq1_name+"_val_1.fq")
        fastq2_name=os.path.basename(fastq2)
        fqout2=os.path.join(workdir,fastq2_name+"_val_2.fq")
        report2=os.path.join(workdir,fastq2_name+"_trimming_report.txt")
        fi.copy_file_add_prefix(fqout2,outdir,out_prefix)  
        fi.copy_file_add_prefix(report2,outdir,out_prefix) 
   
    report1=os.path.join(workdir,fastq1_name+"_trimming_report.txt")
    fi.copy_file_add_prefix(fqout1,outdir,out_prefix)
    fi.copy_file_add_prefix(report1,outdir,out_prefix)


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
   
     
