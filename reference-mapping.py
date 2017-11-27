#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

import argparse
import sys
import os
from utilities import command
from utilities import fileutils
from utilities import properties


def get_args():
    
    global fi
    global properties_file
    global genome_name
    global ref_fasta
    global fastq1
    global fastq2
    global prefix
    global ref_fasta_root_default
    
    fi=fileutils()
    ref_fasta_root_default='ena_ref_fasta'
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script mapping short reads to reference genomes using BWA')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name you are mapping to, 
                        only \'ch\' or \'cp\' permitted.
                        \'ch\' stands for \'Cryptosporidium hominis\' and \'cp\' stands for \'Cryptosporidium parvum\'''',
                        required=True)
    parser.add_argument('-r', '--ref_fasta', type=str, help='Please provide reference genome fasta file', required=False)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', required=True)  
    
    args=parser.parse_args()
    # check args
    fi.check_files_exist([args.properties_file, args.fastq1])   
    if args.genome_name!='ch' and genome_name!='cp':
        print "ERROR: Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted"
        sys.exit(1)    
    if args.ref_fasta is not None:
        fi.check_file_exist(args.ref_fasta)
        
    # define variables        
    properties_file=args.properties_file
    genome_name=args.genome_name 
    if args.ref_fasta is None:
        ref_fasta=ref_fasta_root_default+'_'+genome_name
    fastq1=args.fastq1
    fastq2=args.fastq2  
    prefix=args.prefix    
        
    print "properties_file:",properties_file  
    print "genome_name:",genome_name
    print "refrence_genome fasta file:",ref_fasta
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "prefix:",prefix
  
def run_comm(comm):
    #run command       
    print comm
    returncode, stdout, stderr=command(comm).run(360000)
    if returncode!=0:
        print "\nERROR: {} FAILED!!! \n\nSTDERR: {}\nSTDOUT: {}\n".format(comm,stderr,stdout)    
        sys.exit(1) 
        
def run_mapping(fastqfiles):  
    global ref_fasta          
    fi.change_dir(workdir)     
    bwa_path=getattr(prop,"bwa")
    if ref_fasta == "{}_{}".format(ref_fasta_root_default,genome_name):        
        ref_fasta=getattr(prop,ref_fasta)
    
    #define command         
    fq1=fastqfiles[0]
    fq1_sai=prefix+'.fq1.sai'
    run_comm("{} aln -t 12 {} {} >{}".format(bwa_path,ref_fasta,fq1,fq1_sai))
    if len(fastqfiles)==2:
        fq2=fastqfiles[1]
        fq2_sai=prefix+'.fq2.sai'        
        run_comm("{} aln -t 12 {} {} >{}".format(bwa_path,ref_fasta,fq2,fq2_sai))
        run_comm("{} sampe {} {} {} {} {}>{}.sam".format(bwa_path,ref_fasta,fq1_sai,fq2_sai,fq1,fq2,prefix))
    elif len(fastqfiles)==1:
        run_comm("{} samse {} {} {}>{}.sam".format(bwa_path,ref_fasta,fq1_sai,fq1,prefix))

def initiate():
    print "initiating..."
    global indir
    global outdir
    global fi
    global reads_type     
    global workdir
    global subdir
    
    subdir="reference_mapping"
    workdir=prop.workdir+"/"+subdir
    fi=fileutils()
   
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_into_dest(fastq1,indir)  
    if fastq2 is not None:
        fi.copy_file_into_dest(fastq2,indir)

def execute():
    print "executing..."
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_mapping(fastqfiles)


def post_process():
    print "post_processing..."
    fi.copy_file_into_dest("{}/{}.sam".format(workdir,prefix),outdir)

    
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
