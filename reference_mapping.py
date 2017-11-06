#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

import argparse
import sys
import os
from utilities import command
from utilities import fileutils
from utilities import properties


def get_args():

    global properties_file
    global mapping_software
    global genome_name
    global ref_fasta
    global fastq1
    global fastq2
    global prefix
    global ref_fasta_root_default
    
    ref_fasta_root_default='ena_ref_fasta'
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script mapping short reads to reference genomes using bwa')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted, 
                                                                \'ch\' stands for \'Cryptosporidium% hominis\' and 
                                                                \'cp\' stands for \'Cryptosporidium% parvum\', required=True''', required=True)
    parser.add_argument('-r', '--ref_fasta', type=str, help='Please provide reference genome fasta file', required=False)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
        required=True)  
    
    mapping_software='bwa'
    args=parser.parse_args()
    properties_file=args.properties_file
    genome_name=args.genome_name 
    ref_fasta=args.ref_fasta
    if ref_fasta is None:
        ref_fasta=ref_fasta_root_default+'_'+genome_name
    fastq1=args.fastq1
    fastq2=args.fastq2  
    prefix=args.prefix    
    
    if genome_name!='ch' and genome_name!='cp':
        print "ERROR: Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted"
        sys.exit(1)
        
    print "properties_file:",properties_file  
    print "mapping_software:",mapping_software
    print "genome_name:",genome_name
    print "refrence_genome fasta file:",ref_fasta
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "prefix:",prefix


def get_property(key):
    try:
        return getattr(prop, key) 
    except AttributeError:
        print "ERROR: {} path was not defined in {}. Please define it and re-run.".format(key,properties_file)        
        sys.exit(1)

    
def run_mapping(fastqfiles):            
    fi.change_dir(workdir)     
    bwa_path=get_property(mapping_software)
    ref_fasta_path=get_property(ref_fasta)
        
    #define command         
    fq1=fastqfiles[0]
    fq1_sai=out_prefix+reads_type+'.fq1.sai'
    if len(fastqfiles)==2:
        fq2=fastqfiles[1]
        fq2_sai=out_prefix+reads_type+'.fq2.sai'
        comm="{} aln {} {} >{}; ".format(bwa_path,ref_fasta_path,fq1,fq1_sai)
        comm+=("{} aln {} {} >{}; ".format(bwa_path,ref_fasta_path,fq2,fq2_sai))
        comm+=("{} sampe {} {} {} {} {}>{}{}.sam".format(bwa_path,ref_fasta_path,fq1_sai,fq2_sai,fq1,fq2,out_prefix,reads_type))
    elif len(fastqfiles)==1:
        comm="{} aln {} {}>{}; ".format(bwa_path,ref_fasta_path,fq1,fq1_sai)
        comm+=("{} samse {} {} {}>{}{}.sam".format(bwa_path,ref_fasta_path,fq1_sai,fq1,out_prefix,reads_type)) 
        
    #run command  
    returncode, stdout, stderr=command(comm).run(3600)
    if returncode!=0:
        print "\nERROR: {} failed!!! \n\nSTDERR: {}\nSTDOUT: {}\n".format(comm,stderr,stdout)    
        sys.exit(1) 
    print returncode, stdout, stderr
    print prop.workdir


def initiate():
    print "initiating..."
    global indir
    global outdir
    global fi
    global out_prefix
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
    reads_type='single'
    if fastq2 is not None:
        reads_type='paired'
        fi.copy_file_into_dest(fastq2,indir)
    out_prefix="{}_{}_{}_".format(prefix.upper(),mapping_software,genome_name)
    

def execute():
    print "executing..."
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_mapping(fastqfiles)


def post_process():
    print "post_processing..."
    fi.copy_file_into_dest("{}/{}{}.sam".format(workdir,out_prefix,reads_type),outdir)

    
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
