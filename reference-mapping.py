#! /analysis/xin/parasite/bin/python/bin/python

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
    global rm_dup_sw
    global if_dedup
    
    fi=fileutils()
    ref_fasta_root_default='ena_ref_fasta'
    rm_dup_sw="samtools"

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script mapping short reads to reference genomes using BWA')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file, which including BWA path and workdir', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name you are mapping to, 
                        only \'ch\' or \'cp\' permitted.
                        \'ch\' stands for \'Cryptosporidium hominis\' and \'cp\' stands for \'Cryptosporidium parvum\'''',
                        required=True)
    parser.add_argument('-r', '--ref_fasta', type=str, help='''Please provide reference genome fasta file, 
                        otherwise please provide the path of ena fasta file in the property file''', required=False)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', required=True)  
    parser.add_argument('-de', '--dedup', action='store_true',help='if true, remove duplications after mapping using samtools', default=False)
   
    args=parser.parse_args()
    # check args
    fi.check_files_exist([args.properties_file, args.fastq1])   
    if args.genome_name!='ch' and args.genome_name!='cp':
        print "ERROR: Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted"
        sys.exit(1)    
    if args.ref_fasta is not None:
        fi.check_exist(args.ref_fasta)
        
    # define variables        
    properties_file=args.properties_file
    genome_name=args.genome_name 
    if args.ref_fasta is None:
        ref_fasta=ref_fasta_root_default+'_'+genome_name
    fastq1=args.fastq1
    fastq2=args.fastq2  
    prefix=args.prefix    
    if_dedup=args.dedup
        
    print "properties_file:",properties_file  
    print "genome_name:",genome_name
    print "refrence_genome fasta file:",ref_fasta
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "prefix:",prefix
    print "dedup:",if_dedup

        
def run_mapping(fastqfiles):  
    global ref_fasta          
    global sam_out
    global bam_out
    sam_out="{}/{}.sam".format(workdir,prefix)
    bam_out="{}/{}.bam".format(workdir,prefix)
    bwa=prop.get_attrib("bwa")
    samtools=prop.get_attrib("samtools")
    if ref_fasta == "{}_{}".format(ref_fasta_root_default,genome_name):        
        ref_fasta=prop.get_attrib(ref_fasta)
    
    #define command         
    fq1=fastqfiles[0]
    fq1_sai=prefix+'.fq1.sai'
    command("{} aln -t 12 {} {} >{}".format(bwa,ref_fasta,fq1,fq1_sai)).run_comm(0)
    if len(fastqfiles)==2:
        fq2=fastqfiles[1]
        fq2_sai=prefix+'.fq2.sai'        
        command("{} aln -t 12 {} {} >{}".format(bwa,ref_fasta,fq2,fq2_sai)).run_comm(0)
        command("{} sampe {} {} {} {} {}>{}.sam".format(bwa,ref_fasta,fq1_sai,fq2_sai,fq1,fq2,prefix)).run_comm(0)
    elif len(fastqfiles)==1:
        command("{} samse {} {} {}>{}.sam".format(bwa,ref_fasta,fq1_sai,fq1,prefix)).run_comm(0)
    command("{} sort {}.sam > {}.bam".format(samtools,prefix,prefix)).run_comm(0)

def deduplication():
    global bam_dedup_out
    rm_dup_sw_path=prop.get_attrib(rm_dup_sw)
    bam_dedup_out=bam_out.replace(".bam",".dedup.bam")
    command("{} rmdup -S {} {}".format(rm_dup_sw_path,bam_out,bam_dedup_out)).run_comm(0)

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
    run_mapping(fastqfiles)
    if if_dedup:
        print("removing duplicates...")
        deduplication()

def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir(sam_out,outdir)
    fi.copy_file_to_destdir(bam_out,outdir)
    if if_dedup :
        fi.copy_file_to_destdir(bam_dedup_out,outdir)
    
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
