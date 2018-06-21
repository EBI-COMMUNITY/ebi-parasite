#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

import sys
import StringIO
import pylab
import re
import argparse
import re
import os
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict
import glob
from docutils.nodes import title

def get_args():  
    global fi  
    global prop
    global properties_file
    global genome_name
    global genome_fasta
    global bam_file_pattern
    global bam_request_pattern
    global bam_files
    global prefix
                       
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script creating relocation files for multiple bam files from the same genome 
                                    and automatically open the GUI.''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''if "C. hominis" or "C. parvum" will be used as the genome, 
                                                                 please provide "ch" for "C. hominis" or "cp" for "C. parvum"''', 
                                                                 required=False)  
    parser.add_argument('-f', '--genome_fasta', type=str, help='''Please provide the directory for the genome fasta file, 
                                                                 if "ch" or "cp" is not the genome name.''', 
                                                                 required=False)
    parser.add_argument('-b', '--bam_file_pattern', type=str, help='''Please provide the bam files' pattern with the full path''', required=True)      
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  

    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)  
    if args.genome_name is not None and not args.genome_name=="ch" and not args.genome_name=="cp" :
        print "genome name need to be ch or cp"
        sys.exit(1)      
    if args.genome_fasta is not None:            
        fi.check_exist(args.genome_fasta) 
    bam_file_pattern=args.bam_file_pattern
    bam_files=glob.glob(bam_file_pattern)
    bam_request_pattern="^.*/?(.*?).bam$"
    for bam_file in bam_files:
        if not re.search(bam_request_pattern, bam_file):
            print "bam_file not ended with .bam"
            sys.exit(1)
    fi.check_files_exist(bam_files)  
    
    # define variables         
    properties_file=args.properties_file    
    prop=properties(properties_file)
    if args.genome_name is not None:
        genome_name=args.genome_name
    if args.genome_fasta is not None:    
        genome_fasta=args.genome_fasta   
    else:
        if args.genome_name is None:
            print "If no genome_fasta provided, genome name must be provided as ch or cp."
            sys.exit(1)
        else:
            genome_fasta=prop.get_attrib(genome_name+"_fasta")
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome_name:",genome_name
    print "genome_fasta:",genome_fasta
    print "bam_file_pattern:",bam_file_pattern
    print "prefix:",prefix

def run_relocation():
    global output_fname
    samtools=prop.get_attrib("samtools");
    bcftools=prop.get_attrib("bcftools");
    tabix=prop.get_attrib("tabix");
    progressiveMauve=prop.get_attrib("progressivemauve");
    mauve=prop.get_attrib("mauve");
    all_bam_fasta_fname_str="";
    for bam_fpath in (bam_files):
        bam_fname=os.path.basename(bam_fpath)
        sorted_bam_fname=bam_fname.rstrip(".bam")+".sorted.bam"
        compressed_sorted_bam_fname=sorted_bam_fname+".gz"
        sorted_bam_fasta_fname=sorted_bam_fname+".fa"
        command("{} sort {} > {}".format(samtools, bam_fpath,sorted_bam_fname)).run_comm(0)
        command("{} mpileup -uf {} {} -B | {} call -mv -Oz -o {} ".format(samtools, genome_fasta,sorted_bam_fname,bcftools,compressed_sorted_bam_fname)).run_comm(0)
        command("{} {}".format(tabix, compressed_sorted_bam_fname)).run_comm(0)
        command("cat {} | {} consensus {} > {}".format(genome_fasta,bcftools,compressed_sorted_bam_fname,sorted_bam_fasta_fname)).run_comm(0)
        all_bam_fasta_fname_str=all_bam_fasta_fname_str+" "+sorted_bam_fasta_fname
    output_fname="{}.xmfa".format(prefix)
    command("{} --output={} {}".format(progressiveMauve, output_fname, all_bam_fasta_fname_str)).run_comm(0)
    command("{} {}".format(mauve,output_fname)).run_comm(0)    

def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global type        
    subdir="structure"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.change_dir(workdir)    

def execute():
    print "executing..."
    run_relocation()

def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir(output_fname,outdir)

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run_blast the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   

