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
    global genome_fasta
    global bam_file_pattern
    global bam_request_pattern
    global bam_files
    global map_fpath
    global map_dict
    global prefix
                       
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script creating relocation files for multiple bam files from 
                                    various genomes and automatically open the GUI''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-m', '--map_file', type=str, help='''Please provide the map file, 
                        in which the first column is the full path of the genome fasta file and the second column is
                        the full path of the bam file and the bam files need to ended with .bam''', required=True)   
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  

    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_files_exist([args.properties_file,args.map_file])   
    
    # define variables         
    properties_file=args.properties_file    
    prop=properties(properties_file)
    map_fpath=args.map_file
    fh_map=open(map_fpath, "r")
    map_dict={}
    for line in fh_map:
        line=line.rstrip()        
        (fasta_fpath,bam_fpath)=getVar(line.split(),[0,1])
        fi.check_files_exist([fasta_fpath,bam_fpath])
        map_dict[bam_fpath]=fasta_fpath
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "map_file:",map_fpath
    print "prefix:",prefix

def run_get_relocation():
    global output_fname
    samtools=prop.get_attrib("samtools");
    bcftools=prop.get_attrib("bcftools");
    tabix=prop.get_attrib("tabix");
    progressiveMauve=prop.get_attrib("progressivemauve");
    mauve=prop.get_attrib("mauve");
    all_bam_fasta_fname_str="";
    for bam_fpath in (map_dict.keys()):
        genome_fasta=map_dict[bam_fpath]
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
    command("{} --output={} {}".format(progressiveMauve,output_fname,all_bam_fasta_fname_str)).run_comm(0)
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
    run_get_relocation()

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
   

