#! /usr/bin/env python3

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
import os
import re
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict

def get_args():
    global prop
    global properties_file
    global genome_name
    global gvcf_files
    global gvcf_files_str
    global prefix
    global if_filter
    filter_dict={}
        
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script for getting SNP from bam files using gatk')
    parser.add_argument('-p', '--properties_file', type=str, help='''Please provide the properties file, 
                                                                     which including the paths of workdir''', 
                                                             required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                                                                 which is provided by genome_list.txt''', 
                                                             required=True)
    parser.add_argument('-gv', '--gvcf_files', type=str, help='Please provide gvcf files', required=True)
    parser.add_argument('-f', '--if_filter', default=False, action="store_true", 
                                                         help='whether to filter SNP and INDEL seperately', 
                                                             required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                                                             required=True)  
    
    args=parser.parse_args() 
    # check args
    FI.check_exist(args.properties_file)
    FI.check_files_exist(glob.glob(args.gvcf_files))     
    properties_file=args.properties_file
    prop=properties(properties_file)    
    MISC.check_genome_avl(prop.get_attrib("available_genomes"), args.genome_name)     
    genome_name = args.genome_name
    # define variables  
    properties_file=args.properties_file
    gvcf_files=glob.glob(args.gvcf_files)
    gvcf_files_str=""
    for gvcf_file in gvcf_files:
        gvcf_files_str+=gvcf_file+" "
    gvcf_files_str=gvcf_files_str.rstrip(" ")
    prefix=args.prefix 
    if_filter=args.if_filter
    
    # print args
    print ("properties_file:",str(properties_file))  
    print ("genome_name:", genome_name)
    print ("gvcf_files:",gvcf_files_str)
    print ("if_filter:",if_filter)
    print ("prefix:",prefix)
        
def run():  
    global merged_gvcf_path
    global merged_gvcf_snp_path
    global merged_gvcf_indel_path
    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)
    merged_gvcf_path="{}_merged.snps.indels.g.vcf".format(prefix)
    merged_gvcf_snp_path="{}_merged.snps.g.vcf".format(prefix)
    merged_gvcf_indel_path="{}_merged.indels.g.vcf".format(prefix)
    var_str=""
    for gvcf_file in gvcf_files:
        runID = re.findall("([A-Z]RR\d+)", gvcf_file)[0]
        var_str+=" --variant:{},vcf {}".format(runID, gvcf_file)
    command("gatk -T GenotypeGVCFs -R {} {} -o {}".format(REF_FASTA, var_str, merged_gvcf_path)).run_comm(0)
    command("gatk -T SelectVariants -R {} -V {} -selectType SNP -o {}".format(REF_FASTA,merged_gvcf_path,merged_gvcf_snp_path)).run_comm(0)
    command("gatk -T SelectVariants -R {} -V {} -selectType INDEL -o {}".format(REF_FASTA,merged_gvcf_path,merged_gvcf_indel_path)).run_comm(0)
    if if_filter:
        global filtered_snp_path
        global filtered_indel_path
        tmp_filtered_snp_path="{}_show_filt_snp.g.vcf".format(prefix)
        tmp_filtered_indel_path="{}_show_filt_indel.g.vcf".format(prefix)
        filtered_snp_path="{}_pass_filt.snp.g.vcf".format(prefix)
        filtered_indel_path="{}_pass_filt.indel.g.vcf".format(prefix)
        snp_filter_str="QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        indel_filter_str="QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0"
        command("gatk -T VariantFiltration -R {} --filterExpression \"{}\" --filterName default_snp_filter -V {} -o {}"
                 .format(REF_FASTA,snp_filter_str,merged_gvcf_snp_path,tmp_filtered_snp_path)).run_comm(0)
        command("gatk -T VariantFiltration -R {} --filterExpression \"{}\" --filterName default_indel_filter -V {} -o {}"
                 .format(REF_FASTA,indel_filter_str,merged_gvcf_indel_path,tmp_filtered_indel_path)).run_comm(0)
        command("grep -e \'^#\' -e PASS {} >{}".format(tmp_filtered_snp_path,filtered_snp_path)).run_comm(0)
        command("grep -e \'^#\' -e PASS {} >{}".format(tmp_filtered_indel_path,filtered_indel_path)).run_comm(0)

def initiate():
    print ("initiating...")
    global indir
    global outdir    
    global reads_type     
    global workdir
    global subdir    
    subdir="snp"
    workdir=prop.workdir+"/"+subdir   
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    for gvcf_file in gvcf_files:
        FI.copy_file_to_destdir(gvcf_file,indir)  
    
def execute():
    print ("executing...")
    FI.change_dir(workdir)    
    run()

def post_process():
    print ("post_processing...")
    fpath_array=[merged_gvcf_path, merged_gvcf_snp_path, merged_gvcf_indel_path]
    if if_filter:
        fpath_array.extend([filtered_snp_path,filtered_indel_path])
    for fpath in fpath_array:   
        FI.copy_file_to_destdir(fpath,outdir) 
 
if __name__ == '__main__':
    global FI
    global MISC
    FI = fileutils()
    MISC=misc()    
    get_args()
    getVar = lambda searchList, ind: [searchList[i] for i in ind]
    print ("\n","Properties attributes:")
    print (prop.__dict__)
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")

