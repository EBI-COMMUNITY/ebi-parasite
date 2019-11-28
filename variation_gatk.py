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
    global bam_file
    global prefix_ori
    global if_filter
    global mapping_file
    global runID
        
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script for getting filtered or unfiltered 
                                                    SNP and INDEL vcf and gvcf files from bam 
                                                    files using gatk, and then create statistics 
                                                    summary by using bcf_tools''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', 
                                                    required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                                                                 which is provided by genome_list.txt''', 
                                                    required=True)
    parser.add_argument('-bam', '--bam_file', type=str, help='Please provide one bam file', 
                                                    required=True)
    parser.add_argument('-f', '--if_filter', default=False, action="store_true", help='whether to filter SNP and INDEL', 
                                                    required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                                                    required=True)  
    parser.add_argument('-m', '--mapping_file', type=str, help='''if map file was provided, sample name will be 
                        included in the output file name''', required=False)
    args=parser.parse_args() 
    # check args
    FI.check_exist(args.properties_file)
    FI.check_exist(args.bam_file)     
    properties_file=args.properties_file
    prop=properties(properties_file)    
    MISC.check_genome_avl(prop.get_attrib("available_genomes"), args.genome_name)     
    runID=MISC.get_runID(args.bam_file)
    if args.mapping_file is not None:
        FI.check_exist(args.mapping_file)
    # define variable
    genome_name = args.genome_name
    bam_file=args.bam_file
    properties_file=args.properties_file
    prefix_ori=args.prefix      
    if_filter=args.if_filter
    mapping_file=args.mapping_file
    
    # print args
    print ("properties_file:",str(properties_file))
    print ("genome_name:", genome_name)
    print ("bam_file:",bam_file)
    print ("if_filter:",if_filter)
    print ("prefix:",prefix_ori)
    print ("mapping_file:",mapping_file)
        
def set_prefix():
    global prefix
    global sample_runID
    if mapping_file is None:
        sample_runID=runID
        prefix="{}_{}".format(prefix_ori,runID)
    else:
        sample_runID="{}_{}".format(MISC.get_samples_by_runIDs(mapping_file)[runID],runID)
        prefix="{}_{}".format(prefix_ori,sample_runID)

def call_var():  
    global raw_vcf_path
    global raw_gvcf_path
    global snp_path
    global indel_path
    global picard
    global bcftools
    global rtg
    global snpIndel_bcf_stats

    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)
    bam_file_working=indir+"/"+os.path.basename(bam_file)
    command("picard BuildBamIndex VALIDATION_STRINGENCY=LENIENT I={}".format(bam_file_working)).run_comm(0)
    raw_vcf_path="{}_raw.snps.indels.vcf".format(prefix)
    raw_gvcf_path="{}_raw.snps.indels.g.vcf".format(prefix)
    snp_path="{}_snp.vcf".format(prefix)
    indel_path="{}_indel.vcf".format(prefix)
    ### vcf
    command("gatk -T HaplotypeCaller -R {} -I {} -o {}".format(REF_FASTA,bam_file_working,raw_vcf_path)).run_comm(0)
    command("gatk -T SelectVariants -R {} -V {} -selectType SNP -o {}".format(REF_FASTA,raw_vcf_path,snp_path)).run_comm(0)
    command("gatk -T SelectVariants -R {} -V {} -selectType INDEL -o {}".format(REF_FASTA,raw_vcf_path,indel_path)).run_comm(0)
    ### gvcf
    command("gatk -T HaplotypeCaller -R {} -I {} --emitRefConfidence GVCF -o {}".format(REF_FASTA,bam_file_working,raw_gvcf_path)).run_comm(0)
    if if_filter:
        global filtered_snp_path
        global filtered_indel_path
        tmp_filtered_snp_path="{}_show_filt_snp.vcf".format(prefix)
        tmp_filtered_indel_path="{}_show_filt_indel.vcf".format(prefix)
        filtered_snp_path="{}_pass_filt.snp.vcf".format(prefix)
        filtered_indel_path="{}_pass_filt.indel.vcf".format(prefix)
        snp_filter_str="QD < 2.0 || MQ < 40.0 ||FS > 60.0 || SOR > 3.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0"
        indel_filter_str="QD < 2.0 || ReadPosRankSum < -20.0 || InbreedingCoeff < -0.8 || FS > 200.0 || SOR > 10.0"
        command("gatk -T VariantFiltration -R {} --filterExpression \"{}\" --filterName default_snp_filter -V {} -o {}"
                 .format(REF_FASTA,snp_filter_str,snp_path,tmp_filtered_snp_path)).run_comm(0)
        command("gatk -T VariantFiltration -R {} --filterExpression \"{}\" --filterName default_indel_filter -V {} -o {}"
                 .format(REF_FASTA,indel_filter_str,indel_path,tmp_filtered_indel_path)).run_comm(0)
        command("grep -e \'^#\' -e PASS {} >{}".format(tmp_filtered_snp_path,filtered_snp_path)).run_comm(0)
        command("grep -e \'^#\' -e PASS {} >{}".format(tmp_filtered_indel_path,filtered_indel_path)).run_comm(0)
        final_snp_path=filtered_snp_path
        final_indel_path=filtered_indel_path
    else:
        final_snp_path=snp_path
        final_indel_path=indel_path
    ## qc
    snpIndel_vcf_for_qc="{}.vcf".format(sample_runID)
    snpIndel_bcf_stats="{}.bcf_stats".format(sample_runID)
    command("gatk -T CombineVariants -R {} --variant {} --variant {} -o {} -genotypeMergeOptions UNIQUIFY".format(REF_FASTA,
                                                                   final_snp_path,final_indel_path,snpIndel_vcf_for_qc)).run_comm(0)
    command("bcftools stats {} >{} ".format(snpIndel_vcf_for_qc,snpIndel_bcf_stats)).run_comm(0)

def multiqc():
    global multiqc_outFN_prefix
    multiqc_outFN_prefix=sample_runID+".multiQC"
    command("multiqc -f {} -o {} --filename {} -v".format(snpIndel_bcf_stats,workdir,multiqc_outFN_prefix)).run_comm(0) 

def initiate():
    print ("initiating...")
    global indir
    global outdir    
    global reads_type     
    global workdir
    global subdir    
    global qcdir
    subdir="snp"
    workdir=prop.workdir+"/"+subdir   
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    qcdir=workdir+"/qc/"+prefix_ori
    FI.create_processing_dir(workdir)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    FI.create_processing_dir(qcdir)
    FI.copy_file_to_destdir(bam_file,indir)  
    
def execute():
    print ("executing...")
    FI.change_dir(workdir)
    set_prefix()
    call_var()
    multiqc()

def post_process():
    print ("post_processing...")
    FI.copy_file_to_destdir(raw_vcf_path,outdir)
    out_files=[raw_vcf_path,raw_gvcf_path,snp_path,indel_path]
    if if_filter:
        out_files.extend([filtered_snp_path,filtered_indel_path])
    for fpath in out_files:   
        FI.copy_file_to_destdir(fpath,outdir) 
    for qc_file in [snpIndel_bcf_stats, multiqc_outFN_prefix+".html"]:
        FI.copy_file_to_destdir(qc_file,qcdir)

if __name__ == '__main__':
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()

    get_args()
    print ("\n","Properties attributes:\n",prop.__dict__)
    getVar = lambda searchList, ind: [searchList[i] for i in ind]
 
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")

