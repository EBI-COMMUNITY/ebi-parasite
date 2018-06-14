#! /analysis/xin/parasite/bin/python/bin/python

import argparse
import sys
import os
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict

def get_args():
    global fi
    global properties_file
    global genome_name
    global ref_fasta
    global bam_files
    global bam_files_str
    global prefix
    global filter_dict
    global default_ref_fasta_root
    global if_classify
    fi=fileutils()
    default_ref_fasta_root='ena_ref_fasta'
    filter_dict={}
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script for getting SNP from bam files using samtools and bcftools')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file, which including the paths of samtools and bcftools and workdir', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted, 
                                                                \'ch\' stands for \'Cryptosporidium hominis\' and 
                                                                \'cp\' stands for \'Cryptosporidium parvum\'''', required=True)
    parser.add_argument('-r', '--ref_fasta', type=str, help='Please provide reference genome fasta file', required=False)
    parser.add_argument('-bam', '--bam_files', type=str, nargs='+', help='Please provide one or multiple bam files', required=True)
    parser.add_argument('-qual', '--QUAL_filter', type=str, help='To filter QUAL(phred-scaled quality score) in vcf file, please provide minimum QUAL value', required=False)
    parser.add_argument('-dp', '--DP_filter', type=str, help='To filter DP("Raw read depth") in vcf file, please provide the minimum DP value', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', required=True)  
    parser.add_argument('-c', '--classify', action='store_true',help='set classifying the result VCF file into coding and non-coding file to true', default=False)
    
    args=parser.parse_args() 
    
    # check args
    fi.check_exist(args.properties_file)
    fi.check_files_exist(args.bam_files)     
    if args.genome_name is not "ch" and not "cp":
        print "ERROR: Please provide the genome name you are mapping to, only \'ch\' or \'cp\' permitted"
        sys.exit(1)
    if args.ref_fasta is not None:
        fi.check_exist(args.ref_fasta)         
    
    # define variables  
    properties_file=args.properties_file
    genome_name=args.genome_name 
    if args.ref_fasta is None:
        ref_fasta=default_ref_fasta_root+'_'+genome_name
    else:
        ref_fasta=args.ref_fasta
    bam_files=args.bam_files
    bam_files_str=""
    for bam_file in bam_files:
        bam_files_str += os.path.abspath(bam_file)+" "
    if args.QUAL_filter is not None:
        filter_dict["QUAL"]=args.QUAL_filter
    if args.DP_filter is not None:
        filter_dict["DP"]=args.DP_filter
    prefix=args.prefix      
    if_classify=args.classify

    # print args
    print "properties_file:",str(properties_file)  
    print "refrence_genome fasta file:",ref_fasta
    print "bam_files:",bam_files_str
    print "QUAL_filter:",args.QUAL_filter
    print "DP_filter:",args.DP_filter
    print "prefix:",prefix
    print "if_classify:",if_classify


def run_snp():  
    global var_vcf_path
    global filtered_vcf_path
    samtools=prop.get_attrib("samtools")
    bcftools=prop.get_attrib("bcftools")
    if ref_fasta == "{}_{}".format(default_ref_fasta_root,genome_name):
        ref_fasta_path=prop.get_attrib(ref_fasta)
    var_vcf_path="{}.var.vcf".format(prefix)
    filtered_vcf_path="{}.filtered.vcf".format(prefix)
    #define command    
    command("{} mpileup --skip-indels --BCF --output-tags DP,AD,ADF,ADR,SP -f {} -o {}.raw.bcf {}".format(samtools,ref_fasta_path, prefix, bam_files_str)).run_comm(0)
    command("{} call --skip-variants indels --multiallelic-caller --variants-only -O v {}.raw.bcf -o {}".format(bcftools,prefix,var_vcf_path)).run_comm(0)
      
    filter_str=""
    if len(filter_dict.keys())>0:
        filter_str+="-i\'"
        for filter in filter_dict.keys():
            filter_str+="{}>{} && ".format(filter,filter_dict[filter])
        filter_str=filter_str.rstrip(" && ")+"'"
        command("{} query {} -f'%CHROM %POS %QUAL %DP\\n' {}.var.vcf > {}".format(bcftools,filter_str,prefix,filtered_vcf_path)).run_comm(0)

        
def classify_cds(in_vcf_fpath):      
    global out_coding_fpath   
    global out_noncoding_fpath
    global out_coding_stat_fpath
    cds_map = defaultdict(lambda: defaultdict(str))
    genome_base_map={}
    genome_base_map['ch']=9043938
    genome_base_map['cp']=9102324
    coding_perc_ch_and_cp="75%"
    vcf_coding_snp_num=0
    vcf_noncoding_snp_num=0
    in_ref_cds_fpath = prop.get_attrib(genome_name+"_cds")
    out_coding_fpath=in_vcf_fpath.replace(".vcf",".coding.vcf")
    out_noncoding_fpath=in_vcf_fpath.replace(".vcf",".noncoding.vcf")
    out_coding_stat_fpath=prefix+"_snp_classify.stat"
    try:
        fhin_vcf=open(in_vcf_fpath,"r")
        fhin_ref_cds=open(in_ref_cds_fpath,"r")
        fhout_coding=open(out_coding_fpath,"w")
        fhout_noncoding=open(out_noncoding_fpath,"w")
        fhout_coding_stat=open(out_coding_stat_fpath,"w")
        for line in fhin_ref_cds:
            (acc,start_posi,end_posi)=line.split()
            for i in range(int(start_posi),int(end_posi)+1):
                cds_map[acc][str(i)]=1
        for line in fhin_vcf:
            if not line.startswith("#"):
                acc=line.split()[0]
                posi = line.split()[1]
                if cds_map[acc][posi]:
                    fhout_coding.write(line)
                    vcf_coding_snp_num=vcf_coding_snp_num+1
                else :
                    fhout_noncoding.write(line)
                    vcf_noncoding_snp_num=vcf_noncoding_snp_num+1
                    vcf_total_snp_num=vcf_coding_snp_num+vcf_noncoding_snp_num
            else:
                fhout_coding.write(line)
                fhout_noncoding.write(line)
        fhout_coding_stat.write("There are totally {} bases in {} genome, {} in coding area\n".format(genome_base_map[genome_name], genome_name, coding_perc_ch_and_cp))
        fhout_coding_stat.write("In the result vcf file:\n{} SNPs in coding area, which is {}% among all SNPs\n".format(vcf_coding_snp_num, int(float(vcf_coding_snp_num)/vcf_total_snp_num*100)))
        fhout_coding_stat.write("{} SNPs in non-coding area, which is {}% among all SNPs\n".format(vcf_noncoding_snp_num, int(float(vcf_noncoding_snp_num)/vcf_total_snp_num*100)))
        
    except:
        print ("Unexpected error:", sys.exc_info()[1])
        raise
 
    finally:
        if fhin_vcf is not None:
            fhin_vcf.close()
        if fhin_ref_cds is not None:
            fhin_ref_cds.close()
        if fhout_coding is not None:
            fhout_coding.close()
        if fhout_noncoding is not None:
            fhout_noncoding.close()
        if fhout_coding_stat is not None:
            fhout_coding_stat.close() 
   

def initiate():
    print("initiating...")
    global indir
    global outdir    
    global reads_type     
    global workdir
    global subdir      
    
    subdir="snp" 
    workdir=prop.workdir+"/"+subdir   
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    for bam_file in bam_files:
        fi.copy_file_to_destdir(bam_file,indir)

def execute():        
    print("executing...")
    fi.change_dir(workdir)    
    run_snp()
    if if_classify:
        print("classifing SNPs...\n")  
        classify_cds(var_vcf_path)
        if os.path.isfile(filtered_vcf_path) and os.path.getsize(filtered_vcf_path)>0:
            classify_cds(filtered_vcf_path)


def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir("{}/{}.raw.bcf".format(workdir,prefix),outdir)
    fpath_array=[var_vcf_path,filtered_vcf_path]
    if 'out_coding_fpath' in globals():
	    fpath_array=[var_vcf_path,filtered_vcf_path,out_coding_fpath,out_noncoding_fpath,out_coding_stat_fpath] 
    for fpath in fpath_array :
        if os.path.isfile(fpath):
            if os.path.getsize(fpath) is 0:
                print "{} is empty, so will not be copied into {}".format(fpath,outdir)
            else:
		fi.copy_file_to_destdir(fpath,outdir) 


if __name__ == '__main__':
    global prop
    get_args()
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

