#! /analysis/xin/parasite/bin/python/bin/python

import argparse
import re
import sys
import glob
from utilities import fileutils
from utilities import properties
from collections import defaultdict
from repeats_variation import strv

def get_args():    
    global properties_file
    global genome_name
    global genome_fasta
    global bam_file_pattern
    global genome_gff
    global prefix
    global prop
    global bam_files
    global bam_request_pattern
   
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script creates the STR variation summary file for all bam files, and STR vcf file, 
        annotation vcf file for each paired bam file''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name''', required=False) 
    parser.add_argument('-f', '--genome_fasta', type=str, help='''Please provide the genome fasta file,
        if the genome name is not 'ch' for C. hominis or 'cp' for C. parvum''', required=False)  
    parser.add_argument('-bp', '--bam_file_pattern', type=str, help="Please provide the paired bam files' pattern with the full path, must end with .bam", required=True)     
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  

    
    # check args
    args = parser.parse_args()
    fi.check_exist(args.properties_file)   
    if args.genome_fasta is not None:
        fi.check_exist(args.genome_fasta) 
    bam_files=glob.glob(args.bam_file_pattern)
    bam_request_pattern="^.*/?(.*?).bam$"
    for bam_file in bam_files:
        if not re.search(bam_request_pattern, bam_file):
            print "bam_file not end with .bam"
            sys.exit(1)
    fi.check_files_exist(bam_files) 
    
    # define variables     
    properties_file=args.properties_file    
    prop=properties(properties_file)    
    if ((args.genome_name is not None and args.genome_name!='ch' and args.genome_name!='cp') 
        or args.genome_name is None) and args.genome_fasta is None:
        print '''Please provide the genome fasta file, if the genome name is not 'ch' for C. hominis or 'cp' for C. parvum'''
        sys.exit(1)
    else:
        genome_name=args.genome_name           
    if args.genome_fasta is None:
        genome_fasta=prop.get_attrib(genome_name+"_fasta")
    bam_file_pattern=args.bam_file_pattern
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome_name:",genome_name
    print "genome_fasta:",genome_fasta
    print "bam_file_pattern:",bam_file_pattern
    print "prefix:",prefix 
       
def submit_str_jobs():    
    global new_prefixs
    global subdir
    subdir="repeats"
    new_prefixs=[]
    for bam_file in bam_files:
        bam_file=bam_file.rstrip("\n")
        bam_request_pattern="^(.*/)?(.*?).bam$"
        bam_prefix=re.findall(bam_request_pattern, bam_file)[0][1]
        new_prefix=prefix+"."+bam_prefix
        new_prefixs.append(new_prefix)       
        st=strv(properties_file,genome_name,genome_fasta,bam_file,new_prefix,"true",subdir)
        st.run() 


def summary():
    ann_pattern='ANN=(.*?)\s+'          
    var_dict=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    region_syn_dict=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    region_non_syn_dict=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    trf_dict=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    workdir=prop.workdir+"/"+subdir
    fPath_trf="{}/{}.trf.str".format(workdir,genome_name)
    fPath_summary="{}/out/{}.csv".format(workdir,prefix)
    file_trf=open(fPath_trf, 'r')
    fout_summary=open(fPath_summary, 'w')
    
    fout_summary.write(",".join('''chr_acc start_posi end_posi variant_num synonymous_variant_num missense_variant_num'''.split())+"\n")
    for line in file_trf:
        line=line.rstrip("\n")
        if not line.startswith('#') and re.search("\w+",line):
            (con_acc,region_start,region_end)=get_sub_list(line.split(),[0,1,2])
            trf_dict[con_acc][region_start][region_end]=0
            region_syn_dict[con_acc][region_start][region_end]=0
            region_non_syn_dict[con_acc][region_start][region_end]=0
    for new_prefix in new_prefixs:
        fPath_vcf="{}/out/{}.ann.vcf".format(workdir,new_prefix)
        infile=open(fPath_vcf, 'r')
        for line in infile:
            if not line.startswith( '#' ): 
                (con_acc,posi)=get_sub_list(line.split(),[0,1])
                if re.compile(ann_pattern).search(line) :
                    ann_field=re.findall(ann_pattern, line)[0]
                    annseg_array=ann_field.split(",")
                    for annseg in annseg_array:
                        variant=annseg.split("|")[1]
                        gene_acc=annseg.split("|")[3]
                        var_dict[con_acc][posi][variant]=gene_acc 
                        for region_start in trf_dict[con_acc]:
                            for region_end in trf_dict[con_acc][region_start]:
                                if int(posi)>=int(region_start) and int(posi)<=int(region_end):
                                    trf_dict[con_acc][region_start][region_end]=trf_dict[con_acc][region_start][region_end]+1
                                    if variant=="synonymous_variant":
                                        region_syn_dict[con_acc][region_start][region_end]=region_syn_dict[con_acc][region_start][region_end]+1
                                    if variant=="missense_variant":
                                        region_non_syn_dict[con_acc][region_start][region_end]=region_non_syn_dict[con_acc][region_start][region_end]+1
    for con_acc in sorted(trf_dict.keys()):
        for region_start in sorted(trf_dict[con_acc]):
            for region_end in sorted(trf_dict[con_acc][region_start]):
                fout_summary.write("{},{},{},{},{},{}".format(con_acc,region_start,region_end,
                                                                   str(trf_dict[con_acc][region_start][region_end]),
                                                                   str(region_syn_dict[con_acc][region_start][region_end]),
                                                                   str(region_non_syn_dict[con_acc][region_start][region_end])
                                                                   )
                                   +"\n")
    

if __name__ == '__main__':
    get_sub_list = lambda searchList, ind: [searchList[i] for i in ind]    
    global fi
    fi=fileutils()
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    submit_str_jobs()
    summary()
