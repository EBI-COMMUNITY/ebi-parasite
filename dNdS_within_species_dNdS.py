#! /analysis/xin/parasite/bin/python/bin/python

import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict
import os

def get_args():    
    global properties_file
    global genome
    global gff
    global prefix
    global vcf_file_pattern
    global prop
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script invests genes under selection pressure within species through dNdS')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name, 
                        only "ch" for "C. hominis" or "cp" for "C. parvum" can be used''', required=True) 
    parser.add_argument('-gff', '--genome_gff_file', type=str, help='''Please provide the genome gff file, 
                        only C. hominis or C. parvum gff file can be used''', required=False)  
    parser.add_argument('-f', '--vcf_file_pattern', type=str, help="Please provide vcf files' pattern with full file path", required=True)     
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  
    
    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)   
    if args.genome_gff_file is not None:
        fi.check_exist(args.genome_gff_file) 
        
    # define variables     
    properties_file=args.properties_file    
    prop=properties(properties_file)    
    if args.genome_name!='ch' and args.genome_name!='cp':
        print "only 'ch' or 'cp' can be used as the genome name"
        sys.exit(1)
    else:
        genome=args.genome_name           
    if args.genome_gff_file is None:
        gff=prop.get_attrib(genome+"_gff")
    vcf_file_pattern=args.vcf_file_pattern
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome:",genome
    print "genome_gff:",gff
    print "vcf_file_pattern:",vcf_file_pattern
    print "prefix:",prefix
   
   
def get_genes_from_gff():
    global genome_gene_dict
    genome_gene_dict={}
    pattern="Name=(.*?);"
    file_gff=open (gff,'r')
    for line in file_gff:
        line=line.rstrip("\n")
        if not line.startswith("#") and line.split("\t")[2]=='gene':
            if re.search(pattern,line):
                genome_gene_dict[re.findall(pattern,line)[0]]=line.split("\t")[0] #genome_gene_dict[gene]=chr_acc
            else: 
                print "ERROR: there is no gene name in {}".format(line)
                sys.exit(1)


def run_snpEff():         
    global ann_vcf_fpaths
    annot_sw="snpeff"    
    ann_vcf_fpaths=[]
    snpEff=prop.get_attrib(annot_sw)  
    genome_db={"ch":"c_hominis",
               "cp":"c_parvum"}    
    for vcf_fpath in all_vcf_fpaths:
        vcf_fpath_prefix=os.path.basename(vcf_fpath).rstrip(".vcf")
        ann_vcf_fpath=vcf_fpath_prefix+".ann.vcf"
        ann_vcf_fpaths.append(ann_vcf_fpath)
        command("java -jar {} -c snpEff.config {} {} > {}".format(snpEff,genome_db[genome],vcf_fpath,ann_vcf_fpath)).run_comm(0);
        
        
def analysis_ann_files():    
    global fPath_out     
    no_dNdS_gene_dict={}
    dNdS_gene_dict={}
    snp_hash=defaultdict(lambda: defaultdict(str))
    ann_pattern='ANN=(.*?)\s+'
    fPath_out="{}.{}.ann.csv".format(prefix,genome)
    fileout=open(fPath_out,'w')
    fileout.write("{},{},{},{},\n".format("GENOME_ACC","GENE_NAME","SYNONYMOUS_SNP_NUM","NONSYNONYMOUS_SNP_NUM"))
    for ann_vcf_fpath in ann_vcf_fpaths:
        infile=open(ann_vcf_fpath, 'r')
        for line in infile:
            if not line.startswith( '#' ):
                if re.compile(ann_pattern).search(line) :
                    ann_field=re.findall(ann_pattern, line)[0]
                    annseg_array=ann_field.split(",")
                    for annseg in annseg_array:
                        variant=annseg.split("|")[1]
                        gene_acc=(annseg.split("|")[3])
                        if not re.search('[a-zA-Z]', gene_acc) or not re.search('[a-zA-Z]', variant):
                            print "WARNING: None:"+annseg+"\n"+line
                        if gene_acc in snp_hash.keys() and variant in snp_hash[gene_acc].keys() :
                            snp_hash[gene_acc][variant]=snp_hash[gene_acc][variant]+1
                        else:
                            snp_hash[gene_acc][variant]=1
    get_genes_from_gff()
    for gene in sorted(genome_gene_dict.keys()):
        genome_acc=genome_gene_dict[gene]
        if gene not in snp_hash or ("synonymous_variant" not in snp_hash[gene].keys() and "missense_variant" not in snp_hash[gene].keys()):
            no_dNdS_gene_dict[gene]=genome_acc
        else:
            dNdS_gene_dict[gene]=genome_acc
    for gene in sorted(no_dNdS_gene_dict.keys()):
        fileout.write("{},{},0,0\n".format(no_dNdS_gene_dict[gene],gene))
    for gene in sorted(dNdS_gene_dict.keys()):
        if "synonymous_variant" not in snp_hash[gene].keys():
            fileout.write("{},{},0,{}\n".format(dNdS_gene_dict[gene],gene,str(snp_hash[gene]["missense_variant"])))
        elif "missense_variant" not in snp_hash[gene].keys():
            fileout.write("{},{},{},0\n".format(dNdS_gene_dict[gene],gene,str(snp_hash[gene]["synonymous_variant"])))
        else:
            fileout.write("{},{},{},{}\n".format(dNdS_gene_dict[gene],gene,str(snp_hash[gene]["synonymous_variant"]),str(snp_hash[gene]["missense_variant"])))
    fileout.close()
    #print "no_dNdS_gene "+str(len(no_dNdS_gene_dict.keys()))
    #print "dNdS_gene_dict "+str(len(dNdS_gene_dict.keys()))


def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global type    
    global all_vcf_fpaths
    subdir="dNdS"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    all_vcf_fpaths=glob.glob(vcf_file_pattern)
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    for file_in in (all_vcf_fpaths): ### ? GFF file
        fi.copy_file_to_destdir(file_in,indir)
    fi.change_dir(workdir)    


def execute():
    print "executing..."
    run_snpEff()
    analysis_ann_files()


def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir(fPath_out,outdir)
    

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    fi=fileutils()
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run_blast the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()

