#! /usr/bin/python

import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict
import os

def get_args():    
    global properties_file
    global genome
    global prefix
    global vcf_file_pattern
    global prop
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script invests genes under selection pressure within species through dNdS. 
        Species can be chosen from genome_list.txt''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name, only with those obtained from genome_list.txt''', 
                        required=True) 
    parser.add_argument('-f', '--vcf_file_pattern', type=str, help="Please provide snp vcf files' pattern with full file path", 
                        required=True)     
    parser.add_argument('-pre', '--prefix', type=str,help='Please provide the prefix for the output file.', 
                        required=True)  


    # check args
    args = parser.parse_args() 
    fi=fileutils()    
    fi.check_exist(args.properties_file)  
    properties_file=args.properties_file    
    prop=properties(properties_file)
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()) and args.genome_name!="cryptosporidium_hominis":
        misc.my_exit("{} is not available, please try another genome".format(args.genome_name))     
    if not re.search(".vcf",args.vcf_file_pattern):
        misc.my_exit("vcf_file_pattern need to end up with .vcf")
    genome=args.genome_name
    vcf_file_pattern=args.vcf_file_pattern
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome:",genome
    print "vcf_file_pattern:",vcf_file_pattern
    print "prefix:",prefix
    
   
def get_gff_from_genome_name():
    global genome_gff_fname
    if genome=="cryptosporidium_hominis":        
        command("cp -p {} .".format(prop.get_attrib("ch_gff"))).run_comm(0)
        genome_gff_fname="GCA_002223825.1_C.hominis.v1_genomic.gff"  
    else:
        genome_gff_fname=misc.download("gff3",genome,"")
       
def get_chr_map_from_gff():
    global chr_map
    chr_map={}
    pattern="chromosome:(\d+);Alias=(\w+.\d+)"
    file_fh=open(genome_gff_fname,'r')
    for line in file_fh:
        line=line.rstrip()
        if re.search(pattern,line): 
             result = re.findall(pattern,line)[0]      
             chr_map[result[0]]=result[1]
               
def get_genes_from_gff():
    global genome_gene_dict
    get_chr_map_from_gff()
    genome_gene_dict={}
    pattern_ch="Name=(.*?);"
    pattern="ID=gene:(.*?);"
    file_gff=open(genome_gff_fname,'r')
    for line in file_gff:
        line=line.rstrip("\n")
        if not line.startswith("#") and line.split("\t")[2]=='gene':
            if genome=="cryptosporidium_hominis" and re.search(pattern_ch,line):
                genome_gene_dict[re.findall(pattern_ch,line)[0]]=line.split("\t")[0] #genome_gene_dict[gene]=chr_acc
            elif re.search(pattern,line):
                genome_gene_dict[re.findall(pattern,line)[0]]=chr_map[line.split("\t")[0]]
                #1    ena    gene    29719    31791    .    +    .    ID=gene:cgd1_130;b --cp
                #LN877947.1    EMBL    gene    1367    2680    .    -    .    ID=gene0;Name=CHUDEA1_10;gbkey=Gene;gene_biotype=protein_coding;locus_tag=CHUDEA1_10 --ch
            else: 
                print "ERROR: there is no gene name in {}".format(line)
                sys.exit(1)

def build_snpeff_db():
    snpeff_config=prop.get_attrib("snpeff_config")
    if re.search("\w+",command("grep {} {}".format(genome,snpeff_config)).run_comm(1)):
        return
    snpeff_dir=prop.get_attrib("snpeff_dir")
    db_dir=snpeff_dir+"/data/"+genome
    fi.make_dir(db_dir)
    #Copy the files into snpEff's directory structure
    command("cp {} {}/genes.gff".format(genome_gff_fname, db_dir)).run_comm(0)
    command("cp {} {}/sequences.fa".format(misc.download("fasta",genome,"dna"),db_dir)).run_comm(0)

    #Edit snpEff.config and insert your specific database information:
    command("echo \"{}.genome : {}\" >> {}".format(genome,genome,snpeff_config)).run_comm(0)

    #Build the database
    command("java -jar {} build -gff3 -v {}".format(prop.get_attrib("snpeff"),genome)).run_comm(0)

def run_snpEff():         
    global ann_vcf_fpaths  
    ann_vcf_fpaths=[]
    snpEff=prop.get_attrib("snpeff")    
    build_snpeff_db()
    for vcf_fpath in all_vcf_fpaths:
        vcf_fpath_prefix=os.path.basename(vcf_fpath).rstrip(".vcf")
        ann_vcf_fpath=vcf_fpath_prefix+".ann.vcf"
        ann_vcf_fpaths.append(ann_vcf_fpath)
        command("java -jar {} -c snpEff.config {} {} > {}".format(snpEff,genome,vcf_fpath,ann_vcf_fpath)).run(0);
        
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
    get_gff_from_genome_name()
    run_snpEff()


def post_process():
    print "post_processing..."
    fi.copy_file("snpEff_summary.html","{}/{}_snpEff_summary.html".format(outdir,prefix))
    fi.copy_file("snpEff_genes.txt","{}/{}_snpEff_genes.txt".format(outdir,prefix))
    for ann_vcf_fpath in ann_vcf_fpaths:
        fi.copy_file_to_destdir(ann_vcf_fpath,outdir)
    #fi.copy_file_to_destdir(ann_csv_fname,outdir)


if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    fi=fileutils()
    misc=misc()
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run_blast the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()

