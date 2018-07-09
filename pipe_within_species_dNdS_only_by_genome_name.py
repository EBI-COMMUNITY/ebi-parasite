#! /nfs/production/seqdb/embl/developer/xin/bin/python/bin/python

import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict
import os

def get_gene_list_str():
   return '''Cryptosporidium_hominis
Cryptosporidium_parvum_iowa_ii
Eimeria_tenella
Leishmania_major
Phytophthora_infestans
Phytophthora_ramorum
Plasmodium_berghei
Plasmodium_chabaudi
Plasmodium_falciparum
Plasmodium_knowlesi
Plasmodium_vivax
Plasmodium_yoelii_yoelii
Pythium_ultimum
Theileria_annulata
Theileria_parva
Toxoplasma_gondii
Trypanosoma_cruzi'''
    
def get_args():    
    global properties_file
    global genome
    global prefix
    global vcf_file_pattern
    global prop
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script invests genes under selection pressure within species through dNdS. 
        Species can be chosen from -genome_list, which including 17 genomes. They are the common genomes of protists parasite and existing in snpEff''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-genome_list', '--genome_list', help="This will display the genome name list", action="store_true")
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name, only with those obtained from -genome_list''', 
                        required='-genome_list' not in sys.argv) 
    parser.add_argument('-f', '--vcf_file_pattern', type=str, help="Please provide snp vcf files' pattern with full file path", 
                        required='-genome_list' not in sys.argv)     
    parser.add_argument('-pre', '--prefix', type=str,help='Please provide the prefix for the output file.', 
                        required='-genome_list' not in sys.argv)  


    # check args
    args = parser.parse_args() 
    if args.genome_list:
        print get_gene_list_str()
        sys.exit(0)
    fi=fileutils()    
    fi.check_exist(args.properties_file)  
    properties_file=args.properties_file    
    prop=properties(properties_file)
    genome=args.genome_name
    if genome not in get_gene_list_str().split("\n"):
        print "ERROR: genome_name {} not in the list of -genome_list".format(genome)   
        sys.exit(1)           
    vcf_file_pattern=args.vcf_file_pattern
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome:",genome
    print "vcf_file_pattern:",vcf_file_pattern
    print "prefix:",prefix
    
   
def get_gff_from_genome_name():
    global genome_gff_fname
    if genome=="Cryptosporidium_hominis":        
        genome_gff_fname="GCA_002223825.1_C.hominis.v1_genomic.gff"
        genome_zip_gff_fname=genome_gff_fname+".gz"
        command("cp -p {} .".format(prop.get_attrib("ch_gff"))).run_comm(0)
        command("gunzip -f {}".format(genome_zip_gff_fname)).run_comm(0)
        return
    sub_dir2=""
    gff_ftp=""
    ftp_dir="ftp://ftp.ensemblgenomes.org/pub/current/protists/gff3/"
    sub_dir1_ori=command("curl -s {} | awk '{{print $9}}'".format(ftp_dir)).run_comm(1)
    sub_dirs1=sub_dir1_ori.split();
    for sub_dir1 in sub_dirs1:
        sub_dir2_ori=command("curl -s {}/{}/| awk '{{print $9}}'".format(ftp_dir,sub_dir1)).run_comm(1) 
        sub_dirs2=sub_dir2_ori.split();
        for sub_dir2 in sub_dirs2:
            if re.search("{}.[0-9a-zA-Z]+.[0-9]+.gff3.gz".format(genome), sub_dir2): 
                gff_ftp="{}/{}/{}".format(ftp_dir,sub_dir1,sub_dir2)
            elif not re.search(".gz$", sub_dir2) and sub_dir2!="CHECKSUMS" and sub_dir2!="README":
                file_names_ori=command("curl -s {}/{}/{}/| awk '{{print $9}}'".format(ftp_dir,sub_dir1,sub_dir2)).run_comm(1)
                files_name=file_names_ori.split();
                for file_name in files_name:
                    if re.search("{}.[0-9a-zA-Z]+.[0-9]+.gff3.gz".format(genome),file_name):                    
                        gff_ftp="{}/{}/{}/{}".format(ftp_dir,sub_dir1,sub_dir2,file_name)
    if gff_ftp=="":
        print "ERROR: can not find GFF file for {}".format(genome)
        sys.exit(1)
    else:
        print gff_ftp
    
    #get_gff
    command("wget -N {}".format(gff_ftp)).run_comm(0)
    genome_zip_gff_fname=command("ls {}.*.gff3.gz".format(genome)).run_comm(1)
    print genome_zip_gff_fname+"WWW"
    genome_zip_gff_fname=genome_zip_gff_fname.rstrip()
    command("gunzip -f {}".format(genome_zip_gff_fname)).run_comm(0)
    genome_gff_fname=genome_zip_gff_fname.rstrip(".gz")
       
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
    get_gff_from_genome_name()
    get_chr_map_from_gff()
    genome_gene_dict={}
    pattern_ch="Name=(.*?);"
    pattern="ID=gene:(.*?);"
    file_gff=open(genome_gff_fname,'r')
    for line in file_gff:
        line=line.rstrip("\n")
        if not line.startswith("#") and line.split("\t")[2]=='gene':
            if genome=="Cryptosporidium_hominis" and re.search(pattern_ch,line):
                genome_gene_dict[re.findall(pattern_ch,line)[0]]=line.split("\t")[0] #genome_gene_dict[gene]=chr_acc
            elif re.search(pattern,line):
                genome_gene_dict[re.findall(pattern,line)[0]]=chr_map[line.split("\t")[0]]
                #1    ena    gene    29719    31791    .    +    .    ID=gene:cgd1_130;b --cp
                #LN877947.1    EMBL    gene    1367    2680    .    -    .    ID=gene0;Name=CHUDEA1_10;gbkey=Gene;gene_biotype=protein_coding;locus_tag=CHUDEA1_10 --ch
            else: 
                print "ERROR: there is no gene name in {}".format(line)
                sys.exit(1)


def run_snpEff():         
    global ann_vcf_fpaths  
    ann_vcf_fpaths=[]
    snpEff=prop.get_attrib("snpeff")    
    for vcf_fpath in all_vcf_fpaths:
        vcf_fpath_prefix=os.path.basename(vcf_fpath).rstrip(".vcf")
        ann_vcf_fpath=vcf_fpath_prefix+".ann.vcf"
        ann_vcf_fpaths.append(ann_vcf_fpath)
        command("java -jar {} -c snpEff.config {} {} > {}".format(snpEff,genome,vcf_fpath,ann_vcf_fpath)).run(0);
        
        
def analysis_ann_files():    
    global ann_csv_fname     
    no_dNdS_gene_dict={}
    dNdS_gene_dict={}
    snp_hash=defaultdict(lambda: defaultdict(str))
    ann_pattern='ANN=(.*?)\s+'
    ann_csv_fname="{}.{}.ann.csv".format(prefix,genome)
    fileout=open(ann_csv_fname,'w')
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
    fi.copy_file_to_destdir(ann_csv_fname,outdir)


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

