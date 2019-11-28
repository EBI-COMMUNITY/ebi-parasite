#! /usr/bin/env python3

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
import os
import glob
import re
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict

def get_args():
    global PROP
    global PROPERTIES_FILE
    global genome_name
    global vcf_files
    global vcf_files_str
    global gvcf_file
    global genotype_file
    global prefix
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script for creating images: phylogeny tree, PCA, heatmap, upset, 
                                                 and SNP distribution on chromosome for vcf files''')
    parser.add_argument('-p', '--properties_file', type=str, help='''Please provide the properties file, which including 
                                                                     the paths of samtools and bcftools and workdir''',
                                                             required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                                                                 which is provided by genome_list.txt''', 
                                                             required=True)
    parser.add_argument('-v', '--vcf_files', type=str, help='Please provide one or multiple vcf files', 
                                                             required=True)
    parser.add_argument('-gv', '--gvcf_file', type=str, help='Please provide one  concat gvcf file',
                                                             required=True)
    parser.add_argument('-gt', '--genotype_file', type=str, help='''The file need to contain 3 cols: runID, isolate name, and genotype''', 
                                                             required=True)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                                                             required=True)  
    
    args=parser.parse_args() 
    
    # check args
    FI.check_exist(args.properties_file)
    PROPERTIES_FILE = args.properties_file
    PROP = properties(PROPERTIES_FILE)
    MISC.check_genome_avl(PROP.get_attrib("available_genomes"), args.genome_name)
    FI.check_files_exist(glob.glob(args.vcf_files))     
    FI.check_exist(args.gvcf_file)
    FI.check_exist(args.genotype_file)
    genome_name = args.genome_name
    vcf_files = glob.glob(args.vcf_files)
    if len(vcf_files)==0:
        sys.exit("Input vcf files not exist")
    vcf_files_str=""
    for vcf_file in vcf_files:
        MISC.get_runID(vcf_file)
        vcf_files_str += os.path.abspath(vcf_file)+" "
    gvcf_file=args.gvcf_file    
    genotype_file=args.genotype_file
    prefix=args.prefix      

    # print args
    print ("properties_file:", PROPERTIES_FILE)
    print ("genome name:", genome_name)
    print ("vcf_files:",vcf_files_str)
    print ("gvcf_file:",gvcf_file)
    print ("genotype_file:",genotype_file)
    print ("prefix:",prefix)


def out_files_write(out_str):
    fhout_intersect.write(out_str)
    fhout_jaccard.write(out_str)

def get_intersect():
    global intersect_fPath
    global jaccard_fPath
    global fhout_intersect
    global fhout_jaccard
    vcf_dict={}
    intersect_fPath=workdir+"/intersect.matrix"
    jaccard_fPath=workdir+"/jaccard.matrix"
    fhout_intersect=open(intersect_fPath,'w')
    fhout_jaccard=open(jaccard_fPath,'w')
    out_files_write("name")
    for vcf in vcf_files:
        vcf_dict[vcf]=MISC.get_runID(vcf)
        out_files_write(" "+vcf_dict[vcf])
    out_files_write("\n")
    for vcf1 in vcf_files:
        out_files_write(vcf_dict[vcf1])
        for vcf2 in vcf_files:
            (intersect,jaccard)=getVar(command("bedtools jaccard -a {} -b {} |cut -f1,3|grep -v jaccard".format(vcf1,
                                       vcf2)).run_comm(1).decode("utf-8").rstrip().split(),[0,1]) #overlaps
            fhout_intersect.write(" "+intersect)
            fhout_jaccard.write(" "+jaccard)            
        out_files_write("\n")
    fhout_intersect.close()
    fhout_jaccard.close()
    change_fileName_to_isoName(intersect_fPath)
    change_fileName_to_isoName(jaccard_fPath)

def get_mirror_dict():
    global mirror_dict_runID
    global mirror_dict_vcf
    mirror_dict_runID={}
    mirror_dict_vcf={}
    for vcf_file in vcf_files:
        runID=MISC.get_runID(vcf_file)    
        sample=MISC.get_samples_by_runIDs(genotype_file)[runID]
        mirror_dict_vcf[vcf_file] = sample
        mirror_dict_runID[runID]=sample

def change_fileName_to_isoName(infile):
    with open(infile,'r') as infh:
        infile_str=infh.read()
        for vcf_file in (mirror_dict_vcf.keys()):
            infile_str=infile_str.replace(vcf_file,mirror_dict_vcf[vcf_file])
        for runID in (mirror_dict_runID.keys()):
            infile_str=infile_str.replace(runID,mirror_dict_runID[runID])
    outfile_str=infile_str    
    with open(infile,'w') as outfh:
        outfh.write(outfile_str)

def parse_mat():
    headline=command("head -n 1 {}".format(upset_mat_file)).run_comm(1).decode("utf-8").rstrip()
    eles=headline.split()
    new_headline=""
    for ele in eles[:5]:
        new_headline+=ele+"\t"
    for ele in eles[5:]:
        ele=re.sub("^.*/", "", ele)
        new_headline += ele+"\t"
    new_headline = new_headline.rstrip("\t")
    command("mv {} {}.tmp".format(upset_mat_file, upset_mat_file)).run_comm(0)
    command("sed '1s/.*/{}/' {}.tmp >{}".format(new_headline, upset_mat_file, upset_mat_file)).run_comm(0) # remove dir of each sample from the head line
    change_fileName_to_isoName(upset_mat_file)

def create_upset_matrix():
    global upset_mat_file
    upset_mat_file = "{}_upset.mat".format(prefix)
    command("infoseq -sequence {} -only -name -length -outfile {}/genome.size -nohead -auto".format(REF_FASTA,
                                                                                                    workdir)).run_comm(0)
    command("bedtools multiinter -g {}/genome.size -emtpy -header -i {} >{}/{}".format(workdir, vcf_files_str,
                                                                                      workdir, upset_mat_file)).run_comm(0)
    change_fileName_to_isoName(upset_mat_file)

def call_R():
    gvcf_zip_file = gvcf_file+".gz"
    command("gzip -c {} >{}".format(gvcf_file,gvcf_zip_file)).run_comm(0)
    R_script=PROP.r_script
    cmd="{} -t {} -v {} -f {} -g {} --mUpset {} --mInter {} --mJacc {}".format(R_script, genotype_file, gvcf_zip_file, 
                                                                            REF_FASTA, REF_GFF, 
                                                                            upset_mat_file, intersect_fPath, jaccard_fPath) 
    command(cmd).run_comm(0)    
    phylo_PCA_upset_heatmap_file="{}_phylo_PCA_upset_heatmap.pdf".format(prefix)
    out_pdf_files.append(phylo_PCA_upset_heatmap_file)
    command("mv Rplots.pdf {}".format(phylo_PCA_upset_heatmap_file)).run_comm(0)
    
def get_variant_str():
    outstr=""
    for vcf_file in vcf_files:
        outstr=outstr+" --variant "+vcf_file
    return outstr

def run_vcf_merge():
    global merged_vcf
    variants_str=get_variant_str()
    merged_vcf=workdir+"/"+prefix+"_merged.vcf"
    command("gatk -T CombineVariants -R {} {} -o {} -genotypeMergeOptions UNIQUIFY"
             .format(REF_FASTA,variants_str,merged_vcf)).run_comm(0)

def call_R_dist_on_chr():
    run_vcf_merge()
    chr_vcf_head = command("grep '^#' {}".format(merged_vcf)).run_comm(1).decode("utf-8")
    chr_gff_head = command("grep '^#' {}".format(REF_GFF)).run_comm(1).decode("utf-8").rstrip("\n")
    with open (REF_FASTA,"r") as fh_fasta:
        ref_fasta_str=fh_fasta.read()
    splited_ref_fastas = ref_fasta_str.split(">")
    for splited_ref_fasta in splited_ref_fastas: #get chr from fasta file
        if re.search("\w+",splited_ref_fasta):
            chrom = splited_ref_fasta.split()[0].split()[0]
            chr_vcf="{}_{}.vcf".format(prefix,chrom)
            chr_fasta="{}_{}.fasta".format(prefix,chrom)
            chr_gff="{}_{}.gff".format(prefix,chrom)
            with open (chr_vcf,"w") as fh_chr_vcf:
                chr_vcf_body_str=command("grep '^{}' {}".format(chrom,merged_vcf)).run_comm(1).decode("utf-8")
                fh_chr_vcf.write(chr_vcf_head)
                fh_chr_vcf.write(chr_vcf_body_str)
            with open (chr_fasta,"w") as fh_chr_fasta:
                fh_chr_fasta.write(">{}\n".format(chrom))
                fh_chr_fasta.write("\n".join(splited_ref_fasta.split("\n")[1:]))   
            with open (chr_gff,"w") as fh_chr_gff:
                chr_gff_body_str=command("grep '{}' {} |grep -v '^#'".format(chrom,REF_GFF)).run_comm(1).decode("utf-8").rstrip()
                fh_chr_gff.write(chr_gff_head+"\n")
                fh_chr_gff.write(chr_gff_body_str)
            chr_vcf_gz=chr_vcf+".gz"
            command("gzip -c {} > {}".format(chr_vcf,chr_vcf_gz)).run_comm(0)
            cmd="{} -s {} -v {} -f {} -g {}".format(PROP.r_snp_dist_on_chr, chrom, chr_vcf_gz, chr_fasta, chr_gff)
            command(cmd).run_comm(0)
            snp_dist_on_chr="{}_snp_dist_on_{}.pdf".format(prefix, chrom)
            out_pdf_files.append(snp_dist_on_chr)
            command("mv Rplots.pdf {}".format(snp_dist_on_chr)).run_comm(0)

def initiate():
    print("initiating...")
    global indir
    global outdir
    global workdir
    global subdir
    subdir="visualization"
    workdir=prop.workdir+"/"+subdir
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    FI.create_processing_dir(workdir)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    #for vcf_file in vcf_files:
        #FI.copy_file_to_destdir(vcf_file,indir)
    #FI.copy_file_to_destdir(gvcf_file,indir)

def execute():    
    print("executing...")    
    global REF_FASTA
    global REF_GFF
    global out_pdf_files
    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)
    REF_GFF = "{}/../ref/{}.gff".format(workdir, genome_name)
    out_pdf_files = []
    FI.change_dir(workdir)    
    GFF3 = "{}/../../ref/{}.gff".format(workdir, genome_name)
    get_mirror_dict()
    get_intersect()
    create_upset_matrix()
    call_R()
    call_R_dist_on_chr()
    
def post_process():
    print("post_processing...")
    for out_file in out_pdf_files:
        FI.copy_file_to_destdir(out_file,outdir) 

if __name__ == '__main__':
    global prop
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()

    get_args()
    prop=properties(PROPERTIES_FILE)
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

