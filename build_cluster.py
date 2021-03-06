#! /analysis/xin/parasite/bin/python/bin/python 

import matplotlib
matplotlib.use('Agg')
import pandas as pd
from matplotlib import pyplot as plt
import matplotlib
from scipy.cluster import hierarchy
import numpy as np
from scipy.cluster.hierarchy import dendrogram
import sys
from Bio import Phylo
import StringIO
import pylab
import re
import argparse
import re
import sys
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict
import glob
from docutils.nodes import title

def get_args():  
    global fi  
    global prop
    global properties_file
    global genome_name    
    global vcf_file_pattern
    global vcf_files
    global mapping_file
    global image_title
    global prefix
                       
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script build phylogenetic tree and dendragram for the defined group of vcf files from the same genome')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name, only with those obtained from genome_list.txt''',
                        required=True)
    parser.add_argument('-v', '--vcf_file_pattern', type=str, help='''Please provide the vcf files' pattern with the full path,
                                                                  vcf files must ended with ".vcf" ''', required=True)      
    parser.add_argument('-m', '--mapping_file', type=str, help='''Please provide the mapping file path, which contains one column of 
                                                                read_ID from vcf file and one column of its corresponding label on the tree branch,
                                                                otherwise, the read_ID will be labeled on the tree branch''', required=False)
    parser.add_argument('-t', '--title', type=str, help='''Please provide the title of the image''', required=True)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  

    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)
    properties_file=args.properties_file
    prop=properties(properties_file)
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()) and args.genome_name!="cryptosporidium_hominis":
        misc.my_exit("{} is not available, please try another genome".format(args.genome_name))
    if not re.search(".vcf$",args.vcf_file_pattern):
        misc.my_exit("vcf_file_pattern need to end up with .vcf")
    vcf_file_pattern=args.vcf_file_pattern
    vcf_files=glob.glob(vcf_file_pattern)
    fi.check_files_exist(vcf_files)  
    
    # define variables         
    genome_name=args.genome_name
    mapping_file=args.mapping_file
    image_title=args.title
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome_name:",genome_name
    print "vcf_file_pattern:",vcf_file_pattern
    print "mapping_file:",mapping_file
    print "title:",image_title
    print "prefix:",prefix

def get_variant_str():
    outstr=""
    for vcf_file in vcf_files: 
        outstr=outstr+" --variant "+vcf_file
    return outstr

def get_mapping_dict():
    global mapping_dict
    mapping_dict={}
    fhmap=open(mapping_file,'r')
    for line in fhmap:
        line=line.rstrip()
        (read_id,id_label)=getVar(line.split(),[0,1])
        mapping_dict[read_id]=id_label;
    fhmap.close()

def getNewick(node, newick, parentdist, leaf_names):
    if node.is_leaf():
        return "%s:%.2f%s" % (leaf_names[node.id], parentdist - node.dist, newick)
    else:
        if len(newick) > 0:
            newick = "):%.2f%s" % (parentdist - node.dist, newick)
        else:
            newick = ");"
        newick = getNewick(node.get_left(), newick, node.dist, leaf_names)
        newick = getNewick(node.get_right(), ",%s" % (newick), node.dist, leaf_names)
        newick = "(%s" % (newick)
        return newick

def generate_cluster3():
    global cluster_fname
    genomeAnalysisTK=prop.get_attrib("gatk")       
    command("java -jar {} -T CombineVariants -R {} {} -o {} -genotypeMergeOptions UNIQUIFY".
            format(genomeAnalysisTK,genome_fasta,get_variant_str(),prefix)).run_comm(0)
    plink=prop.get_attrib("plink")
    command("{} --vcf {} -cluster --allow-extra-chr -out {}".format(plink,prefix,prefix)).run_comm(0)
    cluster_fname=prefix+".cluster3"  

def create_image_input_file():
    global dendra_fname
    if mapping_file is not None:
        get_mapping_dict()
    dendra_fname=cluster_fname+'_for_dendra' 
    cluster_f_list=[]
    fh_cluster=open(cluster_fname,'r')
    for line in fh_cluster:
        if re.search("\w+",line):
            cluster_f_list.append(line.rstrip())
    col_num=len(cluster_f_list[0].split(' '))
    fh_dendra=open(dendra_fname,'w')    
    ### write first line in fh_dendra    
    fh_dendra.write(identifier)
    for i in range(col_num-1):
        fh_dendra.write(',n')
    ### write the following lines in fh_dendra
    for i in range(len(cluster_f_list)):
        print "line1="+cluster_f_list[i]
        fh_dendra.write("\n")
        outline=""
        cols=cluster_f_list[i].split()
        for i in range(len(cols)):
            if not i==1:
                if i==0: # get the ID's full name
                    cols[i]=cols[i].split(".")[0]    
                    if mapping_file is not None:
                        cols[i]=mapping_dict[cols[i]]
                outline=outline+cols[i]+","
        outline=outline.rstrip(",")       
        fh_dendra.write(outline) 
    fh_dendra.close()
    fh_cluster.close()
    
def generate_images():
    global out_image_fname
    out_image_fname=prefix+".dendrogram.png"
    df = pd.read_csv(
        filepath_or_buffer=dendra_fname,
    )
    df = df.set_index(identifier)
    del df.index.name
    Z=hierarchy.linkage(df, 'ward')
    fig = plt.figure(figsize=(11.69,8.27),dpi = 100)        
    hierarchy.dendrogram(Z, leaf_font_size=8, orientation="left", labels=df.index,get_leaves=True)
    plt.title(image_title,fontsize=18)     
    plt.tight_layout()   
    plt.savefig(out_image_fname)
        
def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global type        
    global genome_fasta
    subdir="cluster"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    if genome_name=="cryptosporidium_hominis":
        genome_fasta=prop.get_attrib("ch_fasta")
    else:    
        genome_fasta=misc.download("fasta",genome_name,"dna")
    for file_in in ([genome_fasta,mapping_file]):
        if file_in is not None:
            fi.copy_file_to_destdir(file_in,indir)
    fi.change_dir(workdir)    

def execute():
    print "executing..."
    global identifier
    identifier="ID"
    generate_cluster3()
    create_image_input_file()
    generate_images() 
    
def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir(workdir+"/"+out_image_fname,outdir)

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
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
   

