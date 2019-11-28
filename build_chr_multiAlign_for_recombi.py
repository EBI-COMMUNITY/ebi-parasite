#!/usr/bin/env python3  

## version 1, 28 Nov 2019 by Xin Liu

import re
import argparse
import sys
import glob
import os
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from collections import defaultdict

def get_args():  
    global fi  
    global prop
    global properties_file
    global genome_name    
    global bam_file_pattern
    global bam_files
    global mapping_file
    global prefix
    global bam_key_pattern
                       
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script build all individual chromosome multiple alignment for recombination')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', 
                        required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name available 
                                                                 in genome_list.txt only''', 
                        required=True)
    parser.add_argument('-bp', '--bam_file_pattern', type=str, help='''Please provide the bam files' pattern 
                                                                       with the full path, ending with .bam, with runID 
                                                                       in the bam file name''', 
                        required=True)
    parser.add_argument('-m', '--mapping_file', type=str, help='''Please provide the mapping file path, containing one 
                                                                  column of the runID and the other column is the expression
                                                                  displayed in the multiple alignment file description line''', 
                        required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', 
                        required=True)  

    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)
    properties_file=args.properties_file
    prop=properties(properties_file)
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()):
        misc.my_exit("{} is not available, please try another genome".format(args.genome_name))
    if not re.search(".bam$",args.bam_file_pattern):
        misc.my_exit("bam_file_pattern need to end up with .bam")
    bam_file_pattern=args.bam_file_pattern
    bam_files=glob.glob(bam_file_pattern)
    bam_key_pattern="[A-Z]RR\d{6,}"
    for bam_file in bam_files:
        if not re.search(bam_key_pattern, bam_file):
            misc.my_exit("There is no runID in the bam file {}".format(bam_file))
    fi.check_files_exist(bam_files)  
    
    # define variables         
    genome_name=args.genome_name
    mapping_file=args.mapping_file
    prefix=args.prefix   
    
    print ("properties_file:",properties_file)
    print ("genome_name:",genome_name)
    print ("bam_file_pattern:",bam_file_pattern)
    print ("mapping_file:",mapping_file)
    print ("prefix:",prefix)
   
def get_chr_fasta():
    # write fasta file for each chromosome, with '>bam_name' 
    global chr_fastas
    chr_fastas=[]
    chr_bam_fasta_dict = defaultdict(lambda: defaultdict(str))
    map_dict = {}
    # get mapping relationship runID <-> sample name
    if mapping_file is not None:
        fhmap = open(mapping_file, 'r')
        for line in fhmap:
            line = line.rstrip()
            [run, genotype] = line.split()
            map_dict[run] = genotype
        fhmap.close()
    # get all fasta files from pilon
    # chr_bam_fasta_dict will hold the fasta sequence for each chr from each sample (bam)
    fasta_files = glob.glob("{}/{}_*.fasta".format(workdir, prefix))
    for fasta_file in fasta_files:
        bam_key = re.findall("({}).*?.fasta".format(bam_key_pattern), fasta_file)[0]
        file_str = open(fasta_file, 'r').read().strip("\n")
        fasta_segs = file_str.split(">")
        for fasta_seg in fasta_segs:
            if re.search("\w", fasta_seg): # to skip the very beginning, which is nothing
                chrom = fasta_seg.split("\n")[0].rstrip("_pilon")
                chr_bam_fasta_dict[chrom][bam_key] = "\n".join(fasta_seg.split("\n")[1:])
    for chrom in sorted(chr_bam_fasta_dict.keys()):
        out_fasta_file = "{}/{}_{}.fa".format(workdir,prefix,chrom)
        chr_fastas.append(out_fasta_file)
        fh_out = open(out_fasta_file, 'w')
        for bam_key in sorted(chr_bam_fasta_dict[chrom].keys()):
            if bam_key in map_dict.keys():
                seg_desc = map_dict[bam_key]
            else:
                seg_desc = bam_key
            fh_out.write(">{}\n".format(seg_desc))
            fh_out.write(chr_bam_fasta_dict[chrom][bam_key]+"\n")
        fh_out.close()

def rename_xmfa_iso(mAlign_map_file, out_xmfa_file):
    getVar = lambda searchList, ind: [searchList[i] for i in ind]
    tmp_file = "{}/pmauve.tmp".format(workdir)
    replace_dict = {}
    fhmap = open(mAlign_map_file, 'r')
    for line in fhmap:
        line = line.rstrip()
        [col1, isoname] = line.split()
        replace_dict[col1+":"] = isoname+":"
    fhmap.close()
    command("mv {} {}".format(out_xmfa_file, tmp_file)).run_comm(0)
    fhin = open(tmp_file, 'r')
    fhout = open(out_xmfa_file, 'w')
    for line in fhin:
        line = line.rstrip()
        if re.search("^> (\d+:)", line):
            str_to_be_replaced = re.findall("^> (\d+:)", line)[0]
            line = line.replace(str_to_be_replaced, replace_dict[str_to_be_replaced])
        fhout.write(line + "\n")
    fhin.close()
    fhout.close()
    command("rm {}".format(tmp_file)).run_comm(0)

def build_multiAlign_file():
    global mAlign_map_file
    global out_xmfa_files
    for bam_file in bam_files:
        # run pilon, and this will create a fasta file
        #command("samtools index {} {}.bai".format(bam_file,bam_file)).run_comm(0)
        pilon_out_file = "{}/{}_{}".format(workdir, prefix, os.path.basename(bam_file).rstrip(".bam"))
        #command("pilon --genome {} --bam {} --output {} --vcf".format(REF_FASTA, bam_file, pilon_out_file)).run_comm(0)
    get_chr_fasta()
    #chr_fastas=["forSimone_LN877948.fa","forSimone_LN877949.fa","forSimone_LN877950.fa","forSimone_LN877951.fa","forSimone_LN877952.fa","forSimone_LN877953.fa","forSimone_LN877954.fa","forSimone_LN877947.fa"]
    for chr_fasta in chr_fastas:
        chr_fasta_prefix=os.path.basename(chr_fasta).replace(".fa","")
        out_xmfa_file=chr_fasta_prefix+".xmfa"
        command("progressiveMauve --output={} {} ".format(out_xmfa_file,chr_fasta)).run_comm(0)
        mAlign_map_file = "{}/{}_multiAlign.map".format(workdir, chr_fasta_prefix)
        cmd="grep '>' {} | cat -n | sed -e \"s/>//g\" | awk '{{split($0, a ,\" \");print (a[1]\" \"a[2])}}' > {}".format(chr_fasta, mAlign_map_file)
        command(cmd).run_comm(0)
        rename_xmfa_iso(mAlign_map_file, out_xmfa_file)
        out_xmfa_files.append(out_xmfa_file)
  

def initiate():
    print ("initiating...")
    global workdir
    global indir
    global outdir
    global type        
    global REF_FASTA
    subdir="recombination"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)    
    for file_in in ([REF_FASTA,mapping_file]):
        if file_in is not None:
            fi.copy_file_to_destdir(file_in,indir)
    fi.change_dir(workdir)    

def execute():
    print ("executing...")
    build_multiAlign_file()
 
def post_process():
    print ("post_processing...")
    for out_xmfa_file in out_xmfa_files:
        fi.copy_file_to_destdir(workdir+"/"+out_xmfa_file,outdir)


if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    misc=misc()
    get_args()
    print ("\n","Properties attributes:")
    print (prop.__dict__)
    
    #run_blast the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")

