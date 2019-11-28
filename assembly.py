#! /usr/bin/env python3

# version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
import os
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc

def get_args():
    global properties_file
    global prop
    global genome_name
    global fastq1
    global fastq2
    global prefix_ori
    global runID
    global mapping_file
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script assembles short reads by using spades and 
                                                    provides statistics summary based on the result assemblies 
                                                    by using QUAST''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                                                                 which need to be in genome_list.txt''', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the forward fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the reverse fastq file', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                         required=True)  
    parser.add_argument('-m', '--mapping_file', type=str, help='''if map file was provided, whose first column is runID
                                                                  and the second column is sample name, seperated by "\t",
                                                                  sample name will be included in the output files' names''', 
                         required=False)

    # check args
    args = parser.parse_args()
    FI.check_files_exist([args.properties_file,args.fastq1]) 
    properties_file=args.properties_file
    prop=properties(properties_file)
    MISC.check_genome_avl(prop.get_attrib("available_genomes"), args.genome_name)
    runID=MISC.get_runID(args.fastq1)
    if args.fastq2 is not None:
        FI.check_exist(args.fastq2)
        MISC.get_runID(args.fastq2)
    if args.mapping_file is not None:
        FI.check_exist(args.mapping_file)

    # define variables
    genome_name = args.genome_name
    fastq1=args.fastq1
    fastq2=args.fastq2
    prefix_ori=args.prefix    
    mapping_file=args.mapping_file    
    print ("properties_file:",properties_file)
    print ("genome_name:", genome_name)
    print ("fastq1:",fastq1)
    print ("fastq2:",fastq2)
    print ("prefix:",prefix_ori)
    print ("mapping_file:",mapping_file)

def set_prefix():
    global prefix
    global sample_runID
    if mapping_file is None: 
        prefix = "{}_{}".format(prefix_ori, runID)
    else:
        sample_runID="{}_{}".format(MISC.get_samples_by_runIDs(mapping_file)[runID],runID)
        prefix = "{}_{}".format(prefix_ori, sample_runID)

def assembly(fastqfiles):            
    global scaffold_fasta ## for QC
    if mapping_file is None:
        scaffold_fasta = "{}.fasta".format(runID)
    else:
        scaffold_fasta = "{}.fasta".format(sample_runID)
    if len(fastqfiles)==2:
        command("spades.py -1 {} -2 {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],fastqfiles[1],workdir)).run_comm(0)
    elif len(fastqfiles)==1:
        command("spades.py -s {} --careful --cov-cutoff auto -o {}".format(fastqfiles[0],workdir)).run_comm(0)
    command("cp scaffolds.fasta {}".format(scaffold_fasta)).run_comm(0)

def qc_by_quast():
    global qc_files
    qc_files=["report.txt","report.tsv","report.pdf","report.html"]
    REF_FASTA = "{}/../../ref/{}.fasta".format(workdir, genome_name)
    GFF3 = "{}/../../ref/{}.gff".format(workdir, genome_name)
    command("quast {} -r {} -g {} -o {}".format(scaffold_fasta,REF_FASTA,GFF3,workdir)).run_comm(0)

def multiqc():
    command("multiqc -f report.tsv -o {} --filename {}.multiQC -v".format(workdir,prefix)).run_comm(0)

def initiate():
    print ("initiating...")
    global indir
    global outdir
    global qcdir
    global workdir
    subdir="assembly"
    workdir="{}/{}/{}".format(prop.workdir,subdir,runID)
    workdir_root="{}/{}".format(prop.workdir,subdir)
    indir=workdir_root+"/in/"
    outdir=workdir_root+"/out/"
    qcdir=workdir_root+"/qc/"+prefix_ori
    FI.create_processing_dir(workdir)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    FI.create_processing_dir(qcdir)
    FI.copy_file_to_destdir(fastq1,indir)  
    if fastq2 is not None:
        FI.copy_file_to_destdir(fastq2,indir)

def execute():
    print ("executing...")
    FI.change_dir(workdir)
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    set_prefix()
    assembly(fastqfiles)
    qc_by_quast()
    multiqc()

def post_process():
    print ("post_processing...")
    out_files = ["scaffolds.fasta","spades.log"]
    for out_file in out_files:
        FI.copy_file_add_prefix(out_file,outdir,prefix+"_")
    for qcfile in qc_files:
        FI.copy_file_to_destdir(qcfile,qcdir)
    command("cp -p {}.multiQC*.html {}".format(prefix,qcdir)).run_comm(0)
    
if __name__ == '__main__':
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()

    get_args()

    print ("\n","Properties attributes:\n",prop.__dict__)
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")     
