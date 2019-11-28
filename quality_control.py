#! /usr/bin/env python3 

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
import os
import re

def get_args():    
    global properties_file
    global prop
    global fastq1
    global fastq2
    global prefix
    global if_dedup
    global runID
    global mapping_file
 
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script provides quality control on fastq files using trim_galore 
                                                  and remove duplicated reads using clumpify. The quality of the original
                                                  and filtered reads are monitored by fastQC and summarized by multiQC''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file.', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the second fastq file.', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                         required=True)  
    parser.add_argument('-de', '--dedup', action='store_true',help='if remove all the exact read duplications', default=False)
    parser.add_argument('-m', '--mapping_file', type=str, help='''if map file was provided, whose first column is runID
                                                                  and the second column is sample name, seperated by "\t",
                                                                  sample name will be included in fastQC and multiQC
                                                                  output file names''', required=False)

    # check args
    args = parser.parse_args()
    FI.check_files_exist([args.properties_file, args.fastq1])
    runID=MISC.get_runID(args.fastq1)
    if args.fastq2 is not None:
        FI.check_exist(args.fastq2)
        MISC.get_runID(args.fastq2)

    # define variables
    properties_file=args.properties_file
    prop=properties(properties_file)
    fastq1=args.fastq1
    fastq2=args.fastq2
    prefix=args.prefix
    if_dedup=args.dedup
    mapping_file=args.mapping_file
    
    print ("properties_file:",properties_file)
    print ("fastq1:",fastq1)
    print ("fastq2:",fastq2)
    print ("prefix:",prefix)
    print ("dedup:",if_dedup)
    print ("mapping_file:",mapping_file)
   
def fastQC(fq_in):
    if mapping_file is not None:
        fq="{}/{}".format(workdir,os.path.basename(fq_in).replace(runID,sample_runID))
        command("cp {} {}".format(fq_in,fq)).run_comm(0)
    else:
        fq=fq_in
    command("gzip -c {} > {}.gz".format(fq, fq)).run_comm(0)
    command("fastqc -o {} --noextract -f fastq {}.gz".format(workdir, fq)).run_comm(0)
    fastqc_out_file_name="{}_fastqc.zip".format(os.path.basename(fq).rstrip(".fq.gz"))
    return fastqc_out_file_name 

def run_trim_galore(fastqfiles):            
    global fqout1
    global fqout2
    global tg_out_files
    tg_out_files = []
    fastq1_name=os.path.basename(fastq1)
    fastq1_name_base=fastq1_name.rstrip(".fastq")
    report1=os.path.join(workdir,fastq1_name+"_trimming_report.txt")
    if fastq2 is None:
        fqout1=os.path.join(workdir,fastq1_name_base+"_trimmed.fq")
        tg_out_files=[fqout1, report1]        
    else:
        fqout1=os.path.join(workdir,fastq1_name_base+"_val_1.fq")
        fastq2_name=os.path.basename(fastq2)
        fastq2_name_base=fastq2_name.rstrip(".fastq")
        fqout2=os.path.join(workdir,fastq2_name_base+"_val_2.fq")
        report2=os.path.join(workdir,fastq2_name+"_trimming_report.txt")
        tg_out_files=[fqout1,report1,fqout2,report2]
    #pair fastq files    
    if len(fastqfiles)==2:
        command( "trim_galore --paired -q 20 "+fastqfiles[0]+" "+fastqfiles[1]).run_comm(0)
    #single fastq files    
    elif len(fastqfiles)==1:
        command("trim_galore -q 20 "+fastqfiles[0]).run_comm(0)

def deduplication():
    global deduped_fqs
    fq_dedup_out1=fqout1.replace(".fq",".dedup.fq")
    command("clumpify in={} out={} dedupe=t".format(fqout1,fq_dedup_out1)).run_comm(0)
    deduped_fqs=[fq_dedup_out1]
    if fastq2 is not None:
        fq_dedup_out2=fqout2.replace(".fq",".dedup.fq")    
        command("clumpify in={} out={} dedupe=t".format(fqout2,fq_dedup_out2)).run_comm(0)
        deduped_fqs.append(fq_dedup_out2)

def run_QC():
    global fastqc_out_files
    global fastqc_out_files_str
    global multiqc_out_prefix
    global sample_runID
    ## define variables
    if mapping_file is not None:
        sample=MISC.get_samples_by_runIDs(mapping_file)[runID]
        sample_runID="{}_{}".format(sample, runID)
    else:
        sample_runID=runID
    fastqc_out_files = []
    fastqc_out_files_str = ""
    if fastq2 is not None:
        infiles = [fq_ori1,fqout1,fq_ori2,fqout2]
    else:
        infiles = [fq_ori1,fqout1]
    ## fastQC
    for infile in infiles:
        fastqc_out_files.append(fastQC(infile))
    for fastqc_out_file in fastqc_out_files:
        fastqc_out_files_str+=fastqc_out_file+" "
    fastqc_out_files_str=fastqc_out_files_str.rstrip(" ")
    ## multiQC
    multiqc_out_prefix=sample_runID+".multiQC"
    command("multiqc -f {} -o {} --filename {} -v".format(fastqc_out_files_str,workdir,multiqc_out_prefix)).run_comm(0)
  
def initiate():
    print ("initiating...")
    global workdir
    global indir
    global outdir
    global qcdir
    global fq_ori1
    global fq_ori2
        
    subdir="quality"
    workdir=prop.workdir+"/"+subdir
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    qcdir=workdir+"/qc/"+prefix
    FI.create_processing_dir(workdir)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    FI.create_processing_dir(qcdir)
    fq_ori1 = "{}/{}_ori.fq".format(indir,os.path.basename(fastq1.rstrip(".fastq")))
    FI.copy_file(fastq1, fq_ori1)
    if fastq2 is not None:
        fq_ori2 = "{}/{}_ori.fq".format(indir, os.path.basename(fastq2.rstrip(".fastq")))
        FI.copy_file(fastq2, fq_ori2)

def execute():
    print ("executing...")
    FI.change_dir(workdir)
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_trim_galore(fastqfiles)
    if if_dedup:
        print("removing exact duplicates...")  
        deduplication()
    run_QC()

def post_process():
    print ("post_processing...")
    fastqc_postfix = "_fastqc.zip"
    for tg_out_file in tg_out_files:
        FI.copy_file_add_prefix(tg_out_file,outdir,prefix+"_")
    for fastqc_out_file in fastqc_out_files:
        FI.copy_file_to_destdir(fastqc_out_file,qcdir)
        FI.copy_file_to_destdir(fastqc_out_file.replace(".zip", ".html"),qcdir)
    command("cp -p {}.html {}".format(multiqc_out_prefix,qcdir)).run_comm(0)
    if if_dedup:
        for deduped_fq in deduped_fqs:
            FI.copy_file_add_prefix(deduped_fq,outdir,prefix+"_")

if __name__ == '__main__':
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()
    
    get_args()
    
    print ("\n","Properties attributes:\n", prop.__dict__)
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")
     
