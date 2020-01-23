#! /usr/bin/env python3 

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
from utilities import command
from utilities import fileutils
from utilities import properties
import os
import re
import glob

def get_args():    
    global prop
    global properties_file
    global prefix
    global fi
  
    fi=fileutils()
 
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script creates multiQC html file using fastqc, 
                                                    bcftools, snpEff, QUAST, and QualiMap 
                                                    output files''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file, which including the paths of workdir', 
                        required=True)    
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file', 
                        required=True)  
    # check args
    args = parser.parse_args()
    fi.check_exist(args.properties_file)
    # define variables
    properties_file=args.properties_file
    prop=properties(properties_file)
    prefix=args.prefix
    print ("properties_file:",str(properties_file))    
    print ("prefix:",prefix)
  
def multiQC():  
    global multiQC_in_files
    global multiQC_files_str
    multiQC_files_str = "multiqc -f"
    assembly_p = "assembly/qc/{}_*/report.tsv".format(prefix)
    fastQC_p = "reference_mapping/qc/{}/*fastqc.zip".format(prefix)
    #bowtie2_p = "reference_mapping/qc/{}/*.log".format(prefix) 
    picard_p = "reference_mapping/qc/{}/*.bam".format(prefix)
    qualimap_p1 = "reference_mapping/qc/{}/*/genome_results.txt".format(prefix)
    qualimap_p2 = "reference_mapping/qc/{}/*/raw_data_qualimapReport".format(prefix)
    bcf_p = "/snp/qc/{}/*.bcf_stats".format(prefix)
    snpEff_p = "dNdS/qc/{}/*.txt".format(prefix) 
    for qc_p in (assembly_p, fastQC_p, picard_p, qualimap_p1, qualimap_p2, bcf_p, snpEff_p):
        multiQC_files_str = get_multiQC_files_str_and_cpToIn(qc_p, multiQC_files_str)
    if multiQC_files_str != "multiqc -f":
        print ("multiQC_files_str="+multiQC_files_str)
        cmd = "{} -o {} --filename {}.multiQC -v".format(multiQC_files_str, workdir, prefix)
        command(cmd).run_comm(0)
 
def get_multiQC_files_str_and_cpToIn(qc_p, multiQC_files_str):
    qc_fullP_p = "{}/../{}".format(workdir, qc_p)
    qc_files = glob.glob(qc_fullP_p)
    if len(qc_files) != 0:
        multiQC_files_str += " " + qc_fullP_p
    return multiQC_files_str

def initiate():
    global workdir  
    global indir
    global outdir
    print ("initiating...")
    subdir="qc_report"
    workdir=prop.workdir+"/"+subdir
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(workdir)
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)

def execute():
    print ("executing...")
    fi.change_dir(workdir)
    multiQC()

def post_process():
    print ("post_processing...")
    fi.copy_file("{}.multiQC.html".format(prefix),outdir)

if __name__ == '__main__':
    get_args()
    global multiQC_files_str
    print ("\n","Properties attributes:")
    print (prop.__dict__)
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")
     
