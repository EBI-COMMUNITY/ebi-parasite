#! /usr/bin/env python3

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import sys
import os
import re
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc

def get_args():
    global properties_file
    global prop
    global mapping_tool
    global genome_name
    global fastq1
    global fastq2
    global platform
    global dna_library
    global if_dedup
    global if_recom
    global prefix_ori
    global mapping_file
    global runID

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script mapping short reads to reference genomes
                                                    using BWA or bowtie2, with the following statistics 
                                                    summary files created by fastQC, qualiMap, and multipleQC''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-t', '--mapping_tool', type=str, help='Please choose mapping tool, bwa or bowtie2', required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name
                                                                 which need to be in genome_list.txt''', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='''Please provide the paired forward fastq file or 
                                                              the single fastq file, containing runID''', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='''Please provide the paired reverse fastq file, 
                                                              containing runID''', required=False)
    parser.add_argument('-f', '--platform', type=str, help='Please provide the plat_form', required=True)
    parser.add_argument('-l', '--dna_library', type=str, help='Please provide the DNA library', required=True)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the subdir for holding the output file', required=True)  
    parser.add_argument('-de', '--dedup', action='store_true',help='if true, remove duplications after mapping using samtools', default=False)
    parser.add_argument('-recom', '--recombination', action='store_true',help='''if true, mapping with bowtie2 with be use local 
                                                                               alignment, which is for recombination analysis''', default=False)
    parser.add_argument('-m', '--mapping_file', type=str, help='''if map file was provided, sample name will be 
                                                                  included in the output file name''', required=False) 
    args=parser.parse_args()
    # check args
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
    if not args.mapping_tool=="bwa" and not args.mapping_tool=="bowtie2":
        sys.exit("only bwa and bowtie2 are available for the mapping_tool")
    
    # define variables        
    mapping_tool = args.mapping_tool
    genome_name = args.genome_name
    fastq1 = args.fastq1
    fastq2 = args.fastq2  
    platform = args.platform
    dna_library = args.dna_library
    prefix_ori = args.prefix   
    if_dedup = args.dedup
    if_recom = args.recombination
    mapping_file=args.mapping_file    
    print ("properties_file:",properties_file  )
    print ("genome_name:", genome_name)
    print ("mapping_tool:",mapping_tool)
    print ("fastq1:",fastq1)
    print ("fastq2:",fastq2)
    print ("platform:",platform)
    print ("dna_library:",dna_library)
    print ("prefix:",prefix_ori)
    print ("dedup:",if_dedup)
    print ("recombination:",if_recom)
    print ("mapping_file:",mapping_file)

def set_prefix():
    global prefix
    global sample_runID
    if mapping_file is None:
        sample_runID=runID
        prefix = "{}_{}".format(prefix_ori,runID)
    else:
        sample=MISC.get_samples_by_runIDs(mapping_file)[runID]
        sample_runID=sample+"_"+runID
        prefix = "{}_{}".format(prefix_ori,sample_runID)
        
def run_mapping(fastqfiles):  
    global REF_FASTA          
    global sam_out
    global bam_sorted
    global fq1
    global bowtie2_log
    global bowtie2_log_for_multiqc
    sam_out="{}/{}.sam".format(workdir,prefix)
    bam_sorted="{}/{}.bam".format(workdir,prefix)
    bowtie2_log=""
    
    #define command         
    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)
    fq1=fastqfiles[0]
    if mapping_tool=="bwa":
        fq1_sai=prefix+'.fq1.sai'
        command("bwa aln -t 12 {} {} > {}".format(REF_FASTA,fq1,fq1_sai)).run_comm(0)
        if len(fastqfiles)==2:
            fq2=fastqfiles[1]
            fq2_sai=prefix+'.fq2.sai'        
            command("bwa aln -t 12 {} {} > {}".format(REF_FASTA,fq2,fq2_sai)).run_comm(0)
            command("bwa sampe {} {} {} {} {} | samtools sort > {}".format(REF_FASTA,fq1_sai,fq2_sai,fq1,fq2,bam_sorted)).run_comm(0)
        elif len(fastqfiles)==1:
            command("bwa samse {} {} {} | {} sort > {}".format(REF_FASTA,fq1_sai,fq1,samtools,bam_sorted)).run_comm(0)
    elif mapping_tool=="bowtie2":
        ref_fasta_key=re.findall(".*/(.*?).fa",REF_FASTA)[0]
        ref_fasta_prefix=workdir+"/../ref/"+ref_fasta_key        
        bowtie2_log="{}/{}_bowtie2.log".format(workdir,prefix)
        bowtie2_log_for_multiqc=sample_runID+".log"
        local_align_str=""
        if if_recom:
            local_align_str="--local"
        if len(fastqfiles)==1:
            command("bowtie2 {} -p 8 -x {} -U {} -S {} >& {}".format(local_align_str,ref_fasta_prefix,fq1,sam_out,bowtie2_log)).run_comm(0)
        if len(fastqfiles)==2:
            command("bowtie2 {} -p 8 -x {} -1 {} -2 {} -S {} >& {}".format(local_align_str,ref_fasta_prefix,fq1,fastqfiles[1],sam_out,bowtie2_log)).run_comm(0)
        command("samtools sort {} >{}".format(sam_out, bam_sorted)).run_comm(0)   
        FI.copy_file(bowtie2_log, bowtie2_log_for_multiqc)

def add_group(fastqfiles):
    global grouped_bam
    global picard
    grouped_bam="{}/{}_grouped.bam".format(workdir,sample_runID)
    fq_id_pat="^\@(.*?)\."
    first_line_fq1=command("head -n 1 "+fq1).run_comm(1).decode("utf-8").rstrip()
    fq_id=re.findall(fq_id_pat,first_line_fq1)[0]
    cmd_str="picard AddOrReplaceReadGroups I={} O={} RGID={} RGPU={} RGSM={} RGLB={} RGPL={} VALIDATION_STRINGENCY=LENIENT"
    command(cmd_str.format(bam_sorted,grouped_bam,fq_id,"NA",fq_id,dna_library,platform)).run_comm(0) 

def deduplication():
    global dedup_bam
    dedup_bam=bam_sorted.replace(".bam","_grouped_dedup.bam")
    metrics_file=bam_sorted.replace(".bam","_metrix.txt")
    cmd_str="picard MarkDuplicates INPUT={} OUTPUT={} METRICS_FILE={} VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=true"
    command(cmd_str.format(grouped_bam,dedup_bam,metrics_file)).run_comm(0)

def fastqc(bam):
    global fastqc_file
    global bam_for_fastqc
    #sample_runID already exists in the bam file name
    bam_for_fastqc = os.path.basename(bam).replace("_grouped_dedup","").replace("_grouped","").replace(prefix_ori+"_","") 
    FI.copy_file(bam, "{}/{}".format(workdir,bam_for_fastqc))
    command("fastqc -o {} --noextract -f bam_mapped {}".format(workdir,bam_for_fastqc)).run_comm(0)
    fastqc_file=bam_for_fastqc.replace(".bam","_fastqc.zip")

def qualimap():
    global qualimap_dir
    qualimap_dir = "{}/{}".format(qcdir,sample_runID)
    FI.create_processing_dir(qualimap_dir)
    command("qualimap bamqc -bam {} -outdir {} -outformat HTML".format(bam_for_fastqc, qualimap_dir)).run_comm(0)

def multiqc():
    global multiQC_outFN
    global fastqc_html
    if mapping_tool=="bowtie2":
        multiQC_input=fastqc_file+" "+bowtie2_log_for_multiqc
    else:
        multiQC_input=fastqc_file
    multiQC_input += " "+qualimap_dir
    multiQC_outFN="{}.multiQC".format(sample_runID)
    command("multiqc -f {} -o {} --filename {} -v".format(multiQC_input,workdir,multiQC_outFN)).run_comm(0)
    fastqc_html="{}_{}_fastqc.html".format(prefix_ori,sample_runID)

def initiate():
    print ("initiating...")
    global indir
    global outdir
    global workdir
    global subdir
    global qcdir
    subdir="reference_mapping"
    workdir=prop.workdir+"/"+subdir
    fi=fileutils()
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    qcdir=workdir+"/qc/"+prefix_ori
    FI.create_processing_dir(workdir)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    FI.create_processing_dir(qcdir)
    FI.copy_file_to_destdir(fastq1,indir)
    if fastq2 is not None:
        FI.copy_file_to_destdir(fastq2,indir)

def execute():
    print ("executing...")
    global out_bam
    FI.change_dir(workdir)
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    set_prefix()
    run_mapping(fastqfiles)
    add_group(fastqfiles)
    if if_dedup:
        print("removing duplicates...")
        deduplication()
        out_bam = dedup_bam
    else:
        out_bam = grouped_bam
    fastqc(out_bam)
    qualimap()
    multiqc()

def post_process():
    print ("post_processing...")
    FI.copy_file_to_destdir(out_bam, outdir)
    pre_to_rm=prefix_ori+"_"
    FI.copy_file(bam_for_fastqc, "{}/{}".format(qcdir,bam_for_fastqc.replace(pre_to_rm,"")))
    FI.copy_file(fastqc_file, "{}/{}".format(qcdir,fastqc_file.replace(pre_to_rm,"")))
    FI.copy_file(fastqc_html,"{}/{}".format(qcdir,fastqc_html.replace(pre_to_rm,"")))
    FI.copy_file_to_destdir(multiQC_outFN+".html",qcdir)
    if mapping_tool=="bowtie2": 
        ## copy bowtie2_log to outdir and qcdir with name changed (as bowtie2.log not work)
        FI.copy_file_to_destdir(bowtie2_log, outdir)
        FI.copy_file_to_destdir(bowtie2_log_for_multiqc, qcdir)

if __name__ == '__main__':
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()

    get_args()
    print ("\n","Properties attributes:\n",prop.__dict__)
    
    getVar = lambda searchList, ind: [searchList[i] for i in ind]   
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
    print (os.path.realpath(__file__)+" DONE")
