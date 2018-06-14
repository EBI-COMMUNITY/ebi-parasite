#!  /analysis/xin/parasite/bin/python/bin/python

import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict


class strv():
    def __init__(self,properties_file_in,genome_name_in,genome_fasta_in,bam_file_in,prefix_in,if_anno_in,subdir_in):
        global genome_name
        global genome_fasta 
        global prop        
        global bam_file
        global prefix
        global if_anno
        global fi
        global get_sub_list
        global subdir
        prop=properties(properties_file_in)  
        genome_name = genome_name_in
        genome_fasta=genome_fasta_in
        bam_file = bam_file_in
        prefix = prefix_in
        if_anno= if_anno_in        
        subdir=subdir_in 
        fi=fileutils()
        get_sub_list = lambda searchList, ind: [searchList[i] for i in ind] 

    
    def run_strviper(self):
        global fName_str_vcf
        global fName_str
        global fName_str_var
        fName_str_var="{}.strv".format(prefix)
        fName_str="{}.trf.str".format(genome_name)
        fName_str_vcf=prefix+".vcf"
        trf=prop.get_attrib("trf")  
        samtools=prop.get_attrib("samtools")    
        jsat=prop.get_attrib("jsat")     
        genome_fasta_pattern="^(.*/)?(.*?)$"
        genome_fasta_fname=re.findall(genome_fasta_pattern, genome_fasta)[0][1]
        trf_out_name=genome_fasta_fname+".2.7.7.80.10.50.500.dat" 
        command("{} {} 2 7 7 80 10 50 500 -d -h".format(trf,genome_fasta)).run_comm(0)
        command("{} parseTRF --input {} --output {} --format str".format(jsat,trf_out_name,fName_str)).run_comm(0) 
        command("{} sort {} >{}.sort.bam".format(samtools,bam_file,prefix)).run_comm(0)
        command("{} sam2fragment --input {}.sort.bam --output {}.fragment".format(jsat,prefix,prefix)).run_comm(0)
        command("{} sortFragment --input {}.fragment --output {}.sort.fragment".format(jsat,prefix,prefix)).run_comm(0)
        command("{} fragment2var  --trfFile {}.trf.str --output {} {}.sort.fragment".format(jsat,genome_name,fName_str_var,prefix)).run_comm(0)
        command("{} strv2vcf  --input {} --output {}.vcf --reference {}".format(jsat,fName_str_var,prefix,genome_fasta)).run_comm(0)

    def run_snpEff(self):         
        print "run annotation..."
        global fName_ann_vcf
        annot_sw="snpeff"    
        snpEff=prop.get_attrib(annot_sw)  
        genome_db={"ch":"c_hominis",
               "cp":"c_parvum"}   
        fName_ann_vcf=prefix+".ann.vcf"
        command("java -jar {} -c snpEff.config {} {} > {}".format(snpEff,genome_db[genome_name],fName_str_vcf,fName_ann_vcf)).run_comm(0);    
  
    def initiate(self):
        print "initiating..."
        global workdir
        global indir
        global outdir
        workdir=prop.workdir+"/"+subdir 
        indir=workdir+"/in/"
        outdir=workdir+"/out/"
        fi.create_processing_dir(indir)
        fi.create_processing_dir(outdir)
        for file_in in (genome_fasta,bam_file): 
            fi.copy_file_to_destdir(file_in,indir)
            fi.change_dir(workdir)    

    def execute(self):
        print "executing..."
        self.run_strviper()
        if if_anno:
             self.run_snpEff()
 
    def post_process(self):
        print "post_processing..."
        fi.copy_file_to_destdir(fName_str,outdir)
        fi.copy_file_to_destdir(fName_str_var,outdir)
        fi.copy_file_to_destdir(fName_str_vcf,outdir)   
        if 'fName_ann_vcf' in globals() and fi.check_exist("{}/{}".format(workdir,fName_ann_vcf)): #ERR311213_1.xin_test.ann.vcf
            fi.copy_file_to_destdir(fName_ann_vcf,outdir)
    

    def run(self):                
        #run_blast the initiation code
        self.initiate() 
        #execute the main part of the program
        self.execute()
        #post execution code
        self.post_process()



'''
head ch.trf.str
##Short Tandem Repeat obtained from trf for LN877947.1
##Parameters: 2 7 7 80 10 50 500
#H:chr    start    end    period    unitNo    score
LN877947.1    32910    32964    6    9.7    62.0
LN877947.1    39676    39710    6    5.8    61.0
LN877947.1    78917    78944    3    9.3    56.0

'''
