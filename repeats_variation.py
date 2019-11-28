
import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from build_snpEff_db import snpEff_db
from collections import defaultdict

get_sub_list = lambda searchList, ind: [searchList[i] for i in ind]

class strv():
    def __init__(self,properties_file,genome_name,genome_fasta,bam_file,prefix,if_anno,subdir):
        self.properties_file = properties_file
        self.prop = properties(properties_file)  
        self.genome_name = genome_name
        self.genome_fasta = genome_fasta
        self.bam_file = bam_file
        self.prefix = prefix
        self.if_anno = if_anno        
        self.subdir = subdir 
        self.fi = fileutils()
    
    def run_strviper(self):
        ## STRViper: Short Tandem Repeat Variation Indentification from Paired-End Reads
        global fName_str_vcf
        global fName_str
        global fName_str_var
        fName_str_var="{}.strv".format(self.prefix)
        fName_str="{}.trf.str".format(self.genome_name)
        fName_str_vcf=self.prefix+".vcf"
        trf=self.prop.get_attrib("trf")  
        samtools=self.prop.get_attrib("samtools")    
        jsat=self.prop.get_attrib("jsat")     
        genome_fasta_pattern="^(.*/)?(.*?)$"
        genome_fasta_fname=re.findall(genome_fasta_pattern, self.genome_fasta)[0][1]
        trf_out_name=genome_fasta_fname+".2.7.7.80.10.50.500.dat" 
        # tandem repeats finder to detect tandem repeat in genome fasta file,output is trf_out_name 
        # parameter details, see *1, output table explanation, see *2
        command("{} {} 2 7 7 80 10 50 500 -d -h".format(trf,self.genome_fasta)).run_comm(0)
        # convert it to the format that STRViper will read in, output explanation see *2, not only parsing, 
        # but also skip some regions as better quality regions related
        command("{} parseTRF --input {} --output {} --format str".format(jsat,trf_out_name,fName_str)).run_comm(0) 
        command("{} sort {} >{}.sort.bam".format(samtools,self.bam_file,self.prefix)).run_comm(0)
        # Extract fragment size from bamfile 
        command("{} sam2fragment --input {}.sort.bam --output {}.fragment".format(jsat,self.prefix,self.prefix)).run_comm(0)
        # Sort fragment list
        command("{} sortFragment --input {}.fragment --output {}.sort.fragment".format(jsat,self.prefix,self.prefix)).run_comm(0)
        # Make the variation calls
        command("{} fragment2var  --trfFile {}.trf.str --output {} {}.sort.fragment".format(jsat,self.genome_name,fName_str_var,self.prefix)).run_comm(0)
        # Convert variation calls in strv format to vcf format
        command("{} strv2vcf  --input {} --output {}.vcf --reference {}".format(jsat,fName_str_var,self.prefix,self.genome_fasta)).run_comm(0)

    def run_snpEff(self):         
        print ("run annotation...")
        global fName_ann_vcf
        annot_sw="snpeff"    
        snpEff=self.prop.get_attrib("snpeff")
        snpeff_db=snpEff_db(self.properties_file,self.genome_name)
        snpeff_db.build_snpeff_db()
        fName_ann_vcf=self.prefix+".ann.vcf"
        command("java -jar {} -c snpEff.config {} {} > {}".format(snpEff,self.genome_name,fName_str_vcf,fName_ann_vcf)).run_comm(0);    
  
    def initiate(self):
        print ("initiating...")
        global workdir
        global indir
        global outdir
        workdir=self.prop.workdir+"/"+self.subdir 
        indir=workdir+"/in/"
        outdir=workdir+"/out/"
        self.fi.create_processing_dir(indir)
        self.fi.create_processing_dir(outdir)
        for file_in in (self.genome_fasta,self.bam_file): 
            self.fi.copy_file_to_destdir(file_in,indir)
        self.fi.change_dir(workdir)    

    def execute(self):
        print ("executing...")
        self.run_strviper()
        if self.if_anno:
             self.run_snpEff()
 
    def post_process(self):
        print ("post_processing...")
        self.fi.copy_file_to_destdir(fName_str,outdir)
        self.fi.copy_file_to_destdir(fName_str_var,outdir)
        self.fi.copy_file_to_destdir(fName_str_vcf,outdir)   
        if 'fName_ann_vcf' in globals() and self.fi.check_exist("{}/{}".format(workdir,fName_ann_vcf)): #ERR311213_1.xin_test.ann.vcf
            self.fi.copy_file_to_destdir(fName_ann_vcf,outdir)
    

    def run(self):                
        #run_blast the initiation code
        self.initiate() 
        #execute the main part of the program
        self.execute()
        #post execution code
        self.post_process()


'''
*1
trf Parameters:
2 7 7 80 10 50 500 
    Match, Mismatch, and Delta: Weights for match, mismatch and indels. These parameters are for Smith-Waterman style local alignment using wraparound dynamic programming. Lower weights allow alignments with more mismatches and indels. A match weight of 2 has proven effective with mismatch and indel penalties in the range of 3 to 7. Mismatch and indel weights are interpreted as negative numbers. A 3 is more permissive and a 7 less permissive. The recomended values for Match Mismatch and Delta are 2, 7, and 7 respectively.
    PM and PI: Probabilistic data is available for PM values of 80 and 75 and PI values of 10 and 20. The best performance can be achieved with values of PM=80 and PI=10. Values of PM=75 and PI=20 give results which are very similar, but often require as much as ten times the processing time when compared with values of PM=80 and PI=10.
    Minscore: The alignment of a tandem repeat must meet or exceed this alignment score to be reported. For example, if we set the matching weight to 2 and the minimun score to 50, assuming perfect alignment, we will need to align at least 25 characters to meet the minimum score (for example 5 copies with a period of size 5).
    Maxperiod: Period size is the program's best guess at the pattern size of the tandem repeat. The program will find all repeats with period size between 1 and 2000, but the output can be limited to a smaller range.

    -m: This is an optional parameter and when present instructs the program to generate a masked sequence file. The masked sequence file is a FASTA format file containing a copy of the sequence with every location that occurred in a tandem repeat changed to the letter 'N'. The word "masked" is added to the sequence description line just after the '>' character.
    -f: If this option is present, flanking sequence around each repeat is recorded in the alignment file. This may be useful for PCR primer determination. Flanking sequence consists of the 500 nucleotides on each side of a repeat.
    -d: A data file is produced if this option is present. This file is a text file which contains the same information, in the same order, as the summary table file, plus consensus pattern and repeat sequences. This file contains no labeling and is suitable for additional processing, for example with a perl script, outside of the program.
    -h: suppress HTML output (this automatically switches -d to ON) 

*2
1 Indices of the repeat relative to the start of the sequence.
2 Period size of the repeat.
3 Number of copies aligned with the consensus pattern.
4 Size of consensus pattern (may differ slightly from the period size).
5 Percent of matches between adjacent copies overall.
6 Percent of indels between adjacent copies overall.
7 Alignment score.
8 Percent composition for each of the four nucleotides.
9 Entropy measure based on percent composition.

*3 
"score" means "Alignment Minscore" from trf_out, bigger maybe better, here threhold is 50

diff between VNTR and STR:
https://pediaa.com/difference-between-vntr-and-str/
'''
