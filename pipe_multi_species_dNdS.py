#! /usr/bin/python

import argparse
import re
import sys
import time
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict
from Bio import SeqIO

def get_args():    
    global properties_file
    global g_names_str
    global fi
    global prop
    global min_homo

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script invests genes under selection pressure among multiple species through dNdS')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g', '--genome_names', type=str, help='''Please provide the genome names, seperating by "," 
        with the format of XXX, YYY, ZZZ''', required=True)   
    parser.add_argument('-min', '--min_homo', type=int, help='''Please provide the minimum poteintial homologue numbers in one group, 
        if not defined, 4 will be used as the default''', required=False)
    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)

    # define variables     
    properties_file=args.properties_file    
    prop=properties(properties_file)
    g_names_str=args.genome_names
    if args.min_homo is None:
        min_homo=4  
    else:
        min_homo=args.min_homo
    print "properties_file:",properties_file
    print "gnames:",g_names_str
    print "min_homo:",str(min_homo)

def download(genome):
    ftp_dir=""
    cDNA_fasta_fname=""
    ftp_root_dir="ftp://ftp.ensemblgenomes.org/pub/current/protists/fasta/"
    sub_dirs1=command("curl -s {} | awk '{{print $9}}'".format(ftp_root_dir)).run_comm(1).split()
    if genome.lower() in sub_dirs1:
        ftp_dir="{}/{}/cdna/".format(ftp_root_dir,genome.lower())
    else:
        for sub_dir1 in sub_dirs1:
            sub_dirs2=command("curl -s {}/ | awk '{{print $9}}'".format(ftp_root_dir+sub_dir1)).run_comm(1).split()
            if genome.lower() in sub_dirs2:
                ftp_dir="{}/{}/{}/cdna/".format(ftp_root_dir,sub_dir1,genome.lower())
    if ftp_dir=="":
        print "can not find cDNA fasta file for {}".format(genome)
        sys.exit(1);        
    else:
        cDNA_fasta_gz=command("curl -s {} | awk '{{print $9}}' | grep cdna".format(ftp_dir)).run_comm(1);
        cDNA_fasta_gz=cDNA_fasta_gz.rstrip()
        command("curl -o {} {}".format(cDNA_fasta_gz,ftp_dir+cDNA_fasta_gz)).run_comm(0)
        cDNA_fasta_fname=cDNA_fasta_gz.rstrip(".gz")
        command("gunzip -c {} > {}".format(cDNA_fasta_gz,cDNA_fasta_fname)).run_comm_no_exit(1)
        fasta_genome_map[cDNA_fasta_fname]=genome
    return cDNA_fasta_fname 

def get_genome_cds_fasta():
    global g_names
    global fasta_files
    global fasta_genome_map
    fasta_genome_map={}
    fasta_files=[]
    g_names=g_names_str.split(",")
    for g_name in g_names:
        g_name=g_name.rstrip()
        g_name=g_name.lstrip()
        fasta_files.append(download(g_name))

def make_ctl_file(fa_file,tree_file,out_file,output_num_s):
    fhin=open(prop.get_attrib("model_ctl"),'r')
    last_part_ctl=fhin.read()
    fhin.close()
    ctl_file="{}.ctl".format(output_num_s)
    fhout=open (ctl_file, "w")
    fhout.write("seqfile = {}\n".format(fa_file))
    fhout.write("treefile = {}\n".format(tree_file))
    fhout.write("outfile  = {}\n".format(out_file))
    fhout.write(last_part_ctl)
    fhout.close()       
    return ctl_file

def check_seqlen_and_change(alignment_fa):
    fhin=open(alignment_fa,'r')
    whole_fa_str=fhin.read()
    fhin.close()
    blocks=whole_fa_str.split(">")
    (acc,seq)=blocks[1].split('\n',1) # the first block
    seq=seq.rstrip()
    seqlen=len(seq)-seq.count('\n')
    seq_postfix=""
    if seqlen%3!=0:
        for i in range(3-seqlen%3):
            seq_postfix=seq_postfix+'-'
        command("mv {} {}.tmp".format(alignment_fa, alignment_fa)).run_comm(0)
        fhout=open(alignment_fa,'w')
        for block in blocks:
            if re.search("\w+",block):
                (acc,seq)=block.split('\n',1)
                seq=seq.rstrip()+seq_postfix
                fhout.write(">{}\n{}\n".format(acc,seq))
        fhout.close()
        command("rm {}.tmp".format(alignment_fa)).run_comm(0)
        
        
def align_phy_tree_codeml():
    global codeml_outputs
    global acc_fasta_map
    codeml_outputs=[]
    blast_comd=prop.get_attrib("crb-blast")
    clustalo=prop.get_attrib("clustalo")
    seqret=prop.get_attrib("seqret")
    seqret_para="-sformat fasta -sprotein1 -osformat phylip"
    raxmlHPC=prop.get_attrib("raxmlhpc")
    raxmlHPC_para="-p 12345 -n T7"
    codeml=prop.get_attrib("codeml")

    tsv_files=[]    
    acc_matchList={}
    acc_seq={}
    acc_fasta_map={}
    orthologs_prefix="clustalO_input"

    for i in range(len(g_names)):
        for j in range(len(g_names)):
            if j>i:
                fasta1=fasta_files[i]
                fasta2=fasta_files[j]
                command("{} --query {} --target {} --threads 8 --evalue 1e-5 --output {}_{}.tsv".format(blast_comd,fasta1,fasta2,fasta1,fasta2)).run_comm(0)
                tsv_files.append("{}_{}.tsv".format(fasta1,fasta2))

    for tsv_file in tsv_files:
        fh_tsv=open(tsv_file,'r')
        print "working on ..."+tsv_file
        for line in fh_tsv:
            line=line.rstrip()
            (acc1,acc2)=sorted(getVar(line.split("\t"),[0,1]))
            # if acc1 is key, append acc2 to the list of value,
            # else if acc2 is key, add acc1 as key and add acc2 value list to acc1 value list
            # else check whether acc1 or acc2 is the value of acc_matchList,
            ### if yes, add acc1 and acc2 to the value list (duplicates will be removed in the next step
            ### if no, add acc1 to key and acc2 to its value list
            if acc1 in acc_matchList:
                acc_matchList[acc1].append(acc2)
                print "acc_matchList..."+acc1+" "+acc2
            elif acc2 in acc_matchList:
                acc_matchList[acc1]=[]
                acc_matchList[acc1].extend(acc_matchList[acc2])
                print "acc_matchList..."+acc1+" "+acc2
            else: #check whether acc1 or acc2 is the value of acc_matchList
                for key, value in sorted(acc_matchList.items()):
                    if value == acc1 or value == acc2:
                        acc_matchList[key].extend([acc1,acc2])
                        break
                # if acc1 and acc2 are not key and are not value
                acc_matchList[acc1]=[]
                acc_matchList[acc1].append(acc2)
          
        fh_tsv.close()

    # make the value list for each key unique
    for key in sorted(acc_matchList):
        value_list=acc_matchList[key]
        tmp_dict={}
        for value in value_list:
            tmp_dict[value]=1
        acc_matchList[key]=sorted(tmp_dict)

    print "Matching acc with seq and fasta file name..."
    for fasta_file in (fasta_files):
        print "here start reading fasta file..."+fasta_file
        tmp_fasta_map=SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
        for acc in tmp_fasta_map:
            acc_seq[acc]=str(tmp_fasta_map[acc].seq)
            acc_fasta_map[acc]=fasta_file

    print "Writing fasta file and phy file for each matched acc list..."
    output_num=0
    rMap={}

    command("find . -name \"RAxML_*\" -delete").run_comm(0)
    for key in sorted(acc_matchList):
        if len(acc_matchList[key])>=min_homo-1:
            output_num=output_num+1
            ortholog_fa_fname=orthologs_prefix+"_"+str(output_num)
            print "writing to..."+ortholog_fa_fname
            fh_ortholog_fa=open(ortholog_fa_fname,'w')
            # fh_ortholog_fa.write(">{}_{}\n{}".format(_homo-1key,acc_fasta_map[key],acc_seq[key])) #acc, fasta file index, seq
            fh_ortholog_fa.write(">{}\n{}".format(key,acc_seq[key])) #acc, fasta file index, seq
            #fasta_files[acc_fasta_map[key]]=1
            for value in sorted(acc_matchList[key]):
            #    if acc_fasta_map[value] not in fasta_files:
             #       fasta_files[acc_fasta_map[value]]=1
            #    fh_ortholog_fa.write("\n>{}_{}\n{}".format(value,acc_fasta_map[value],acc_seq[value])) #acc, seq
                fh_ortholog_fa.write("\n>{}\n{}".format(value,acc_seq[value])) #acc,seq
            fh_ortholog_fa.close()
            # Alignment of Orthologs
            output_num_str=str(output_num)
            command("{} -i {} -t DNA --force -o clustalO_output_{}.fa -v".format(clustalo,ortholog_fa_fname,output_num_str)).run_comm(0)
            # get seqlen for the first sequence, if not the multiple of 3, add -- to the end of each sequence
            check_seqlen_and_change("clustalO_output_{}.fa".format(output_num_str))
            # convert fasta format to phylip format for raxmlHPC
            time.sleep(3)
            command("{}  clustalO_output_{}.fa clustalO_output_{}.phy {}".format(seqret,output_num_str,output_num_str,seqret_para)).run_comm(0)
            time.sleep(3)
            # build tree
            command("{} -m GTRGAMMA {} -s clustalO_output_{}.phy -n {}".format(raxmlHPC,raxmlHPC_para,output_num_str,output_num_str)).run_comm(0)
            time.sleep(3)
            # run codeml
            ctl_file=make_ctl_file("clustalO_output_{}.fa".format(output_num_str),"RAxML_bestTree.{}".format(output_num_str), 
                                   "codeml_out_{}".format(output_num_str),output_num_str)
            time.sleep(3)
            command("yes "" |{} {}".format(codeml,ctl_file)).run_comm_no_exit(1)
            codeml_outputs.append("codeml_out_{}".format(output_num_str))
    

def report():
    print "report..."
    fh_report=open("report",'w')
    for codeml_output in codeml_outputs:
        omega=command("grep omega {}".format(codeml_output)).run_comm(1) 
        if re.search("\d+",omega):
            fh_report.write(omega)
            accs_str=command("grep '^#' {} | awk '{{print $2}}'".format(codeml_output)).run_comm(1)
            accs=accs_str.rstrip().split()
            for acc in accs:
                fh_report.write("{}_{}\n".format(acc,fasta_genome_map[acc_fasta_map[acc]]))
            fh_report.write("\n")        
    fh_report.close()       


def initiate():
    print "initiating..."
    global workdir
    global indir
    global outdir
    global type        
    subdir="dNdS"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    #for file_in in (fasta_files):
     #   fi.copy_file_to_destdir(file_in,indir)
    fi.change_dir(workdir)    

def execute():
    print "executing..."
    get_genome_cds_fasta()
    align_phy_tree_codeml()
    report()   

def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir("report",outdir)

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    get_args()
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run_blast the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
     
