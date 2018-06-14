#! /analysis/xin/parasite/bin/python/bin/python

import argparse
import re
import sys
from utilities import command
from utilities import fileutils
from utilities import properties
from collections import defaultdict

def get_args():    
    global properties_file
    global cds_fna1
    global cds_faa1
    global cds_fna2
    global cds_faa2
    global genome1
    global genome2
    global map_file
    global filter_eval
    global filter_identity
    global prefix
    global fi
    global makeblastdb_sw
    global blastn_sw
    global prop
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script invests genes under selection pressure between two species through dNdS')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', required=True)
    parser.add_argument('-g1', '--genome_name1', type=str, help='''Please provide the first genome name, otherwise, 
                        "ch" for "C. hominis" will be used''', required=False)   
    parser.add_argument('-g2', '--genome_name2', type=str, help='''Please provide the second genome name, otherwise, 
                        "cp" for "C. parvum" will be used''', required=False)     
    parser.add_argument('-fn1', '--cds_fna1', type=str, help='Please provide the first cds fna file, otherwise, ch fna file will be used.', required=False)
    parser.add_argument('-fn2', '--cds_fna2', type=str, help='Please provide the second cds fna file, otherwise, cp fna file will be used.', required=False)
    parser.add_argument('-fa1', '--cds_faa1', type=str, help='Please provide the first cds faa file, otherwise, ch faa file will be used.', required=False)
    parser.add_argument('-fa2', '--cds_faa2', type=str, help='Please provide the second cds faa file, otherwise, cp faa file will be used.', required=False)
    parser.add_argument('-m', '--map', type=str, help='''Please provide the file for mapping the chromosome accessions, 
                                                      one pair in each line and separated by tab, 
                                                      otherwise, no chromosome information will be provided in the output file''', required=False)    
    parser.add_argument('-fi', '--filter_identity', type=str, help='the identity percentage for filtering the blast hits.', required=False)
    parser.add_argument('-fe', '--filter_eval', type=str, help='the eval for filtering the blast hits.', required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', required=True)  
    
    # check args
    args = parser.parse_args()
    fi=fileutils()
    fi.check_exist(args.properties_file)
    for opt_arg_fpath in (args.cds_fna1,args.cds_fna2,args.cds_faa1,args.cds_faa2,args.map):
        if opt_arg_fpath is not None:            
            fi.check_exist(opt_arg_fpath)     

    # define variables     
    makeblastdb_sw="makeblastdb"
    blastn_sw="blastn"
    default_gname1="ch"
    default_gname2="cp"
    map_file="None"
    filter_eval="0"
    filter_identity="0"
    filter_length="0"
        
    properties_file=args.properties_file    
    prop=properties(properties_file)
    if args.genome_name1 is not None:
        genome1=args.genome_name1
    else:
        genome1=default_gname1        
    if args.genome_name2 is not None:
        genome2=args.genome_name2
    else:
        genome2=default_gname2        
    if args.cds_fna1 is not None:
        cds_fna1=args.cds_fna1
    else:
        cds_fna1=prop.get_attrib(genome1+"_cds_fna")        
    if args.cds_faa1 is not None:
        cds_faa1=args.cds_faa1
    else:
        cds_faa1=prop.get_attrib(genome1+"_cds_faa")  
    if args.cds_fna2 is not None:
        cds_fna2=args.cds_fna2
    else:
        cds_fna2=prop.get_attrib(genome2+"_cds_fna")
    if args.cds_faa2 is not None:
        cds_faa2=args.cds_faa2
    else:
        cds_faa2=prop.get_attrib(genome2+"_cds_faa")
    if args.map is not None:
        map_file=args.map
    if args.filter_eval is not None:
        filter_eval=args.filter_eval
    if args.filter_identity is not None:
        filter_identity=args.filter_identity
    prefix=args.prefix   
    
    print "properties_file:",properties_file
    print "genome1:",genome1
    print "genome2:",genome2
    print "cds_fna1:",cds_fna1
    print "cds_faa1:",cds_faa1
    print "cds_fna2:",cds_fna2
    print "cds_faa2:",cds_faa2   
    print "filter_eval:",filter_eval
    print "filter_identity_perc:",filter_identity   
    print "prefix:",prefix
   

def build_seq_dict(fpath):
    return_dict={}
    pattern="^\|([A-Z0-9.]+)_.*?_([A-Z0-9.]+).*\[gbkey\=CDS\](\w+)"
    str_one_file=open(fpath, 'r').read().replace('\n', '')    
    block_array=str_one_file.split(">lcl")
    block_array.pop(0)
    for block in block_array:
        if re.search(pattern,block): 
            (con_acc,protein_acc,seq)=re.findall(pattern,block)[0]            
            return_dict[con_acc+"_"+protein_acc]=seq
        else:
            print "error: "+block
            sys.exit(1)
    return return_dict

            
def get_map_info():
    global map_dict
    map_dict={}
    file_in=open (map_file,'r')
    for line in file_in:
        line=line.rstrip("\n")
        (acc1,acc2)=line.split(",")
        map_dict[acc1]=acc2
        
        
def get_acc_combi(pattern, acc_ori):
    if re.search(pattern,acc_ori):
        chr_acc=re.findall(pattern,acc_ori)[0][0]
        acc=re.findall(pattern,acc_ori)[0][0]+re.findall(pattern,acc_ori)[0][1]
    else:
        print "WARNING: "+acc_ori+" has wrong pattern"
        sys.exit(1)
    return (chr_acc,acc)


def run_blast():   
    global blast_filt1    
    global mapped_dict
    mapped_dict={}
    makeblastdb=prop.get_attrib(makeblastdb_sw)
    blastn=prop.get_attrib(blastn_sw)    
    blast_filt1=defaultdict(lambda: defaultdict(str))
    fname_blast_out="{}_{}_{}_cds_blast.out".format(prefix,genome1,genome2)
    acc_pattern=".*?([A-Z0-9.]+)_.*?(_[A-Z0-9.]+)(.*)"
    command("{} -in {} -parse_seqids -dbtype nucl -out {}_databaseBLAST".format(makeblastdb,cds_fna1,genome1)).run_comm(0)
    command('''{} -query {} -db {}_databaseBLAST -task blastn -dust no -outfmt 7 -max_target_seqs 1 -out {}'''.
            format(blastn,cds_fna2,genome1,fname_blast_out)).run_comm(0)   
    
    if map_file!="None":
        get_map_info()
    
    file_blast_out=open(fname_blast_out,'r')
    for line in file_blast_out:
        line=line.rstrip("\n")
        if line.startswith("lcl"):
            (acc2_ori, acc1_ori)=getVar(line.split(),[0,1])   
            (chr_acc1, acc1)=get_acc_combi(acc_pattern, acc1_ori)
            (chr_acc2, acc2)=get_acc_combi(acc_pattern, acc2_ori)
            if acc1 not in blast_filt1.keys() and acc2 not in blast_filt1[acc1].keys():
                (identity_perc, evalue)=getVar(line.split(),[2,10])
                if identity_perc>=float(filter_identity) and evalue>=float(filter_eval):
                    blast_filt1[acc1][acc2]=",".join(getVar(line.split(),[2,3,4,5,6,7,8,9,10,11])) #[2] is the iden_perc
                    if map_file!="None" and map_dict[chr_acc1]==chr_acc2:
                        mapped_dict[acc1+","+acc2]=1
                        
    file_blast_out.close()
    

def get_snp_anno_num(dna_seq1, protein_seq1,dna_seq2, protein_seq2):
    syn_nuclear_posi_list=[]
    syn_num=0
    nonsyn_num=0
    for i in range(0,len(protein_seq1)-1):
        if protein_seq1[i]!=protein_seq2[i]:
            tmp=(i+1)*3
            for j in range(1,3):
                syn_nuclear_posi_list.append(tmp-j)            
    for i in range(0,len(dna_seq1)-1):
        if dna_seq1[i]!=dna_seq2[i]:
            if i in syn_nuclear_posi_list:
                syn_num=syn_num+1
            else:
                nonsyn_num=nonsyn_num+1  
    return (syn_num,nonsyn_num)
            
        
def analysis_blast_result():   
    global file_out_name
    global same_len
    global file_out
    out_dict=defaultdict(lambda: defaultdict(str))
    cds_fna1_dict=build_seq_dict(cds_fna1)
    cds_faa1_dict=build_seq_dict(cds_faa1)
    cds_fna2_dict=build_seq_dict(cds_fna2)
    cds_faa2_dict=build_seq_dict(cds_faa2)
    same_len={}
    file_out_name="{}_{}_{}_dNdS.csv".format(prefix,genome1,genome2)
    file_out=open(file_out_name,'w')
    key1="NO synonymous or nonsynonymous SNPs"
    key2="Same protein seqs But with nonsynonymous SNPs"
    key3="Diff protein seqs But full CDS hits"
    key4="Not full CDS hits"
    titles=",".join('''query_CONacc_PROTEINacc subject_CONacc_PROTEINacc synonymous_SNP_num nonsynonymous_SNP_num %identity alignment_length mismatches 
            gap_opens query_start query_end subject_start subject_end evalue bit_score'''.split())       
        
    for acc1 in blast_filt1.keys():
        for acc2 in blast_filt1[acc1].keys():
            if (len(cds_fna1_dict[acc1])==len(cds_fna2_dict[acc2]) and len(cds_fna1_dict[acc1])==int(blast_filt1[acc1][acc2].split(",")[1])) or cds_faa1_dict[acc1]==cds_faa2_dict[acc2]:                
                (syn_num,nonsyn_num)=get_snp_anno_num(cds_fna1_dict[acc1],cds_faa1_dict[acc1],cds_fna2_dict[acc2],cds_faa2_dict[acc2])
                same_len[acc1+","+acc2]=str(syn_num)+","+str(nonsyn_num)
            if cds_faa1_dict[acc1]==cds_faa2_dict[acc2]:                 
                if cds_fna1_dict[acc1]==cds_fna2_dict[acc2]:    
                    out_dict[key1][acc1+","+acc2]=blast_filt1[acc1][acc2]                   
                else:
                    out_dict[key2][acc1+","+acc2]=blast_filt1[acc1][acc2]  
            else:
                if acc1+","+acc2 in same_len.keys():
                    out_dict[key3][acc1+","+acc2]=blast_filt1[acc1][acc2]  
                else:
                    out_dict[key4][acc1+","+acc2]=blast_filt1[acc1][acc2]
                
    file_out.write(titles+"\n")
    if map_file!="None":
        file_out.write("\nHIT PAIRS ON THE MAPPED CHROMOSOMES:\n\n")
        for key in (key1,key2,key3):#for key in (key1,key2,key3,key4):
            file_out.write("\n{}:\n".format(key))
            for acc_pair in out_dict[key].keys():
                if acc_pair in mapped_dict.keys():
                    if acc_pair in same_len.keys():
                        file_out.write("{},{},{}\n".format(acc_pair,same_len[acc_pair],out_dict[key][acc_pair]))
                    else:
                        file_out.write("{},{},{}\n".format(acc_pair,"Null,Null",out_dict[key][acc_pair]))
            
        file_out.write("\nHIT PAIRS ON THE UNMAPPED CHROMOSOMES:\n\n")       
        for key in (key1,key2,key3):#for key in (key1,key2,key3,key4)
            file_out.write("{}:\n".format(key))
            for acc_pair in out_dict[key].keys():
                if acc_pair not in mapped_dict.keys():
                    if acc_pair in same_len.keys():
                        file_out.write("{},{},{}\n".format(acc_pair,same_len[acc_pair],out_dict[key][acc_pair]))
                    else:
                        file_out.write("{},{},{}\n".format(acc_pair,"Null,Null",out_dict[key][acc_pair])) 
    else:
        for key in (key1,key2,key3):#for key in (key1,key2,key3,key4):
            file_out.write("{}:\n".format(key))
            for acc_pair in out_dict[key].keys():
                if acc_pair in same_len.keys():
                    file_out.write("{},{},{}\n".format(acc_pair,same_len[acc_pair],out_dict[key][acc_pair]))
                else:
                    file_out.write("{},{},{}\n".format(acc_pair,"Null,Null",out_dict[key][acc_pair])) 
    file_out.close()

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
    for file_in in (cds_fna1,cds_fna2,cds_faa1,cds_faa2):
        fi.copy_file_to_destdir(file_in,indir)
    fi.change_dir(workdir)    

def execute():
    print "executing..."
    run_blast()
    analysis_blast_result()

def post_process():
    print "post_processing..."
    fi.copy_file_to_destdir("{}/{}".format(workdir,file_out_name),outdir)

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
   
     
