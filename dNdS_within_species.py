#! /usr/bin/env python3

## version1, 28 Nov 2019 by Xin Liu

import argparse
import re
import sys
import glob
from utilities import command
from utilities import fileutils
from utilities import properties
from utilities import misc
from build_snpEff_db import snpEff_db
from collections import defaultdict
import os

def get_args():    
    global properties_file
    global genome_name
    global prefix
    global vcf_file_pattern
    global go_file
    global mapping_file
    global prop

    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script invests genes under selection pressure within 
                                                    species through dNdS. It creates variation annotation 
                                                    file for each SNP vcf file, and gene variation annotation 
                                                    summary file based on all vcf files by using snpEff''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', 
                        required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name, only with those obtained from genome_list.txt''', 
                        required=True) 
    parser.add_argument('-vp', '--vcf_file_pattern', type=str, help="Please provide snp vcf files' pattern with full file path", 
                        required=True) 
    parser.add_argument('-go', '--go_file', type=str, help="Please provide the full path of the gene ontology file",
                        required=True)
    parser.add_argument('-m', '--mapping_file', type=str, help='''Please provide the mapping file path, which contains one column of 
                                                                read_ID from vcf file and one column of its corresponding sample_name''', 
                        required=False)
    parser.add_argument('-pre', '--prefix', type=str,help='Please provide the prefix for the output file.', 
                        required=True)  
    
    # check args
    args = parser.parse_args() 
    FI.check_exist(args.properties_file)  
    properties_file=args.properties_file    
    prop=properties(properties_file)
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()):
        MISC.my_exit("{} is not available, please try another genome".format(args.genome_name))     
    if not re.search(".vcf$",args.vcf_file_pattern):
        MISC.my_exit("vcf_file_pattern need to end up with .vcf")
    FI.check_exist(args.go_file)
    go_file=args.go_file
    genome_name=args.genome_name
    vcf_file_pattern=args.vcf_file_pattern
    prefix=args.prefix   
    mapping_file=args.mapping_file
    
    print ("properties_file:",properties_file)
    print ("genome_name:",genome_name)
    print ("vcf_file_pattern:",vcf_file_pattern)
    print ("go_file:",go_file)
    print ("mapping_file:",mapping_file)
    print ("prefix:",prefix)
   
def get_go():
    global go_col_names;
    global go_dict;
    go_dict={}
    fhin=open(go_file,'r')
    lines=fhin.readlines()
    go_col_names_ori = lines[0].rstrip().split(",")[2:]
    go_col_names = []
    for go_col_name in go_col_names_ori:
        go_col_name = go_col_name.rstrip("\"").lstrip("\"")
        go_col_names.append(go_col_name)
    print (go_col_names)
    for line in lines[1:]:
        line=line.rstrip()
        go_cols_all=[]
        go_cols_all_ori = line.split('\",\"')
        for go_col in go_cols_all_ori:
            go_col=go_col.rstrip("\"").lstrip("\"")
            go_cols_all.append(go_col)
        gene = go_cols_all[0]
        go_cols_all[4]="N" if go_cols_all[4] == "0" else "Y"
        go_cols_all[5]="N" if go_cols_all[5] == "null" else "Y"
        go_dict[gene]=go_cols_all[2:]
    fhin.close()

def run_snpEff():         
    global ann_stat_fpath
    global ann_stat_fpaths
    global ann_vcf_fpaths
    ann_vcf_fpaths=[]
    ann_stat_fpaths=[]
    snpeff_db=snpEff_db(properties_file,genome_name)    
    snpeff_db.build_snpeff_db()
    for vcf_fpath in all_vcf_fpaths:
        runID=MISC.get_runID(vcf_fpath)
        vcf_prefix="{}_{}".format(prefix, os.path.basename(vcf_fpath).rstrip(".vcf"))
        ann_vcf_fpath=vcf_prefix+"_ann.vcf"
        ann_vcf_fpaths.append(ann_vcf_fpath)
        if mapping_file is not None:
            ann_stat_fpath="{}_{}.ann_stats".format(mirror[runID],runID)
        else:           
            ann_stat_fpath=runID+".ann_stats"
        ann_stat_fpaths.append(ann_stat_fpath)
        command("snpEff -c snpEff.config {} {} -csvStats {} > {}".format(genome_name, vcf_fpath, ann_stat_fpath, ann_vcf_fpath)).run(0);

def get_genes_from_gff(): 
    global genome_gene_dict
    genome_gene_dict={}
    pattern="ID=(.*?);"
    REF_GFF="{}/../ref/{}.gff".format(workdir, genome_name)
    with open (REF_GFF,"r") as gff_fhin:
        for line in gff_fhin:
            line=line.rstrip("\n")
            if re.search("\w+",line) and not line.startswith("#") and line.split("\t")[2]=='gene':
                if re.search(pattern,line):
                    gene=re.findall(pattern,line)[0]
                    (chro,start,end)=getVar(line.split("\t"),[0,3,4])
                    length=str(int(end)-int(start))
                    genome_gene_dict[gene]=[chro,start,end,length]
                else:
                    print ("ERROR: there is no gene name in {}".format(line))
                    sys.exit(1)

def get_sum(my_map):
    sum=0
    for key in my_map.keys():
        sum=sum+int(my_map[key])
    return str(sum)

def get_top_sum(my_map):
    sum=0
    for key in my_map.keys():
        sum=sum+int(get_sum(my_map[key]))
    return str(sum)

def get_sum_on_variant(my_map,variant):
    sum=0;
    for key in my_map.keys():
        sum=sum+my_map[key][variant]
    return str(sum)
 
def get_sample_runID_dict():
    global mirror
    mirror=MISC.get_samples_by_runIDs(mapping_file)

def analysis_ann_files():
    global fPath_out
    global snp_hash
    global sample_names
    global posi_snp_hash
    sample_names=[]
    value_cols_num=3 # totals
    no_dNdS_gene_dict={}
    dNdS_gene_dict={}
    sample_dict=defaultdict(lambda: defaultdict(str))
    posi_snp_hash=defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(str)))))))
    snp_hash=defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(str)))))))
    ann_pattern='ANN=(.*?)[\s;]'
    for ann_vcf_fpath in sorted(ann_vcf_fpaths):
        runID=MISC.get_runID(ann_vcf_fpath)
        if mapping_file is not None :        
            sample_name=mirror[runID]
        else:
            sample_name=runID
        print ("here working on "+ann_vcf_fpath)
        sample_names.append(sample_name)
        infile=open(ann_vcf_fpath, 'r')
        for line in infile:
            if not line.startswith('#') and re.search("ANN=",line):
                (chro,posi,ref_base)=getVar(line.split(),[0,1,3])
                if re.compile(ann_pattern).search(line) :
                    ann_field=re.findall(ann_pattern, line)[0]
                    annseg_array=ann_field.split(",")
                    for annseg in annseg_array:
                        alt=annseg.split("|")[0]
                        variant=annseg.split("|")[1]
                        gene=(annseg.split("|")[3])
                        if not re.search('[a-zA-Z]', gene) or not re.search('[a-zA-Z]', variant):
                            print ("WARNING: None:"+"annseg="+annseg+" gene="+gene+" variant="+variant+"\n"+line)                            
                        elif not re.search("upstream",variant) and not re.search("downstream",variant):
                            snp_hash[gene][chro][posi][ref_base][sample_name][variant]=alt
                            posi_snp_hash[gene][chro][posi][ref_base][variant][alt][sample_name]=1
        infile.close()

def write_gene_matrix():
    global gene_variant_dict
    global gene_locus
    global gene_matrix_file
    variant_dict={}
    gene_locus=defaultdict(lambda: defaultdict(lambda: defaultdict(lambda:(str))))
    gene_variant_dict=defaultdict(lambda: defaultdict(lambda:(str)))
    gene_matrix_file="{}_gene_matrix.csv".format(prefix)
    fileout=open(gene_matrix_file,'w')
    fileout.write("{}\t{}\t{}\t{}".format("GENE","CHROMOSOME","ref_posi","ref_base"))
    for sample_name in sorted(sample_names):
        fileout.write("\t{}_dN_var\t{}_dS_var".format(sample_name,sample_name,sample_name))
    fileout.write("\t{}\t{}\n".format("support_dN_isolates","support_dS_isolates"))
    for gene in sorted(genome_gene_dict.keys()):
        genome_acc=genome_gene_dict[gene][0]
        if gene not in snp_hash.keys():
            fileout.write("{}\t{}".format(gene,genome_acc))
            for i in range(0,2*len(sample_names)+4):
                fileout.write("\tNA")
            fileout.write("\n")
        else:
            gene_variant_dict[gene]["dN"]=0
            gene_variant_dict[gene]["dS"]=0
            for chro in sorted(snp_hash[gene].keys()):
                for posi in sorted(snp_hash[gene][chro].keys()):
                    gene_posi_samples_dN_count=0
                    gene_posi_samples_dS_count=0
                    locus_posi="{}_{}".format(chro,posi)                    
                    gene_locus[gene]["dN"][locus_posi]=0
                    gene_locus[gene]["dS"][locus_posi]=0
                    for ref_base in sorted(snp_hash[gene][chro][posi].keys()):
                        fileout.write("{}\t{}\t{}\t{}".format(gene,chro,posi,ref_base))
                        var_samples_dN=[]
                        var_samples_dS=[]
                        for sample_name in sorted(sample_names):
                            tmp_list=["NA","NA"]
                            if sample_name in snp_hash[gene][chro][posi][ref_base].keys():
                                for variant in snp_hash[gene][chro][posi][ref_base][sample_name].keys():
                                    ## dS
                                    if re.search("synonymous_variant",variant) or re.search("intron_variant",variant): 
                                        tmp_list[1]=snp_hash[gene][chro][posi][ref_base][sample_name][variant]
                                        gene_posi_samples_dS_count=gene_posi_samples_dS_count+1
                                        gene_variant_dict[gene]["dS"]+=1
                                        gene_locus[gene]["dS"][locus_posi]+=1                                         
                                        variant_dict[variant]=1
                                        var_samples_dS.append(sample_name)
                                    ## dN
                                    else: 
                                        tmp_list[0]=snp_hash[gene][chro][posi][ref_base][sample_name][variant]
                                        gene_posi_samples_dN_count=gene_posi_samples_dN_count+1
                                        gene_variant_dict[gene]["dN"]+=1
                                        gene_locus[gene]["dN"][locus_posi]+=1
                                        variant_dict[variant]=1
                                        var_samples_dN.append(sample_name)
                            fileout.write("\t{}\t{}".format(tmp_list[0],tmp_list[1]))
                        fileout.write("\t{}\t{}".format(str(gene_posi_samples_dN_count),str(gene_posi_samples_dS_count)))
                        fileout.write("\n")
    fileout.close()
    for variant in sorted(variant_dict.keys()):
        print ("variants="+variant)
       
def get_str(gene_locus_geneKey_dict,gene_len): 
    # gene_locus[gene]["dN"][locus_posi]=0 gene_locus_geneKey_dict=gene_locus[gene]
    outstr=""
    count_dict=defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    for num in alt_isolates_nums:
        count_dict["dN"][num]=0
        count_dict["dS"][num]=0
    for dN_dS in gene_locus_geneKey_dict.keys():
        for locus_posi in gene_locus_geneKey_dict[dN_dS].keys():
            iso_num=gene_locus_geneKey_dict[dN_dS][locus_posi]
            for num in alt_isolates_nums:
                if iso_num>=int(num):
                    count_dict[dN_dS][num]+=1
    for num in alt_isolates_nums:
        dS_num=count_dict["dS"][num]
        dN_num=count_dict["dN"][num]
        dN_weight="{:.1f}".format(float(1000*dN_num)/float(gene_len))
        dS_weight="{:.1f}".format(float(1000*dS_num)/float(gene_len))
        if dN_num!=0 and dS_num!=0:
            dNdS="{:.1f}".format(float(dN_num)/float(dS_num))
        else:
            dNdS="{}/{}".format(str(dN_num),str(dS_num))
        outstr+="{}\t{}\t{}\t{}\t{}\t".format(str(dN_num),dN_weight,str(dS_num),dS_weight,dNdS)
    return outstr.rstrip("\t")

def write_gene_matrix_summary():
    global alt_isolates_nums
    global sample_num
    global gene_matrix_summary_file
    gene_matrix_summary_file = "{}_gene_summary.csv".format(prefix)
    sample_num=len(sample_names)
    q25=str(int(0.25*sample_num))
    q50=str(int(0.5*sample_num))
    q75=str(int(0.75*sample_num))
    alt_isolates_nums=[q25,q50,q75]
    ann_types=["dN","dS"]
    fPath_out="{}_gene_summary_ori.csv".format(prefix)
    fileout=open(fPath_out,'w')
    # write the column names
    fileout.write("{}\t{}\t{}\t{}\t{}".format("GENE","CHROMOSOME","STRAT","END","LENGTH"))
    for go_col in go_col_names:
        fileout.write("\t{}".format(go_col))
    for num_range in alt_isolates_nums:    
        for ann_type in ann_types:
            fileout.write("\t{}_locusVar_num(>={}_isolate(s))\t{}_locusVar_num(>={}_isolate(s))/kb"
                          .format(ann_type,num_range,ann_type,num_range))
        fileout.write("\tdN/dS(>={}_isolate(s))".format(num_range))
    fileout.write("\n")
    # write the content 
    for gene in sorted(gene_variant_dict.keys()):
        fileout.write("{}".format(gene))
        fileout.write("\t{}\t{}\t{}\t{}".format(genome_gene_dict[gene][0],genome_gene_dict[gene][1],
                                                genome_gene_dict[gene][2],genome_gene_dict[gene][3]
                                               ))
        for each_go in go_dict[gene]:
            fileout.write("\t{}".format(each_go))
        fileout.write("\t{}\n".format(get_str(gene_locus[gene],genome_gene_dict[gene][3])))
    fileout.close()
    cmd = "cat {} | awk '{{split($0,a,\"\\t\");if(a[14]!=0 || a[16]!=0) print $0}}' >{}"
    command(cmd.format(fPath_out, gene_matrix_summary_file)).run_comm(0)

def get_max_alt_num():
    global max_alt_num
    max_alt_num=0
    for gene in sorted(posi_snp_hash.keys()):
        if gene in genome_gene_dict.keys():
            for chro in sorted(posi_snp_hash[gene].keys()):
                for posi in sorted(posi_snp_hash[gene][chro].keys()):
                    for ref_base in sorted(posi_snp_hash[gene][chro][posi].keys()):
                        for variant in sorted (posi_snp_hash[gene][chro][posi][ref_base].keys()):
                            alt_num=len(posi_snp_hash[gene][chro][posi][ref_base][variant].keys())
                            if max_alt_num<alt_num:
                                max_alt_num=alt_num

def get_ref_base_isolates(alt_samples):
    out_samples=[]
    out_sample_str=""
    for sample in sample_names:
        if sample not in alt_samples:
            out_samples.append(sample)
            out_sample_str=out_sample_str+","+sample
    out_sample_str=out_sample_str[1:]
    return (str(len(out_samples)),out_sample_str)

def get_alt_ref_diff(ref_base,alt):
    tmp_dict={}
    sorted_list=[]
    tmp_dict[len(ref_base)]=ref_base
    tmp_dict[len(alt)]=alt
    for my_len in sorted(tmp_dict.keys(), key=int):
        sorted_list.append(tmp_dict[my_len])
    return "diff:"+sorted_list[1].lstrip(sorted_list[0])
 
def write_posi_file():
    global gene_posi_file
    gene_posi_file = "{}_posi.csv".format(prefix)
    get_max_alt_num()
    fhout=open(gene_posi_file,'w')
    fhout.write("{}\t{}\t{}\t{}\t{}\t".format("gene","chro","posi","isolates_ref_base","isolate_name_ref_base"))
    for i in range(0,max_alt_num):
        fhout.write("{}\t{}\t".format("isolates_alt_base","isolate_name"))
    fhout.write("{}\t{}\t{}".format("singleton","var_type","gene_length"))
    fhout.write("\n")
    for gene in sorted(posi_snp_hash.keys()):
        if gene in genome_gene_dict.keys():
            for chro in sorted(posi_snp_hash[gene].keys()):
                for posi in sorted(posi_snp_hash[gene][chro].keys()):
                    for ref_base in sorted(posi_snp_hash[gene][chro][posi].keys()):
                        for variant in sorted (posi_snp_hash[gene][chro][posi][ref_base].keys()):
                            if re.search("synonymous_variant",variant) or re.search("intron_variant",variant):
                                 dNordS="dS"
                            else:
                                 dNordS="dN"
                            alt_num=len(posi_snp_hash[gene][chro][posi][ref_base][variant].keys())
                            alt_str=""
                            alt_samples=[]
                            per_var_sample_num=0
                            for alt in sorted(posi_snp_hash[gene][chro][posi][ref_base][variant].keys()):
                                 sample_num_alt=len(posi_snp_hash[gene][chro][posi][ref_base][variant][alt])
                                 per_var_sample_num+=sample_num_alt
                                 if sample_num_alt==1:
                                     singleton='Y'
                                 else:
                                     singleton='N'
                                 sample_names_str=""
                                 alt_samples.extend(sorted(posi_snp_hash[gene][chro][posi][ref_base][variant][alt].keys()))
                                 for sample_name in sorted(posi_snp_hash[gene][chro][posi][ref_base][variant][alt].keys()):
                                     sample_names_str=sample_names_str+","+sample_name
                                 sample_names_str=sample_names_str[1:]
                                 alt_str=alt_str+"{}_{}\t{}\t".format(str(sample_num_alt),alt,sample_names_str)
                            if per_var_sample_num<sample_num and per_var_sample_num>0:
                                (ref_base_isolate_num,ref_base_isolates_str)=get_ref_base_isolates(alt_samples)
                                alt_num_diff=max_alt_num-alt_num
                                for i in range(0,alt_num_diff):
                                    alt_str=alt_str+"N/A\tN/A\t"
                                alt_str=alt_str.rstrip()
                                fhout.write("{}\t{}\t{}\t".format(gene,chro,posi))
                                fhout.write("{}_{}\t{}\t{}\t".format(ref_base_isolate_num,ref_base,ref_base_isolates_str,alt_str))
                                fhout.write("{}\t{}\t{}".format(singleton,dNordS,genome_gene_dict[gene][3])) #gene_len
                                fhout.write("\n")
    fhout.close()
    print ("sample number="+str(sample_num))

def initiate():
    print ("initiating...")
    global workdir
    global indir
    global outdir
    global qcdir
    global all_vcf_fpaths
    subdir="dNdS"
    workdir=prop.workdir+"/"+subdir 
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    qcdir=workdir+"/qc/"+prefix
    all_vcf_fpaths=glob.glob(vcf_file_pattern)
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    for file_in in (all_vcf_fpaths): ### ? GFF file
        FI.copy_file_to_destdir(file_in,indir)
    FI.change_dir(workdir)    

def execute():
    print ("executing...")
    if mapping_file is not None:
        get_sample_runID_dict()
    run_snpEff()
    get_go()
    analysis_ann_files()
    get_genes_from_gff()
    write_gene_matrix()
    write_gene_matrix_summary()
    write_posi_file()

def post_process():
    print ("post_processing...")
    for summary_file in (gene_matrix_file, gene_matrix_summary_file, gene_posi_file):
        FI.copy_file_to_destdir(summary_file, outdir)        
    for ann_stat_fpath in ann_stat_fpaths:
        multiQC_ann_stat_fpath = ann_stat_fpath.replace(".ann_stats",".txt")
        FI.copy_file(ann_stat_fpath,"{}/{}".format(qcdir,multiQC_ann_stat_fpath))

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()
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

