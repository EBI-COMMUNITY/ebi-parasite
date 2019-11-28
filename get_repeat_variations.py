#! /usr/bin/env python3

## version 1, 28 Nov 2019 by Xin Liu

import argparse
import re
import sys
import os
import glob
from Bio.Seq import Seq
from utilities import fileutils
from utilities import properties
from utilities import misc
from utilities import command
from collections import defaultdict
from repeats_variation import strv
import numpy as np

def get_args():    
    global properties_file
    global genome_name
    global vcf_file_pattern
    global mapping_file
    global prefix
    global prop
    global vcf_files
    global vcf_request_pattern
   
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='''Script creates the genome Short Tandem Repeat (STR) 
                                                    variation summary file based on all vcf files and 
                                                    multiple alignment files for all repeat regions''')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file.', 
                                                             required=True)
    parser.add_argument('-g', '--genome_name', type=str, help='''Please provide the genome name available in 
                                                                 genome_list.txt only or cryptosporidium_hominis''', 
                                                             required=False) 
    parser.add_argument('-vp', '--vcf_file_pattern', type=str, help='''Please provide the paired vcf files' pattern 
                        with the full path, must end with .vcf''', required=True)     
    parser.add_argument('-m', '--mapping_file', type=str, help='''Please provide the mapping file path, which 
                                                                  contains one column of read_ID from vcf file 
                                                                  and one column of its corresponding sample_name''',
                                                             required=False)
    parser.add_argument('-pre', '--prefix', type=str, help='Please provide the prefix for the output file.', 
                                                             required=True)  

    # check args
    args = parser.parse_args()
    FI.check_exist(args.properties_file)
    properties_file=args.properties_file
    prop=properties(properties_file)    
    if args.genome_name not in (line.rstrip() for line in open(prop.get_attrib("available_genomes")).readlines()):
        MISC.my_exit("{} is not available, please try another genome".format(args.genome_name)) 
    vcf_files=glob.glob(args.vcf_file_pattern)
    if len(vcf_files) == 0:
        sys.exit("ERROR: no vcf files provided")
    vcf_request_pattern="^.*/?(.*?).vcf$"
    for vcf_file in vcf_files:
        if not re.search(vcf_request_pattern, vcf_file):
            sys.exit("vcf_file not end with .vcf")
    FI.check_files_exist(vcf_files) 
    
    # define variables     
    properties_file = args.properties_file    
    prop = properties(properties_file)    
    genome_name = args.genome_name           
    vcf_file_pattern = args.vcf_file_pattern
    prefix = args.prefix   
    mapping_file = args.mapping_file
    
    print ("properties_file:",properties_file)
    print ("genome_name:",genome_name)
    print ("vcf_file_pattern:",vcf_file_pattern)
    print ("mapping_file:",mapping_file)
    print ("prefix:",prefix) 
    
def get_chrom_len():
    global chrom_len_dict
    chrom_len_dict = {}
    chro_ori_lines = command("grep '>' {}".format(REF_FASTA)).run_comm(1).decode("utf-8").rstrip().split("\n")
    for chro_ori_line in chro_ori_lines:
        (chrom, chr_len) = re.findall(">(.*?) .*?length=(\d+) ", chro_ori_line)[0]
        chrom_len_dict[chrom] = chr_len

def mapping_sample():
    global mapping_dict
    mapping_dict = {}
    for vcf_file in vcf_files:
        runID=MISC.get_runID(vcf_file)
        sample=MISC.get_samples_by_runIDs(mapping_file)[runID]
        mapping_dict[runID]=sample

def get_info_from_gff3():
    global gff3_dict
    global gene_bed_fPath
    gff3_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    gene_bed_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(str))))
    REF_GFF = "{}/../ref/{}.gff".format(workdir, genome_name)
    fh_gff3 = open(REF_GFF, "r") 
    gene_bed_fPath = "{}/gene.bed".format(workdir)
    tmp_fPath = "{}/tmp.bed".format(workdir)
    fh_tmp = open(tmp_fPath, "w")
    for line in fh_gff3:
        if not re.search('^#',line):
            (chrom,mol_type,start,end,strand,phase,gene_part_ori)=getVar(line.split(),[0,2,3,4,6,7,8])
            if mol_type == 'gene':
                gene_name = re.findall("ID=(.*?);",gene_part_ori)[0]
                fh_tmp.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(chrom, start, end, strand, phase, gene_name))
            if mol_type == 'CDS':
                gene_name = re.findall("Parent=(.*?)-",gene_part_ori)[0]
            if mol_type == 'gene' or mol_type == 'CDS':
                gff3_dict[chrom][mol_type][(start,end)] = gene_name
    fh_gff3.close()
    fh_tmp.close()
    command("sort -k1,1 -k2,2n {} > {}".format(tmp_fPath, gene_bed_fPath)).run_comm(0)
    command("rm {}".format(tmp_fPath)).run_comm(0)

def get_60Nstr():
    out_str=""
    for i in range(0,60):
        out_str+="N"
    return out_str

def get_up_start(ori_start):
    up_start = int(ori_start) - UP_DOWN_SEQ_LEN
    if up_start <= 0:
        up_start = 1
    return str(up_start)

def get_down_end(ori_end, chomo_len):
    down_end = int(ori_end) + UP_DOWN_SEQ_LEN
    if down_end > int(chomo_len):
        down_end = chomo_len
    return str(down_end)

def check_seq_hunN(region):
    command("echo {} >{}".format(region, tmp_region_fname)).run_comm(0)
    region_fasta_str_ori = command("samtools faidx -r {} {}".format(tmp_region_fname, REF_FASTA)).run_comm(1).decode("utf-8").rstrip()
    region_fasta_str = region_fasta_str_ori.replace("\n","")
    print ("region_fasta_str="+region+"\n"+region_fasta_str)
    if re.search(sixtyNstr, region_fasta_str):
        print ("here 60N:" + region)
        return 1
    else:
        return 0

def get_exp_ref_region(ref_chrom, ref_start, ref_end, ref_chrom_len):
    solid_exp_start = get_up_start(ref_start)
    solid_exp_end = get_down_end(ref_end, ref_chrom_len)
    solid_upstream_region = "{}:{}-{}".format(ref_chrom, solid_exp_start, ref_start)
    solid_downstream_region = "{}:{}-{}".format(ref_chrom, ref_end, solid_exp_end)
    if check_seq_hunN(solid_upstream_region): ## if solid_upstream are 100N
        exp_ref_region_start = ref_start
        print ("upstream 100N: "+ref_chrom+", "+ref_start+", "+ref_end)
    else:
        exp_ref_region_start = solid_exp_start
    if check_seq_hunN(solid_downstream_region): ## if solid_upstream are 100N
        exp_ref_region_end = ref_end
        print ("100N: "+ref_chrom+", "+ref_start+", "+ref_end)
    else:
        exp_ref_region_end = solid_exp_end
    exp_ref_region = "{}:{}-{}".format(ref_chrom, exp_ref_region_start, exp_ref_region_end)
    exp_ref_region_len = str(int(exp_ref_region_start) - int(exp_ref_region_end) + 1)
    return (exp_ref_region, exp_ref_region_len)

def get_region_ann(region):
    [str_chro,str_start,str_end] = re.findall("^(.*?):(.*?)-(.*?)$", region)[0]
    print ("region: {} {} {}".format(str_chro,str_start,str_end))
    out_geneName = "N/A"
    out_mol_type = "N/A"
    if_found = 0
    ## check all CDSs in this chromosome to see whether the str fall into any CDS region
    for cds_region in gff3_dict[str_chro]['CDS'].keys():
        (cds_start,cds_end) = cds_region
        if int(str_start) >= int(cds_start) and int(str_end) <= int(cds_end):
            out_geneName = gff3_dict[str_chro]['CDS'][cds_start,cds_end]
            out_mol_type = "CDS({}-{})".format(cds_start, cds_end)
            if_found = 1
            break
    ## if str not found in and CDS, check all genes in this chromosome
    if not if_found:
        for gene_region in gff3_dict[str_chro]['gene'].keys():
            (gene_start,gene_end) = gene_region
            if int(str_start) >= int(gene_start) and int(str_end) <= int(gene_end):
                out_geneName = "{}({}-{})".format(gff3_dict[str_chro]["gene"][gene_start,gene_end], gene_start,gene_end)
                out_mol_type = "non_CDS"
                break
    return (out_geneName,out_mol_type)

def get_closest_region_and_gene(tmp_bed_fPath):
    global closest_dict
    closest_dict = {}
    ref_regions_bed_fPath = "{}/ref_regions1.sorted.bed".format(workdir)
    ref_regions_bed_cp_fPath = "{}/ref_regions2.sorted.bed".format(workdir)
    command("sort -k 1,1 -k2,2n {} > {}".format(tmp_bed_fPath, ref_regions_bed_fPath)).run_comm(0)
    command("rm {}".format(tmp_bed_fPath)).run_comm(0)
    command("cp {} {}".format(ref_regions_bed_fPath, ref_regions_bed_cp_fPath)).run_comm(0)
    cmd = "bedtools closest -a {} -b {} -io -D a | bedtools closest -a - -b {} -io -D b | cut -f1-3,7-10,13,14-16,19,20"
    closest_out_lines = command(cmd.format(ref_regions_bed_fPath,
                                           ref_regions_bed_cp_fPath,
                                           gene_bed_fPath)).run_comm(1).decode("utf-8").rstrip().split("\n")
    for line in closest_out_lines:
        (ref_chrom, ref_start, ref_end) = getVar(line.split("\t"), [0,1,2])
        ref_region = "{}:{}-{}".format(ref_chrom, ref_start, ref_end)
        closest_region_and_gene = re.findall(".*?\t.*?\t.*?\t(.*?)$", line)[0]
        closest_dict[ref_region] = closest_region_and_gene
    
def ref_region_analsis():
    global ref_exp_region_dict
    global ref_exp_region_fasta
    global ann_dict
    global sixtyNstr
    sixtyNstr = get_60Nstr()
    ref_exp_region_dict = {}
    ann_dict = {}
    # tandem repeats finder to detect tandem repeat in genome fasta file,output is TRF_out_name 
    # parameter details, see *1, output table explanation, see *2
    command("trf {} 2 7 7 80 10 50 500 -f -d -m".format(REF_FASTA)).run_comm(0)
    # convert it to the format that STRViper will read in, output explanation see *2, not only parsing, 
    # but also skip some regions as better quality regions related
    ## define file names
    ref_fasta_pattern = "^(./)?(.*?)$"
    ref_fasta_fname = re.findall(ref_fasta_pattern, REF_FASTA)[0][1]
    trf_out_name = ref_fasta_fname+".2.7.7.80.10.50.500.dat"
    fName_str = "{}.trf.str".format(genome_name)
    exp_ref_region_fname = "{}.ref.exp.regions".format(prefix)
    fhout_exp_ref_region = open (exp_ref_region_fname, "w")
    ref_exp_region_fasta = "{}.ref.exp.fasta".format(prefix)
    tmp_bed_fPath = "{}/tmp.bed".format(workdir)
    fh_bed_tmp = open(tmp_bed_fPath, "w")
    ## run command    
    command("jsat parseTRF --input {} --output {} --format str".format(trf_out_name,fName_str)).run_comm(0)
    cmd_str = "grep -v '^#' {} | grep -v '^$' | awk '{{print $1\" \"$2\" \"$3\" \"$4\" \"$5}}'"
    ref_regions_str = command(cmd_str.format(fName_str)).run_comm(1).decode("utf-8").rstrip().split("\n")
    for ref_region_str in ref_regions_str:
        (ref_chrom, ref_start, ref_end, ref_str_unit_len, ref_unit_num_ori) = getVar(ref_region_str.split(),[0,1,2,3,4])
        ref_unit_num = str(int(float(ref_unit_num_ori)))
        (exp_ref_region, exp_ref_region_len) = get_exp_ref_region(ref_chrom, ref_start, ref_end, chrom_len_dict[ref_chrom])
        ref_region = "{}:{}-{}".format(ref_chrom, ref_start, ref_end)
        ref_len = str(abs(int(ref_start)-int(ref_end)) + 1)
        cmd = "grep '^[0-9]' {} | awk '{{if ($1=={} && $2=={} && $3=={}) print $14}}'"
        ref_str_seq = command(cmd.format(trf_out_name,ref_start,ref_end,ref_str_unit_len)).run_comm(1).decode("utf-8").rstrip()
        ref_exp_region_dict[exp_ref_region] = (exp_ref_region_len, ref_region, ref_str_seq, ref_str_unit_len, ref_len, ref_unit_num)
        fhout_exp_ref_region.write(exp_ref_region + "\n")
        ann_dict[ref_region] = get_region_ann(ref_region)
        print ("ann_dict1={} {}".format(ann_dict[ref_region][0], ann_dict[ref_region][1]))
        fh_bed_tmp.write("{}\t{}\t{}\t{}:{}:{}\t.\t.\n".format(ref_chrom, ref_start, ref_end,
                                                               ref_str_seq, ref_str_unit_len,ref_unit_num))
    fhout_exp_ref_region.close()
    fh_bed_tmp.close()
    command("samtools faidx -r {} {} -o {}".format(exp_ref_region_fname, REF_FASTA, ref_exp_region_fasta)).run_comm(0)
    get_closest_region_and_gene(tmp_bed_fPath)

def get_iso_region(ref_e_region,ref_region, iso_e_region):
    [ref_region_chr, ref_e_start, ref_e_end] = re.split(':|-', ref_e_region)
    [ref_region_chr, ref_region_start, ref_region_end] = re.split(':|-', ref_region)
    [iso_chrom,iso_e_start,iso_e_end] = re.split(':|-', iso_e_region)
    upstream_len=int(ref_region_start)-int(ref_e_start)
    downstream_len=int(ref_e_end)-int(ref_region_end)
    if upstream_len <0 or downstream_len<0:
        sys.exit("upstream_len="+str(upstream_len)+"downstream_len="+str(downstream_len))
    return "{}:{}-{}".format(iso_chrom, 
                             str(int(iso_e_start)+upstream_len), 
                             str(int(iso_e_end)-upstream_len))


def exp_region_analysis(vcf_fname_key, blast_out_fname, consensus_fasta_name):
    print ("blast_out_fname="+blast_out_fname)
    for ref_exp_region in ref_exp_region_dict.keys():
        if_multi_best_hit = "N"
        ## get best hit line
        cmd = "grep {} {} |grep -v '#' ".format(ref_exp_region, blast_out_fname)
        blast_out_line_ori = command(cmd).run_comm(1).decode("utf-8").rstrip()
        blast_out_lines = blast_out_line_ori.split("\n")
        if len(blast_out_lines) > 1:
            if_multi_best_hit = "Y"
        blast_out_line = blast_out_lines[0]
        print (blast_out_line)
        ## get expanded regions and isolates' best hit info 
        (ref_e_region, iso_chrom, iso_e_start, iso_e_end, iso_identity, align_len) = getVar(blast_out_line.split(), [0, 1, 8, 9, 2, 3])
        [ref_region_chr, ref_e_start, ref_e_end] = re.split(':|-', ref_e_region)
        
        ## ref str info
        (ref_region, ref_str_seq, str_unit_len, ref_len, ref_unit_num) = getVar(ref_exp_region_dict[ref_e_region], [1,2,3,4,5])
        ref_e_region_len = abs(int(ref_e_start) - int(ref_e_end)) + 1
        expanded_len = ref_e_region_len - int(ref_len)
        iso_repeat_region_len = abs(int(iso_e_start) - int(iso_e_end)) + 1 - expanded_len
        iso_unit_num = int(int(iso_repeat_region_len)/int(str_unit_len))
        if iso_unit_num < 0:
            iso_unit_num = "0"
        else:
            iso_unit_num_str = str(iso_unit_num)
        iso_e_region = "{}:{}-{}".format(iso_chrom, iso_e_start, iso_e_end)
        ## get iso_region
        iso_region = get_iso_region(ref_e_region,ref_region, iso_e_region)
        ## annotation
        (gene_name, cds) = ann_dict[ref_region]
        ref_key = "{}\t{}\t{}\t{}\t{}\t{}\t{}".format(ref_region, ref_str_seq, str_unit_len, 
                                                      ref_unit_num, ref_len, gene_name, cds)
        [ref_region_fasta_seq,ref_region_seq]=get_region_seq(ref_region,REF_FASTA)
        iso_region_seq = get_region_seq(iso_region,consensus_fasta_name)[1]
        print ("iso_region_seq="+iso_region_seq+"\nref_region_seq="+ref_region_seq)
        if ref_region_seq != iso_region_seq:
            ## diff_out_indi_dict holding all individaul repeat regions info, 
            ## including all ref regions and iso unit number and vcf file name
            diff_out_indi_dict[ref_key][vcf_fname_key] = (iso_unit_num_str, if_multi_best_hit)
        multi_align_indi_dict[(ref_region, ref_exp_region)][vcf_fname_key] = (iso_e_region, consensus_fasta_name)

def vcf_analysis():
    global multi_align_indi_dict
    iso_name = "NA"
    multi_align_indi_dict = defaultdict(lambda: defaultdict(str))
    for vcf_fpath in (vcf_files):
        vcf_fname_key = re.findall("(.*)\.", os.path.basename(vcf_fpath))[0]
        compressed_vcf_name = vcf_fname_key+".bgzip"
        consensus_fasta_name = vcf_fname_key+".consensus.fasta"
        command("bgzip -c {} > {}".format(vcf_fpath, compressed_vcf_name)).run_comm(0)
        command("tabix {} -f".format(compressed_vcf_name)).run_comm(0)
        command("cat {} | bcftools consensus {} > {}".format(REF_FASTA, compressed_vcf_name,
                                                       consensus_fasta_name)).run_comm(0)
        command("makeblastdb -in {}  -parse_seqids -dbtype nucl -out {}.DBblast".format(consensus_fasta_name,
                                                                               vcf_fname_key)).run_comm(0)
        blast_out_fname = vcf_fname_key+".ref.exp.blast_out"
        command("blastn -query {} -db {}.DBblast -dust no -outfmt 7 -max_target_seqs 1 -out {}".format(ref_exp_region_fasta, vcf_fname_key, blast_out_fname)).run_comm(0)
        for map_key in mapping_dict.keys():
            if re.search(map_key, consensus_fasta_name):
                iso_name = mapping_dict[map_key]
        if iso_name == "NA":
            iso_name = re.findall("(.*?)_IbA10G2", consensus_fasta_name)[0]
        exp_region_analysis(vcf_fname_key, blast_out_fname, consensus_fasta_name)

def get_diff_out_dict():
    global diff_out_dict
    global multi_align_dict
    diff_out_dict = defaultdict(lambda: defaultdict(str))
    multi_align_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(str)))
    ## get average diff
    ## "variation_num" is how many different repeat region lenghthes were found
    for key1 in diff_out_indi_dict.keys():
        diff_sum = 0.0
        out_isos = ""
        if_multi_best_hit_str = ""
        if_multi_best_hit_dict = {}
        multi_best_hit_N_num = 0
        multi_best_hit_Y_num = 0
        num_of_variation_dict = {}
        for vcf_file_key in diff_out_indi_dict[key1].keys():
            [iso_unit_num_str, multi_best_hit] = diff_out_indi_dict[key1][vcf_file_key]
            num_of_variation_dict[iso_unit_num_str] = 1
            if  multi_best_hit == 'N':
                multi_best_hit_N_num += 1
            else: 
                multi_best_hit_Y_num += 1
        variation_num = len(num_of_variation_dict.keys())
        variation_num_str = str(variation_num)
        multi_best_hit_N_str = str(multi_best_hit_N_num)
        multi_best_hit_Y_str = str(multi_best_hit_Y_num)
        closest_info = closest_dict[key1.split("\t")[0]] ## closest_dict[ref_region]
        for vcf_file_key in sorted(diff_out_indi_dict[key1].keys()):
            (repeat_len, multi_best_hit) = diff_out_indi_dict[key1][vcf_file_key]
            diff_out_dict[key1 + "\t" + closest_info][vcf_file_key] = (repeat_len, variation_num_str, 
                                                                       multi_best_hit_N_str, 
                                                                       multi_best_hit_Y_str)

def write_output():
    global output_summary_fpath
    output_summary_fpath = "{}/{}.summary".format(workdir, prefix)
    head_line_format = ""
    for i in range (0, 18):
        head_line_format += "{}\t"
    head_line_format = head_line_format.rstrip("\t")
    head_line = head_line_format.format("ref_str_region", "str_unit_seq", "ref_str_unit_len", 
                                        "ref_str_unit_num", "ref_str_region_len",
                                    	"gene_name", "cds",
                                        "closest_repeat_chromosome", "closest_repeat_start", "closest_repeat_end",
                                        "closest_repeat_unit_seq:unit_length:unit_num", "closest_repeat_distrance",
                                        "closest_gene_chromosome", "closest_gene_start", "closest_gene_end",
                                        "closest_gene_name", "closest_gene_distrance",
                                        "var_iso_num" 
                                        )
    fhout = open(output_summary_fpath, "w")
    fhout.write(head_line+"\n")
    for key1 in diff_out_dict.keys():
        var_iso_num = len(diff_out_dict[key1].keys())
        fhout.write("{}\t{}\n".format(key1, str(var_iso_num)))
    fhout.close()

def replace_the_first_line(ori_str, new_first_line):
    lines = ori_str.split("\n")
    new_lines = lines[1:]
    out_str = new_first_line
    for line in new_lines:
        out_str += "\n"+line
    return out_str

def get_seq_from_fasta(fasta_seq):
    return "\n".join(str(fasta_seq).split("\n")[1:])

def get_region_seq(region, seq_file):
    [chrom, region_start, region_end]=re.findall("(.*):(.*)-(.*)",region)[0]
    if int(region_start)>int(region_end):
        region="{}:{}-{}".format(chrom,region_end,region_start)
    command("echo {} >{}".format(region, tmp_region_fname)).run_comm(0)
    region_fasta_seq=command("samtools faidx -r {} {}".format(tmp_region_fname,
                          seq_file)).run_comm(1).decode("utf-8").rstrip()
    if int(region_start)<int(region_end):
        return [region_fasta_seq,get_seq_from_fasta(region_fasta_seq)]
    else:
        region_fasta_revCom_seq=Seq(region_fasta_seq).reverse_complement()
        print ("region_fasta_revCom_seq="+region_fasta_revCom_seq)
        return [region_fasta_revCom_seq,get_seq_from_fasta(region_fasta_revCom_seq)]

def multi_align():
    #eon_info_key = "{}_{}N{}Y".format(ref_exp_region, multi_best_hit_N_str, multi_best_hit_Y_str)
    # multi_align_indi_dict[ref_region][vcf_fname_key] = [ref_exp_region, iso_region, consensus_fasta_name]
    for ref_info in multi_align_indi_dict.keys():
        (ref_region, ref_e_region) = ref_info
        print ("region2: "+ ref_region +" " + ref_e_region)
        out_fasta_fname = "{}.fasta".format(ref_region)
        multi_align_outF = "{}_{}.clustalo_num".format(prefix, ref_region)
        ## write ref to fasta file
        command("echo {} >{}".format(ref_e_region, tmp_region_fname)).run_comm(0)
        command("samtools faidx -r {} {} >{}".format(tmp_region_fname, REF_FASTA, out_fasta_fname)).run_comm(0)
        command("sed -i \"1s/.*/>ref_{}/\" {}".format(ref_e_region, out_fasta_fname)).run_comm(0)
        for vcf_file_key in multi_align_indi_dict[ref_info].keys():
            (iso_e_region, iso_concensus_fasta) = multi_align_indi_dict[ref_info][vcf_file_key]
            ## write isolates to fasta file
            iso_fasta_lines_str = get_region_seq(iso_e_region,iso_concensus_fasta)[0]
            iso_fasta_lines = replace_the_first_line(iso_fasta_lines_str, ">{}_{}".format(vcf_file_key, iso_e_region))
            command("echo \"{}\" >> {}".format(iso_fasta_lines, out_fasta_fname)).run_comm(0)
        ## do alignment for each expanded region
        command("clustalo --infile {} --threads 8 --verbose --outfmt clustal --resno --outfile {} --output-order input-order --seqtype dna --force"
                .format(out_fasta_fname, multi_align_outF)).run_comm(0)
 
def run(): 
    global UP_DOWN_SEQ_LEN
    global REF_FASTA 
    global tmp_region_fname
    global diff_out_indi_dict
    diff_out_indi_dict = defaultdict(lambda: defaultdict(str))
    tmp_region_fname = prefix+".tmp.region"
    UP_DOWN_SEQ_LEN = 100
    REF_FASTA = "{}/../ref/{}.fasta".format(workdir, genome_name)
    get_chrom_len() # get all individual chr len
    mapping_sample() 
    get_info_from_gff3() # get gff3_dict[chrom][mol_type][(start,end)] = gene_name
    ref_region_analsis()
    vcf_analysis()
    get_diff_out_dict()
    write_output()
    multi_align()

def initiate():
    print ("initiating...")
    global workdir
    global indir
    global outdir
    subdir="repeats"
    workdir=prop.workdir+"/"+subdir
    indir=workdir+"/in/"
    outdir=workdir+"/out/"
    FI.create_processing_dir(indir)
    FI.create_processing_dir(outdir)
    FI.change_dir(workdir)

def post_process():
    print ("post_processing...")
    # create multi alignment to the working dir and the copy to out/multi_align
    FI.copy_file_to_destdir(output_summary_fpath,outdir)
    clustalo_num_files = glob.glob("{}_*.clustalo_num".format(prefix))
    multi_align_dir = "{}/multi_align".format(outdir)
    FI.make_dir(multi_align_dir)
    for clustalo_num_file in clustalo_num_files:
        FI.copy_file_to_destdir("{}/{}".format(workdir, clustalo_num_file), multi_align_dir)

if __name__ == '__main__':
    getVar = lambda searchList, ind: [searchList[i] for i in ind]    
    global FI
    global MISC
    FI=fileutils()
    MISC=misc()
    get_args()
    print ("\n","Properties attributes:")
    print (prop.__dict__)
    #submit_str_jobs()
    #get_info_from_gff3()
    #summary2()
    #summary1()
    initiate()
    run()
    post_process()
