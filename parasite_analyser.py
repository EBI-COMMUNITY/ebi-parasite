import sys
import argparse
import os
import subprocess
from lxml import etree
import xml.etree.ElementTree as ET
from utilities import properties

#Parameters:
#taxonomic division, -d virus/bacteria
#Assembler: spade, velvet
#read files 
#properties in case
#
#
#
#
#
#
#



def get_args():

    global properties_file
    global fastq1
    global fastq2
    global assemblers
    global division
 
    # Assign description to the help doc
    parser = argparse.ArgumentParser(
        description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=False)
    parser.add_argument('-d', '--division', type=str, help='Please provide the taxonomic division, virus or bacteria', required=False)
    parser.add_argument('-a', '--assembly_programs', type=str, help='Please provide the assembler, velvet or spade', required=False)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the first fastq file', required=False)
    
    
    args = parser.parse_args()
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    division=args.division
    if args.assembly_programs==None:
        if division==None:
            assemblers=['velvet','spade']
        else:
            if division=='virus':
                assemblers=['spade']
            elif division=='bactria':
                assemblers=['velvet']
    else:
        assemblers = args.assembly_programs.split(",")
        print "TODO: split the asssemblies and put it in assemblers"



def get_properties(property_file):
		with open(property_file) as f:
			 lines = f.readlines()
		
		global path_to_quasr
		global path_to_usearch 
		global path_to_spades 
		global data_base_list_file 
		global illumina_core 
		global biopython
                global teilenq
                global velvetoptimiser	
	
		path_to_quasr_provided = False
		path_to_usearch_provided = False
		path_to_spades_provided = False
		data_base_list_file_provided = False
		illumina_core_provided = False
		biopython_provided=False		
                teilenq_provided=False
                velvetoptimiser_provided=False	
		
		for l in lines:
			pair = l.strip().split(":")
			if pair[0].lower() == 'path_to_quasr':
				path_to_quasr = pair[1]
				path_to_quasr_provided = True
			elif pair[0].lower() == 'path_to_usearch':
				path_to_usearch = pair[1]
				path_to_usearch_provided = True
			elif pair[0].lower() == 'path_to_spades':
				path_to_spades = pair[1]
				path_to_spades_provided = True
			elif pair[0].lower() == 'data_base_list_file':
				data_base_list_file = pair[1]
				data_base_list_file_provided = True
			elif pair[0].lower() == 'illumina_core':
				illumina_core = pair[1]
				illumina_core_provided = True
			elif pair[0].lower() == 'biopython':
				biopython = pair[1]
				biopython_provided = True
                        elif pair[0].lower() == 'teilenq':
                                teilenq = pair[1]
                                teilenq_provided = True
					
		if path_to_quasr_provided == False:
			print "ERROR: path_to_quasr missing"
			sys.exit(1)
		if path_to_usearch_provided == False:
		   print "ERROR: path_to_usearch missing"
		   sys.exit(1)
		if path_to_spades_provided == False:
		   print "ERROR: path_to_spades missing"
		   sys.exit(1)
		if data_base_list_file_provided == False:
		   print "ERROR: data_base_list_file missing"
		   sys.exit(1)
		if illumina_core_provided == False:
		   print "ERROR: llumina_core missing"
		   sys.exit(1)
		if biopython_provided == False:
		   print "ERROR: biopython missing"
		   sys.exit(1)
                if teilenq_provided == False:
                   print "ERROR: teilenq missing"
                   sys.exit(1)


def run(command):
		print "running the command"
		print command
		sp = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
		
		out, err = sp.communicate()
		if out:
			print "standard output of subprocess:"
			print out
			data=out.split('\n')
			i=0
			for line in data:
				if 'error' in line.lower(): 
					message=data[i-1]+'\n'+data[i]
					self.error_list.append(message.replace("'",""))
				i=i+1

		if err:
			print "standard error of subprocess:"
			print err
			data=err.split('\n')
			i=0
			for line in data:
				if 'error' in line.lower():
					message=data[i-1]+'\n'+data[i]
					self.error_list.append(message.replace("'",""))
				i=i+1
		if sp.returncode!=0:
			self.error_list.append(err.replace("'",""))
		print >> sys.stderr, err



if __name__ == '__main__':
    
    get_args()
    print assemblers
    print fastq1
    print fastq2
    #command=''
    #print command 
    #run(command)



