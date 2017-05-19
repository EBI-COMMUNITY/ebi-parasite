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
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
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


def initiate():
	

def execute_stage1():
	print "stage 1 has started!"


def execute_stage2():
    print "stage 2 has started!"

def execute_stage3():
    print "stage 3 has started!"

def execute_stage4():
    print "stage 4 has started!"

def execute_stage5():
    print "stage 5 has started!"

def execute_stage6():
    print "stage 6 has started!"

def execute_stage7():
    print "stage 7 has started!"

def execute_stage8():
	print "stage 8 has started!"




if __name__ == '__main__':
    
    get_args()

    execute_stage1()
    execute_stage2()
    execute_stage3()
    execute_stage4()
    execute_stage5()
    execute_stage6()
    execute_stage7()
    execute_stage8()

    print assemblers
    print fastq1
    print fastq2
    #command=''
    #print command 
    #run(command)



