import argparse
import sys
from utilities import command
from utilities import fileutils
from utilities import properties


def get_args():    
    global properties_file
    global fastq1
    global fastq2
    global qcs
    default_qcs="trim_galore"
    
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the first fastq file', required=False)
    parser.add_argument('-qcs', '--qc_software', type=str, help='Please provide the quality control software', required=False)
    
    args = parser.parse_args()
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2
    qcs=args.qc_software    
    if qcs is None:
        qcs = default_qcs

    print "properties_file:",properties_file
    print "fastq1:",fastq1
    print "fastq2:",fastq2
    print "qcs:",qcs
   
   
def run_trim_galore(fastqfiles):            
    #check whether qc_software path was defined in the property file
    try:
        qcs_path=getattr(prop, qcs) 
    except AttributeError:
        print "ERROR: {} path was not defined in {}. Please define it and re-run.".format(qcs,properties_file)        
        sys.exit(1)        
    #pair fastq files    
    if len(fastqfiles)==2:
        comm=qcs_path+" --paired -q 20 "+fastqfiles[0]+" "+fastqfiles[1]
    #single fastq files    
    elif len(fastqfiles)==1:
        comm=qcs_path+" -q 20 "+fastqfiles[0]
    else:
        print "unkown number of fq files."
    comm_obj=command(comm)
    returncode, stdout, stderr=comm_obj.run(3600)
    print returncode, stdout, stderr
    print prop.workdir
   
    
def initiate():
    print "initiating..."
    global indir
    global outdir
    fi=fileutils()
    indir=prop.workdir+"/quality/in/"
    outdir=prop.workdir+"/quality/out/"
    fi.create_processing_dir(indir)
    fi.create_processing_dir(outdir)
    fi.copy_file_into_dest(fastq1,indir)
    if fastq2 is not None:
        fi.copy_file_into_dest(fastq2,indir)


def execute():
    print "executing..."
    fastqfiles=[]
    fastqfiles.append(fastq1)
    if fastq2 is not None:
        fastqfiles.append(fastq2)
    run_trim_galore(fastqfiles)


def post_process():
    print "post_processing..."
    fi=fileutils()
    fqout1=fastq1+"_val_1.fq"
    if fastq2 is not None:
        fqout2=fastq2+"_val_2.fq"
    report=fastq1+"_trimming_report.txt"
    fi.copy_file_into_dest(fqout1,outdir)
    if fastq2 is not None:
        fi.copy_file_into_dest(fqout2,outdir)
    fi.copy_file_into_dest(report,outdir)
    


if __name__ == '__main__':
    get_args()
    global prop
    prop=properties(properties_file)
    print "\n","Properties attributes:"
    print prop.__dict__
    
    #run the initiation code
    initiate()
 
    #execute the main part of the program
    execute()

    #post execution code
    post_process()
   
     
