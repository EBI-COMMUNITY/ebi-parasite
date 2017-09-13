import argparse
from utilities import command
from utilities import fileutils



def get_args():

    global properties_file
    global fastq1
    global fastq2
 
    # Assign description to the help doc
    parser = argparse.ArgumentParser(description='Script assembles short reads based on some criteria')
    parser.add_argument('-p', '--properties_file', type=str, help='Please provide the properties file', required=True)
    parser.add_argument('-fq1', '--fastq1', type=str, help='Please provide the first fastq file', required=True)
    parser.add_argument('-fq2', '--fastq2', type=str, help='Please provide the first fastq file', required=False)
    
    
    args = parser.parse_args()
    properties_file=args.properties_file
    fastq1=args.fastq1
    fastq2=args.fastq2

    print "properties_file:",properties_file
    print "fastq1:",fastq1
    print "fastq2:",fastq2

   
def run_trim_galore(fastqfiles):
        
        #pair fastq files
		if len(fastqfiles)==2:
			comm=self.prop.trim_galore+" --paired -q 20 "+ fastqfiles[0]+" "+fastqfiles[1]
        #single fastq files    
		elif len(fastqfiles)==2:
			comm=self.prop.trim_galore+" -q 20 "+ fastqfiles[0]
        else
            print "unkown number of fq files."
		print comm
		comm_obj=command(comm)
		returncode, stdout, stderr=comm_obj.run(3600)
		print returncode, stdout, stderr
		print self.prop.workdir
		


def initiate():
    #Here please provide the initiation code
     print "Here please do the initiation code"
     #fi=fileutils()
	 #indir=self.prop.workdir+"/quality/in/"
     #outdir=self.prop.workdir+"/quality/out/"
	 #fi.create_processing_dir(indir)
	 #fi.create_processing_dir(outdir)
	 #fi.copy_src_into_dest(self.fq1,indir)
	 #fi.copy_src_into_dest(self.fq2,indir)


def execute():
    #Here please provide the execucation code
    print "Here please do the execucation  code"
    run_trim_galore()
    
    #Execution of the external softwares, you need to use command class from utilities
    #Example: from utilities import command
    #Then you can create an object of command class and then use a function:
    #Example:
    #comm = command("python /usr/bin/spades.py -1 8605_7_1_test.fastq -2 8605_7_2_test.fastq -k 21,33,55,77,99,127 --careful --only-assembler -o sp    ades_small_dir")
    #returncode,stdout,stderr=command.run(timeout=5)

def post_process():
    #Here please provide the post process code
    print "Here please do the post process  code"
    	fi=fileutils()
		indir=self.prop.workdir+"/quality/in/"
		outdir=self.prop.workdir+"/quality/out/"
		fi.create_processing_dir(indir)
		fi.create_processing_dir(outdir)
		fi.copy_src_into_dest(self.fq1,indir)
		fi.copy_src_into_dest(self.fq2,indir)
		fqout1=self.fq1+"_val_1.fq"
		fqout2=self.fq2+"_val_2.fq"
		report=self.fq1+"_trimming_report.txt"
		fi.copy_src_into_dest(fqout1,indir)
		fi.copy_src_into_dest(fqout2,indir)
		fi.copy_src_into_dest(report,indir)



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
   
     
