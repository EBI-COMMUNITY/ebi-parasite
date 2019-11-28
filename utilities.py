from os import kill
import os
import sys
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen
import shutil
import re
import time

'''
How to run it:
You need to first import it this Class into your code
Example: from utilities import command
Then you can create an object of command class and then use a function:
Example:
command = Command("python /usr/bin/spades.py -1 8605_7_1_test.fastq -2 8605_7_2_test.fastq -k 21,33,55,77,99,127 --careful --only-assembler -o spades_small_dir")
returncode,stdout,stderr=command.run(timeout=5)

'''

class command():
    def __init__(self, command):
        self.command = command
        print (self.command)

    def run(self, timeout = -1):
        env=None
        kill_tree = True
        class Alarm(Exception):
            pass    
        def alarm_handler(signum, frame):
            raise Alarm
        p=Popen(self.command, shell = True, stdout = PIPE, stderr = PIPE, env = env)
        if timeout != -1:
            signal(SIGALRM, alarm_handler)
            alarm(timeout)
        try:
            stdout, stderr = p.communicate()	
            if timeout != -1:
                alarm(0)
        except Alarm as a:
            pids = [p.pid]
            if kill_tree:
                pids.extend(self.get_process_children(p.pid))
            for pid in pids:
                # This is to avoid OSError: no such process in case process dies before getting to this line
                try: 
                    kill(pid, SIGKILL)
                except OSError:
                    pass
            return -1,'',''
        return p.returncode, stdout, stderr

    
    def run_comm(self,if_out_return):
        returncode, stdout, stderr=self.run(360000)
        if returncode and stderr:
            sys.stderr.write("\nERROR: {} FAILED!!! \n\nSTDERR: {}\nSTDOUT: {}\n".format(self.command,stderr,stdout))
            sys.exit(1)    
        if if_out_return:
            return stdout

    def run_comm_no_exit(self,if_out_return):
        returncode, stdout, stderr=self.run(360000)
        if returncode and stderr:
            sys.stderr.write("\nERROR: FAILED {} \n\nSTDERR: {}\nSTDOUT: {}\n".format(self.command,stderr,stdout))
        if if_out_return:
            return stdout


    def get_process_children(self,pid):
        p = Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
	         stdout = PIPE, stderr = PIPE)
        stdout, stderr = p.communicate()
        return [int(p) for p in stdout.split()]


class fileutils():
    def change_dir(self,path):
        try: 
            os.chdir(path)
        except OSError as err:
            sys.stderr.write("OS error: {}".format(err))
            raise

    def make_dir(self,path):        
        if not os.path.isdir(path):
            try:
                os.makedirs(path)
            except OSError as err:
                sys.stderr.write("OS error: {}".format(err))
                raise

    def file_to_string(self,path):
        try:
            self.check_exist(path)
            with open(path, 'r') as myfile:
                data = myfile.read().rstrip()                
            myfile.close()
            return data
        except IOError as e:
            sys.stderr.write("I/O error({0}): {1}".format(e.errno, e.strerror))
            raise
        except: 
            sys.stderr.write("Unexpected error:", sys.exc_info()[0])
            raise
        
    def create_processing_dir(self,directory):
        try: 
            if not os.path.exists(directory):
                os.makedirs(directory)
        except OSError as err:
            sys.stderr.write("OS error: {}".format(err))
            raise

    def rm_file(self,filein):
        try:
            os.system("rm "+filein)
        except OSError as err:
            sys.stderr.write("OS error: {}".format(err))
            raise

    def rename(self,filein,fileout):
        try:
            os.rename(filein, fileout)
        except OSError as err:
            sys.stderr.write("OS error: {}".format(err))
            raise

    def copy_file(self,srcfile,destfile):
        try: 
            shutil.copy(srcfile, destfile)
            print ("srcfile={}\ndestfile={}".format(srcfile, destfile))
            time.sleep(15)
        except OSError as err:
            sys.stderr.write('OS error: {}'.format(err))
            raise
			
    def copy_file_to_destdir(self,srcfile,destdir):
        name=os.path.basename(srcfile)
        dest_file=os.path.join(destdir, name)
        self.copy_file(srcfile,dest_file)	

    def add_file_prefix(self,source_fpath,prefix):    	
        try:
            dest_fpath=os.path.join(os.path.dirname(source_fpath),prefix+os.path.basename(source_fpath))
            shutil.move(source_fpath, dest_fpath)
        except OSError as err:
            sys.stderr.write('OS error: {}'.format(err)) 
            raise

    
    def copy_file_add_prefix(self,source_fpath,outdir,prefix):
        self.copy_file_to_destdir(source_fpath,outdir)
        source_file2=os.path.join(outdir,os.path.basename(source_fpath))		
        self.add_file_prefix(source_file2,prefix)
		
    
    def check_exist(self,path):
        if not os.path.exists(path):
            sys.stderr.write("ERROR: {} does not exist!\n".format(path))
            sys.exit(1)
        else:
            return 1


    def check_files_exist(self,fpaths):
        for fpath in fpaths:
            self.check_exist(fpath)
        return 1

    def unzip(self,fpath):
        unzip_fpath=fpath.rstrip(".gz")
        command("gunzip -c {} > {}".format(fpath,unzip_fpath)).run_comm(0)
        return unzip_fpath

class properties():
    def __init__(self, property_file):
        self.propertyf=property_file
        with open(property_file) as f:
            lines = f.readlines()
		
        for line in lines:
            if line:
                pair=line.strip().split(":")
                property_name=pair[0].lower()
                property_content=pair[1].strip("\n")				
            setattr(self, property_name, property_content)


    def get_attrib(self, attrib):
        if not hasattr(self, attrib):
           misc.my_exit("ERROR: {} path was not defined in {}. Please define it and re-run.\n".format(attrib,self.propertyf))
        else:
            fileutils().check_exist(getattr(self,attrib))
            return getattr(self,attrib)

class misc():

    def my_exit(self,error_message):
        sys.stderr.write(error_message+"\n")
        sys.exit(1)

    def download(self,download_dir,file_type,genome,genome_type): #file_type: gff or fasta; genome_type:DNA or cDNA
        fileutils().change_dir(download_dir)
        ftp_dir=""
        if genome.lower()!="cryptosporidium_hominis":
            ftp_root_dir="ftp://ftp.ensemblgenomes.org/pub/current/protists/{}/".format(file_type)
            sub_dirs1=command("curl -s {} | awk '{{print $9}}'".format(ftp_root_dir)).run_comm(1).split()
            if genome.lower() in sub_dirs1:
                if file_type=="gff3":
                    ftp_dir="{}/{}/".format(ftp_root_dir,genome.lower())
                else: #fasta
                    ftp_dir="{}/{}/{}/".format(ftp_root_dir,genome.lower(),genome_type)
            else:
                for sub_dir1 in sub_dirs1:
                    sub_dirs2=command("curl -s {}/ | awk '{{print $9}}'".format(ftp_root_dir+sub_dir1)).run_comm(1).split()
                    if genome.lower() in sub_dirs2:
                        if file_type=="gff3":
                            ftp_dir="{}/{}/{}/".format(ftp_root_dir,sub_dir1,genome.lower())
                        else:  #fasta
                            ftp_dir="{}/{}/{}/{}/".format(ftp_root_dir,sub_dir1,genome.lower(),genome_type)
            if ftp_dir=="":
                misc.my_exit("can not find DNA fasta file for {}".format(genome))
            else:
                if file_type=="gff3":
                    out_gz=command("curl -s {} | awk '{{print $9}}' | grep gff3.gz |grep -v chromosome".format(ftp_dir)).run_comm(1);
                else:#fasta
                    if genome_type=="dna":
                        out_gz=command("curl -s {} | awk '{{print $9}}' | grep dna.toplevel.fa.gz".format(ftp_dir)).run_comm(1);
                    else:
                        out_gz=command("curl -s {} | awk '{{print $9}}' | grep cdna.all.fa.gz".format(ftp_dir)).run_comm(1);
                out_gz=out_gz.rstrip()
                command("curl -o {} {}".format(out_gz,ftp_dir+out_gz)).run_comm(0)
                out_fname=out_gz.rstrip(".gz")
                command("gunzip -c {} > {}".format(out_gz,out_fname)).run_comm_no_exit(1)
        else:
            ftp_root_dir="https://cryptodb.org/common/downloads/Current_Release/ChominisUdeA01/"
            out_fname=genome.lower()+"."+file_type;
            cmd="curl -s {}/gff/data/ | grep buildNumber | awk \'{{split($0,a,\"\\\"\");print a[2]}}\'".format(ftp_root_dir)
            print ("cmd="+cmd)
            version=str(command(cmd).run_comm(1).decode("utf-8").rstrip())
            print ("version="+version)
            if file_type=="gff3":
                ftp_last_part="gff/data/CryptoDB-{}_ChominisUdeA01.gff".format(version)
            if file_type=="fasta":
                ftp_last_part="fasta/data/CryptoDB-{}_ChominisUdeA01_Genome.fasta".format(version)
            if file_type=="cds":             
                ftp_last_part="fasta/data/CryptoDB-{}_ChominisUdeA01_AnnotatedCDSs.fasta".format(version)
            command("curl {}/{} -o {}/{}".format(ftp_root_dir,ftp_last_part,download_dir,out_fname)).run_comm(0)    
        return download_dir+"/"+out_fname

    def clean_fasta(self,fasta_file):
        command("mv {} {}_old".format(fasta_file,fasta_file)).run_comm(0)
        command("sed -e \'/^[^>]/s/[^ATGCatgc]/N/g\' {}_old > {}".format(fasta_file,fasta_file)).run_comm(0)
        command("rm {}_old".format(fasta_file)).run_comm(0)
        return fasta_file

    def check_genome_avl(self, available_genomes_f, genome_name):
        genome_list_fh = open(available_genomes_f)
        if genome_name not in (line.rstrip() for line in genome_list_fh.readlines()):
            misc.my_exit("{} is not available, please try another genome".format(genome_name))
        genome_list_fh.close()

    def get_runID(self, infile):
        runID_pat = "([A-Z]RR\d+)"
        if not re.search(runID_pat, infile):
            misc.my_exit("{} not containing runID".format(infile))
        else:
            return re.findall(runID_pat,infile)[0]

    def get_samples_by_runIDs(self,mapfile):
        run_sample_dict = {}
        fhmap = open(mapfile,'r')
        for line in fhmap:
            line=line.rstrip()
            run = line.split()[0]
            sample=line.split()[1]
            run_sample_dict[run]=sample
        fhmap.close()
        return run_sample_dict

