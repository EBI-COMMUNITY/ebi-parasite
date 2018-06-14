from os import kill
import os
import sys
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen
import shutil
import re


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
        print self.command

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
            print "\nERROR: {} FAILED!!! \n\nSTDERR: {}\nSTDOUT: {}\n".format(self.command,stderr,stdout)
            sys.exit(1)    
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
	    print("OS error: {}".format(err))
            raise

    def create_processing_dir(self,directory):
        try: 
    	    if not os.path.exists(directory):
	        os.makedirs(directory)
	except OSError as err:
	    print("OS error: {}".format(err))
            raise


    def copy_file(self,srcfile,destfile):
	try: 
		#shutil.copy (srcfile, destfile)
  	    shutil.copy(srcfile, destfile)
	except OSError as err:
	    print('OS error: {}'.format(err))
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
            print ('OS error: {}'.format(err)) 
            raise

    
    def copy_file_add_prefix(self,source_fpath,outdir,prefix):
	self.copy_file_to_destdir(source_fpath,outdir)
	source_file2=os.path.join(outdir,os.path.basename(source_fpath))		
	self.add_file_prefix(source_file2,prefix)
		
    
    def check_exist(self,path):
        if not os.path.exists(path):
            print("ERROR: {} does not exist!".format(path))
            sys.exit(1)
        else:
            return 1


    def check_files_exist(self,fpaths):
        for fpath in fpaths:
            self.check_exist(fpath)
        return 1

    

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
           print "ERROR: {} path was not defined in {}. Please define it and re-run.".format(attrib,self.propertyf)
           sys.exit(1)     
        else:
            fileutils().check_exist(getattr(self,attrib))
            return getattr(self,attrib)
