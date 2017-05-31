from os import kill
import os
import sys
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen
import shutil

'''
How to run it:
You need to first import it this Class into your code
Example: from utilities import Command
Then you can create an object of Command class and then use a function:
Example:
command = Command("python /usr/bin/spades.py -1 8605_7_1_test.fastq -2 8605_7_2_test.fastq -k 21,33,55,77,99,127 --careful --only-assembler -o spades_small_dir")
returncode,stdout,stderr=command.run(timeout=5)

'''

class command(object):
	def __init__(self, command):
		self.command = command
		

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
		except Alarm:
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

	def get_process_children(self,pid):
		p = Popen('ps --no-headers -o pid --ppid %d' % pid, shell = True,
			  stdout = PIPE, stderr = PIPE)
		stdout, stderr = p.communicate()
		return [int(p) for p in stdout.split()]






class fileutils(object):


	def delete(self,dir,error_list):
		try:
			shutil.rmtree(dir)
		except shutil.Error as e:
			message='Directory not copied. Error: %s' %e
			error_list.append(message.replace("'",""))
			print(message)

	def create_processing_dir(self,directory):
		
		try: 
			if not os.path.exists(directory):
				os.makedirs(directory)
		except OSError as err:
			print("OS error: {0}".format(err))


	def move_file(self,file):

		try: 
			if not os.path.exists(file):
				os.makedirs(directory)
		except OSError as err:
			print("OS error: {0}".format(err))

	def copy_file(self,srcfile,destfile):

		try: 
			#shutil.copy (srcfile, destfile)
			shutil.copytree(srcfile, destfile)
		except shutil.Error as e:
			message='Directory not copied. Error: %s' %e
			print(message)
		except OSError as e:
			message='Directory not copied. Error: %s' %e
			print(message)
			

	def copy_file_into_dest(self,src, dest):
			name=os.path.basename(src)
			dest_file=os.path.join(dest, name)
			try:
				shutil.copy(src, dest_file)
			except shutil.Error as e:
				message='Directory not copied. Error: %s' %e
				print(message)
			except OSError as e:
				message='Directory not copied. Error: %s' %e
				print(message)

	def copy_dir_into_dest(self,src, dest):
			name=os.path.basename(src)
			dest_file=os.path.join(dest, name)
			try:
				shutil.copytree(src, dest_file)
			except shutil.Error as e:
				message='Directory not copied. Error: %s' %e
				print(message)
			except OSError as e:
				message='Directory not copied. Error: %s' %e
				print(message)




class db(object):

	def get_connection(db_user,db_password,db_host,db_database):
		conn = MySQLdb.connect(user=db_user, passwd=db_password, host=db_host,db=db_database)
		return conn


class properties(object):
	
	def __init__(self, property_file):
		with open(property_file) as f:
			 lines = f.readlines()
			 #print lines
		
		workdir_provided=False
		trim_galore_provided=False
		#workdir_input_provided=False
		#archivedir_provided=False
		#dbuser_provided=False
		#dbpassword_provided=False
		#dbhost_provided=False
		#dbname_provided=False
		#max_core_job_provided=False
		#path_to_quasr_provided = False
		#path_to_usearch_provided = False
		#path_to_spades_provided = False
		#data_base_list_file_provided = False
		#illumina_core_provided = False
		#biopython_provided=False        
		#teilenq_provided=False
		#velvetoptimiser_provided=False  
		
			
		
		for l in lines:
			#print l
			pair=l.strip().split(":")
			if pair[0].lower()=='workdir':
				self.workdir=pair[1].strip("\n")
				workdir_provided=True
			elif pair[0].lower()=='trim_galore': 
				self.trim_galore=pair[1].strip("\n")
				trim_galore_provided=True
			# elif pair[0].lower()=='max_core_job':
			#     self.max_core_job=pair[1].strip("\n")
			#     max_core_job_provided=True
			# elif pair[0].lower()=='workdir_input':
			#     self.workdir_input=pair[1].strip("\n")
			#     workdir_input_provided=True
			# elif pair[0].lower()=='archivedir':
			#     self.archivedir=pair[1].strip("\n")
			#     archivedir_provided=True
			# elif pair[0].lower()=='dbuser':
			#     self.dbuser=pair[1].strip("\n")
			#     dbuser_provided=True
			# elif pair[0].lower()=='dbpassword':
			#     self.dbpassword=pair[1].strip("\n")
			#     dbpassword_provided=True
			# elif pair[0].lower()=='dbhost':
			#     self.dbhost=pair[1].strip("\n")
			#     dbhost_provided=True
			# elif pair[0].lower()=='dbname':
			#     self.dbname=pair[1].strip("\n")
			#     dbname_provided=True
			# elif pair[0].lower() == 'path_to_quasr':
			#     path_to_quasr = pair[1].strip("\n")
			#     path_to_quasr_provided = True
			# elif pair[0].lower() == 'path_to_usearch':
			#     path_to_usearch = pair[1].strip("\n")
			#     path_to_usearch_provided = True
			# elif pair[0].lower() == 'path_to_spades':
			#     path_to_spades = pair[1].strip("\n")
			#     path_to_spades_provided = True
			# elif pair[0].lower() == 'data_base_list_file':
			#     data_base_list_file = pair[1].strip("\n")
			#     data_base_list_file_provided = True
			# elif pair[0].lower() == 'illumina_core':
			#     illumina_core = pair[1].strip("\n")
			#     illumina_core_provided = True
			# elif pair[0].lower() == 'biopython':
			#     biopython = pair[1].strip("\n")
			#     biopython_provided = True
			# elif pair[0].lower() == 'teilenq':
			#     teilenq = pair[1].strip("\n")
			#     teilenq_provided = True
				
		
		if workdir_provided==False:
		   self.workdir=os.getcwd()
		if trim_galore_provided==False:
		   print "ERROR: path_to_trim_galore missing"
		   sys.exit(1)
		# if max_core_job_provided==False:
		#     self.max_core_job=10
		# if workdir_input_provided==False:
		#     self.workdir_input=''
		# if archivedir_provided==False:
		#    self.archivedir=''
		# if dbuser_provided==False:
		#    self.dbuser_provided=''
		# if dbpassword_provided==False:
		#    self.dbpassword=''
		# if dbhost_provided==False:
		#    self.dbhost=''
		# if dbname_provided==False:
		#    self.dbname=''
		# if path_to_quasr_provided == False:
		#     print "ERROR: path_to_quasr missing"
		#     sys.exit(1)
		# if path_to_usearch_provided == False:
		#    print "ERROR: path_to_usearch missing"
		#    sys.exit(1)
		# if path_to_spades_provided == False:
		#    print "ERROR: path_to_spades missing"
		#    sys.exit(1)
		# if data_base_list_file_provided == False:
		#    print "ERROR: data_base_list_file missing"
		#    sys.exit(1)
		# if illumina_core_provided == False:
		#    print "ERROR: llumina_core missing"
		#    sys.exit(1)
		# if biopython_provided == False:
		#    print "ERROR: biopython missing"
		#    sys.exit(1)
		# if teilenq_provided == False:
		#    teilenq=''

		   


