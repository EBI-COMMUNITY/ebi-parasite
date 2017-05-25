from utilities import command



class trim_galore:



	def __init__(self,fq1,fq2,prop, pair):
		self.fq1=fq1
		self.fq2=fq2
		self.prop=prop
		self.pair=pair
		error_list=list()
		self.error_list=error_list


	def run_trim_galore(self):
                print self.pair
                print "3:",self.prop.path_to_trim_galore
		if self.pair=='True':
			comm=self.prop.path_to_trim_galore+" --paired -q 20 "+ self.fq1+" "+self.fq2
		else:
			comm=self.prop.path_to_trim_galore+" -q 20 "+ self.fq1
		comm_obj=command(comm)
		returncode, stdout, stderr=comm_obj.run(3600)
		print returncode, stdout, stderr


	def post_process():
	    print "TODO:"


	def execute(self):
            self.run_trim_galore()

