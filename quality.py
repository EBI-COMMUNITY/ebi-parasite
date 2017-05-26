from utilities import command
from utilities import file


class trim_galore:



	def __init__(self,fq1,fq2,prop, pair):
		self.fq1=fq1
		self.fq2=fq2
		self.prop=prop
		self.pair=pair
		error_list=list()
		self.error_list=error_list


	def run_trim_galore(self):
		if self.pair==True:
			comm=self.prop.trim_galore+" --paired -q 20 "+ self.fq1+" "+self.fq2
		else:
			comm=self.prop.trim_galore+" -q 20 "+ self.fq1
		print comm
		comm_obj=command(comm)
		returncode, stdout, stderr=comm_obj.run(3600)
		print returncode, stdout, stderr


	def post_process():
		fi=file()
        indir=prop.workdir+"/quality/in/"
        outdir=prop.workdir+"/quality/out/"
        fi.create_processing_dir(indir)
        fi.create_processing_dir(outdir)

	    


	def execute(self):
            self.run_trim_galore()

