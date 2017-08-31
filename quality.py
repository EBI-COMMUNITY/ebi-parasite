from utilities import command
from utilities import fileutils


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
		print self.prop.workdir
		self.post_process()


	def post_process(self):
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
		


	def execute(self):
			self.run_trim_galore()
			

