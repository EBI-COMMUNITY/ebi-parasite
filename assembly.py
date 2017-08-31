import argparse



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



def initiate():
    #Here please provide the initiation code
     print "Here please do the initiation code"


def execute():
    #Here please provide the execucation code
    print "Here please do the execucation  code"
    
    #Execution of the external softwares, you need to use command class from utilities
    #Example: from utilities import command
    #Then you can create an object of command class and then use a function:
    #Example:
    #command = Command("python /usr/bin/spades.py -1 8605_7_1_test.fastq -2 8605_7_2_test.fastq -k 21,33,55,77,99,127 --careful --only-assembler -o sp    ades_small_dir")
    #returncode,stdout,stderr=command.run(timeout=5)

def post_process():
    #Here please provide the post process code
    print "Here please do the post process  code"



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
   
     
