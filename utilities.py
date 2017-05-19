from os import kill
from signal import alarm, signal, SIGALRM, SIGKILL
from subprocess import PIPE, Popen

'''
How to run it:
You need to first import it this Class into your code
Example: from utilities import Command
Then you can create an object of Command class and then use a function:
Example:
command = Command("python /usr/bin/spades.py -1 8605_7_1_test.fastq -2 8605_7_2_test.fastq -k 21,33,55,77,99,127 --careful --only-assembler -o spades_small_dir")
returncode,stdout,stderr=command.run(timeout=5)

'''

class Command(object):
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








def delete(dir):
    try:
        shutil.rmtree(dir)
    except shutil.Error as e:
        message='Directory not copied. Error: %s' %e
	error_list.append(message.replace("'",""))
        print(message)
