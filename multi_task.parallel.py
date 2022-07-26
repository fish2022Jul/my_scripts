import multiprocessing
import time
import subprocess
import re
import sys

def do_shell(line):
    subprocess.call(line,shell=True)

def do_multi_step(sample_name):
    subprocess.call("bowtie2 --no-unal --score-min L,16,1 --local -L 16  -p 2 -x  /newlustre/home/fubeide/2b_rad_2016/grass_carp_lww0403/00.data/all.74406.BcgI -f -U test.{}.all.fa -S {}.sam".format(sample_name,sample_name),shell=True)
#    subprocess.call("SAMFilter.pl -m 25 -i {}.sam -o filtered.{}.sam".format(sample_name,sample_name),shell=True)    
    subprocess.call("SAMBaseCounts.BcgI.pl -i {}.sam -r /newlustre/home/fubeide/2b_rad_2016/grass_carp_lww0403/00.data/all.74406.BcgI.fa -c 5 -o {}.tab".format(sample_name,sample_name),shell=True)
p=multiprocessing.Pool(6)
res_l=[]
print("The time is {0}".format(time.ctime()))
with open(sys.argv[1]) as fp:
    for line in fp:
        m=re.search('([CXkg]\d+)\.',line)
#        res=p.apply_async(do_shell,args=(line,))

        res=p.apply_async(do_multi_step,args=(m.group(1),))
        res_l.append(res)

p.close()
p.join()

for res in res_l:
    print(res.get())

print("The time is {0}".format(time.ctime()))
