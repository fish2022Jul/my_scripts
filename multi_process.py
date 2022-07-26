import multiprocessing 
import time
import subprocess

import sys

def do_shell(line):
    subprocess.call(line,shell=True)
p=multiprocessing.Pool(int(sys.argv[2]))
res_l=[]
print("The time is {0}".format(time.ctime()))
with open(sys.argv[1]) as fp:
    for line in fp:
        res=p.apply_async(do_shell,args=(line,))
        res_l.append(res)
        time.sleep(1)

p.close()
p.join()

for res in res_l:
    print(res.get())

print("The time is {0}".format(time.ctime()))
