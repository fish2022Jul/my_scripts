import re
import sys
marker=dict()
#marker.setdefault(0)

with open(sys.argv[1])as fh:
    for line in fh.readlines():
        if(re.search(r'^##|INDEL',line)):
            continue
        t = list(line.split())
        name=t[0]+"_"+t[1]
        if(marker.setdefault(name,0)==1):
            continue
        if(re.match(r'^#CH',line)):
            print("id,",end="")
        else:
            print(name,",1",sep="",end="")
        for i in range(9,len(t)):
            if(re.match(r'^0',t[i])):
                print(",AA",end="")
            elif(re.match(r'^1',t[i])):
                print(",BB",end="")
            elif(re.match(r'\.:\.',t[i])):
                print(",--",end="")
            elif(re.search(r'(.*)\.sor',t[i])):
                print(",",re.search(r'^(.*).sor',t[i]).group(1),sep="",end="")
        print()
        marker[name]=1
