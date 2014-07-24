import os, time

f = open('failed.txt')
for l in f:
    os.environ['job']=l.strip()
    print os.environ['job']
    os.system('qsub -koe -me -M rclsa@fnal.gov -l nodes=1 -v job -q medium@d0cabsrv3 cabpythia.sh')
    time.sleep(1)
f.close()
