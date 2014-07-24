#!/usr/bin/env python

import os, re, subprocess, sys

cabnum = 3
proc = subprocess.Popen(['qstat', '@d0cabsrv%d' % cabnum],stdout=subprocess.PIPE)
thelist = []
for line in proc.stdout:
    mine = re.search('rclsa', line)
    if mine is None:
        continue
    num = int(line.rstrip().split('.')[0])
    if num>2394900:
        proc2 = subprocess.Popen(['qstat', '-f', '%d.d0cabsrv%d.fnal.gov' % (num, cabnum)],stdout=subprocess.PIPE)    
        for line2 in proc2.stdout:
            job = re.search('job=(\d+),', line2)
            if job is None:
                continue
            jobnum = int(job.group(1))
            thelist.append(jobnum)

thelist.sort()
lastjob = thelist[0]
missing = []
for num in thelist:
    if num!=lastjob:
        missing.append(lastjob)
    lastjob = num + 1

print missing
print lastjob
