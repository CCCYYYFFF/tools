#!/usr/bin/python
#-*-coding=utf-8-*-

import sys
import os
import glob

indir,outdir,resultdir = sys.argv[1:]
indir = os.path.abspath(indir)
absname = os.path.dirname(indir)
outdir = os.path.abspath(outdir)
resultdir = os.path.abspath(resultdir)
ab1files = glob.glob('%s/*.ab1' %indir)
all = {}
for i in ab1files:
	name = (os.path.basename(i)).split(".")
	sample = name[0]
	key = "BRAF_01_r"
	if sample in all:
		all[sample].update({key:i})
	else:
		all.update({sample:{key:i}})

if(not os.path.exists(outdir)):
	os.makedirs(outdir)

for i in all:
	lis_name =outdir+"/"+str(i)+".lis"
	F = open(lis_name,'w')
	for j in all[i]:
		F.write(j+"="+all[i][j]+"\n")
	F.write("output="+resultdir+"/"+i+"/\n")
	F.close()

