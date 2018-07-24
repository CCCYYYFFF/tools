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
	name = (os.path.basename(i)).split("_")[1]
	f = name.split("-")
	if(f[0] == "MTHFR"):
		sample = f[1]+"_"+f[2]
		genotyping = "M"+f[3][0:2]
		order = f[3][2]
	else:
		sample = f[0]+"_"+f[1]
		genotyping = f[2]
		order = f[3]

	if(order == "R"):
		key = "MTHFR_"+genotyping+"_f"
	elif(order == "F"):
		key = "MTHFR_"+genotyping+"_r"
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

