#!/usr/bin/python
# Copyright (C) 2018-04-02 zhaoerying All Rights Reserved

import sys
import os
import subprocess
import time
import re
import ConfigParser

from A_Load import load_input 
from A_Load import ab1_phred_v2
from A_Load import ab1_trace_dump
from B_Align import normal_sw_align
from D_Anno import anno_genome
from D_Anno import anno_cosmic
from F_Report import draw_pic
from optparse import OptionParser

def printtime(message, *args):
	if args:
		message = message % args
	print "[ " + time.strftime('%X') + " ] " + message
	sys.stdout.flush()
	sys.stderr.flush()

def RunCommand(command,description):
	printtime(' ')
	printtime('Task    : ' + description)
	printtime('Command : ' + command)
	printtime(' ')
	stat = subprocess.call(command,shell=True)
	if stat != 0:
		printtime('ERROR: command failed with status %d' % stat)
		sys.exit(1)

class ref_info:
	def __init__(self, cf, real_key):

		self.path = "%s/O_Database/%s.fa" % (os.path.split(os.path.realpath(__file__))[0], cf.get("Key2db", real_key))

		self.chr = 'na'
		self.start = 0
		self.end = 0
		self.region = 'na'
		self.seq = ''

	def getinfo(self):

		input = open(self.path, 'r')
		for line in input:
			line = line.strip()
			if not line or line[0] == '#':
				continue
			if line[0] == '>':
				name = re.sub(r'>','',line)
				arr_name = re.split(r'_',name)
				self.chr = arr_name[0]
				self.start = arr_name[1]
				self.end = arr_name[2]
				self.region = arr_name[3]
			else:
				self.seq += line

def get_r_q(new_query, new_ref):
	dict_r_q = {}
	pos_r_q = {}

	len_r = len(new_ref)
	r_pos = -1
	q_pos = -1
	for i in range(0, len_r):
		if new_query[i] != '-' :
			q_pos += 1
		if new_ref[i] != '-' :
			r_pos += 1
			dict_r_q[r_pos] = new_query[i]
			pos_r_q[r_pos] = q_pos

	return (dict_r_q,pos_r_q)

def get_pos(obj_ref,obj_dump,ref_pos,strand):
	q_pos_dict = {}
	obj_ref_dump = normal_sw_align.Sanger_align(obj_ref.seq, obj_dump.seq)
	obj_ref_dump.score()
	obj_ref_dump.mergeFa()
	obj_ref_dump.get_pos_relation()

	q_pos = int(ref_pos)-int(obj_ref.start)
	r_pos = obj_ref_dump.query_ref_pos[q_pos]
	if strand == 'r':
		r_pos = obj_dump.seq_len - r_pos - 1
	q_pos_dict[r_pos] = "%s,%s,%s" % (obj_ref.chr,ref_pos,obj_ref.seq[q_pos]+","+obj_ref.seq[q_pos])
	return (q_pos_dict)

def correctRun(raw_poly1, raw_poly2, ref_seq, out_log):

	obj_cor_poly1_ref = normal_sw_align.Sanger_align(raw_poly1, ref_seq)
	obj_cor_poly1_ref.score() 
	obj_cor_poly1_ref.mergeFa()

	obj_cor_poly2_ref = normal_sw_align.Sanger_align(raw_poly2, ref_seq)
	obj_cor_poly2_ref.score() 
	obj_cor_poly2_ref.mergeFa() 

	len_q1 = len(raw_poly1)
	len_q2 = len(raw_poly2)
	min_len_q = len_q1
	if min_len_q > len_q2:
		min_len_q = len_q2

	out_log.write("raw_poly1: %s\n" % raw_poly1)
	out_log.write("raw_poly2: %s\n" % raw_poly2)
	out_log.write("loop q1: %s\n" % obj_cor_poly1_ref.new_query)
	out_log.write("loop r1: %s\n" % obj_cor_poly1_ref.new_ref)
	out_log.write("loop q2: %s\n" % obj_cor_poly2_ref.new_query)
	out_log.write("loop r2: %s\n" % obj_cor_poly2_ref.new_ref)
	'''
	print "raw_poly1: %s" % raw_poly1
	print "raw_poly2: %s" % raw_poly2
	print "loop q1  : %s" % obj_cor_poly1_ref.new_query
	print "loop r1  : %s" % obj_cor_poly1_ref.new_ref
	print "loop q2  : %s" % obj_cor_poly2_ref.new_query
	print "loop r2  : %s" % obj_cor_poly2_ref.new_ref
	'''
	len_ref = len(ref_seq)

	len_cor_poly1 = len(obj_cor_poly1_ref.new_ref)

	dict_r_q2 = {}
	pos_r_q = {}
	dict_q2_r = {}
	pos_q_r = {}

	(dict_r_q2,pos_r_q) = get_r_q( obj_cor_poly2_ref.new_query, obj_cor_poly2_ref.new_ref)
	(dict_q2_r,pos_q_r) = get_r_q( obj_cor_poly2_ref.new_ref, obj_cor_poly2_ref.new_query)

	re_flag = 0

	q_pos = -1
	r_pos = -1

	new_poly1 = ""
	new_poly2 = ""

	# obj_cor_poly1_ref.new_query, obj_cor_poly1_ref.new_ref
	for i in range(0,len_cor_poly1):

		q1_base = obj_cor_poly1_ref.new_query[i]
		r1_base = obj_cor_poly1_ref.new_ref[i]

		if q1_base != "-":
			q_pos += 1
		if r1_base != '-':
			r_pos += 1

		'''
		if q_pos < tmp_q_pos :
			new_poly1 += raw_poly1[q_pos]
			new_poly2 += raw_poly2[q_pos]
			continue
		'''

		if q_pos == -1 or r_pos == -1:
			if q_pos != -1:
				new_poly1 += raw_poly1[q_pos]
				new_poly2 += raw_poly2[q_pos]				
			continue

		if r_pos >= len_ref:
			print "cor ref all done "
			return (new_poly1, new_poly2)
		if q_pos >= min_len_q:
			print "cor query all done "
			return (new_poly1, new_poly2)
		if (dict_r_q2[r_pos] == '-' and r1_base == '-') or (q1_base == '-' and dict_q2_r[q_pos] == '-'):# 
			print "cor stop "
			return (new_poly1, new_poly2)

		if q1_base != '-' and r1_base != '-':
			#print "%s %s %s" % (r1_base, q1_base, raw_poly2[q_pos])
			if raw_poly2[q_pos] == r1_base:

				new_poly1 += raw_poly2[q_pos]
				new_poly2 += q1_base
			else:
				new_poly1 += q1_base
				new_poly2 += raw_poly2[q_pos]
		elif q1_base == '-' and dict_r_q2[r_pos] == '-':
			continue
		elif r1_base == '-'	and dict_q2_r[q_pos] == '-':
			new_poly1 += q1_base
			new_poly2 += raw_poly2[q_pos]			
			continue
		elif q1_base != '-' and r1_base == '-': # poly1 insert
			r_pos += 1
			#print "poly1 insert start"
			if r_pos >= len_ref:
				return (new_poly1,new_poly2)
			q2_pos = pos_r_q[r_pos]

			tmp_n1 = new_poly1
			tmp_n2 = new_poly2
			
			if q_pos<q2_pos:
				new_poly1 = tmp_n1 + raw_poly2[q2_pos:]
				new_poly2 = tmp_n2 + raw_poly2[q_pos:q2_pos] + raw_poly1[q_pos:]
			elif q_pos > q2_pos:
				new_poly1 = tmp_n1 + tmp_n2[q2_pos:q_pos] + raw_poly2[q_pos:]
				new_poly2 = tmp_n2[:q2_pos] + raw_poly1[q_pos:]
			else:
				new_poly1 = tmp_n1 + raw_poly2[q_pos:]
				new_poly2 = tmp_n2 + raw_poly1[q_pos:]

			re_flag = 1

		elif q1_base == '-' and r1_base != '-': # poly1 del	
			q_pos += 1
			#print "poly1 del start"
			q2_pos = pos_r_q[r_pos]

			#print "new1: %s" % new_poly1
			#print "new2: %s" % new_poly2

			tmp_n1 = new_poly1
			tmp_n2 = new_poly2
			
			if q_pos<q2_pos:
				new_poly1 = tmp_n1 + raw_poly2[q2_pos:]
				new_poly2 = tmp_n2 + raw_poly2[q_pos:q2_pos] + raw_poly1[q_pos:]
			elif q_pos > q2_pos:
				new_poly1 = tmp_n1 + tmp_n2[q2_pos:q_pos] + raw_poly2[q_pos:]
				new_poly2 = tmp_n2[:q2_pos] + raw_poly1[q_pos:]
			else:
				new_poly1 = tmp_n1 + raw_poly2[q_pos:]
				new_poly2 = tmp_n2 + raw_poly1[q_pos:]

			re_flag = 1
			
		else:
			print "cor #### ?"


		if re_flag == 1:

			(new_poly1,new_poly2) = correctRun(new_poly1, new_poly2, ref_seq, out_log)
			break

	return (new_poly1,new_poly2)

def Master(inputfile, confini, javapath):
	print "\n\n\n"
	#01 load conf.ini
	cf = ConfigParser.ConfigParser()
	cf.read(confini)

	#02 load input_file
	obj_ab1s = load_input.Input(inputfile)
	obj_ab1s.getinfo()
	
	outdir_report_root = obj_ab1s.outdir + '/Report'
	if not os.path.exists(outdir_report_root):
		os.makedirs(outdir_report_root)	

	report_QC  = open(outdir_report_root+'/report_QC.txt','w')
	report_Mut = open(outdir_report_root+'/report_Mut.txt','w')

	cmd_cp = "cp %s %s " % (inputfile, outdir_report_root)
	RunCommand(cmd_cp,"cp inputfile")
	cmd_cp = "cp %s %s " % (confini, outdir_report_root)
	RunCommand(cmd_cp,"cp confini")

	tmp_file = outdir_report_root+"/raw_anno.vcf"
	raw_vcf = open(tmp_file, 'w')

	Mut_keys = {}
	Mut_info = {}

	#03 for ab1 to analysis
	for tmp_key in sorted(obj_ab1s.ab1s):
		print "\n\n\n"
		print tmp_key
		print obj_ab1s.ab1s[tmp_key]

		#03.1 creat dirs
		outdir_analysis = obj_ab1s.outdir + '/Analysis/' + tmp_key
		outdir_report = obj_ab1s.outdir + '/Report/' + tmp_key
		if not os.path.exists(outdir_analysis):
			os.makedirs(outdir_analysis)
			os.makedirs(outdir_report)

		out_log = open(outdir_analysis + '/analysis.log','w')
		out_mask = open(outdir_report + '/'+tmp_key+'.mask','w')
		out_dump_seq = open(outdir_report + '/'+tmp_key+'.txt','w')

		#cp
		cmd_cp = "cp %s %s " % (obj_ab1s.ab1s[tmp_key], obj_ab1s.outdir+ '/Report')
		RunCommand(cmd_cp,"cp ab1")

		out_log.write("\n"+ tmp_key + "\n" + "ab1 : "+obj_ab1s.ab1s[tmp_key])

		#03.2 get ref.fa
		arr_keys = re.split(r'_',tmp_key) 
		real_key = "%s_%s" % (arr_keys[0],arr_keys[1])
		real_strand = arr_keys[2]

		obj_ref = ref_info(cf, real_key)
		obj_ref.getinfo()
		
		#03.3 load phred 
		# cf.get("Key2db", real_key)
		#(self, ab1, outdir, key, strand, phred_E=1e-2, phred_num=20, l_cutoff=30, r_cutoff=400, indel_ratio=8.5, qual_nearby_num=10)
		obj_pred = ab1_phred_v2.phred(obj_ab1s.ab1s[tmp_key], outdir_analysis, tmp_key, real_strand, cf.get("Phred", "phred_E"), cf.get("Phred", "phred_num"), cf.get("Phred", "%s-cutoff_s" % tmp_key), cf.get("Phred", "%s-cutoff_e" % tmp_key), cf.get("Phred", "indel_ratio"), cf.get("Phred", "qual_nearby_num"), float(cf.get("Phred", "auto_cutoff_h3_ratio")))
		obj_pred.run_phred()
		obj_pred.format_poly_file()
		obj_pred.load_phred_file()
		obj_pred.filter_h3()
		if int(cf.get("Filter", "ruijin_model")) == 1:
			obj_pred.ruijin_cor()
		elif int(cf.get("Filter", "ruijin_model")) == 2:
			obj_pred.ruijin_cor2()
		obj_pred.cor_strand()



		#03.4 load trace_dump
		obj_dump = ab1_trace_dump.dump(obj_ab1s.ab1s[tmp_key], outdir_analysis, tmp_key)
		obj_dump.trace_dump()
		obj_dump.load_dump_file()

		out_dump_seq.write(obj_dump.out_seq)
		
		out_log.write("\n\nobj_pred\nN_num: " +str(obj_pred.N_num)+"\nN_ratio: " +str(obj_pred.N_ratio)+"\npoly1: "+obj_pred.normal_poly1_seq+"\npoly2: "+obj_pred.normal_poly2_seq)
		out_log.write("\n\nobj_dump\nN_num: " +str(obj_dump.N_num)+"\nN_ratio: " +str(obj_dump.N_ratio)+"\nseq: "+obj_dump.seq)

		qc_1 = "phred-N_ratio:%s(<%s)" % (int(obj_pred.N_ratio), int(cf.get("Filter", "obj_phred-N_ratio")))
		qc_2 = "dump-N_ratio:%s(<%s)" % (int(obj_dump.N_ratio), int(cf.get("Filter", "obj_dump-N_ratio")))
		qc_3 = "dump_seq_len:%s(<%s)" % (len(obj_dump.seq), int(cf.get("Filter", "obj_dump_seq_len")))
		out_log.write("%s\n%s\n%s\n" % (qc_1, qc_2, qc_3))

		if int(obj_pred.N_ratio) >= int(cf.get("Filter", "obj_phred-N_ratio")) or int(obj_dump.N_ratio) > int(cf.get("Filter", "obj_dump-N_ratio"))  or len(obj_dump.seq)<int(cf.get("Filter", "obj_dump_seq_len")):
			qc_no = "?"
			if int(obj_pred.N_ratio) >= int(cf.get("Filter", "obj_phred-N_ratio")) :
				qc_no = "phred-N_ratio:%s(<%s)" % (int(obj_pred.N_ratio), int(cf.get("Filter", "obj_phred-N_ratio")))
			elif int(obj_dump.N_ratio) > int(cf.get("Filter", "obj_dump-N_ratio")):
				qc_no = "dump-N_ratio:%s(<%s)" % (int(obj_dump.N_ratio), int(cf.get("Filter", "obj_dump-N_ratio")))
			elif len(obj_dump.seq)<int(cf.get("Filter", "obj_dump_seq_len")) :
				qc_no = "dump_seq_len:%s(<%s)" % (len(obj_dump.seq), int(cf.get("Filter", "obj_dump_seq_len")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "null_report"
			continue

		qc_1 = "poly1_poly2_s_ratio:%s(>=%s)" % (int(obj_pred.p1_p2_s_ratio), int(cf.get("Filter", "poly1_poly2_s_ratio")))
		out_log.write("%s\n" % (qc_1))

		# phred filter
		if int(obj_pred.p1_p2_s_ratio) < int(cf.get("Filter", "poly1_poly2_s_ratio")) :
			qc_no = "poly1_poly2_s_ratio:%s(>=%s)" % (int(obj_pred.p1_p2_s_ratio), int(cf.get("Filter", "poly1_poly2_s_ratio")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )			
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue			

		#03.5 get dump poly position
		obj_poly1_dump = normal_sw_align.Sanger_align(obj_pred.normal_poly1_seq, obj_dump.seq)
		obj_poly1_dump.score() 
		obj_poly1_dump.mergeFa() 

		qc_1 = "poly1_dump_align: no_align"
		out_log.write("%s\n" % (qc_1))

		if obj_poly1_dump.no_align == 1:
			qc_no = "poly1_dump_align: no_align"

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "null_report"
			continue			
		obj_poly1_dump.get_pos_relation()

		out_log.write("\n\nobj_poly1_dump\nalign_qe - align_qs : "+ str(obj_poly1_dump.align_qe) +" - "+ str(obj_poly1_dump.align_qs))
		out_log.write("\nq: "+obj_poly1_dump.new_query+"\nr: "+obj_poly1_dump.new_ref+"\n")

		qc_1 =  "poly1_dump_align_len:%s(>=%s)" % (obj_poly1_dump.align_qe - obj_poly1_dump.align_qs, int(cf.get("Filter", "poly1_dump_align_len")))
		out_log.write("%s\n" % (qc_1))

		#03.6 QC dump_poly
		if obj_poly1_dump.align_qe - obj_poly1_dump.align_qs < int(cf.get("Filter", "poly1_dump_align_len")) :
			qc_no = "poly1_dump_align_len:%s(>=%s)" % (obj_poly1_dump.align_qe - obj_poly1_dump.align_qs, int(cf.get("Filter", "poly1_dump_align_len")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )

			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue
		
		# for test snp indel model
		'''
		obj_ref.seq = 				'AAACCCATGTATGAAGTACAGTGGAAGGTTGTTGAGGAGATAAATGGAAACAATTATGTTTACATAGACCCAACACAACT'

		obj_pred.normal_poly1_seq = 'AAACCCATGTATGAAGTACAATCGGTGGAAGGTTGTTGAGGAGATAAATGGAAACAATTATGTTTACATAGACCCAACACAACT'
		obj_pred.normal_poly2_seq = 'AAACCCATGTATGAAGTACAGTGGAAGGTTGTTGAGGAGATAAATGGAAACATCGAATTATGTTTACATAGACCCAACACAACT'		
		'''

		#03.7 pre correct Poly1
		obj_cor_poly1_ref = normal_sw_align.Sanger_align(obj_pred.normal_poly1_seq, obj_ref.seq)
		obj_cor_poly1_ref.score() 
		obj_cor_poly1_ref.mergeFa() 

		qc_1 = "poly1_ref_align_len:%s(>=%s)" % (obj_cor_poly1_ref.align_qe - obj_cor_poly1_ref.align_qs, int(cf.get("Filter", "poly1_ref_align_len")))
		out_log.write("%s\n" % (qc_1))

		if obj_cor_poly1_ref.align_qe - obj_cor_poly1_ref.align_qs < int(cf.get("Filter", "poly1_ref_align_len")) :
			qc_no = "poly1_ref_align_len:%s(>=%s)" % (obj_cor_poly1_ref.align_qe - obj_cor_poly1_ref.align_qs, int(cf.get("Filter", "poly1_ref_align_len")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )

			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue

		#03.9 pre correct Poly2
		obj_cor_poly2_ref = normal_sw_align.Sanger_align(obj_pred.normal_poly2_seq, obj_ref.seq)
		obj_cor_poly2_ref.score() 
		obj_cor_poly2_ref.mergeFa() 

		qc_1 = "poly2_ref_align_len:%s(>=%s)" % (obj_cor_poly2_ref.align_qe - obj_cor_poly2_ref.align_qs, int(cf.get("Filter", "poly2_ref_align_len")))
		out_log.write("%s\n" % (qc_1))

		if obj_cor_poly2_ref.align_qe - obj_cor_poly2_ref.align_qs < int(cf.get("Filter", "poly2_ref_align_len")) :
			qc_no = "poly2_ref_align_len:%s(>=%s)" % (obj_cor_poly2_ref.align_qe - obj_cor_poly2_ref.align_qs, int(cf.get("Filter", "poly2_ref_align_len")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )			
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue

		(cor_poly1, cor_poly2) = correctRun(obj_pred.normal_poly1_seq, obj_pred.normal_poly2_seq, obj_ref.seq, out_log)

		# get poly1 ref muts
		obj_poly1_ref = normal_sw_align.Sanger_align(cor_poly1, obj_ref.seq)
		obj_poly1_ref.score() 
		obj_poly1_ref.mergeFa() 
		obj_poly1_ref.get_SnpIndel()

		qc_1 = "poly1_ref_align_len:%s(>=%s)" % (obj_poly1_ref.align_qe - obj_poly1_ref.align_qs,int(cf.get("Filter", "poly1_ref_align_len")) )
		out_log.write("%s\n" % (qc_1))

		#03.8 QC cor_poly1_ref
		if obj_poly1_ref.align_qe - obj_poly1_ref.align_qs < int(cf.get("Filter", "poly1_ref_align_len")) :

			qc_no = "?"
			if obj_poly1_ref.align_qe - obj_poly1_ref.align_qs < int(cf.get("Filter", "poly1_ref_align_len")) :
				qc_no = "poly1_ref_align_len:%s(>=%s)" % (obj_poly1_ref.align_qe - obj_poly1_ref.align_qs,int(cf.get("Filter", "poly1_ref_align_len")) )

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )		
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no		
			print "draw_pic"
			print "null_report"
			continue

		# get poly2 ref muts
		obj_poly2_ref = normal_sw_align.Sanger_align(cor_poly2, obj_ref.seq)
		obj_poly2_ref.score() 
		obj_poly2_ref.mergeFa() 
		obj_poly2_ref.get_SnpIndel()

		qc_1 = "poly2_ref_align_len:%s(>=%s)" % (obj_poly2_ref.align_qe - obj_poly2_ref.align_qs,int(cf.get("Filter", "poly2_ref_align_len")) )
		out_log.write("%s\n" % (qc_1))

		#03.10 QC cor_poly2_ref
		if obj_poly2_ref.align_qe - obj_poly2_ref.align_qs < int(cf.get("Filter", "poly2_ref_align_len")) :
			qc_no = "?"
			if obj_poly2_ref.align_qe - obj_poly2_ref.align_qs < int(cf.get("Filter", "poly2_ref_align_len")) :
				qc_no = "poly2_ref_align_len:%s(>=%s)" % (obj_poly2_ref.align_qe - obj_poly2_ref.align_qs,int(cf.get("Filter", "poly2_ref_align_len")) )
			#elif int(obj_poly2_ref.s_ratio) < int(cf.get("Filter", "poly2_ref_s_ratio")) :
				#qc_no = "poly2_ref_s_ratio:%s(>=%s)" % (int(obj_poly2_ref.s_ratio) , int(cf.get("Filter", "poly2_ref_s_ratio")))
			#elif int(obj_poly2_ref.nosame_num) > int(cf.get("Filter", "poly2_ref_nosame_num")):
				#qc_no = "poly2_ref_nosame_num:%s(<%s)" % (int(obj_poly2_ref.nosame_num), int(cf.get("Filter", "poly2_ref_nosame_num")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )

			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue

		print obj_poly1_ref.SnpIndel
		print obj_poly2_ref.SnpIndel
		if obj_poly1_ref.SnpIndel or obj_poly2_ref.SnpIndel:
			(obj_poly1_ref.SnpIndel, obj_poly2_ref.SnpIndel, hete_indel_pos) = filter_raw_mut(obj_poly1_ref.SnpIndel, obj_poly2_ref.SnpIndel, real_strand)
			if real_strand == 'f':
				if hete_indel_pos ==0:
					obj_poly1_ref.same_ratio(obj_poly1_ref.align_qs, obj_poly1_ref.align_qe)
					obj_poly2_ref.same_ratio(obj_poly2_ref.align_qs, obj_poly2_ref.align_qe)
				else:

					obj_poly1_ref.same_ratio(obj_poly1_ref.align_qs, hete_indel_pos)
					obj_poly2_ref.same_ratio(obj_poly2_ref.align_qs, hete_indel_pos)

			else:
				if hete_indel_pos ==0:
					obj_poly1_ref.same_ratio(obj_poly1_ref.align_qs, obj_poly1_ref.align_qe)
					obj_poly2_ref.same_ratio(obj_poly2_ref.align_qs, obj_poly2_ref.align_qe)
				else:

					obj_poly1_ref.same_ratio(hete_indel_pos, obj_poly1_ref.align_qe)
					obj_poly2_ref.same_ratio(hete_indel_pos, obj_poly2_ref.align_qe)
		else:
			obj_poly1_ref.same_ratio(obj_poly1_ref.align_qs, obj_poly1_ref.align_qe)
			obj_poly2_ref.same_ratio(obj_poly2_ref.align_qs, obj_poly2_ref.align_qe)
		print obj_poly1_ref.SnpIndel
		print obj_poly2_ref.SnpIndel

		#print "obj_poly2_ref.s_ratio : %s" % obj_poly2_ref.s_ratio

		qc_1 = "poly1_ref_s_ratio:%s(>=%s)" % (int(obj_poly1_ref.s_ratio) , int(cf.get("Filter", "poly1_ref_s_ratio")))
		qc_2 = "poly2_ref_s_ratio:%s(>=%s)" % (int(obj_poly2_ref.s_ratio) , int(cf.get("Filter", "poly2_ref_s_ratio")))
		out_log.write("%s\n%s\n" % (qc_1, qc_2))

		if int(obj_poly2_ref.s_ratio) < int(cf.get("Filter", "poly2_ref_s_ratio")) or int(obj_poly1_ref.s_ratio) < int(cf.get("Filter", "poly1_ref_s_ratio")) :
			qc_no = "?"
			if int(obj_poly1_ref.s_ratio) < int(cf.get("Filter", "poly1_ref_s_ratio")) :
				qc_no = "poly1_ref_s_ratio:%s(>=%s)" % (int(obj_poly1_ref.s_ratio) , int(cf.get("Filter", "poly1_ref_s_ratio")))

			elif int(obj_poly2_ref.s_ratio) < int(cf.get("Filter", "poly2_ref_s_ratio")) :
				qc_no = "poly2_ref_s_ratio:%s(>=%s)" % (int(obj_poly2_ref.s_ratio) , int(cf.get("Filter", "poly2_ref_s_ratio")))

			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )

			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			out_log.write("\n\nMut_result\nNO muts\n")
			out_log.write("low quality : %s \n" % qc_no)
			out_mask.write("mask(start:%s,end:%s)" % (0, 0))

			print "low quality : %s " % qc_no
			print "draw_pic"
			print "null_report"
			continue

		muts_info = {}
		muts_sum = {}

		# get mask info 
		mask_s = obj_poly1_ref.align_qs
		if mask_s < obj_poly2_ref.align_qs:
			mask_s = obj_poly2_ref.align_qs

		mask_e = obj_poly1_ref.align_qe
		if mask_e > obj_poly2_ref.align_qe:
			mask_e = obj_poly2_ref.align_qe
		mask_e = mask_e -1
		print "\nraw_mask: %s -- %s \n" % (mask_s, mask_e)

		if mask_s in obj_poly1_dump.query_ref_pos:
			mask_s = obj_poly1_dump.query_ref_pos[mask_s]
			if real_strand == 'r':
				mask_s = obj_dump.seq_len - mask_s - 1
		else:
			mask_s = 0

		if mask_e in obj_poly1_dump.query_ref_pos:
			mask_e = obj_poly1_dump.query_ref_pos[mask_e]
			if real_strand == 'r':
				mask_e = obj_dump.seq_len - mask_e - 1
		else:
			mask_e = 0

		print "\nref_mask: %s -- %s \n" % (mask_s, mask_e)

		if mask_s > mask_e:
			tmp = mask_s
			mask_s = mask_e
			mask_e = tmp
		if mask_s != 0:
			mask_s = mask_s -1
		else:
			mask_s = 30
		if mask_e != 0:
			mask_e = mask_e +1

		

		out_log.write("\n\npoly1_align\n1q: "+obj_poly1_ref.new_query+"\n    "+obj_poly1_ref.new_tag+"\n1r: "+obj_poly1_ref.new_ref+"\nmuts\n")
		out_log.write("\n\npoly2_align\n2q: "+obj_poly2_ref.new_query+"\n    "+obj_poly2_ref.new_tag+"\n2r: "+obj_poly2_ref.new_ref+"\nmuts\n")

		#print obj_poly2_ref.new_query
		#print obj_poly2_ref.new_ref


		if obj_poly2_ref.SnpIndel or obj_poly1_ref.SnpIndel:
			out_log.write(" poly1_mut : %s\n" % str(obj_poly1_ref.SnpIndel))
			out_log.write(" poly2_mut : %s\n" % str(obj_poly2_ref.SnpIndel))
			print obj_poly1_ref.SnpIndel
			print obj_poly2_ref.SnpIndel

			region_start = re.split(r',',cf.get("Region", "%s_s" % tmp_key ))
			region_end   = re.split(r',',cf.get("Region", "%s_e" % tmp_key ))

			print "region_start: %s, region_end: %s \n" % (region_start, region_end)

			out_log.write(" region : s:%s e:%s \n" % (region_start, region_end))
			max_ref_alt_len = 0
			alt_num = {'A':0, 'T':0, 'C':0, 'G':0}
			# get obj_poly1_ref SnpIndel
			for p1_num in obj_poly1_ref.SnpIndel:
				#mut_q_pos,mut_s,mut_e,mut_r,mut_q
				arr = re.split(r',',obj_poly1_ref.SnpIndel[p1_num])
				start = str(int(arr[1]) + int(obj_ref.start) - 1)

				if max_ref_alt_len < len(arr[2]):
					max_ref_alt_len = len(arr[2])
				if max_ref_alt_len < len(arr[3]):
					max_ref_alt_len = len(arr[3])

				if arr[3] == "G" :
					alt_num['G'] += 1
				elif arr[3] == "C":
					alt_num['C'] += 1
				elif arr[3] == "T":
					alt_num['T'] += 1
				elif arr[3] == "A":
					alt_num['A'] += 1

				#print "s: %s, e: %s, n: %s" % (region_start, region_end, start)
				flag_region = 0
				tmp_i = 0
				for region_s in region_start:
					region_e = region_end[tmp_i]
					if int(start) >= int(region_s) and int(start) <= int(region_e) :
						flag_region = 1
						break
					tmp_i += 1

				if flag_region == 0:
					out_log.write(" poly1_out_of_region : %s n: %s\n" % (obj_poly1_ref.SnpIndel[p1_num], start))
					print " poly1_out_of_region : %s n: %s\n" % (obj_poly1_ref.SnpIndel[p1_num], start)
					continue

				mut_key = "%s,%s,%s,%s" % (obj_ref.chr, start, arr[2], arr[3])

				if muts_sum.has_key(mut_key):
					muts_sum[mut_key] += 1
				else:
					muts_sum[mut_key] = 1
				muts_info[mut_key] = "%s,%s" % (arr[0], obj_ref.region)

			# get obj_poly2_ref SnpIndel
			for p2_num in obj_poly2_ref.SnpIndel:
				#mut_q_pos,mut_s,mut_e,mut_r,mut_q
				arr = re.split(r',',obj_poly2_ref.SnpIndel[p2_num])
				start = str(int(arr[1]) + int(obj_ref.start) - 1)

				if max_ref_alt_len < len(arr[2]):
					max_ref_alt_len = len(arr[2])
				if max_ref_alt_len < len(arr[3]):
					max_ref_alt_len = len(arr[3])

				if arr[3] == "G" :
					alt_num['G'] += 1
				elif arr[3] == "C":
					alt_num['C'] += 1
				elif arr[3] == "T":
					alt_num['T'] += 1
				elif arr[3] == "A":
					alt_num['A'] += 1

				flag_region = 0
				tmp_i = 0
				for region_s in region_start:
					region_e = region_end[tmp_i]
					if int(start) >= int(region_s) and int(start) <= int(region_e) :
						flag_region = 1
						break
					tmp_i += 1

				if flag_region == 0:
					out_log.write(" poly2_out_of_region : %s n: %s\n" % (obj_poly2_ref.SnpIndel[p2_num], start))
					print " poly2_out_of_region : %s n: %s\n" % (obj_poly2_ref.SnpIndel[p2_num], start)
					continue

				mut_key = "%s,%s,%s,%s" % (obj_ref.chr, start, arr[2], arr[3])

				if muts_sum.has_key(mut_key):
					muts_sum[mut_key] += 1
				else:
					muts_sum[mut_key] = 1
				muts_info[mut_key] = "%s,%s" % (arr[0], obj_ref.region)

			max_alt_num = 0
			max_alt_base = ""
			for base in alt_num:
				if max_alt_num < alt_num[base]:
					max_alt_num = alt_num[base]
					max_alt_base = base

			if max_alt_num >= int(cf.get("Filter", "max_alt_num")) or max_ref_alt_len >= int(cf.get("Filter", "max_ref_alt_len")):
				qc_no = ""
				if max_alt_num >= int(cf.get("Filter", "max_alt_num")):
					qc_no = "max_alt_num:%s:%s(<%s)" % (max_alt_base, max_alt_num, int(cf.get("Filter", "max_alt_num")))
				if max_ref_alt_len > int(cf.get("Filter", "max_ref_alt_len")):
					qc_no = "max_ref_alt_len:%s(<%s)" % (max_ref_alt_len, int(cf.get("Filter", "max_ref_alt_len")))

				report_QC.write("%s\t%s\t%s\n" % (tmp_key, "low quality", qc_no ) )			
				draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), 0, 0, int(cf.get("Report", "draw_part_len")), {}, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
				draw.draw_master()
				out_log.write("\n\nMut_result\nNO muts\n")
				out_log.write("low quality : %s \n" % qc_no)
				out_mask.write("mask(start:%s,end:%s)" % (0, 0))

				print "low quality : %s " % qc_no 
				print "draw_pic"
				print "null_report"
				continue	

			# merge one ab1

			hete_sum = 0
			hete_report = 0
			hete_indel_pos = 1000
			for mut_key in list(muts_info.keys()):

				# filter N
				obj_N = re.search(r'N',mut_key,re.M|re.I)
				if obj_N:
					del muts_info[mut_key]
					out_log.write(" mut : %s\n" % mut_key)
					out_log.write(" filter : N in %s\n" % mut_key)
					print " mut : %s" % mut_key
					print " filter : N in %s" % mut_key
					continue	

				arr_key = re.split(r',', mut_key)
				tag = 'INDEL'
				if len(arr_key[2])==1 and len(arr_key[3])==1 and arr_key[2] != '-' and arr_key[3] != '-':
					tag = 'SNP'
				if muts_sum[mut_key] == 2:
					muts_info[mut_key] += ',homo'
				else:
					muts_info[mut_key] += ',hete'
					hete_sum += 1

				
				arr_info = re.split(r',',muts_info[mut_key])
				# filter ---> anno

				print "  mut_KEY : %s " % mut_key	
				out_log.write("  mut_KEY : %s \n" % mut_key)		
				######################################################################### filter SNP
				if tag == 'SNP' : # snp
					##### filter ---> ratio #(p1, p2, h1, h2, arr[1])
					
					if muts_sum[mut_key] == 1 :
						now_info = re.split(r',', obj_pred.base_cor_info[int(arr_info[0])])		

						out_log.write("	%s snp_raw_ratio : %s\n" % (mut_key, float(now_info[2])/float(now_info[3])))
						#print int(arr_info[0])
						#print now_info
						if float(now_info[2])/float(now_info[3]) > float(cf.get("Filter", "snp_ratio")) :
							del muts_info[mut_key]
							out_log.write("	filter raw_ratio : %s\n" % (float(now_info[2])/float(now_info[3])))
							print "	filter raw_ratio : %s" % (float(now_info[2])/float(now_info[3]))
					##### filter ---> qual	
					out_log.write("	%s snp_qual : %s\n" % (mut_key, obj_pred.qual[int(arr_info[0])] ))
					out_log.write("	%s snp_per_qual : %s\n" % (mut_key, obj_pred.qual_nearby[int(arr_info[0])] ))

					if mut_key in muts_info and ( obj_pred.qual[int(arr_info[0])] < int(cf.get("Filter", "snp_qual")) or obj_pred.qual_nearby[int(arr_info[0])] < int(cf.get("Filter", "snp_per_qual"))):
						del muts_info[mut_key]
						out_log.write("	filter qual : %s %s\n" % (obj_pred.qual[int(arr_info[0])], obj_pred.qual_nearby[int(arr_info[0])]))
						print "	filter qual : %s %s" % (obj_pred.qual[int(arr_info[0])], obj_pred.qual_nearby[int(arr_info[0])])

				##### filter INDEL					
				elif tag == 'INDEL': 
					
					if muts_sum[mut_key] == 1 :
						if hete_indel_pos > int(arr_info[0]) and len(arr_key[2]) != len(arr_key[3]): 
							hete_indel_pos = int(arr_info[0])

						# tmp add ?  ruijin 
						if int(cf.get("Filter", "ruijin_model")) != 0 and len(arr_key[2]) == len(arr_key[3]):
							out_log.write("	filter ruijin_INDEL : %s \n" % (muts_sum[mut_key]))
							print "	filter ruijin_INDEL : %s " % (muts_sum[mut_key])
							del muts_info[mut_key]

						'''
						if int(obj_pred.p1_p2_s_ratio) > int(cf.get("Filter", "hete_indel_s_ratio")) :
							out_log.write("	filter INDEL : %s %s\n" % (muts_sum[mut_key], obj_pred.p1_p2_s_ratio))
							print "	filter INDEL : %s %s" % (muts_sum[mut_key], obj_cor.similar)
							del muts_info[mut_key]
						'''
				else:
					print "	no_filter"

			
			if muts_info : # 
			
				q_pos_dict = {}
				for mut_tmp in list(muts_info.keys()):
					hh = re.split(r',', muts_info[mut_tmp])[2]
					

					# filter after indel --------
					q_pos = int(re.split(r',', muts_info[mut_tmp])[0])
					if (q_pos > hete_indel_pos and real_strand == 'f') or (q_pos < hete_indel_pos and real_strand == 'r' and hete_indel_pos != 1000) :
						del muts_info[mut_tmp]
						out_log.write(" mut : %s\n" % mut_tmp)
						out_log.write(" filter : after hete_indel\n")
						print " mut : %s" % mut_tmp
						print " filter : after hete_indel"
						continue					
					# filter target region

					#++++++++++++++

					if q_pos in obj_poly1_dump.query_ref_pos:
						r_pos = obj_poly1_dump.query_ref_pos[q_pos]
						#print "mut: %s | rpos_f : %s | rpos_r : %s" % (mut_tmp, r_pos, obj_dump.seq_len - r_pos - 1)
						if real_strand == 'r':
							r_pos = obj_dump.seq_len - r_pos - 1
						q_pos_dict[r_pos] = "%s,%s" % (mut_tmp, hh)
					# obj_poly1_dump.query_ref_pos[int(arr_info[0])]

				#print "q_pos_dict : %s" % str(q_pos_dict)
				draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), mask_s, mask_e, int(cf.get("Report", "draw_part_len")), q_pos_dict, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
				draw.draw_master()

			if muts_info :
				#anno snpeff
				report_QC.write("%s\t%s\t%s\n" % (tmp_key, "good quality", "have muts") )
				for mut_tmp in list(muts_info.keys()):
					if Mut_keys.has_key(mut_tmp):
						Mut_keys[mut_tmp] += ";%s:%s" %  (tmp_key, hh)

					else:
						Mut_keys[mut_tmp] = "%s:%s" %  (tmp_key, hh)

					arr= re.split(r',',mut_tmp)
					print arr

					raw_vcf.write("%s\t%s\t%s\t%s\t%s\n" % (arr[0], arr[1],".", arr[2], arr[3]))

			else:
				report_QC.write("%s\t%s\t%s\n" % (tmp_key, "good quality", "no muts") )
				q_pos_dict = {}
				if(cf.get("Report", "draw_pos") == "1"):
					q_pos_dict = get_pos(obj_ref,obj_dump,cf.get("Ref_Pos",real_key),real_strand)
				draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), mask_s, mask_e, int(cf.get("Report", "draw_part_len")), q_pos_dict, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
				draw.draw_master()
	
			print "qual_is_ok"	
			print "final_result : %s " % muts_info
			out_log.write("\n\nMut_result\nqual_is_ok\n" + str(muts_info))

		else:
			report_QC.write("%s\t%s\t%s\n" % (tmp_key, "good quality", "no muts") )
			q_pos_dict = {}
			if(cf.get("Report", "draw_pos") == "1"):
				q_pos_dict = get_pos(obj_ref,obj_dump,cf.get("Ref_Pos",real_key),real_strand)
			draw = draw_pic.Draw(obj_dump, "%s/%s" % (outdir_report,tmp_key), mask_s, mask_e, int(cf.get("Report", "draw_part_len")), q_pos_dict, "header_tag", "tail_tag", int(cf.get("Report", "draw_complete")), int(cf.get("Report", "draw_part")))
			draw.draw_master()
			print "NO muts\n"
			out_log.write("\n\nMut_result\nNO muts\nqual_is_ok\n")

		out_mask.write("mask(start:%s,end:%s)" % (mask_s, mask_e))
		out_log.close()
		out_mask.close()
		out_dump_seq.close()

	# snpeff anno
	raw_vcf.close()

	dir_pwd = (os.path.split(os.path.realpath(__file__))[0])

	cosmic_dict = {}
	if int(cf.get("Report", "cosmic")) == 1:
		cosmic_dict = anno_cosmic.LoadAnno("%s/O_Database/hg38_cosmic85.sort.GIST.txt" % dir_pwd )

	if Mut_keys:
		
		snpeff = "%s/D_Anno/snpEff_latest_core/snpEff/snpEff.jar" % dir_pwd


		hgvs1LetterAa = " "
		if int(cf.get("anno", "hgvs1LetterAa")) == 1:
			hgvs1LetterAa = " -hgvs1LetterAa "

		#/share/biosoft/Software/jdk8/bin/java

		genome = cf.get("anno", "genome")
		#com_snpeff = "%s -Xmx4g -jar %s %s %s -noStats %s > %s" % (javapath, snpeff, genome, outdir_report_root+"/raw_anno.vcf", hgvs1LetterAa, outdir_report_root+"/result_anno.vcf")
				
		#RunCommand(com_snpeff,"anno snpeff")

		Mut_info = {}
		#Mut_info = get_snpeff_info(outdir_report_root+"/result_anno.vcf")

		#04 merge all report
		if Mut_info:
			for tmp_mut in Mut_info:
				if int(cf.get("Report", "cosmic")) == 1:
					cosmic_id = anno_cosmic.GetCosmic(tmp_mut,cosmic_dict)
					report_Mut.write("%s\t%s\t%s\t%s\n" % (Mut_keys[tmp_mut], tmp_mut, Mut_info[tmp_mut], cosmic_id))
				else:
					report_Mut.write("%s\t%s\t%s\n" % (Mut_keys[tmp_mut], tmp_mut, Mut_info[tmp_mut]))

	report_QC.close()
	report_Mut.close()

def filter_raw_mut(SnpIndel1, SnpIndel2, strand):

	new_SnpIndel1 = {}
	new_SnpIndel2 = {}
	hete_indel_pos1 = 0
	hete_indel_pos2 = 0

	if SnpIndel1:
		for num in SnpIndel1:
			#mut_q_pos,mut_s,mut_r,mut_q
			arr = re.split(r',',SnpIndel1[num])
			if len(arr[2]) != len(arr[3]):
				tmp_pos = int(arr[0])
				if strand == 'f' and (hete_indel_pos1 > tmp_pos or hete_indel_pos1 == 0)  :
					hete_indel_pos1 = tmp_pos
				if strand == 'r' and hete_indel_pos1 < tmp_pos:
					hete_indel_pos1 = tmp_pos
	if SnpIndel2:
		for num in SnpIndel2:
			#mut_q_pos,mut_s,mut_r,mut_q
			arr = re.split(r',',SnpIndel2[num])
			if len(arr[2]) != len(arr[3]):
				tmp_pos = int(arr[0])
				if strand == 'f' and (hete_indel_pos2 > tmp_pos or hete_indel_pos2 == 0) :
					hete_indel_pos2 = tmp_pos
				if strand == 'r' and hete_indel_pos2 < tmp_pos:
					hete_indel_pos2 = tmp_pos

	hete_indel_pos = 0
	if hete_indel_pos1 != 0:
		hete_indel_pos = hete_indel_pos1
	elif hete_indel_pos2 != 0:
		hete_indel_pos = hete_indel_pos2

	print "\nhete1: %s  hete2: %s  hete: %s\n" % (hete_indel_pos1, hete_indel_pos2, hete_indel_pos)

	if hete_indel_pos2 != 0 and hete_indel_pos1!= 0 :
		if strand == 'f':
			if hete_indel_pos2 > hete_indel_pos1:
				hete_indel_pos = hete_indel_pos1
			else:
				hete_indel_pos = hete_indel_pos2
		else:
			if hete_indel_pos2 < hete_indel_pos1:
				hete_indel_pos = hete_indel_pos1
			else:
				hete_indel_pos = hete_indel_pos2

	new_num = 0
	if hete_indel_pos != 0 :
		if strand == 'f':
			if SnpIndel1:
				for num in SnpIndel1:
					arr = re.split(r',',SnpIndel1[num])
					tmp_pos = int(arr[0])
					if tmp_pos <= hete_indel_pos:
						print "%s <= %s" % (tmp_pos, hete_indel_pos) 
						new_SnpIndel1[new_num] = SnpIndel1[num]
						new_num += 1
			if SnpIndel2:
				for num in SnpIndel2:
					arr = re.split(r',',SnpIndel2[num])
					tmp_pos = int(arr[0])
					if tmp_pos <= hete_indel_pos:
						print "%s <= %s" % (tmp_pos, hete_indel_pos) 
						new_SnpIndel2[new_num] = SnpIndel2[num]
						new_num += 1
			hete_indel_pos += 10
		else:
			if SnpIndel1:
				for num in SnpIndel1:
					arr = re.split(r',',SnpIndel1[num])
					tmp_pos = int(arr[0])
					if tmp_pos >= hete_indel_pos:
						new_SnpIndel1[new_num] = SnpIndel1[num]
						new_num += 1
			if SnpIndel2:
				for num in SnpIndel2:
					arr = re.split(r',',SnpIndel2[num])
					tmp_pos = int(arr[0])
					if tmp_pos >= hete_indel_pos:
						new_SnpIndel2[new_num] = SnpIndel2[num]
						new_num += 1
			hete_indel_pos -= 10
	else:
		new_SnpIndel1 = SnpIndel1
		new_SnpIndel2 = SnpIndel2		


	return (new_SnpIndel1, new_SnpIndel2, hete_indel_pos)

def get_snpeff_info(file):
	snpeff = open(file,'r')
	gene_ok = {'UGT1A1':'UGT1A1',}
	Mut_info = {}
	for line in snpeff:
		line = line.strip()
		if not line or line[0] == '#':
			continue
		arr = re.split(r'\t', line)

		mut_key = "%s,%s,%s,%s" % (arr[0],arr[1],arr[3],arr[4])

		# select first
		info = (re.split(r',', arr[7]))

		print_info = ""

		flag = 0
		first_anno = ""
		flag_pro = 0
		first_line = 0
		for per_anno in info:
			gene_info = re.split(r'\|',per_anno)
			i = 0
			for tmp in gene_info:
				if not tmp:
					gene_info[i] = '-'
				i += 1

			if gene_info[0] != "ANN=%s" % arr[4] and gene_info[0] != arr[4] :
				continue

			obj_nm = re.search(r'NM',gene_info[6],re.M|re.I)
			if not obj_nm:
				continue

			first_line += 1
			if first_line == 1 :
				first_anno = "%s\t%s\t%s\t%s\t%s\t%s" % (gene_info[3],gene_info[6],gene_info[7],gene_info[1],gene_info[9],gene_info[10])
			if gene_info[10] != '-' and flag_pro == 0:
				flag_pro = 1

				first_anno = "%s\t%s\t%s\t%s\t%s\t%s" % (gene_info[3],gene_info[6],gene_info[7],gene_info[1],gene_info[9],gene_info[10])
				
			if gene_info[3] in gene_ok:
				flag = 1
				print_info = "%s\t%s\t%s\t%s\t%s\t%s" % (gene_info[3],gene_info[6],gene_info[7],gene_info[1],gene_info[9],gene_info[10])

		if flag == 0:
			print_info = first_anno
		# 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'
		# ANN=A|synonymous_variant|LOW|SLC25A18|SLC25A18|transcript|XM_011546149.2|protein_coding|3/11|c.102G>A|p.Ala34Ala|1284
		# 2	233760233	.	C	CAT	.	.	
		# CAT|upstream_gene_variant|MODIFIER|LOC100286922|LOC100286922|transcript|NR_037695.1|pseudogene||n.-4890_-4889dupAT|||||4889|,
		Mut_info[mut_key] = print_info

	return Mut_info

def main():
	usage = '''
	
	Copyright (c) BLD 2018
	Description:	Sanger_Pipline
	Author:		zhaoerying@celloud.cn
	Date:		2018/05/21
	Version:	$ver

	usage: %prog [options] arg1 arg2	''' 
	parser = OptionParser(usage=usage)
	parser.add_option('-i', '--input',	help='input input file (required)', dest='input')
	parser.add_option('-p', '--conf',	help='input conf parameters file (required)', dest='conf')
	parser.add_option('-j', '--java',	help='input java path', dest='java')

	(options, args) = parser.parse_args()

	if not options.input or not options.conf :
		parser.print_help()
		exit(1)

	if not os.path.exists(options.input):
		printtime('ERROR: No input file found at: ' + options.input)
		sys.exit(1)

	if not os.path.exists(options.conf): 
		printtime('ERROR: No conf file found at: ' + options.conf)
		sys.exit(1)
	java = 'java'
	if options.java and os.path.exists(options.java):
		java = options.java


	Master(options.input, options.conf, java)


if __name__ == '__main__':
    main()

