#!/usr/bin/python
# Copyright (C) 2018-04-02 zhaoerying All Rights Reserved

import sys
import os
import subprocess
import time
import re

try:
	import json
except:
	import simplejson as json
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

class Draw:
	def __init__( self, obj_dump, out_file, mask_s, mask_e, extend_len, q_pos_dict, pic_header, pic_tail,big_pic = 1,small_pic=1):

		self.big_pic = big_pic
		self.small_pic = small_pic

		self.obj_dump = obj_dump
		self.out_file = out_file
		self.mask_s = mask_s
		self.mask_e = mask_e

		self.dump_mask_s_X = 0
		self.dump_mask_e_X = 0

		self.extend_len = extend_len
		self.q_pos_dict = q_pos_dict

		self.pic_header = pic_header
		self.pic_tail = pic_tail

		self.x_extend = 5
		self.pic_s = 0
		self.pic_e = self.obj_dump.NPoints
		self.hight = -8

		self.head_l = 0
		self.pic_region = 600
		self.tail_l = 0

	def draw_master(self):

		cut_pos = 1
		if self.big_pic == 1:
			cut_pos = 0
			self.draw_Pic(self.q_pos_dict, ".svg", cut_pos)

		if self.small_pic == 1:
			cut_pos = 1
			for q_pos in self.q_pos_dict: # mut_tag   \s \n is no
				self.pic_s = 0
				self.pic_e = self.obj_dump.NPoints
				self.q_pos = int(q_pos)
				print "self.pic_s: %s  self.pic_e: %s" % (self.pic_s, self.pic_e)
				self.draw_Pic(self.q_pos_dict, "." + self.q_pos_dict[q_pos] + ".svg", cut_pos)

	def draw_Pic(self, mut_pos_tag, name, cut_pos):

		# get pic_start pic_end
		self.get_draw_start_end(cut_pos)

		if cut_pos == 1:
			self.hight = -1   #it can change
			self.get_y_h(cut_pos)
		elif cut_pos == 0:
			self.hight = -8
			self.get_y_h(cut_pos)

		out = open(self.out_file+name, 'w')

		X = (int(self.pic_e) - int(self.pic_s)) * self.x_extend
		Y = self.pic_region + self.head_l + self.tail_l

		y1 = self.head_l
		y2 = self.pic_region + self.head_l

		line = svg_paper(X, Y)
		line += svg_line(0, y1, X, y1,"black",4)
		line += svg_line(0, y1, 0, y2,"black",4)
		line += svg_line(X, y2, 0, y2,"black",4)
		line += svg_line(X, y1, X, y2,"black",4)

		line += self.draw_line()
		line += self.draw_mut_tag(mut_pos_tag)
		if cut_pos == 0:
			line += self.draw_mask(X)

		line += svg_end()

		out.write(line)
		out.close()

		# svg2png
		pwd_dir = (os.path.split(os.path.realpath(__file__))[0])
		cmd = "java -jar %s/util/svg2png-0.0.1.jar %s %s" % (pwd_dir, self.out_file+name, self.out_file+name+".png")
		#RunCommand(cmd, "svg2png")

	def draw_mask(self,X):
		line = ""
		y1 = self.head_l
		y2 = self.pic_region + self.head_l
		if int(self.mask_s) == 0 and int(self.mask_e) == 0:
			line += svg_rect(0,y1,X,y2,"grey",0.3)
		else:
			line += svg_rect(0,y1,self.dump_mask_s_X +20,y2,"grey",0.3)
			line += svg_rect(self.dump_mask_e_X -20,y1,X,y2,"grey",0.3)

		return line

	def get_draw_start_end(self, cut_pos):

		#print "cut_pos == 1 (%s)" % cut_pos
		if int(cut_pos) == 1 :			
			for i in range(int(self.extend_len),1,-1):
				num_s = int(self.q_pos) - i - 1
				if str(num_s) in self.obj_dump.num_pos :
					self.pic_s = self.obj_dump.num_pos[str(num_s)]
					break
			for i in range(int(self.extend_len),1,-1):
				num_e = int(self.q_pos) + int(self.extend_len) + 1
				if str(num_e) in self.obj_dump.num_pos :
					self.pic_e = self.obj_dump.num_pos[str(num_e)]
					break

	def get_y_h(self, cut_pos):
		arr = []
		if int(cut_pos) == 0:
			for tmp_key in self.obj_dump.line_pos_base_high:
				h = self.obj_dump.line_pos_base_high[tmp_key]
				arr.append(h)
				ratio = 1
		elif int(cut_pos) == 1:
			for tmp_key in self.obj_dump.line_pos_base_high:
				art = re.split(r',',tmp_key)
				if int(art[1]) >= int(self.pic_s) and int(art[1]) <= int(self.pic_e) :
					h = self.obj_dump.line_pos_base_high[tmp_key]
					arr.append(h)
					ratio = 1.2
		arr1 = sorted(map(int,arr))
		#self.arr = arr1
		self.hight = arr1[int(self.hight)]*ratio

	def draw_base(self):

		pass

	def draw_line(self):
		#self.line_pos_base_high[base_key] = arr[0]
		line = ""

		for base_key in self.obj_dump.line_pos_base_high :
			arr = re.split(r',',base_key)
			if int(arr[1]) < int(self.pic_s) or int(arr[1]) > int(self.pic_e) : 
				continue
			color = "yellow"
			if arr[0] == "A":
				color = "green"
			elif arr[0] == "C":
				color = "blue"
			elif arr[0] == "G":
				color = "black"
			elif arr[0] == "T":
				color = "red"

			if self.obj_dump.line_pos_base_high.has_key("%s,%s" % (arr[0], int(arr[1])+1)) and int(arr[1])+1 < self.pic_e :
				new_x1 = ( int(arr[1]) - int(self.pic_s) ) * self.x_extend
				new_x2 = ( int(arr[1])+1 - int(self.pic_s) ) * self.x_extend
				ratio1 = float(self.obj_dump.line_pos_base_high[base_key])/float(self.hight)
				if ratio1 > 1:
					ratio1 = 1
				ratio2 =  float(self.obj_dump.line_pos_base_high["%s,%s" % (arr[0],int(arr[1])+1)])/float(self.hight)
				if ratio2 > 1:
					ratio2 = 1
				new_y1 = self.head_l + self.pic_region*(1-ratio1)
				new_y2 = self.head_l + self.pic_region*(1-ratio2)
				line += svg_line(new_x1, new_y1, new_x2, new_y2, color,2)

				#print "svg_line(new_x1: %s, new_y1: %s, new_x2: %s, new_y2: %s, h : %s,2)" % (new_x1, new_y1, new_x2, new_y2,self.obj_dump.line_pos_base_high[base_key])

		return line

	def draw_mut_tag(self, mut_pos_tag):

		line = ""
		
		for num in self.obj_dump.num_base:
			base = self.obj_dump.num_base[num]
			pos  = self.obj_dump.num_pos[num]

			if int(pos) <= int(self.pic_s) or int(pos) >= int(self.pic_e) : 
				continue

			color = "sandybrown"
			if base == "A":
				color = "green"
			elif base == "C":
				color = "blue"
			elif base == "G":
				color = "black"
			elif base == "T":
				color = "red"

			x = ( int(pos) - int(self.pic_s) ) * self.x_extend

			if int(self.mask_s) == int(num) :
				self.dump_mask_s_X = x 
			if int(self.mask_e) == int(num) :
				self.dump_mask_e_X = x

			y_num = self.head_l + 25
			y_base = self.head_l + 50

			line += svg_mid_txt(x, y_base, 20, color, base)
			line += svg_mid_txt(x, y_num, 20, 'black', num)

			if int(num) in mut_pos_tag:

				line += svg_circle(10, x,self.head_l+70,"darkslateblue");
				line += svg_mid_txt(x,105+self.head_l,20,"darkslateblue",mut_pos_tag[int(num)],0);

		return line

	def draw_ref_alt_seq(self):
		pass

	def pic_header(self):	
		pass

def svg_paper(width, height, fill = "white"):
	
	line = "<?xml version=\"1.0\" standalone=\"no\"?>\n"
	line += "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \n\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd \">\n\n"
	line += "<svg width=\"%s\" height=\"%s\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n" % (width, height)
	line += "<rect x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" fill=\"%s\" />\n" % (0, 0, width, height, fill)
	return line

def svg_polygon(fill, stroke, opacity, points):  #colorfill,colorstroke,coloropacity,point1,point2,...

	merge_point = " ".join(points)
	line = "<polygon fill=\"%s\" stroke=\"%s\" opacity=\"%s\" points=\"%s\" />\n" % (fill, stroke, opacity, merge_point)
	return line

def svg_circle(r, x, y, fill): #&svg_circle(x,y,r,color,[info])

	line = "<circle r=\"%s\" cx=\"%s\" cy=\"%s\" fill=\"%s\" />\n" % (r, x, y, fill)
	return line

def svg_txt(x, y, size, color, text, vertical = 0):  #&svg_txt(x,y,size,color,text,[vertical,0/1/2/3]);

	svg_matrix=''
	if vertical==0 :
		svg_matrix="1 0 0 1"	
	if vertical==1 :
		svg_matrix="0 1 -1 0"	
	if vertical==2 :
		svg_matrix="-1 0 0 -1"
	if vertical==3 :
		svg_matrix="0 -1 1 0"	
	line = "<text fill=\"%s\" transform=\"matrix(%s %s %s)\" font-family=\"ArialNarrow-Bold\" font-size=\"%s\">%s</text>\n" % (color, svg_matrix, x, y, size, text)
	return line

def svg_mid_txt(x, y, size, color, text, vertical = 0): #&svg_mid_txt(x,y,size,color,text,[vertical,0/1/2/3]);

	svg_matrix=''
	if vertical==0 :
		svg_matrix="1 0 0 1"	
	if vertical==1 :
		svg_matrix="0 1 -1 0"	
	if vertical==2 :
		svg_matrix="-1 0 0 -1"
	if vertical==3 :
		svg_matrix="0 -1 1 0"	
	line = "<text fill=\"%s\" transform=\"matrix(%s %s %s)\" text-anchor=\"middle\" font-family=\"ArialNarrow-Bold\" font-size=\"%s\">%s</text>\n" % (color, svg_matrix, x, y, size, text)
	return line

def svg_dashed(x1,y1,x2,y2,color,dash="10 5", width = 0):  #&svg_line(x1,y1,x2,y2,color,"10 5",[width])

	line = "<line x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\" style=\"stroke-dasharray:%s;fill:none;stroke:%s\"/>\n" % (x1,y1,x2,y2,dash,color)
	if width != 0 :
		line = "<line x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\" style=\"stroke-dasharray:%s;fill:none;stroke:%s;stroke-width:%s\"/>\n" % (x1,y1,x2,y2,dash,color,width)
	return line

def svg_line(x1,y1,x2,y2,color,width = 0):  #&svg_line(x1,y1,x2,y2,color,[width])

	line = "<line fill=\"%s\" stroke=\"%s\" x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\"/>\n" % (color,color,x1,y1,x2,y2)
	if width != 0 : 
		line = "<line fill=\"%s\" stroke=\"%s\" stroke-width=\"%s\" x1=\"%s\" y1=\"%s\" x2=\"%s\" y2=\"%s\"/>\n" % (color,color,width,x1,y1,x2,y2)
	return line


def svg_rect(x,y,width,height,color,opacity = 1):  #&svg_rect(x,y,width,height,color,[opacity])

	line = "<rect x=\"%s\" y=\"%s\" width=\"%s\" height=\"%s\" fill=\"%s\" opacity=\"%s\"/>\n" % (x,y,width,height,color,opacity)
	return line;


def svg_end():  #end

	return "</svg>\n"




def main():
	usage = '''
	
Copyright (c) BLD 2018
Description:	Draw_svg_for_ab1
Author:		zhaoerying@celloud.cn
Date:		2018/04/02
Version:	$ver

usage: %prog [options] arg1 arg2	''' 
	parser = OptionParser(usage=usage)
	parser.add_option('-v', '--vcf',	help='input VCF file (required)', dest='vcf')
	parser.add_option('-t', '--vcf-type',	help='VCF type (required)', dest='type')
	parser.add_option('-i', '--anno-info',	help='VCF annovar file (required)', dest='anno')
	parser.add_option('-k', '--file-name',	help='output file name (required)', dest='key')
	parser.add_option('-u', '--url-funtion',	help='get anno url (required)', dest='url')
	parser.add_option('-o', '--output-dir',	help='outputdir (default: current)', dest='outdir', default='.')

	(options, args) = parser.parse_args()

	if not options.vcf or not options.type or not options.key or not options.url or not options.outdir or not options.anno:
		parser.print_help()
		exit(1)

	if not os.path.exists(options.vcf):
		printtime('ERROR: No vcf file found at: ' + options.reference)
		sys.exit(1)

	if not os.path.exists(options.anno): 
		printtime('ERROR: No annovar file found at: ' + options.anno)
		sys.exit(1)

	if not os.path.exists(options.outdir):
		os.makedirs(options.outdir)


if __name__ == '__main__':
    main()
