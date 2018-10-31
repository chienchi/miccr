#!/usr/bin/env python3
__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = "0.0.1"
__date__      = "2018/10/29"
__copyright__ = "BSD-3"

import argparse as ap
import textwrap as tw
import sys
import os
import time
import subprocess
import pandas as pd
from re import search, findall
from multiprocessing import Pool
from itertools import chain
import taxonomy as t
import copy

def parse_params(ver):
	class SmartFormatter(ap.HelpFormatter):
		def _split_lines(self, text, width):
			if text.startswith('R|'):
				return text[2:].splitlines()  
			# this is the RawTextHelpFormatter._split_lines
			return ap.HelpFormatter._split_lines(self, text, width)

	p = ap.ArgumentParser(prog='pangia.py', description="""PanGIA Bioinformatics %s""" % ver, formatter_class=SmartFormatter)

	eg = p.add_mutually_exclusive_group(required=True)

	eg.add_argument('-i', '--input', 
			metavar='[FASTQ]', nargs='+', type=str,
	  				help="Input one or multiple FASTQ file(s). Use space to separate multiple input files.")

	p.add_argument('-d', '--database',
			metavar='[FASTA/MMI]', type=str, nargs=1,
					help="Name/path of readmapper's index [default: None]")

	p.add_argument('-dp', '--dbPath',
			metavar='[PATH]', type=str, default=None,
					help="""Path of databases. If dbPath isn't specified but a path is provided in "--database" option, this path of database will also be used in dbPath. 
					Otherwise, the program will search "database/" in program directory.
					[default: database/]""")

	p.add_argument('-x','--platform',
			type=str, default='asm10',
					choices=['asm5', 'asm10', 'map-pb', 'map-ont'],
					help="""R|You can specify one of the following platform:\n"""
						 """"asm5"    : compare mapping results with the background;\n"""
						 """"asm10"   : score based on uniqueness;\n"""
						 """"map-pb"  : bg * standalone;\n"""
						 """"map-ont" : bg * standalone;\n"""
						 """[default: 'asm10']""" )

	p.add_argument( '-c','--cpu', metavar='<INT>', type=int, default=1,
					help="Number of cpus [default: 1]")

	p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
					help="Output directory [default: .]")

	p.add_argument( '--debug', action="store_true",
					help="Debug mode. Provide verbose running messages and keep all temporary files.")

	args_parsed = p.parse_args()

	"""
	Checking options
	"""
	if args_parsed.input and not args_parsed.database:
		p.error( '--database option is missing.' )

	if args_parsed.input and args_parsed.sam:
		p.error( '--input and --same are incompatible options.' )

	if not args_parsed.dbPath:
		if args_parsed.database and "/" in args_parsed.database[0] and os.path.isfile( args_parsed.database[0] + ".amb" ):
			db_dir = search( '^(.*?)[^\/]+$', args_parsed.database[0] )
			args_parsed.dbPath = db_dir.group(1)
		else:
			bin_dir = os.path.dirname(os.path.realpath(__file__))
			args_parsed.dbPath = bin_dir + "/database"

	# glob database path
	if args_parsed.database:
		for db in args_parsed.database:
			if args_parsed.dbPath and not "/" in db:
				db = args_parsed.dbPath+"/"+db
			if not os.path.isfile( db ):
				p.error( 'Database not found: %s' % db )
			if os.path.isfile( db + ".mmi" ):
				db = db + ".mmi"
		args_parsed.database = db

	if not args_parsed.prefix:
		if args_parsed.input:
			name = search('([^\/\.]+)\..*$', args_parsed.input[0] )
			args_parsed.prefix = name.group(1)
		elif args_parsed.sam:
			name = search('([^\/]+).\w+.\w+$', args_parsed.sam[0].name )
			args_parsed.prefix = name.group(1)
		else:
			args_parsed.prefix = "pangia"

	if not args_parsed.tempdir:
		args_parsed.tempdir = args_parsed.outdir+"/"+args_parsed.prefix+"_tmp"

	if not args_parsed.singleEnd:
		args_parsed.singleEnd = "auto"

	return args_parsed

def dependency_check(cmd):
	proc = subprocess.Popen("which " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	outs, errs = proc.communicate()
	return outs.decode().rstrip() if proc.returncode == 0 else False

def isDescendant( taxid, taxid_ant ):
	fullLineage = t.taxid2fullLineage( taxid )
	if "|%s|" % taxid_ant in fullLineage:
		return True
	else:
		return False

def lineageLCR(taxids):
	""" lineageLCR
	
	Find late common rank (LCR) from two lineage

	Arguments:
	  arg1 <LIST>: taxids

	Returns:
	  level <STR>: LCR name of current and quering lineage
	  name <STR> : LCR level of current and quering lineage
	"""

	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']

	merged_dict = t._autoVivification()
	for tid in taxids:
		lng = t.taxid2lineageDICT(tid, 1, 1)
		for r in lng:
			ttid = lng[r]['taxid']
			if ttid in merged_dict[r]:
				merged_dict[r][ttid] += 1
			else:	
				merged_dict[r][ttid] = 1

	for r in ranks:
		if len(merged_dict[r]) == 1:
			for ttid in merged_dict[r]:
				# skip if no tid in this rank
				if ttid==0:
					continue
				tname = t.taxid2name(ttid)
				return r, tname, merged_dict

	return "root", "root", merged_dict

def flagInHeader( flag, header ):
	if flag:
		h = header.split('|')
		for f in flag:
			if f in h[3]:
				return True

	return False

def processPAF(paf, cpus):
	df = pd.from_

def readMapping( fa, db, cpus, platform, paf, logfile ):
	"""
	mapping fa to database
	"""
	sam_list = []
	num_input_contigs = 0

	print_message( "Mapping to %s..." % db, argvs.silent, begin_t, logfile )

	bash_cmd = "set -o pipefail; set -x;"
	mp_cmd   = "minimap2 -x %s -t%s %s %s" % (platform, cpus, db, fa)
	cmd      = "%s %s 2>> %s > %s" % (bash_cmd, mp_cmd, logfile, paf)

	if argvs.debug: print_message( "[DEBUG] CMD: %s"%cmd, argvs.silent, begin_t, logfile )

	proc = subprocess.Popen( cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
	outs, errs = proc.communicate()
	exitcode = proc.poll()

	if exitcode!=0:
		print_message( "[ERROR] error occurred while running read mapping (code: %s, message: %s)."%(exitcode, errs), argvs.silent, begin_t, logfile, True )

	return True

def print_message(msg, silent, start, logfile, errorout=0):
	message = "[%s] %s\n" % (timeSpend(start), msg)
	#loging
	with open( logfile, "a" ) as f:
		f.write( message )
		f.close()
	
	if errorout:
		sys.exit( message )
	elif not silent:
		sys.stderr.write( message )

def timeSpend( start ):
	done = time.time()
	elapsed = done - start
	return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

if __name__ == '__main__':
	argvs    = parse_params( __version__ )
	begin_t  = time.time()
	paf      = "%s/%s.paf" % (argvs.outdir, argvs.prefix) if not os.path.exists(argvs.paf) else argvs.paf
	logfile  = "%s/%s.log" % (argvs.outdir, argvs.prefix)

	#create output directory if not exists
	if not os.path.exists(argvs.outdir):
		os.makedirs(argvs.outdir)

	#load taxonomy
	print_message( "Loading taxonomy information...", argvs.silent, begin_t, logfile )
	t.loadTaxonomy( argvs.dbPath )
	print_message( "Done.", argvs.silent, begin_t, logfile )

	# if reads provided
	if argvs.input:
		print_message( "Running minimap2...", argvs.silent, begin_t, logfile )
		readMapping( argvs.input, argvs.database, argvs.cpus, paf, logfile )

	print_message( "Processing PAF file... ", argvs.silent, begin_t, logfile )
	processPAF(paf, argvs.cpus)
	print_message( "Done processing SAM file, %s alignment(s)."% (mapped_r_cnt), argvs.silent, begin_t, logfile )