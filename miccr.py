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
import gc
import re
import subprocess
import pandas as pd
import taxonomy as t
import numpy as np
from bitstring import BitArray, BitStream
from multiprocessing import Pool

def parse_params(ver):
    class SmartFormatter(ap.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines()  
            # this is the RawTextHelpFormatter._split_lines
            return ap.HelpFormatter._split_lines(self, text, width)

    p = ap.ArgumentParser(prog=sys.argv[0], description="""MInimap2 Contig ClassifieR (MICCR) %s""" % ver, formatter_class=SmartFormatter)

    eg = p.add_mutually_exclusive_group(required=True)

    eg.add_argument('-i', '--input', 
            metavar='[FASTQ]', nargs='*', type=str,
                    help="Input one or multiple FASTQ file(s). Use space to separate multiple input files.")

    eg.add_argument('-f', '--paf', 
            metavar='[PAF]', type=str,
                    help="Input one PAF file.")

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

    p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
                    help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

    p.add_argument( '-c','--cpu', metavar='<INT>', type=int, default=1,
                    help="Number of cpus [default: 1]")

    p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
                    help="Output directory [default: .]")

    p.add_argument( '--silent', action="store_true",
                    help="Disable all messages.")

    p.add_argument( '--debug', action="store_true",
                    help="Debug mode. Provide verbose running messages and keep all temporary files.")

    args_parsed = p.parse_args()

    """
    Checking options
    """
    if args_parsed.input and not args_parsed.database:
        p.error( '--database option is missing.' )

    if args_parsed.input and args_parsed.paf:
        p.error( '--input and --paf are incompatible options.' )

    if not args_parsed.dbPath:
        if args_parsed.database and "/" in args_parsed.database[0] and os.path.isfile( args_parsed.database[0] + ".amb" ):
            db_dir = re.search( '^(.*?)[^\/]+$', args_parsed.database[0] )
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
            name = re.search('([^\/\.]+)\..*$', args_parsed.input[0] )
            args_parsed.prefix = name.group(1)
        elif args_parsed.paf:
            name = re.search('([^\/]+).\w+.\w+$', args_parsed.paf )
            args_parsed.prefix = name.group(1)
        else:
            args_parsed.prefix = "miccr"

    return args_parsed

def dependency_check(cmd):
    proc = subprocess.Popen("which " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outs, errs = proc.communicate()
    return outs.decode().rstrip() if proc.returncode == 0 else False

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

def get_extra_regions(mask, qstart, qend, qlen, add=None):
    # init bitstring "add" if it's empty
    if not add:
        add = BitArray( "0b%s%s%s"%("0"*qstart, "1"*(qend-qstart), "0"*(qlen-qend)) )
    add_mask = add&(mask^add)
    mask = mask|add
    
    p = re.compile('1+')
    iterator = p.finditer(add_mask.bin)
    
    return (mask, [match.span() for match in iterator])

def aggregate_ctg(ctg_df):
    """
    aggregate alignments
    """
    qlen = ctg_df.iloc[0].qlen.item()
    ctg_mask = BitArray("0b"+"0"*qlen)

    # aggregate region and find LCA taxonomy
    df1 = ctg_df.groupby(['ctg','qstart','qend'])['taxid'].aggregate(
        {'lca_taxid': t.lca_taxid, 'hit_cnt':'count'}
    )

    df1['lca_rank'] = df1['lca_taxid'].apply(t.taxid2rank)
    df1['lca_name'] = df1['lca_taxid'].apply(t.taxid2name)

    df2 = ctg_df.groupby(['ctg','qstart','qend']).aggregate(
        {'match_bp': sum, 'mapping_bp': sum, 'score': 'first'}
    )

    ctg_df_agg = pd.concat([df1, df2], axis=1)
    ctg_df_agg.sort_values(by=['score','qstart','qend'], ascending=False, inplace=True)
    ctg_df_agg.reset_index(level=['ctg','qstart','qend'], inplace=True)
    
    ctg_df_agg['qlen'] = qlen
    ctg_df_agg['region'] = np.nan
    ctg_df_agg['agg_len'] = int
    
    # calculate
    for i in range(0, len(ctg_df_agg)):
        qstart = ctg_df_agg.loc[i,'qstart']
        qend   = ctg_df_agg.loc[i,'qend']
        (ctg_mask, regions) = get_extra_regions(ctg_mask, qstart, qend, qlen)

        if regions:
            rlen=0
            for r in regions: rlen += r[1]-r[0]
            ctg_df_agg.loc[i,'agg_len'] = rlen
            ctg_df_agg.loc[i,'region'] = str(regions)
    
    return ctg_df_agg.dropna()

def processPAF(paf, cpus):
    df = pd.read_csv(
        paf,
        sep='\t',
        header=None,
        names=['qname','qlen','qstart','qend','strand','tname','tlen','tstart','tend','match_bp','mapping_bp','mqua','tp','cm','score','opt1','opt2'],
        index_col=['qname'],
        usecols=['qname','qlen','qstart','qend','strand','tname','tlen','tstart','tend','match_bp','mapping_bp','score'],
    )
    df['ctg'] = df.index

    print_message( "Done.", argvs.silent, begin_t, logfile )

    # sorting by contig, score, then qstart and qend
    df.sort_values(by=['ctg', 'score', 'qstart', 'qend'], ascending=False, inplace=True)

    # only keep rows with max score for the same mapped regions
    df['score_max'] = df.groupby(['ctg','qstart','qend'])['score'].transform(max)
    df = df[ df['score']==df['score_max'] ]
    df['taxid'] = df['tname'].apply(t.acc2taxid)
    df['score'] = df['score'].str.replace('s1:i:','').astype(int)

    #clean memory
    gc.collect()

    print_message( "Arggregating alignments using %s subprocesses..."%cpus, argvs.silent, begin_t, logfile )
    pool = Pool(processes=cpus)
    jobs = []
    results = []

    cnt=0
    tol_jobs=len(df)
    for ctgname in df.index.unique().tolist():
        print(ctgname)
        results.append( aggregate_ctg(df.loc[[ctgname]]) )
        cnt+=1
        if argvs.debug: print_message( "[DEBUG] Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

    # for ctgname in df.index.unique().tolist():
    #     jobs.append( pool.apply_async(aggregate_ctg, df.loc[[ctgname]]) )

    # #wait for all jobs to finish
    # tol_jobs = len(jobs)
    # cnt=0
    # for job in jobs:
    #     results.append( job.get() )
    #     cnt+=1
    #     if argvs.debug: print_message( "[DEBUG] Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

    #clean up
    pool.close()

    return pd.concat(results)

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
        readMapping( argvs.input, argvs.database, argvs.cpu, paf, logfile )
        print_message( "Done mapping reads.", argvs.silent, begin_t, logfile )

    print_message( "Processing PAF file... ", argvs.silent, begin_t, logfile )
    dfctg = processPAF(paf, argvs.cpu)
    print_message( "Done processing PAF file, %s alignment(s)."% (len(dfctg['ctg'].unique().tolist())), argvs.silent, begin_t, logfile )

    print(dfctg)
