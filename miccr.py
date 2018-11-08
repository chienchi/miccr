#!/usr/bin/env python3
__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__version__   = "0.0.1"
__date__      = "2018/10/29"
__copyright__ = "BSD-3"

import sys
import os
import time
import gc
import re
import argparse as ap
import subprocess
import taxonomy as t
import numpy as np
from multiprocessing import Pool
import pkg_resources
pkg_resources.require("pandas>=0.23.0")
import pandas as pd

def parse_params(ver):
    class SmartFormatter(ap.HelpFormatter):
        def _split_lines(self, text, width):
            if text.startswith('R|'):
                return text[2:].splitlines() 
            # this is the RawTextHelpFormatter._split_lines
            return ap.HelpFormatter._split_lines(self, text, width)

    p = ap.ArgumentParser(prog=sys.argv[0], description="""MInimap2 Contig ClassifieR (MICCR) %s"""%ver, formatter_class=SmartFormatter)

    eg = p.add_mutually_exclusive_group(required=True)

    eg.add_argument('-i', '--input',
            metavar='[FASTA]', type=str,
                    help="Input one or multiple contig files in FASTA format. Use space to separate multiple input files.")

    eg.add_argument('-f', '--paf',
            metavar='[PAF]', type=str,
                    help="Input a PAF alignment file.")

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
                         """"asm5"    : Long assembly to reference mapping (avg divergence < 5%%);\n"""
                         """"asm10"   : Long assembly to reference mapping (avg divergence < 10%%);\n"""
                         """"asm20"   : Long assembly to reference mapping (avg divergence < 20%%);\n"""
                         """"map-pb"  : PacBio/Oxford Nanopore read to reference mapping;\n"""
                         """"map-ont" : Slightly more sensitive for Oxford Nanopore to reference mapping;\n"""
                         """[default: 'asm10']""")

    p.add_argument( '-t','--numthreads', metavar='<INT>', type=int, default=1,
                    help="Number of cpus [default: 1]")

    p.add_argument( '-c','--stdout', action="store_true",
                    help="Output to STDOUT [default: False]")

    p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
                    help="Output directory [default: .]")

    p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
                    help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

    p.add_argument( '-m','--minLcaProp', metavar='<FLOAT>', type=float, default=0.1,
                    help="LCA classified segments more than a proportion of contig length [default: 0.1]")

    p.add_argument( '--silent', action="store_true",
                    help="Disable all messages.")

    p.add_argument( '-v','--verbose', action="store_true",
                    help="Provide verbose running messages.")

    args_parsed = p.parse_args()

    """
    Checking options
    """
    if args_parsed.input and not args_parsed.database:
        p.error( '--database option is missing.' )

    if args_parsed.input and args_parsed.paf:
        p.error( '--input and --paf are incompatible options.' )

    if not args_parsed.dbPath:
        if args_parsed.database and "/" in args_parsed.database[0]:
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

def contig_mapping( fa, db, cpus, platform, paf, logfile ):
    """
    mapping contig sequences to database using minimap2
    """
    sam_list = []
    num_input_contigs = 0

    print_message( "Mapping to %s..." % db, argvs.silent, begin_t, logfile )

    bash_cmd = "set -o pipefail; set -x;"
    mp_cmd   = "minimap2 -x %s -t%s %s %s" % (platform, cpus, db, fa)
    cmd      = "%s %s 2>> %s > %s" % (bash_cmd, mp_cmd, logfile, paf)

    if argvs.verbose: print_message( "[DEBUG] CMD: %s"%cmd, argvs.silent, begin_t, logfile )

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
    """
    Take a bitmask of a contig, a given alignment region (or in bitmask). This function
    will return combined bitmask and additional covered regions.
    """
    if not add:
        add = int( "%s%s%s"%("0"*qstart, "1"*(qend-qstart), "0"*(qlen-qend)), 2)
    
    add_mask = add&(mask^add)
    mask = mask|add
    
    p = re.compile('1+')
    bitstr = bin(add_mask).replace('0b','')
    bitstr = "0"*(qlen-len(bitstr))+bitstr
    iterator = p.finditer(bitstr)
    
    return (mask, [match.span() for match in iterator])

def aggregate_ctg(df, cnames):
    """
    aggregate alignments to taxonomic annotate contigs
    argvs: df <pd.DataFrame>, cnames <list>: list of contig indexes
    return: return a <pd.DataFrame> of annotation results
    """
    ctg_df_list=[]

    for cname in cnames:
        ctg_df = df.loc[[cname]]
        qlen = ctg_df.iloc[0].qlen.item()
        #ctg_mask = bitarray("0"*qlen)
        ctg_mask = int(0)

        # aggregate region and find LCA taxonomy
        ctg_df_agg = ctg_df.groupby(['ctg','qstart','qend']).aggregate(
            {'taxid': t.lca_taxid, 'match_bp': sum, 'mapping_bp': sum, 'score': 'first', 'tname':'count'}
        )

        ctg_df_agg = ctg_df_agg.rename(columns={'tname':'hit_count','taxid':'lca_taxid'})

        ctg_df_agg['lca_taxid'] = ctg_df_agg['lca_taxid'].astype(str)
        ctg_df_agg = ctg_df_agg[ctg_df_agg.lca_taxid != '0']

        # sort by score first because we want to process segments with best score first
        ctg_df_agg.sort_values(by=['score','qstart','qend'], ascending=False, inplace=True)
        ctg_df_agg.reset_index(level=['ctg','qstart','qend'], inplace=True)

        # get lca ranks and taxas
        ctg_df_agg['lca_rank'] = ctg_df_agg['lca_taxid'].apply(t.taxid2rank)
        ctg_df_agg['lca_name'] = ctg_df_agg['lca_taxid'].apply(t.taxid2name)

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
        
        ctg_df_list.append(ctg_df_agg.dropna())
    
    return pd.concat(ctg_df_list)

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
    print_message( "Done loading PAF file.", argvs.silent, begin_t, logfile )

    # only keep rows with max score for the same mapped regions
    if argvs.verbose: print_message( "Filtering out secondary alignments for each mapped segment...", argvs.silent, begin_t, logfile )
    df['score'] = df['score'].str.replace('s1:i:','').astype(int)
    df['score_max'] = df.groupby(['ctg','qstart','qend'])['score'].transform(max)
    df = df[ df['score']==df['score_max'] ]
    if argvs.verbose: print_message( "Done.", argvs.silent, begin_t, logfile )

    if argvs.verbose: print_message( "Converting acc# of mapped reference to taxid...", argvs.silent, begin_t, logfile )
    df['taxid'] = df['tname'].apply(t.acc2taxid)
    df = df[df.taxid != 'None'] # dropping alignments with no taxid
    if argvs.verbose: print_message( "Done.", argvs.silent, begin_t, logfile )

    #clean memory
    gc.collect()

    print_message( "Aggregating alignments using %s subprocesses..."%cpus, argvs.silent, begin_t, logfile )
    pool = Pool(processes=cpus)
    jobs = []
    results = []

    ctgnames = df.index.unique().tolist()
    n = 200 if len(ctgnames)/1500 < cpus else 1500
    chunks = [ctgnames[i:i + n] for i in range(0, len(ctgnames), n)]

    for chunk in chunks:
        jobs.append( pool.apply_async(aggregate_ctg, (df,chunk) ) )

    tol_jobs = len(jobs)
    cnt=0
    for job in jobs:
        results.append( job.get() )
        cnt+=1
        if argvs.verbose: print_message( "Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

    #clean up
    pool.close()

    return pd.concat(results)

def lca_aggregate_ctg(dfctg, minLenProp):
    #['CONTIG','LENGTH','START','END','LCA_TAXID','LCA_RANK','LCA_NAME','HIT_COUNT','SCORE','AGG_LENGTH','AVG_IDENTITY','REGION']
    lca_dfctg = dfctg[ dfctg['AGG_LENGTH'] > dfctg['LENGTH']*minLenProp ].groupby(['CONTIG']).aggregate({
        'LENGTH': 'first',
        'START': min,
        'END': max,
        'AGG_LENGTH': sum,
        'LCA_TAXID': t.lca_taxid,
        'MATCH_BP': sum,
        'MAPPING_BP': sum,
        'HIT_COUNT': sum,
        'SCORE': max,
        'REGION': lambda x: "[%s]"%', '.join(x.str.strip('[]'))
    })

    lca_dfctg['AVG_IDENTITY'] = lca_dfctg['MATCH_BP']/lca_dfctg['MAPPING_BP']
    lca_dfctg['LCA_RANK'] = lca_dfctg['LCA_TAXID'].apply(t.taxid2rank)
    lca_dfctg['LCA_NAME'] = lca_dfctg['LCA_TAXID'].apply(t.taxid2name)
    
    return lca_dfctg

if __name__ == '__main__':
    argvs    = parse_params(__version__)
    begin_t  = time.time()
    paf      = "%s/%s.paf" % (argvs.outdir, argvs.prefix) if not argvs.paf else argvs.paf
    logfile  = "%s/%s.log" % (argvs.outdir, argvs.prefix)
    outfile_ctg = sys.stdout if argvs.stdout else "%s/%s.ctg.tsv"%(argvs.outdir, argvs.prefix)
    outfile_lca = sys.stdout if argvs.stdout else "%s/%s.lca_ctg.tsv"%(argvs.outdir, argvs.prefix)

    print_message( "MInimap2 Contig ClassifieR (MICCR) v%s"%__version__, argvs.silent, begin_t, logfile )
    
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
        contig_mapping( argvs.input, argvs.database, argvs.numthreads, argvs.platform, paf, logfile )
        print_message( "Done mapping reads.", argvs.silent, begin_t, logfile )

    # processing PAF
    print_message( "Processing PAF file... ", argvs.silent, begin_t, logfile )
    dfctg = processPAF(paf, argvs.numthreads)
    print_message( "Done processing PAF file.", argvs.silent, begin_t, logfile )

    # output annotation of contig segments
    print_message( "Writing contig classification results...", argvs.silent, begin_t, logfile )
    dfctg['avg_identity'] = dfctg['match_bp']/dfctg['mapping_bp']
    dfctg = dfctg.rename(columns={'ctg':'contig', 'qstart':'start', 'qend':'end', 'qlen':'length', 'agg_len':'agg_length'})
    dfctg.columns = dfctg.columns.str.upper()

    display_cols=['CONTIG','LENGTH','START','END','LCA_TAXID','LCA_RANK','LCA_NAME','HIT_COUNT','SCORE','AGG_LENGTH','AVG_IDENTITY','REGION']

    dfctg.to_csv(
        outfile_ctg,
        sep='\t',
        header=True,
        index=False,
        columns=display_cols
    )
    print_message( "Done.", argvs.silent, begin_t, logfile )
    
    # output LCA of contig
    print_message( "Writing contig LCA classification results...", argvs.silent, begin_t, logfile )
    lca_dfctg = lca_aggregate_ctg(dfctg, argvs.minLcaProp)
    lca_dfctg.to_csv(
        outfile_lca,
        sep='\t',
        header=True,
        index=True,
        columns=display_cols[1:]
    )
    print_message( "Done.", argvs.silent, begin_t, logfile )