#!/usr/bin/env python3

import taxonomy as t

import sys;
t.loadTaxonomy( sys.argv[1] )

# ctg.tsv
# ./add_lineage.py taxonomy_db test/HMPeven_2.ctg.tsv  > test/HMPeven.ctg.tsv.lineage
file=sys.argv[2];
f =  open (file, "r")
rank_l=['superkingdom','phylum','class','order','family','genus','species']
header = f.readline().strip()
print("%s\t%s" % (header, "\t".join(rank_l)))

for line in f:
    if not line.strip():continue
    rank_n=[]
    temp = [ x.strip() for x in line.split('\t')]
    for i in range(len(rank_l)):
        name = t.taxid2nameOnRank(temp[4],rank_l[i])
        rank_n.append(name)

    print("%s\t%s" % ("\t".join(temp),"\t".join(rank_n)))



