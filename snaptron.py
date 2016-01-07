#!/usr/bin/env python

import sys
import os
import subprocess
import re
import shlex
from collections import namedtuple
import urllib2
import operator
import time

operators={'>=':operator.ge,'<=':operator.le,'>':operator.gt,'<':operator.lt,'=':operator.eq,'!=':operator.ne}
DEBUG_MODE=False
TABIX="tabix"
TABIX_INTERVAL_DB='all_SRA_introns_ids_stats.tsv.gz'
TABIX_DB_PATH='/data2/gigatron2'
TABIX_DBS={'chromosome':TABIX_INTERVAL_DB,'samples_count':'sample_count.gz','coverage_sum':'by_coverage_sum.gz','coverage_avg':'by_coverage_avg.gz','coverage_median':'by_coverage_median.gz'}
SAMPLE_MD_FILE='/data2/gigatron2/all_illumina_sra_for_human_ids.tsv'
SAMPLE_IDS_COL=7
SAMPLE_ID_COL=0
INTRON_ID_COL=0

INTRON_URL='http://localhost:8090/solr/gigatron/select?q='
SAMPLE_URL='http://localhost:8090/solr/sra_samples/select?q='

INTRON_HEADER='gigatron_id	chromosome	start	end	strand	donor	acceptor	samples	read_coverage_by_sample	samples_count	coverage_count	coverage_sum	coverage_avg	coverage_median'
SAMPLE_HEADER=""
INTRON_HEADER_FIELDS=INTRON_HEADER.split('\t')
INTRON_HEADER_FIELDS_MAP={}
i=0
for field in INTRON_HEADER_FIELDS:
   INTRON_HEADER_FIELDS_MAP[field]=i 
   i+=1


def run_tabix(qargs,rquerys,tabix_db,filter_set=None,sample_set=None,filtering=False,debug=True):
    tabix_db = "%s/%s" % (TABIX_DB_PATH,tabix_db)
    if debug:
        sys.stderr.write("running %s %s %s\n" % (TABIX,tabix_db,qargs))
    if not filtering:
        sys.stdout.write("Type\t%s\n" % (INTRON_HEADER))
    ids_found=set()
    tabixp = subprocess.Popen("%s %s %s" % (TABIX,tabix_db,qargs),stdout=subprocess.PIPE,shell=True)
    for line in tabixp.stdout:
        fields=[]
        #build either filter set or sample set or both
        if sample_set != None or filter_set != None:
             fields=line.rstrip().split("\t")
             if filter_set != None and fields[INTRON_ID_COL] not in filter_set:
                 continue
             if filtering:
                 ids_found.add(fields[INTRON_ID_COL])
                 continue
             if sample_set != None:
                 sample_ids=set(fields[SAMPLE_IDS_COL].split(","))
                 sample_set.update(sample_ids)
        #filter return stream based on range queries (if any)
        if rquerys:
            if len(fields) == 0:
                fields=line.rstrip().split("\t")
            skip=False
            for rfield in rquerys.keys():
                (op,val)=rquerys[rfield]
                if rfield not in INTRON_HEADER_FIELDS_MAP:
                    sys.stderr.write("bad field %s in range query,exiting\n" % (rfield))
                    sys.exit(-1)
                fidx = INTRON_HEADER_FIELDS_MAP[rfield]
                if not op(float(fields[fidx]),val):
                    skip=True
                    break
            if skip:
                continue
        #now just stream back the result
        if not filtering:
            sys.stdout.write("I\t%s" % (line))
    exitc=tabixp.wait() 
    if exitc != 0:
        raise RuntimeError("%s %s %s returned non-0 exit code\n" % (TABIX,TABIX_DB,args))
    if filtering:
        return ids_found

def stream_samples(sample_set,sample_map):
    sys.stdout.write("Type\t%s\n" % (SAMPLE_HEADER))
    for sample_id in sample_set:
        sys.stdout.write("S\t%s\n" % (sample_map[sample_id]))

#use to load samples metadata to be returned
#ONLY when the user requests by overlap coords OR
#coords and/or non-sample metadata thresholds (e.g. coverage)
#otherwise we'll return the whole thing from SOLR instead
def load_sample_metadata(file_):
    fmd={}
    #dont need the hash-on-column headers just yet
    with open(file_) as f:
       for line in f:
           fields=line.rstrip().split("\t")
           fmd[fields[0]]=line
    return fmd

import urllib
MAX_SOLR_ROWS=1000000000
def stream_solr(solr_query,filter_set=None,sample_set=None,debug_mode=False):
    header_just_id='intropolis_sample_id_i'
    solr_url = "%s%s&wt=csv&csv.separator=%%09&rows=%d&fl=%s" % (SAMPLE_URL,urllib.quote_plus(solr_query),MAX_SOLR_ROWS,header_just_id)
    if debug_mode:
        sys.stderr.write("opening %s\n" % (solr_url))
    solrR = urllib2.urlopen(solr_url)
    if debug_mode:
        sys.stderr.write("streaming solr results now (querying done)\n")
    line=solrR.readline()
    while(line):
        line=line.rstrip()
        fields=line.split("\t")
        if sample_set != None:
            sample_set.add(fields[SAMPLE_ID_COL]) 
        line=solrR.readline()

comp_op_pattern=re.compile(r'([=><!]+)')
def range_query_parser(rangeq):
    rquery={}
    if rangeq is None or len(rangeq) < 1:
        return (None,None,rquery)
    fields = rangeq.split(',')
    (first_tdb,first_rquery)=(None,None) #first_col,first_val1,first_val2)
    #for (field_,tdb) in TABIX_DBS.iteritems():
    for field in fields:
        m=comp_op_pattern.search(field)
        (col,op_,val)=re.split(comp_op_pattern,field)
        val=float(val)
        if not m or not col or col not in TABIX_DBS:
            continue
        op=m.group(1)
        if op not in operators:
            sys.stderr.write("bad operator %s in range query,exiting\n" % (str(op)))
            sys.exit(-1)
        rquery[col]=(operators[op],val)
        if first_tdb:
            continue 
        #only do the following for the first range query
        tdb=TABIX_DBS[col]
        first_tdb=tdb
        extension=""
        if op == '=':
            extension="-%d" % (val)
        if op == '<=':
            extension="-%d" % (val)
            val=1
        if op == '<':
            extension="-%d" % (val-1)
            val=1
        if op == '>':
            val=val+1
        first_rquery="1:%d%s" % (val,extension)
    return (first_tdb,first_rquery,rquery)

#this does the reverse: given a set of sample ids,
#return all the introns associated with each sample
import cPickle
def stream_introns_from_samples(sample_set):
    start = time.time()
    sample2introns={}
    if os.path.exists("./sample2introns.pkl"):
        #f=open("./introns.pkl","rb")
        #introns=cPickle.load(f)
        #f.close()
        f=open("./sample2introns.pkl","rb")
        sample2introns=cPickle.load(f)
        f.close()
    else:
        f=open("/data2/gigatron2/all_SRA_introns_ids_stats.tsv.new","r")
        for line in f:
            if "gigatron_id" in line:
                continue
            fields=line.rstrip().split('\t')
            sample_ids=fields[SAMPLE_IDS_COL].split(',')
            #introns[fields[0]]=line
            for sample_id in sample_ids:
                if sample_id not in sample2introns:
                    sample2introns[sample_id]=set()
                sample2introns[sample_id].add(int(fields[0]))
        f.close()
    if not os.path.exists("./sample2introns.pkl"):
        #f=open("./introns.pkl","wb")
        #cPickle.dump(introns,f,cPickle.HIGHEST_PROTOCOL)
        #f.close()
        f=open("./sample2introns.pkl","rb")
        cPickle.dump(sample2introns,f,cPickle.HIGHEST_PROTOCOL)
        f.close()
    end = time.time()
    taken=end-start
    print("loaded introns in %d" % (taken))
    introns_seen=set()
    for sample_id in sample_set:
        #for intron_id in sample2introns[sample_id]:
        #    if intron_id in seen:
        #        continue
        introns_seen.update(sample2introns[sample_id])
    sys.stdout.write("s2I\t%s\n" % (str(len(introns_seen))))
    #for intron_id in introns_seen:
    #    sys.stdout.write("I\t%s\n" % (introns[intron_id]))
    
#cases:
#1) just interval (one function call)
#2) interval + range query(s) (one tabix function call + field filter(s))
#3) one or more range queries (one tabix range function call + field filter(s))
#4) interval + sample (2 function calls: 1 solr for sample filter + 1 tabix using sample filter)
#5) sample (1 solr call -> use interval ids to return all intervals)
#6) sample + range query(s) (2 function calls: 1 solr for sample filter + 1 tabix using sample filter + field filter)
def main():
    input_ = sys.argv[1]
    DEBUG_MODE_=DEBUG_MODE
    if len(sys.argv) > 2:
        DEBUG_MODE_=True
    (intervalq,rangeq,sampleq) = input_.split('|')
    samples = load_sample_metadata(SAMPLE_MD_FILE)
    if DEBUG_MODE_:
        sys.stderr.write("loaded %d samples metadata\n" % (len(samples)))
    sample_set=set()
    if len(sampleq) >= 1:
       stream_solr(sampleq,sample_set=sample_set)
       if DEBUG_MODE_:
           sys.stderr.write("found %d introns in solr\n" % (len(filter_set)))
       if len(intervalq) == 0 and len(rangeq) == 0:
           stream_introns_from_samples(sample_set)
    #whether or not we use the interval query as a filter set or the whole query
    (first_tdb,first_rquery,rquery) = range_query_parser(rangeq)
    if len(intervalq) >= 1:
        run_tabix(intervalq,rquery,TABIX_INTERVAL_DB,sample_set=sample_set,debug=DEBUG_MODE_)
    elif len(rangeq) >= 1:
        run_tabix(first_rquery,rquery,first_tdb,sample_set=sample_set,debug=DEBUG_MODE_)
    if DEBUG_MODE_:
        sys.stderr.write("found %d samples\n" % (len(sample_set)))
    #stream_samples(sample_set,samples)

if __name__ == '__main__':
    main()
