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

import gzip

import lucene
from java.io import File
from org.apache.lucene.analysis.standard import StandardAnalyzer
from org.apache.lucene.document import Document, Field
from org.apache.lucene.search import IndexSearcher
from org.apache.lucene.index import IndexReader
from org.apache.lucene.queryparser.classic import QueryParser
from org.apache.lucene.queryparser.classic import MultiFieldQueryParser
from org.apache.lucene.search import BooleanClause
from org.apache.lucene.store import SimpleFSDirectory
from org.apache.lucene.util import Version

import snaputil 

operators={'>=':operator.ge,'<=':operator.le,'>':operator.gt,'<':operator.lt,'=':operator.eq,'!=':operator.ne}
DEBUG_MODE=True
TABIX="tabix"
#TABIX_INTERVAL_DB='all_SRA_introns_ids_stats.tsv.gz'
TABIX_INTERVAL_DB='all_SRA_introns_ids_stats.tsv.new2_w_sourcedb2.gz'
TABIX_DB_PATH='/data2/gigatron2'
TABIX_DBS={'chromosome':TABIX_INTERVAL_DB,'length':'by_length.gz','snaptron_id':'by_id.gz','samples_count':'sample_count.gz','coverage_sum':'by_coverage_sum.gz','coverage_avg':'by_coverage_avg.gz','coverage_median':'by_coverage_median.gz'}
SAMPLE_MD_FILE='/data2/gigatron2/all_illumina_sra_for_human_ids.tsv'
SAMPLE_IDS_COL=12
SAMPLE_ID_COL=0
INTRON_ID_COL=0
LUCENE_MAX_HITS=1000
SAMPLE_QUERY_DELIMITER='==='
SAMPLE_QUERY_FIELD_DELIMITER='::'

FLOAT_FIELDS=set(['coverage_avg','coverage_median'])

DATA_SOURCE='SRA'
#may have to adjust this parameter for performance (# of tabix calls varies inversely with this number)
MAX_DISTANCE_BETWEEN_IDS=1000
#INTRON_URL='http://localhost:8090/solr/gigatron/select?q='
#SAMPLE_URL='http://localhost:8090/solr/sra_samples/select?q='

#setup headers for both the original intron list and the sample metadata list
INTRON_HEADER='snaptron_id	chromosome	start	end	length	strand	annotated?	left_motif	right_motif	left_annotated?	right_annotated?	samples	read_coverage_by_sample	samples_count	coverage_sum	coverage_avg	coverage_median	source_dataset_id'

SAMPLE_HEADER='intropolis_sample_id_i	run_accession_s	sample_accession_s	experiment_accession_s	study_accession_s	submission_accession_s	sra_ID_s	run_ID_s	run_alias_t	run_date_t	updated_date_t	spots_s	bases_s	run_center_t	experiment_name_t	run_attribute_t	experiment_ID_s	experiment_alias_t	experiment_title_t	study_name_t	sample_name_t	design_description_t	library_name_t	library_strategy_s	library_source_t	library_selection_t	library_layout_t	library_construction_protocol_t	read_spec_t	platform_t	instrument_model_t	platform_parameters_t	experiment_url_link_s	experiment_attribute_t	sample_ID_s	sample_alias_t	taxon_id_s	common_name_t	description_t	sample_url_link_s	sample_attribute_t	study_ID_s	study_alias_t	study_title_t	study_type_t	study_abstract_t	center_project_name_t	study_description_t	study_url_link_s	study_attribute_t	related_studies_t	primary_study_s	submission_ID_s	submission_comment_t	submission_center_t	submission_lab_t	submission_date_t	sradb_updated.x_s	fastq_ID_s	file_name_t	md5_s	bytes_s	audit_time_s	sradb_updated_file_s	date_download_t	URL_s	layout_s	URL_R_s	md5_R_s	cell_type_t	tissue_t	cell_line_t	strain_t	age_s	disease_t	population_t	sex_s	source_name_s'

INTRON_HEADER_FIELDS=INTRON_HEADER.split('\t')
INTRON_HEADER_FIELDS_MAP={}
for (i,field) in enumerate(INTRON_HEADER_FIELDS):
   INTRON_HEADER_FIELDS_MAP[field]=i 

SAMPLE_HEADER_FIELDS=SAMPLE_HEADER.split('\t')
SAMPLE_HEADER_FIELDS_MAP={}
for (i,field) in enumerate(SAMPLE_HEADER_FIELDS):
   SAMPLE_HEADER_FIELDS_MAP[field]=i 

#setup lucene reader for sample related searches
lucene.initVM()
analyzer = StandardAnalyzer(Version.LUCENE_4_10_1)
reader = IndexReader.open(SimpleFSDirectory(File("lucene/")))
searcher = IndexSearcher(reader)

def run_tabix(qargs,rquerys,tabix_db,filter_set=None,sample_set=None,filtering=False,print_header=True,debug=True):
    tabix_db = "%s/%s" % (TABIX_DB_PATH,tabix_db)
    if debug:
        sys.stderr.write("running %s %s %s\n" % (TABIX,tabix_db,qargs))
    if not filtering and print_header:
        sys.stdout.write("DataSource:Type\t%s\n" % (INTRON_HEADER))
    ids_found=set()
    tabixp = subprocess.Popen("%s %s %s | cut -f 2-" % (TABIX,tabix_db,qargs),stdout=subprocess.PIPE,shell=True)
    for line in tabixp.stdout:
        fields=[]
        #build either filter set or sample set or both
        if (sample_set != None and len(sample_set) > 0) or (filter_set != None and len(filter_set) > 0):
             fields=line.rstrip().split("\t")
             if filter_set != None and fields[INTRON_ID_COL] not in filter_set:
                 #sys.stderr.write("field %s not in filter_set\n" % (fields[INTRON_ID_COL]))
                 continue
             if sample_set != None:
                 sample_set.update(set(fields[SAMPLE_IDS_COL].split(",")))
        #filter return stream based on range queries (if any)
        if rquerys:
            if len(fields) == 0:
                fields=line.rstrip().split("\t")
            skip=False
            for rfield in rquerys.keys():
                (op,rval)=rquerys[rfield]
                if rfield not in INTRON_HEADER_FIELDS_MAP:
                    sys.stderr.write("bad field %s in range query,exiting\n" % (rfield))
                    sys.exit(-1)
                fidx = INTRON_HEADER_FIELDS_MAP[rfield]
                val = float(fields[fidx])
                if rfield not in FLOAT_FIELDS:
                    val = int(val)
                if not op(val,rval):
                    skip=True
                    break
            if skip:
                continue
        if filtering:
            ids_found.add(fields[INTRON_ID_COL])
            continue
        #now just stream back the result
        else:
            sys.stdout.write("%s:I\t%s" % (DATA_SOURCE,line))
    exitc=tabixp.wait() 
    if exitc != 0:
        raise RuntimeError("%s %s %s returned non-0 exit code\n" % (TABIX,TABIX_DB,args))
    if filtering:
        return ids_found

#based on the example code at
#http://graus.nu/blog/pylucene-4-0-in-60-seconds-tutorial/
def search_samples_lucene(sample_map,sampleq,sample_set,stream_sample_metadata=False):
    (fields,queries,booleans) = lucene_query_parse(sampleq)
    query = MultiFieldQueryParser.parse(Version.LUCENE_4_10_1, queries, fields, booleans, analyzer)
    #query = MultiFieldQueryParser.parse(Version.LUCENE_4_10_1, ['human AND adult AND brain'], ['description_t'], [BooleanClause.Occur.MUST], analyzer)
    hits = searcher.search(query, LUCENE_MAX_HITS)
    if DEBUG_MODE: 
        sys.stderr.write("Found %d document(s) that matched query '%s':\n" % (hits.totalHits, sampleq))
    if stream_sample_metadata:
        sys.stdout.write("DataSource:Type\t%s\n" % (SAMPLE_HEADER))
    for hit in hits.scoreDocs:
        doc = searcher.doc(hit.doc)
        sid = doc.get("intropolis_sample_id_i")
        #track the sample ids if asked to
        if sample_set != None:
            sample_set.add(sid)
        #stream back the full sample metadata record from the in-memory dictionary
        if stream_sample_metadata:
            sys.stdout.write("%s:S\t%s\n" % (DATA_SOURCE,sample_map[sid]))
            
def stream_samples(sample_set,sample_map):
    sys.stdout.write("DataSource:Type\t%s\n" % (SAMPLE_HEADER))
    for sample_id in sample_set:
        sys.stdout.write("%s:S\t%s\n" % (DATA_SOURCE,sample_map[sample_id]))

def lucene_query_parse(sampleq):
    queries_ = sampleq.split(SAMPLE_QUERY_DELIMITER)
    fields = []
    queries = []
    booleans = []
    for query_tuple in queries_:
        (field,query) = query_tuple.split(SAMPLE_QUERY_FIELD_DELIMITER)
        fields.append(field)
        query = query.replace('AND',' AND ')
        #sys.stderr.write("query + fields: %s %s\n" % (query,field))
        queries.append(query)
        booleans.append(BooleanClause.Occur.MUST)
    return (fields,queries,booleans)

#use to load samples metadata to be returned
#ONLY when the user requests by overlap coords OR
#coords and/or non-sample metadata thresholds (e.g. coverage)
#otherwise we'll return the whole thing from SOLR instead
def load_sample_metadata(file_):
    start = time.time()
    fmd=snaputil.load_cpickle_file("%s.pkl" % (file_))
    if fmd:
        end = time.time()
        taken = end-start
        #sys.stderr.write("time taken to load samples from pickle: %d\n" % taken)
        return fmd
    start = time.time()
    fmd={}
    #dont need the hash-on-column headers just yet
    with open(file_,"r") as f:
       for line in f:
           line = line.rstrip()
           fields=line.split("\t")
           fmd[fields[0]]=line
    end = time.time()
    taken = end-start
    #sys.stderr.write("time taken to load samples from normal: %d\n" % taken)
    snaputil.store_cpickle_file("%s.pkl" % (file_),fmd)
    return fmd

#do multiple searches by a set of ids
def search_introns_by_ids(snaptron_ids):
    sid_queries = []
    start_sid = 1    
    end_sid = 1
    #coalesce the snaptron_ids into ranges
    #to avoid making too many queries (n+1) to Tabix
    for sid in sorted(snaptron_ids):
        sid = int(sid)
        dist = abs(sid - start_sid)
        if dist <= MAX_DISTANCE_BETWEEN_IDS:
           end_sid = sid
        else:
           sid_queries.append("1:%d-%d" % (start_sid,end_sid))
           start_sid = sid
           end_sid = sid
    sid_queries.append("1:%d-%d" % (start_sid,end_sid))
    print_header = True
    for query in sid_queries:
        if DEBUG_MODE:
            sys.stderr.write("query %s\n" % (query))
        run_tabix(query,None,TABIX_DBS['snaptron_id'],filter_set=snaptron_ids,print_header=print_header,debug=DEBUG_MODE)
        print_header = False


comp_op_pattern=re.compile(r'([=><!]+)')
def range_query_parser(rangeq,snaptron_ids,rquery_will_be_index=False):
    rquery={}
    #snaptron_ids = []
    if rangeq is None or len(rangeq) < 1:
        return (None,None,rquery)
    fields = rangeq.split(',')
    (first_tdb,first_rquery)=(None,None)
    for field in fields:
        m=comp_op_pattern.search(field)
        (col,op_,val)=re.split(comp_op_pattern,field)
        if not m or not col or col not in TABIX_DBS:
            continue
        op=m.group(1)
        if op not in operators:
            sys.stderr.write("bad operator %s in range query,exiting\n" % (str(op)))
            sys.exit(-1)
        #queries by id are a different type of query, we simply record the id
        #and then move on, if there is only a id query that will be caught higher up
        if col == 'snaptron_id':
            snaptron_ids.update(set(val.split('-')))
            continue
        val=float(val)
        if col not in FLOAT_FIELDS:
            val=int(val)
        if first_tdb:
            rquery[col]=(operators[op],val)
            continue 
            #return (None,None,None,only_ids)
        #add first rquery to the rquery hash if we're not going to be
        #used as an index 
        #OR the case where it's floating point and we need to work around
        #Tabix's lack of support for that
        if not rquery_will_be_index or 'avg' in col or 'median' in col:
            rquery[col]=(operators[op],val)
        #if we are used for the index,
        #then for 2nd pass columns where the value could be 0 (GTEx)
        #we need to add another predicate to avoid the 0
        #since tabix doesn't handle 0's and will return them
        #even for a >=1 query
        #elif col[-1] == '2' and val > 0.0:
        #    rquery[col]=(operators['>'],0.0)
        #only do the following for the first range query
        tdb=TABIX_DBS[col]
        first_tdb=tdb
        extension=""
        #since tabix only takes integers, round to nearest integer
        val = int(round(val))
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
def intron_ids_from_samples(sample_set):
    start = time.time()
    sample2introns=snaputil.load_cpickle_file("./sample2introns.pkl")
    #print("setting up sample2intron map")
    if not sample2introns:
        sample2introns={}
        #f=open("/data2/gigatron2/all_SRA_introns_ids_stats.tsv.new","r")
        f=gzip.open("%s/%s" % (TABIX_DB_PATH,TABIX_INTERVAL_DB),"r")
        print("opened gzip file for introns")
        num_lines = 0
        for line in f:
            if "gigatron_id" in line or "snaptron_id" in line:
                continue
            fields=line.rstrip().split('\t')
            sample_ids=fields[SAMPLE_IDS_COL].split(',')
            #just map the intron id
            for sample_id in sample_ids:
                if sample_id not in sample2introns:
                    sample2introns[sample_id]=set()
                sample2introns[sample_id].add(int(fields[0]))
            num_lines+=1
            #if num_lines % 10000 == 0:
            #    print("loaded %d introns" % (num_lines))
        f.close()
    #print("about to write pkl file")
    snaputil.store_cpickle_file("./sample2introns.pkl",sample2introns)
    #print("pkl file written")
    end = time.time()
    taken=end-start
    if DEBUG_MODE:
        sys.stderr.write("loaded %d samples2introns in %d\n" % (len(sample2introns),taken))
    introns_seen=set()
    for sample_id in sample_set:
        introns_seen.update(sample2introns[sample_id])
    if DEBUG_MODE:
        sys.stderr.write("s2I\t%s\n" % (str(len(introns_seen))))
    return introns_seen
    
#cases:
#1) just interval (one function call)
#2) interval + range query(s) (one tabix function call + field filter(s))
#3) one or more range queries (one tabix range function call + field filter(s))
#4) interval + sample (2 function calls: 1 lucene for sample filter + 1 tabix using snaptron_id filter)
#5old) sample (1 lucene call -> use interval ids to return all intervals)
#5) sample (1 lucene call -> use snaptron_ids to do a by_ids search (multiple tabix calls))
#6) sample + range query(s) (2 function calls: 1 lucene for sample filter + 1 tabix using snaptron_id filter + field filter)
def main():
    input_ = sys.argv[1]
    DEBUG_MODE_=DEBUG_MODE
    if len(sys.argv) > 2:
        DEBUG_MODE_=True
    (intervalq,rangeq,sampleq) = input_.split('|')
    sample_map = load_sample_metadata(SAMPLE_MD_FILE)
    if DEBUG_MODE_:
        sys.stderr.write("loaded %d samples metadata\n" % (len(sample_map)))
    snaptron_ids = set()
    #if we have any sample related queries, do them first to get sample id set
    if len(sampleq) >= 1:
        #stream_solr(sampleq,sample_set=sample_set)
        sample_set = set()
        search_samples_lucene(sample_map,sampleq,sample_set,stream_sample_metadata=False)
        if DEBUG_MODE_:
            sys.stderr.write("found %d samples matching sample metadata fields/query\n" % (len(sample_set)))
        snaptron_ids = intron_ids_from_samples(sample_set)
        #if no snaptron_ids were found we're done in keeping with the strict AND policy (currently)
        if len(snaptron_ids) == 0:
            return
    #whether or not we use the interval query as a filter set or the whole query
    rquery_index = len(intervalq) < 1 and len(rangeq) >= 1
    #now parse the range querie(s) (this can include one by_id hack search which would be a list of one or more ids)
    (first_tdb,first_rquery,rquery) = range_query_parser(rangeq,snaptron_ids,rquery_will_be_index=rquery_index)
    #UPDATE if someone wants both by_ids and other queries we'll give them both by passing snaptron_ids later to the tabix method as a filter
    #otherwise we do by_id queries only
    #sample_set is used to track samples we pick up from introns queried to stream them later (if desired)
    sample_set = set()
    if len(snaptron_ids) > 0 and len(intervalq) == 0 and (len(rangeq) == 0 or not first_tdb):
        search_introns_by_ids(snaptron_ids)
    #back to usual processing, interval queries come first possibly with filters from the point range queries and/or ids
    elif len(intervalq) >= 1:
        run_tabix(intervalq,rquery,TABIX_INTERVAL_DB,filter_set=snaptron_ids,sample_set=sample_set,debug=DEBUG_MODE_)
    #if there's no interval query to use with tabix, use a point range query (first_rquery) with additional filters from the following point range queries and/or ids
    elif len(rangeq) >= 1:
        run_tabix(first_rquery,rquery,first_tdb,filter_set=snaptron_ids,sample_set=sample_set,debug=DEBUG_MODE_)

if __name__ == '__main__':
    main()

