#!/usr/bin/env python2.7
import sys
import os
import urllib
import urllib2
import argparse
import gzip
import csv
import re

import clsnapconf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP
from SnaptronIteratorLocal import SnaptronIteratorLocal

GTEX_TISSUE_COL=65

fmap = {'thresholds':'rfilter','filters':'sfilter','region':'regions'}
breakpoint_patt = re.compile(r'(^[^:]+-[^:]+$)|(^COSF\d+$)|(^\d+$)')
def parse_query_argument(record, fieldnames, groups, header=True):
    endpoint = 'snaptron'
    query=[]
    for field in fieldnames:
        if len(record[field]) > 0:
            if field == 'thresholds' or field == 'filters':
                predicates = re.sub("=",":",record[field])
                predicates = predicates.split('&')
                query.append("&".join(["%s=%s" % (fmap[field],x) for x in predicates]))
            elif field == 'group':
                groups.append(record[field])
            else:
                mapped_field = field
                if field in fmap:
                    mapped_field = fmap[field]
                query.append("%s=%s" % (mapped_field,record[field]))
            if field == 'region' and breakpoint_patt.search(record[field]) is not None:
                endpoint = 'breakpoint'
    if not header:
        query.append("header=0")
    return (query,endpoint)


def parse_command_line_args(args):
    fieldnames = []
    endpoint = 'snaptron'
    for field in clsnapconf.FIELD_ARGS.keys():
        if field in vars(args) and vars(args)[field] is not None:
            fieldnames.append(field)
    groups = []
    (query,endpoint) = parse_query_argument(vars(args), fieldnames, groups, header=args.function is None and not args.noheader)
    return (["&".join(query)], groups, endpoint)


def parse_query_params(args):
    if args.query_file is None:
        return parse_command_line_args(args)
    endpoint = 'snaptron'
    queries = []
    groups = []
    with open(args.query_file,"r") as cfin:
        creader = csv.DictReader(cfin,dialect=csv.excel_tab)
        get_header = args.function is None and not args.noheader
        for (i,record) in enumerate(creader):
            (query, endpoint) = parse_query_argument(record, creader.fieldnames, groups, header=get_header)
            queries.append("&".join(query))
            get_header = False
    #assume the endpoint will be the same for all lines in the file
    return (queries,groups,endpoint)

def calc_jir(a, b):
    numer = b - a
    denom = a + b + 1
    return numer/float(denom)


def junction_inclusion_ratio_bp(args, results, group_list, sample_records):
    sample_stats = results['samples']
    (group_a_g1, group_a_g2, group_b_g1, group_b_g2) = group_list
    group_a = group_a_g1[:-2]
    group_b = group_b_g1[:-2]
    g1 = group_a_g1[-1:]
    g2 = group_a_g2[-1:]
    
    sample_scores = {}
    for sample in sample_stats:
        if group_a_g1 not in sample_stats[sample]:
            sample_stats[sample][group_a_g1]=0
        if group_b_g1 not in sample_stats[sample]:
            sample_stats[sample][group_b_g1]=0
        if group_a_g2 not in sample_stats[sample]:
            sample_stats[sample][group_a_g2]=0
        if group_b_g2 not in sample_stats[sample]:
            sample_stats[sample][group_b_g2]=0

        sample_stats[sample][g1] = [calc_jir(sample_stats[sample][group_a_g1], sample_stats[sample][group_b_g1]), sample_stats[sample][group_a_g1], sample_stats[sample][group_b_g1]]
        sample_stats[sample][g2] = [calc_jir(sample_stats[sample][group_a_g2], sample_stats[sample][group_b_g2]), sample_stats[sample][group_a_g2], sample_stats[sample][group_b_g2]]
        sample_scores[sample] = min(sample_stats[sample][g1][0], sample_stats[sample][g2][0])

    missing_sample_ids = set()
    counter = 0
    if not args.noheader:
        sys.stdout.write("analysis_score\t%s jir/a/b raw counts\t%s jir/a/b raw counts\tsample metadata\n" % (g1,g2))
    for sample in sorted(sample_scores.keys(),key=sample_scores.__getitem__,reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        score = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        sample_record = sample_records[sample]
        sys.stdout.write("%s\t%s\t%s\t%s\n" % (str(score),":".join([str(x) for x in sample_stats[sample][g1]]),":".join([str(x) for x in sample_stats[sample][g2]]),sample_record))


def junction_inclusion_ratio(args, results, group_list, sample_records):
    sample_stats = results['samples']
    group_a = group_list[0]
    group_b = group_list[1]
    sample_scores = {}
    for sample in sample_stats:
        if group_a not in sample_stats[sample]:
            sample_stats[sample][group_a]=0
        if group_b not in sample_stats[sample]:
            sample_stats[sample][group_b]=0
        #numer = sample_stats[sample][group_b] - sample_stats[sample][group_a]
        #denom = sample_stats[sample][group_b] + sample_stats[sample][group_a] + 1
        #sample_scores[sample]=numer/float(denom)
        sample_scores[sample] = calc_jir(sample_stats[sample][group_a], sample_stats[sample][group_b])

    missing_sample_ids = set()
    counter = 0
    if not args.noheader:
        sys.stdout.write("analysis_score\t%s raw count\t%s raw count\tsample metadata\n" % (group_a,group_b))
    for sample in sorted(sample_scores.keys(),key=sample_scores.__getitem__,reverse=True):
        counter += 1
        if args.limit > -1 and counter > args.limit:
            break
        score = sample_scores[sample]
        if sample not in sample_records:
            missing_sample_ids.add(sample)
            continue
        sample_record = sample_records[sample]
        sys.stdout.write("%s\t%d\t%d\t%s\n" % (str(score),sample_stats[sample][group_a],sample_stats[sample][group_b],sample_record))

def track_exons(args, results, record, group):
    exons = results['exons']
    fields = record.split('\t')
    snid = fields[clsnapconf.INTRON_ID_COL]
    for (type_,col) in {'start':clsnapconf.INTERVAL_START_COL,'end':clsnapconf.INTERVAL_END_COL}.iteritems():
        coord = int(fields[col])
        if coord not in exons:
            exons[coord]={type_:set()}
        if type_ not in exons[coord]:
            exons[coord][type_]=set()
        exons[coord][type_].add(snid)

def filter_exons(args, results, group_list, sample_records):
    #only used if filtering by length range
    (rlen1,rlen2) = (None,None)
    if args.exon_length is not None:
        (rlen1,rlen2) = map(lambda x: int(x), args.exon_length.split('-'))
    #filter the joint list of both intron start and end coordinates
    coords = sorted(results['exons'].keys())
    end = None
    end_ids = None
    sys.stdout.write("Type\tLeft End Snaptron IDs\tRight End Snaptron IDs\tstart\tend\tlength\n")
    for (i,coord) in enumerate(coords):
        fields = results['exons'][coord]
        #end here is intron end (exon start)
        if 'end' in fields:
            end = coord
            end_ids = fields['end']
            #look for joining exon ends after this entry
            for coord2 in coords[i+1:]:
                fields2 = results['exons'][coord2]
                #start here is intron start (exon end)
                if 'start' in fields2:
                    start_ids = fields2['start']
                    dist = coord2 - 1 - end
                    if dist >= clsnapconf.MIN_EXON_SIZE and (rlen1 is None or (dist >= int(rlen1) and dist <= int(rlen2))):
                        sys.stdout.write("exon\t%s\t%s\t%d\t%d\t%d\n" % (",".join(end_ids),",".join(start_ids),end+1,coord2-1,dist))

TISSUE_SPECIFICITY_FUNC='ts'
SHARED_SAMPLE_COUNT_FUNC='shared'
JIR_FUNC='jir'
TRACK_EXONS_FUNC='exon'
FUNCTION_TO_TYPE={TRACK_EXONS_FUNC:'not-shared', JIR_FUNC:'not-shared', None:None, TISSUE_SPECIFICITY_FUNC:'shared',SHARED_SAMPLE_COUNT_FUNC:'shared'}
def count_samples_per_group(args, results, record, group):
    sample_stats = results['samples']
    fields = record.split('\t')
    samples = fields[clsnapconf.SAMPLE_IDS_COL].split(',')
    sample_covs = fields[clsnapconf.SAMPLE_IDS_COL+1].split(',')
    start_value = 0
    if args.function == TISSUE_SPECIFICITY_FUNC:
        start_value = sys.maxint
    #track samples shared across the flanking junctions here
    if 'shared' in results and group not in results['shared']:
        results['shared'][group] = set()
    for (i,sample_id) in enumerate(samples):
        #this can happen with GTEx
        if int(sample_covs[i]) == 0:
            continue
        #if we're doing tissue spec. then make sure we get shared samples, otherwise skip
        if FUNCTION_TO_TYPE[args.function] == 'shared' and group in results['groups_seen']:
            #haven't seen this sample before, so must not be shared
            if sample_id not in sample_stats or group not in sample_stats[sample_id]:
                continue
            else:
                results['shared'][group].add(sample_id)
        if sample_id not in sample_stats:
            sample_stats[sample_id]={}
        if group not in sample_stats[sample_id]:
            sample_stats[sample_id][group]=start_value
        if FUNCTION_TO_TYPE[args.function] != TISSUE_SPECIFICITY_FUNC:
            sample_stats[sample_id][group]+=int(sample_covs[i])
        else:
            #initially we used coverage of just those shared samples, but now we just do present or not
            #sample_stats[sample_id][group]=min(sample_stats[sample_id][group],int(sample_covs[i]))
            sample_stats[sample_id][group]=1

def tissue_specificity(args, results, group_list, sample_records):
    sample_stats = results['samples']
    sys.stdout.write("group\tsample_id\tshared\ttissue\n")
    for group in group_list:
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            sys.stderr.write("No shared samples between splice junctions for %s\n" % (group))
        #for sample_id in results['shared'][group]:
        for sample_id in sample_records.keys():
            present = 0
            #if sample_id in sample_stats and group in sample_stats[sample_id]:
            if sample_id in results['shared'][group]:
                present = sample_stats[sample_id][group]
            sfields = sample_records[sample_id].split("\t")
            tissue = sfields[GTEX_TISSUE_COL]
            sys.stdout.write("%s\t%s\t%d\t%s\n" % (group, sample_id, present, tissue))

def report_shared_sample_counts(args, results, group_list, sample_records):
    sys.stdout.write("group\tshared_sample_counts\n")
    for group in group_list:
        if group not in results['shared'] or len(results['shared'][group]) == 0:
            sys.stderr.write("No shared samples between splice junctions for %s, skipping\n" % (group))
            #sys.stdout.write("%s\t0\n" % (group))
            continue
        count = len(results['shared'][group])
        sys.stdout.write("%s\t%d\n" % (group, count))

def download_sample_metadata(args):
    sample_records = {}
    gfout = None
    if clsnapconf.CACHE_SAMPLE_METADTA:
        parts = clsnapconf.SERVICE_URL.split('/')
        cache_file = os.path.join(args.tmpdir,"snaptron_sample_metadata_cache.%s.tsv.gz" % parts[-1])
        if os.path.exists(cache_file):
            with gzip.open(cache_file,"r") as gfin:
                for line in gfin:
                    line = line.rstrip()
                    fields = line.split('\t')
                    sample_records[fields[0]]=line
            del sample_records['']
            return sample_records
        else:
            gfout = gzip.open(cache_file,"w")
    response = urllib2.urlopen("%s/samples?all=1" % (clsnapconf.SERVICE_URL))
    all_records = response.read()
    all_records = all_records.split('\n')
    for line in all_records:
        fields = line.split('\t')
        sample_records[fields[0]]=line
        if gfout is not None:
            gfout.write("%s\n" % (line))
    if gfout is not None:
        gfout.close()
    del sample_records['']
    return sample_records

iterator_map = {True:SnaptronIteratorLocal, False:SnaptronIteratorHTTP}
def process_queries(args, query_params_per_region, groups, endpoint, function=None, local=False):
    results = {'samples':{},'queries':[],'exons':{'start':{},'end':{}}}
    if FUNCTION_TO_TYPE[args.function] == 'shared':
        results['groups_seen']=set()
        results['shared']={}
    first = True
    for (group_idx, query_param_string) in enumerate(query_params_per_region):
        sIT = iterator_map[local](query_param_string, endpoint)
        group = None
        if group_idx < len(groups):
            group = groups[group_idx]
        #assume we get a header in this case and don't count it against the args.limit
        counter = -1
        if args.noheader:
            counter = 0
        for record in sIT:
            if function is not None:
                function(args, results, record, group)
            elif endpoint == 'breakpoint':
                (region, contains, group_) = record.split('\t')
                if region == 'region':
                    continue
                (query, _) = parse_query_argument({'region':region,'contains':contains,'group':group_}, ['region', 'contains', 'group'], groups, header=False)
                results['queries'].append("&".join(query))
            else:
                counter += 1
                if args.limit > -1 and counter > args.limit:
                    break
                group_label = ''
                if len(groups) > 0 and len(groups[group_idx]) > 0:
                    group_label = "%s\t" % (groups[group_idx])
                    if not args.noheader and first and counter == 0:
                        group_label = 'group\t'
                sys.stdout.write("%s%s\n" % (group_label, record))
        if group is not None and 'groups_seen' in results:
            results['groups_seen'].add(group)
        first = False
    return results


compute_functions={JIR_FUNC:(count_samples_per_group,junction_inclusion_ratio),'jirbp':(count_samples_per_group,junction_inclusion_ratio_bp),TRACK_EXONS_FUNC:(track_exons,filter_exons),TISSUE_SPECIFICITY_FUNC:(count_samples_per_group,tissue_specificity),SHARED_SAMPLE_COUNT_FUNC:(count_samples_per_group,report_shared_sample_counts),None:(None,None)}
def main(args):
    #parse original set of queries
    (query_params_per_region, groups, endpoint) = parse_query_params(args)
    #get original functions (if passed in)
    (count_function, summary_function) = compute_functions[args.function]
    #process original queries
    results = process_queries(args, query_params_per_region, groups, endpoint, function=count_function, local=args.local)
    #print "results length %d" % len(results['samples'])
    #we have to do a double process if doing a breakpoint query since we get the coordinates
    #in the first query round and then query them in the second (here)
    if endpoint == 'breakpoint':
        #update original None functions with the JIR
        (count_function, summary_function) = compute_functions['jirbp']
        #now process the coordinates that came back from the first breakpoint query using the normal query + JIR
        results = process_queries(args, results['queries'], groups, 'snaptron', function=count_function, local=args.local)
    #if either the user wanted the JIR to start with on some coordinate groups OR they asked for a breakpoint, do the JIR now
    if args.function or endpoint == 'breakpoint':
        sample_records = {}
        if args.function != TRACK_EXONS_FUNC:
            sample_records = download_sample_metadata(args)
        group_list = set()
        map(lambda x: group_list.add(x), groups)
        group_list = sorted(group_list)
        summary_function(args,results,group_list,sample_records)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Snaptron command line client')
    for (field,settings) in clsnapconf.FIELD_ARGS.iteritems():
        parser.add_argument("--%s" % field, metavar=settings[0], type=settings[1], default=settings[2], help=settings[3])

    parser.add_argument('--query-file', metavar='/path/to/file_with_queries', type=str, default=None, help='path to a file with one query per line where a query is one or more of a region (HUGO genename or genomic interval) optionally with one or more thresholds and/or filters specified and/or contained flag turned on')

    parser.add_argument('--function', metavar='jir', type=str, default=None, help='function to compute between specified groups of junctions ranked across samples; currently only supports Junction Inclusion Ratio "jir" and exon finding "exon"')

    parser.add_argument('--tmpdir', metavar='/path/to/tmpdir', type=str, default=clsnapconf.TMPDIR, help='path to temporary storage for downloading and manipulating junction and sample records')
    
    parser.add_argument('--limit', metavar='1', type=int, default=-1, help='# of records to return, defaults to all (-1)')

    parser.add_argument('--exon-length', metavar='50-60', type=str, default=None, help='length range of exons to look for within queried region when function is set to "exon", defaults to None (print out all exons in queried region when function is "exon")')
    
    
    parser.add_argument('--local', action='store_const', const=True, default=False, help='if running Snaptron modeules locally (skipping WSI)')
    parser.add_argument('--noheader', action='store_const', const=True, default=False, help='turn off printing header in ouput')
    

    #returned format (UCSC, and/or subselection of fields) option?
    #intersection or union of intervals option?

    args = parser.parse_args()
    if args.region is None and args.thresholds is None and args.filters is None and args.query_file is None:
        sys.stderr.write("Error: no discernible arguments passed in, exiting\n")
        parser.print_help()
        sys.exit(-1)
    if args.function == TISSUE_SPECIFICITY_FUNC and 'gtex' not in clsnapconf.SERVICE_URL:
        sys.stderr.write("Error: attempting to do tissue specificity (ts) on a non-GTEx Snaptron instance. Please change the SERVICE_URL setting in clsnapconf.py file to be the GTEx Snaptron instance before running this function; exiting\n")
        sys.exit(-1)
    main(args)
