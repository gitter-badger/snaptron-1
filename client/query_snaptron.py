#!/usr/bin/env python2.7
import sys
import urllib
import urllib2
import time
import json
import argparse

import clsnapconf
from SnaptronIteratorHTTP import SnaptronIteratorHTTP

#TODO use python tmp
TMPDIR='/tmp'

def parse_regions(args):
    regions = []
    if '/' in args.regions:
        with open(args.regions,"r") as fin:
            for line in fin:
                line = line.rstrip()
                regions.append(line)
    else:
        regions = args.regions.split('|')
    return regions

import csv
import re
def parse_query_params(args):
    queries = []
    with open(args.query_file,"r") as cfin:
        creader = csv.DictReader(cfin,['region','contains','thresholds','filters'],dialect=csv.excel_tab)
        for record in creader:
            query=[]
            #have to have a region
            #TODO format/name check here
            query.append("regions=%s" % record['region'])
            if len(record['contains']) > 0:
                query.append('contains=1')
            if record['thresholds'] is not None:
                thresholds = record['thresholds']
                thresholds = re.sub("=",":",thresholds)
                thresholds = thresholds.split('&')
                query.append("&".join(["rfilter=%s" % x for x in thresholds]))
            if record['filters'] is not None:
                filters = record['filters']
                filters = filters.split('&')
                query.append("&".join(["sfilter=%s" % x for x in filters]))
            queries.append("&".join(query))
    return queries


def main(args):
    query_params_per_region=parse_query_params(args)
    endpoint = 'snaptron'
    for query_param_string in query_params_per_region:
        sIT = SnaptronIteratorHTTP(query_param_string,endpoint)
        for record in sIT:
            sys.stdout.write("%s\n" % (record))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Snaptron command line client')

    parser.add_argument('--query-file', metavar='/path/to/file_with_queries', type=str, required=True, help='path to a file with one query per line where a query is one or more of a region (HUGO genename or genomic interval) optionally with one or more thresholds and/or filters specified and/or contained flag turned on')

    parser.add_argument('--groups', metavar='a=chr#:#-#,chr#:#-#;b=chr#:#-#;...', type=str, default='', help='comma separated list of junction groups defined by genomic intervals; records are prefixed with the group name on output')
    parser.add_argument('--function', metavar='jir', type=str, default='', help='function to compute between specified groups of junctions ranked across samples; currently only supports Junction Inclusion Ratio (JIR)')
    parser.add_argument('--tmpdir', metavar='/path/to/tmpdir', type=str, default=TMPDIR, help='path to temporary storage for downloading and manipulating junction and sample records')

    parser.add_argument('--regions', metavar='chr#:#-#|chr#:#-#|CD99|...', type=str, default='BRCA1', help='one or more genomic regions (gene names and/or genomic intervals) seperated by \'|\' (logical OR) or the name of a file with a list of one or more regions (include one or more \'/\'s to denote a file path)')
    parser.add_argument('--thresholds', metavar='samples_count>=#[&|]coverage_avg=#[&|]...', type=str, default='samples_count>=1000&coverage_avg>100', help='one or more cutoffs for filtering returned junctions; separator is either \'&\' (logical AND) or \'|\' (logical OR)')
    parser.add_argument('--filters', metavar='design_description:cortex[&|]library_layout=PAIRED[&|]...', type=str, default=None, help='one or more filter predicates on the metadata of the samples from which the junctions are derived; separator is either \'&\' (logical AND) or \'|\' (logical OR)')
    parser.add_argument('--contained', action='store_const', const=True, default=False,
                        help='only return junctions whose coordinate span are within one or more of the regions specified')
    #returned format (UCSC, and/or subselection of fields) option?
    #intersection or union of intervals option?

    main(parser.parse_args())
