#!/usr/bin/env python2.7
import sys
import os
import unittest

import snapconf
import snaputil
import snaptron

#set of test interval queries
IQs=['chr1:10160-10161']
RQs=['1:100000-100000','1:5-5']
#holds the set of intropolis ids for each specific query for the original SRA set of inropolis junctions
EXPECTED_IIDS={
               IQs[0]:set(['0','1','2','3','4','5','6','9','10','11','13','14','15','16','17','18','20','21','22','23','24','26','27','28','29','30','31','32','33','34','35','36','37','38']),
               RQs[0]:set(["1042183","4655249","8228366","8228479","8613428","8813669","8992620","9191837","9191838","12439314","16533047","16585013","17895417","19061596","19061599","19659394","20360595","20360638","20360704","20360763","20360865","23543324","23821668","25640534","25640535","26783371","28153140","28463313","31138498","31576775","32888166","32888167","35125664","37187404","37876025","37876029","38244820","38418223","42053224"])
              }

def setUpModule():
    pass

def tearDownModule():
    pass

#shortcuts for snaptron methods used in tests
tc = snaptron.run_tabix
rqp = snaptron.range_query_parser
sbi = snaptron.search_introns_by_ids

tdbs = snapconf.TABIX_DBS

class TestTabixCalls(unittest.TestCase):
    '''
    check tabix for basic queries (both interval and range including ids)
    def run_tabix(qargs,rquerys,tabix_db,filter_set=None,sample_set=None,filtering=False,print_header=True,debug=True):
    returns an id set if filtering==True
    can also populate sample_set if defined
    '''

    def setUp(self):
        pass
   
    def itc(self,interval_query,filter_set=None,sample_set=None,filtering=False):
        '''wrap the normal run_tabix call to hardcode defaults for interval querying'''
        range_query = None
        return tc(interval_query,range_query,tdbs['chromosome'],filter_set=filter_set,sample_set=sample_set,filtering=filtering)
 
    def test_basic_interval(self):
        '''make sure we're getting back an expected set of intropolis ids'''
        chosen_query = 0
        #get intropolis ids
        iids = self.itc(IQs[chosen_query], filtering=True)
        self.assertEqual(iids, EXPECTED_IIDS[IQs[chosen_query]])
    
    def ltc(self,range_query,filter_set=None,sample_set=None,filtering=False):
        '''wrap the normal run_tabix call to hardcode defaults for LENGTH range querying'''
        extra_ranges = None
        return tc(range_query,extra_ranges,tdbs['length'],filter_set=filter_set,sample_set=sample_set,filtering=filtering)
    
    def test_basic_length(self):
        '''make sure we're getting back an expected set of intropolis ids'''
        chosen_query = 0
        #get intropolis ids
        iids = self.ltc(RQs[chosen_query], filtering=True)
        self.assertEqual(iids, EXPECTED_IIDS[RQs[chosen_query]])

if __name__ == '__main__':
    unittest.main()
