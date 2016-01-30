#!/usr/bin/env python2.7
import sys
import os
import unittest

import snapconf
import snaputil
import snaptron

#set of test interval queries
IQs=['chr1:10160-10161']
#holds the set of intropolis ids for each specific query for the original SRA set of inropolis junctions
EXPECTED_IIDS={IQs[0]:set(['0','1','2','3','4','5','6','9','10','11','13','14','15','16','17','18','20','21','22','23','24','26','27','28','29','30','31','32','33','34','35','36','37','38'])}

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
   
    def tc(self,interval_query,filter_set=None,sample_set=None,filtering=False):
        '''wrap the normal run_tabix call to hardcode defaults for interval querying'''
        range_query = None
        return tc(interval_query,range_query,tdbs['chromosome'],filter_set=filter_set,sample_set=sample_set,filtering=filtering)
 
    def test_basic_interval(self):
        '''make sure we're getting back an expected set of intropolis ids'''
        chosen_query = 0
        #get intropolis ids
        iids = self.tc(IQs[chosen_query], filtering=True)
        self.assertEqual(iids, EXPECTED_IIDS[IQs[chosen_query]])

if __name__ == '__main__':
    unittest.main()
