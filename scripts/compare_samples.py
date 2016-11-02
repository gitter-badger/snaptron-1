#!/usr/bin/env python2.7
import sys
import snaputil as su
from bitarray import bitarray
MAX_ID=81063662


#ba1 = su.load_cpickle_file("./sample_ids1_1.pkl", compressed=True)
#ba2 = su.load_cpickle_file("./sample_ids1_2.pkl", compressed=True)

def orthem(ba1,ba2):
    ba_final = ba1 | ba2
    return ba_final

def setthem(ba_final):
    i = 0
    s1=set()
    [s1.add(i) for (i,x) in enumerate(ba_final) if x]
    #for bit in ba_final:
    #    if bit:
    #        s1.add(i)
    #    i+=1
    return s1

ba_final = su.load_cpickle_file("/data3/snaptron/sample_ids/0.pkl.gz", compressed=True)
for i in xrange(1,131):
    ba2 = su.load_cpickle_file("/data3/snaptron/sample_ids/%s.pkl.gz" % str(i), compressed=True)
    if ba2 != None:
        ba_final = orthem(ba2,ba_final)

#ba_final = orthem(ba1,ba2)
s1 = setthem(ba_final)
print len(s1)
#sys.stdout.write(",".join([str(x) for x in s1]))

