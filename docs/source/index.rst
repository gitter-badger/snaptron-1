.. Snaptron documentation master file, created by
   sphinx-quickstart on Mon Apr 25 16:03:47 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

===================
Snaptron User Guide
===================

.. toctree::
   :maxdepth: 2

Quickstart
----------

First, we will present an example query and then break it down to allow the impatient users to get on with their research and skip the longer explanation of the details: ::

  curl "http://stingray.cs.jhu.edu:8090/srav1/snaptron?region=chr6:1-514015&rfilter=samples_count:100"

The above command uses cURL to query the Snaptron web service for all junctions that overlap the coordinate range of ``1-514015`` on chromosome 6 and that have 1 or more reads coverage in exactly 100 samples (for CGI parsing reasons the .:. is used instead of .=. as a range operator).  The return format is a TAB delimited text stream of junction records, one per line including a header as the first line to explain the columns returned.

Gene symbols (exact HUGO gene symbols) can also be used instead of chromosome coordinates: ::

  curl "http://stingray.cs.jhu.edu:8090/srav1/snaptron?regions=CD99&rfilter=samples_count:20"

A Snaptron query is a set of predicates logically AND'ed together from three different query types, each with their own format. Table 1 displays the four different query types.  The three main query types are: region, range over a summary statistics column (.range.), a freetext/field search of the associated sample metadata (.metadata.), and an exact ID retrieval for both exon-exon junctions (Snaptron assigned IDs) and samples (Intropolis-assigned IDs).

A snaptron query may contain only one of the three types of queries or may contain all three, or some combination of two types.  In the example above the region and range query types are present as ``chr6:1-514015`` for the region type and ``samples_count:100`` for the range type.

Table 1. Query Types
--------------------

=============== ================================================================ ============ =============================================== ==================
Query Type      Description                                                      Multiplicity Format                                          Example
--------------- ---------------------------------------------------------------- ------------ ----------------------------------------------- ------------------
Region          chromosome based coordinates range (1-based); HUGO gene name     0-1          chr(1-22,X,Y,M):1-size of chromosome; gene_name chr21:1-500; CD99
Range           range over summary statistic column values                       0 or more    column_name(>:,<:,:)number (integer or float)   coverage_avg>:10
Sample Metadata is-equal-to/contains text (keywords) search over sample metadata 0 or more    fieldname:keyword                               description:cortex
Snaptron IDs    one or more snaptron_ids                                         0 or more    snaptron_id=\d+[,\d+]*                          snaptron_id=5,7,8
=============== ================================================================ ============ =============================================== ==================

Table 2 shows the queryable fields and their query type.  Often the Range query type columns will be used as a way to reduce the number of false positive junctions.  This can be done easily with the two columns: samples_count and coverage_sum.  For example some suggested values from our own research are presented in Table 3.

Table 2. Query Fields
---------------------

================================== ======================= =========================================================== =======================
Field                              Query Type              Range of Values                                             Example
---------------------------------- ----------------------- ----------------------------------------------------------- -----------------------
coordinate* (chromosome:start-end) Interval                chr(1-22,X,Y,M):1-size of chromosome                        chr1:4-100
gene symbol*                       string exact match      HUGO gene symbols                                           CD99
length                             Range                   1-600K                                                      intron_length<:5000
samples_count                      Range                   1-Inf                                                       samples_count>:5  
coverage_sum**                     Range                   1-Inf                                                       coverage_sum>:10
coverage_avg**                     Range                   1.0-Inf                                                     coverage_avg>:5.0
coverage_median**                  Range                   1.0-Inf                                                     coverage_median>:6.0
snaptron_id                        Set                     unique, stable snaptron specific IDs, one per intron record snaptron_id=5
================================== ======================= =========================================================== =======================

[TODO: add sample metadata fields]

\*you can either pass a coordinate string or a gene symbol in the interval query segment, but not both
\*\*for the GTEx dataset append a "2" to the end of this column to get query off the 2nd pass alignment read coverage statistics

Exact HUGO gene symbols can be searched in snaptron SRA instance en lieu of actual coordinates.   If the gene symbol had multiple coordinate ranges that were either on different chromosomes or more than 10,000 bases apart, Snaptron will do multiple tabix lookups and will stream them back in coordinate order per chromosome (the chromosome order itself is sorted via default python sorted which is not 1-22,X,Y,M).

The gene symbol to coordinate mapping is provided from the UCSC RefSeq .Flat. dataset: http://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/refFlat.txt.gz

That dataset maps HUGO gene symbols to RefSeq gene IDs and transcript coordinates.

Table 3.  Suggested Quality Threshold for Selected Range Columns
----------------------------------------------------------------

==============  ======================  ===================
Selected Field  Quality Threshold Type  Threshold Predicate
--------------  ----------------------  -------------------
samples_count   baseline                >:5
samples_count   higher confidence       >:1000
coverage_sum    baseline                >:10?
coverage_sum    higher confidence       >:50?
intron_length   baseline                <:10000?
intron_length   higher confidence       <:3000?
==============  ======================  ===================

The return format is a TAB-delimited series of fields where each line represents a unique intron call.  Table 4 displays the complete list of fields in the return format of the Snaptron web service.  Fields marked with an "*" are queryable as they are specifically indexed.  The ``chromosome``, ``start``, and, ``end`` fields are a special case where the index is a combination of all three of them together.

Table 4. Complete list of Snaptron Fields In Return Format
----------------------------------------------------------

=========== ======================= ======================================= ============================================================================================================== =======================
Field Index Field Name              Type                                    Description                                                                                                    Example
----------- ----------------------- --------------------------------------- -------------------------------------------------------------------------------------------------------------- -----------------------
0           DataSource:Type         Abbrev:Single Character                 Differentiates between a return line of type Intron (I), Sample (S), or Gene (G).                              SRAv1:I
1*          snaptron_id             Integer                                 stable, unique ID for Snaptron junctions                                                                       5
2           chromosome              String                                  Reference ID for genomics coordinates                                                                          chr7
3           start                   Integer                                 beginning (left) coordinate of intron                                                                          10113
4           end                     Integer                                 last (right) coordinate of intron                                                                              10244
5           length                  Integer                                 Length of intron coordinate span                                                                               132
6           strand                  Single Character                        Orientation of intron (Watson or Crick)                                                                        "+" or "-"
7**         annotated?              Boolean Integer                         If both ends of the intron are annotated as *a* splice site in some annotation                                 1 or 0
8           left_motif              String                                  Splice site sequence bases at the left end of the intron                                                       GT
9           right_motif             String                                  Splice site sequence bases at the right end of the intron                                                      AG
10**        left_annotated?         String                                  If the left end splice site is annotated or not and which annotations it appears in (maybe more than once)     aC19,cG19,cG38:1;0
11**        right_annotated?        String                                  If the right end splice site is in an annotated or not, same asleft_annotated?                                 aC19,cG19,cG38:1;0
12          samples                 Comma separated list of Integers IDs    The list of samples which had one or more reads covering the intron(?). IDs are from the IntropolisDB.         5,10,14
13          read_coverage_by_sample Comma separated list of Integers        Coverage of the intron per sample (matches "samples" column position)                                          1,6,20
16*         coverage_avg            Float                                   Average coverage across all samples which had at least 1 read covering the intron in the first pass alignment  8.667
17*         coverage_median         Float                                   Median coverage across all samples which had at least 1 read covering the intron in the first pass alignment   6
18          source_dataset_id       Integer                                 Snaptron ID for the original dataset used (SRA, GTEx, TCGA)                                                    SRAv1=0,SRAv2=1,GTEx=2
=========== ======================= ======================================= ============================================================================================================== =======================

\*\*These fields are not present in the GTEx version of the Snaptron webservice at this time.  They are not queryable in the SRA version.


* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
