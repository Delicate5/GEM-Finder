# GEM-Finder
Framework briefing. Also contain permutation test scripts.

Note
----
Annotating Hi-C interaction was processed under python2, and handling them was processed under python3. We will modify python2 scripts. 

script accessibility
--------------------
download the scripts. Overall directory structure should be fixed when you use GEM-Finder.

requirement
-----------
1) DEG for individual targets
2) CREs for individual targets or merged set
3) Hi-C contact information for individual targets

processing
----------
1. Hi-C proecssing
Bam file - bam2inter.py -> interaction file
interaction map

   chr1    chr1.0.20000    10000   0
   chr1    chr1.20000.40000        30000   0
   chr1    chr1.40000.60000        50000   0
   chr1    chr1.60000.80000        70000   0
   chr1    chr1.80000.100000       90000   0

interaction hit

   chr3    175850000       chr3    176150000       13
   chr3    136490000       chr3    136970000       1
   chr3    184350000       chr3    184870000       4
   chr3    55990000        chr3    56890000        1
   chr3    15850000        chr3    16230000        24

interaction file - fithic -> significant interaction file

significant interaction file

   chr1    fragmentMid1    chr2    fragmentMid2    contactCount    p-value q-value bias1   bias2
   chr9    111350000       chr9    111470000       29      3.938140e-02    2.585804e-01    1.000000e+00    1.000000e+00
   chr9    38670000        chr9    39570000        1       9.363949e-01    1.000000e+00    1.000000e+00    1.000000e+00
   chr9    109170000       chr9    109990000       4       4.131927e-01    1.000000e+00    1.000000e+00    1.000000e+00
   chr9    89950000        chr9    90390000        13      1.979620e-02    1.550282e-01    1.000000e+00    1.000000e+00
   chr9    115490000       chr9    115750000       11      5.531412e-01    1.000000e+00    1.000000e+00    1.000000e+00
   
reannotate & rebuild & sort significant interaction file

   frag1chr        frag1mid        frag2chr        frag2mid        contact_count   p_value q_value bias1   bias2
   chr1    29650000        chr1    30430000        37      1.964721e-12    3.781693e-11    1.000000e+00    1.000000e+00
   chr1    32650000        chr1    32750000        96      5.690006e-04    3.575127e-03    1.000000e+00    1.000000e+00
   chr1    60530000        chr1    61230000        22      9.931711e-04    5.921797e-03    1.000000e+00    1.000000e+00
   chr1    199930000       chr1    200370000       44      2.775619e-07    3.042722e-06    1.000000e+00    1.000000e+00
   chr1    10310000        chr1    10490000        73      1.894558e-05    1.571043e-04    1.000000e+00    1.000000e+00

+ filter-out static interaction

2. CRE & TSS processing
CRE annotation

   chr1:700000-720000	chr1:713306-714971
   chr1:720000-740000	none
   chr1:740000-760000	chr1:742993-743667
   chr1:760000-780000	chr1:761466-763351
   chr1:780000-800000	none

TSS annotation

   chr1:0-20000	none
   chr1:20000-40000	none
   chr1:40000-60000	none
   chr1:60000-80000	OR4F5
   chr1:80000-100000	none
   chr1:100000-120000	none

3. Intergrate information
   inter_id        inter_dir       frag1peak       peak_list       frag2prom       prom_list       dist    inter_type      inter_class     qval
   chr4.190480000.190500000-chr4.190500000.190520000       f       chr4.190480000.190500000        none    chr4.190500000.190520000        none    20000   static  none_none       0.0
   chr4.190480000.190500000-chr4.190500000.190520000       r       chr4.190500000.190520000        none    chr4.190480000.190500000        none    20000   static  none_none       0.0
   chr1.16920000.16940000-chr1.16940000.16960000   f       chr1.16920000.16940000  chr1:16939224-16940627  chr1.16940000.16960000  NBPF1   20000   static  peak_prom       0.0
   chr1.16920000.16940000-chr1.16940000.16960000   r       chr1.16940000.16960000  chr1:16946332-16947616;chr1:16950571-16950811;chr1:16952274-16952712;chr1:16956525-16957941     chr1.16920000.16940000  none    20000   static  peak_none       0.0
   chr1.121260000.121280000-chr1.121480000.121500000       f       chr1.121260000.121280000        chr1:121260615-121261829        chr1.121480000.121500000        none    220000  static  peak_none       0.0

It works in function-level.

1. Annotating
   code for CRE annotation
2. Annotating
   code for DEG TSS annotation
3. Interaction filtering
   code for DEG TSS containing extracting
4. CRE set extracting
   code for CRE forming

