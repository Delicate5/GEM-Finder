# GEM-Finder
Framework briefing. Also contain permutation test scripts, and other preocessing scirpts.

Note
----
Annotating Hi-C interaction was processed under python2, and handling them was processed under python3. We will modify python2 scripts. 

script accessibility
--------------------
download the scripts. Overall directory structure should be fixed when you use GEM-Finder.

requirement-interaction annotation
----------------------------------
1) merged CRE sets
```txt
chr1    713306  714971
chr1    742993  743667
chr1    761466  763351
chr1    839785  840980
chr1    851417  852087
chr1    852433  852827
```
4) GTF files for TSS annotation
```txt
ensembleID      geneID  chrID   tss     tes     strand  geneType
ENSG00000186092.7       OR4F5   chr1    65419   71585   +       protein_coding
ENSG00000235249.1       OR4F29  chr1    367640  368634  +       protein_coding
ENSG00000284662.1       OR4F16  chr1    622053  621059  -       protein_coding
ENSG00000187634.13      SAMD11  chr1    859303  879955  +       protein_coding
ENSG00000188976.11      NOC2L   chr1    894689  879583  -       protein_coding
```
3) Deduplicated&Sorted Hi-C bam files of individual targets
```txt
[ID]  115     chr1    10030   27      100M    chr4    99705   0       CTAACCCTAACCCTAACCCTACCCCTACCCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAAC    FDBFED>=9FFFF=FG8EGF9D?F>EFGDFEGF5AGFEG<FB8GFEEFBDGGFFDGAFBF=CFGGBFAF<DGAG>FGB7A@FBFE@F:FGEF?9DFEEE;  XA:Z:chr1,-10362,26M1D6M1D47M2D21M,5;   MD:Z:21A5A5A66  PG:Z:MarkDuplicates     NM:i:3  AS:i:85 XS:i:73
[ID]   163     chr1    10034   29      100M    =       10059   25      CCCTAACCCTAACCCTACCCCTACCCCTACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCAACCCTAACCCTAACCCTAACCCTA    GGFEEFFGFEFGGGFBGFGFEGGFGGFFCEFFG@EFCFE?FGFGG:FGBDG9DFFGG5BGCDE<@GECB6AE>GE3E4FE=<EEDD8<CFGG@<FFFE4A  XA:Z:chr1,+10366,22M1D6M1D47M2D25M,6;   MD:Z:17A5A5A70  PG:Z:MarkDuplicates     NM:i:3  AS:i:85 XS:i:72
```
Overall workflow-interaction annotation
---------------------------------------
1. Hi-C proecssing
Bam file - bam2inter.py -> interaction file
```txt
chr1    chr1.0.20000    10000   0
chr1    chr1.20000.40000        30000   0
chr1    chr1.40000.60000        50000   0
chr1    chr1.60000.80000        70000   0
chr1    chr1.80000.100000       90000   0
```
interaction hit
```txt
chr3    175850000       chr3    176150000       13
chr3    136490000       chr3    136970000       1
chr3    184350000       chr3    184870000       4
chr3    55990000        chr3    56890000        1
chr3    15850000        chr3    16230000        24
```
interaction file - fithic -> significant interaction file
significant interaction file
```txt
chr1    fragmentMid1    chr2    fragmentMid2    contactCount    p-value q-value bias1   bias2
chr9    111350000       chr9    111470000       29      3.938140e-02    2.585804e-01    1.000000e+00    1.000000e+00
chr9    38670000        chr9    39570000        1       9.363949e-01    1.000000e+00    1.000000e+00    1.000000e+00
chr9    109170000       chr9    109990000       4       4.131927e-01    1.000000e+00    1.000000e+00    1.000000e+00
chr9    89950000        chr9    90390000        13      1.979620e-02    1.550282e-01    1.000000e+00    1.000000e+00
chr9    115490000       chr9    115750000       11      5.531412e-01    1.000000e+00    1.000000e+00    1.000000e+00
```
reannotate & rebuild & sort significant interaction file
```txt
frag1chr        frag1mid        frag2chr        frag2mid        contact_count   p_value q_value bias1   bias2
chr1    29650000        chr1    30430000        37      1.964721e-12    3.781693e-11    1.000000e+00    1.000000e+00
chr1    32650000        chr1    32750000        96      5.690006e-04    3.575127e-03    1.000000e+00    1.000000e+00
chr1    60530000        chr1    61230000        22      9.931711e-04    5.921797e-03    1.000000e+00    1.000000e+00
chr1    199930000       chr1    200370000       44      2.775619e-07    3.042722e-06    1.000000e+00    1.000000e+00
chr1    10310000        chr1    10490000        73      1.894558e-05    1.571043e-04    1.000000e+00    1.000000e+00
```
+ filter-out static interaction
2. CRE & TSS processing
CRE annotation
```txt
chr1:700000-720000	chr1:713306-714971
chr1:720000-740000	none
chr1:740000-760000	chr1:742993-743667
chr1:760000-780000	chr1:761466-763351
chr1:780000-800000	none
```
TSS annotation
```txt
chr1:0-20000	none
chr1:20000-40000	none
chr1:40000-60000	none
chr1:60000-80000	OR4F5
chr1:80000-100000	none
chr1:100000-120000	none
```
3. Intergrate information
```txt
inter_id        inter_dir       frag1peak       peak_list       frag2prom       prom_list       dist    inter_type      inter_class     qval
chr4.190480000.190500000-chr4.190500000.190520000       f       chr4.190480000.190500000        none    chr4.190500000.190520000        none    20000   static  none_none       0.0
chr4.190480000.190500000-chr4.190500000.190520000       r       chr4.190500000.190520000        none    chr4.190480000.190500000        none    20000   static  none_none       0.0
chr1.16920000.16940000-chr1.16940000.16960000   f       chr1.16920000.16940000  chr1:16939224-16940627  chr1.16940000.16960000  NBPF1   20000   static  peak_prom       0.0
chr1.16920000.16940000-chr1.16940000.16960000   r       chr1.16940000.16960000  chr1:16946332-16947616;chr1:16950571-16950811;chr1:16952274-16952712;chr1:16956525-16957941     chr1.16920000.16940000  none    20000   static  peak_none       0.0
chr1.121260000.121280000-chr1.121480000.121500000       f       chr1.121260000.121280000        chr1:121260615-121261829        chr1.121480000.121500000        none    220000  static  peak_none       0.0
```

requirement-DEG linked CRE selection
------------------------------------
Overall workflow-DEG linked CRE selection
-----------------------------------------

References
----------
1. Arya Kaul, Sourya Bhattacharyya & Ferhat Ay 2020. "Identifying statistically significant chromatin contacts from Hi-C data with FitHiC2." Nature Protocols. 15:991-1012, 2020. doi: 10.1038/s41596-019-0273-0.
