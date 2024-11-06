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
It works in function-level.

1. Annotating
   code for CRE annotation
2. Annotating
   code for DEG TSS annotation
3. Interaction filtering
   code for DEG TSS containing extracting
4. CRE set extracting
   code for CRE forming

