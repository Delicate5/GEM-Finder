#!/home/pgs8070/miniconda3/bin/python
from datetime import datetime
import os, sys, time

import a_data_handler as dh
import a_assembly_liftover as al
import b_data_processor as dp
import b_pvalue_maker as pm
import a_GEM_finder as gf

def main(argv):
	for output_d in ['data_processor', 'pvalue_maker']: path_checker(output_d)
	chrlist = sorted(['chr'+str(i+1) for i in range(22)] + ['chrX'])
	ldlist	= ['AFR', 'AMR', 'EAS', 'EUR', 'SAS']

	print(datetime.now().strftime('%H:%M:%S'), "[main] Start Gwas analyzer")

	#GEM-finder_analysis_controlmaking
	# query = []
	# with open('../data/GWAS_Shrinked/study_well_known.txt', 'r') as fp:
	# 	header = next(fp).strip().split('\t')
	# 	for i in fp: query.append(i.strip().split('\t')[0])
	# pm.p_value_maker(chrlist, query, argv[1], 2000, argv[2], argv[3])
	# pm.p_value_merger(query, 40, 'Early', 'ChIP_CRE')
	# pm.p_value_merger(query, 40, 'Mid', 'ChIP_CRE')
	# pm.p_value_merger(query, 40, 'Late', 'ChIP_CRE')
	# pm.p_value_merger(query, 40, 'Full', 'ChIP_CRE')
	print(datetime.now().strftime('%H:%M:%S'), "[main] Finish Gwas analyzer")

def path_checker(output_d):
	if not os.path.exists('../output/'+output_d):
		print(datetime.now().strftime('%H:%M:%S'), "[main] cannot find an output directory for "+output_d+", create it as ../output/"+output_d)
		os.makedirs('../output/'+output_d)

def simple_ldexpander(chrlist, name, query):
	print(datetime.now().strftime('%H:%M:%S'), "[main] Expanding with",name)
	dp.associations_extractor(name, query, assembly = 'hg19')
	dp.Ld_expander(chrlist, name)

main(sys.argv)
