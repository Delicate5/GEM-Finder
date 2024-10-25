from datetime import datetime
import os

output_d 	= '../data/GWAS_Shrinked/'
LD_block_d 	= '../data/LD_block_merged_0.8/'

def associations_shrinker(input_file = '../data/GWAS_Catalog/associations_all.txt'):
	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] association shrinker start")
	with open(input_file, 'r') as fp:
		shrinked = []
		header = next(fp).strip().split('\t')
		for i in fp:
			line = i.strip().split('\t')
			for j in range(len(line)): line[j] = line[j].replace(' ,',',').replace(', ', ',').replace(' - ', '-').replace('\x3b', ';').replace('; ', ';').replace(' ', '_').replace('\xa0', '_').replace('/','_')
			GWAS_Catalog, snp_ID, context, pvalue = line[7], line[21], line[24], line[27]
			chr_ID, chr_Pos, reported_gene, mapped_gene, upstream_gene, downstream_gene = line[11:17]
			if(chr_ID):
				if(';' in snp_ID):
					if(len(snp_ID.split(';')) == len(chr_ID.split(';')) == len(chr_Pos.split(';')) == len(context.split(';'))):
						separator = [[snps, ids, poss, cons] for snps, ids, poss, cons in zip(snp_ID.split(';'), chr_ID.split(';'), chr_Pos.split(';'), context.split(';'))]
						for eachsnp in separator: shrinked.append([eachsnp[1], eachsnp[2], GWAS_Catalog, eachsnp[0], pvalue, reported_gene, mapped_gene, upstream_gene, downstream_gene, eachsnp[3]+'_pgs8070-manualseparation'])
				elif('_x_' in snp_ID):
					if (('_x_' in chr_ID) or ('_x_' in chr_Pos)):
						separator = [[snps, ids, poss] for snps, ids, poss in zip(snp_ID.split('_x_'), chr_ID.split('_x_'), chr_Pos.split('_x_'))]
						shrinked.append([separator[0][1], separator[0][2], GWAS_Catalog, separator[0][0], pvalue, reported_gene, mapped_gene, upstream_gene, downstream_gene, context+'_pgs8070-'+separator[1][0]])
						shrinked.append([separator[1][1], separator[1][2], GWAS_Catalog, separator[1][0], pvalue, reported_gene, mapped_gene, upstream_gene, downstream_gene, context+'_pgs8070-'+separator[0][0]])
				else: shrinked.append([chr_ID, chr_Pos, GWAS_Catalog, snp_ID, pvalue, reported_gene, mapped_gene, upstream_gene, downstream_gene, context])

	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] write shrinked as "+output_d+'association_shrinked.txt')
	with open(output_d+'association_shrinked.txt', 'w') as wp:
		wp.write("#"+'\t'.join(['chr_ID', 'chr_Pos', 'GWAS_Catalog', 'snp_ID', 'p-value', 'reported_gene', 'mapped_gene', 'upstream_gene', 'downstream_gene', 'context'])+'\n')
		for i in shrinked: wp.write('\t'.join(i)+'\n')

	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] write trait ID as "+output_d+"association_list.txt")
	with open(output_d+"association_list.txt", 'w') as wp: wp.write(os.popen("awk 'NR > 1 {print $3}' "+output_d+"association_shrinked.txt | sort | uniq -c | awk '{$1=$1; print}' | tr ' ' '\t'").read())
	
	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] association shrinker finish")

def super_phenotype_extractor(input_file = '../data/GWAS_Shrinked/association_list.txt'):
	all_list = []
	with open(input_file, 'r') as fp:
		header = next(fp).strip().split('\t')
		for i in fp: all_list.append(i.strip().split('\t')[1])
	super_list = []
	for i in all_list:
		checker = 1
		for j in all_list:
			if(j in i and j != i):
				checker = 0
				break
		if(checker): super_list.append(i)
	with open('../data/GWAS_Shrinked/association_super_list.txt', 'w') as wp: wp.write('#GWAS_catalog\n' + '\n'.join(super_list)+'\n')

def well_known_phenotype_extractor(cutoff = 5e-8, candidates = '../data/GWAS_Shrinked/association_list.txt', traits = '../data/GWAS_Shrinked/association_shrinked.txt'):
	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] get phenotypes which have good p-values")
	with open(candidates, 'r') as fp:
		targets = {}
		header = next(fp).strip().split('\t')
		for i in fp:
			line = i.strip().split('\t')
			if(int(line[0]) >= 100): targets[line[1]] = 0

	with open(traits, 'r') as fp:
		header = next(fp).strip().split('\t')
		for i in fp:
			line = i.strip().split('\t')
			if(float(line[4]) <= cutoff):
				for j in targets.keys():
					if j == line[2]:
						targets[j] += 1
						break
	with open('../data/GWAS_Shrinked/study_well_known.txt', 'w') as wp:
		wp.write('#GWAS_Catalog\tnumber_of_association_over_'+str(cutoff)+'\n')
		for k, v in targets.items():
			if(v >= 100) : wp.write(k+'\t'+str(v)+'\n')

	# with open(traits, 'r') as fp:
	# 	header = next(fp).strip().split('\t')
	# 	for i in fp:
	# 		line = i.strip().split('\t')
	# 		if(float(line[4]) <= cutoff):
	# 			for j in targets.keys():
	# 				if j in line[2]: targets[j] += 1
	# finals = {}
	# for i in targets.keys():
	# 	checker = 1
	# 	for j in targets.keys():
	# 		if(j in i and i != j): checker = 0
	# 	if checker: finals[i] = targets[i]
	# with open('../data/GWAS_Shrinked/well_known_phenotypes.txt', 'w') as wp:
	# 	wp.write('#GWAS_Catalog\tnumber_of_association_over_'+str(cutoff)+'\n')
	# 	for k, v in finals.items():
	# 		if(v >= 100) : wp.write(k+'\t'+str(v)+'\n')			

def LD_block_bulider(chrlist, ldlist = ['AFR','AMR','EAS','EUR','SAS'], input_directory = '/home/ajl1213/genome.info/LD/LD_naoki/'):
	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] Build LD block")
	for chrID in chrlist:
		print(datetime.now().strftime('%H:%M:%S'), "[data_handler] Start " + chrID)
		LD_block, LD_count = {}, {}
		for ldID in ldlist:
			print(datetime.now().strftime('%H:%M:%S'), "[data_handler] " + ldID)
			with open(input_directory+chrID+'_'+ldID+'_ld_all.hap.ld', 'r') as fp:
				next(fp)
				for i in fp:
					line = i.split()
					pos1, pos2, r_sqr = line[1], line[2], line[4]
					if (float(r_sqr) > 0.8):
						interaction = pos1+'-'+pos2
						if interaction in LD_count: LD_count[interaction] += 1
						else:
							LD_block[pos1] = []
							LD_block[pos2] = []
							LD_count[interaction] = 1
		
		for k, v in LD_count.items():
			pos1, pos2 	= k.split('-')
			count 		= str(v)
			LD_block[pos1].append(pos2+':'+count)
			LD_block[pos2].append(pos1+':'+count)

		with open(LD_block_d+chrID+'.txt', 'w') as wp:	
			for k, v in LD_block.items(): wp.write(k+'\t'+','.join(v)+'\n')
	print(datetime.now().strftime('%H:%M:%S'), "[data_handler] Done")

def stage_annotated_bedmaker(input_directory = '../data/Project_EC_CRC/EC_peaks/', stages = ['Early', 'Mid', 'Late', 'Full']):
	with open(input_directory+'Merged_EC.bed', 'w') as wp:
		wp.write('\t'.join(['#chr', 'start_pos', 'end_pos', 'stage'])+'\n')
		for stage in stages:
			with open(input_directory+stage+'_EC.bed', 'r') as fp:
				header = next(fp).strip().split('\t')
				for i in fp:
					line = i.strip().split('\t')
					wp.write('\t'.join(line[:3] + [stage])+'\n')
