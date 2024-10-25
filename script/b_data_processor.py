from datetime import datetime
import os, fileinput

output_d 	= '../output/data_processor/'
LD_block_d	= '/home/ajl1213/genome.info/LD/LD_naoki/'

def associations_extractor(name, query, assembly = 'hg38', strict = True):
	print(datetime.now().strftime('%H:%M:%S'), "[data_processor] extract target traits")
	input_file = '../data/GWAS_Shrinked/association_shrinked.txt' if assembly == 'hg38' else '../data/GWAS_hg19_liftover/association_shrinked_hg19.txt'
	if query == None: query = [name]
	targets = []
	with open(output_d+name+'.tsv', 'w') as wp:
		with open(input_file, 'r') as fp:
			wp.write(next(fp))
			if(strict):
				for i in fp:
					line = i.strip().split('\t')
					for trait in query:
						if(trait.lower() == line[2].lower()):
							wp.write(i)
							break
			else:
				for i in fp:
					line = i.strip().split('\t')
					for trait in query:
						if(trait.lower() in line[2].lower()):
							if not (line[2] in targets) : targets.append(line[2])
							wp.write(i)
							break
	if(not strict):
		with open(output_d+name+'_traits.txt', 'w') as wp: wp.write('#target_traits_name\n'+'\n'.join(targets)+'\n')

def Ld_expander(chrlist, name, cutoff = 1):
	print(datetime.now().strftime('%H:%M:%S'), "[data_processor] expand SNP with LD information")
	input_directory = '../data/LD_block_merged_0.8/'

	print(datetime.now().strftime('%H:%M:%S'), "[data_processor] Load SNP information")
	with open(output_d+name+'.tsv', 'r') as fp:
			header = next(fp).strip().split('\t')
			targets = {}
			for chrID in chrlist: targets[chrID] = {}
			for i in fp:
				line = i.strip().split('\t')
				if(line[0]) == 'Y' : continue
				chrID, pos, snpID = 'chr'+line[0], line[1], line[3]
				targets[chrID][pos] = snpID

	print(datetime.now().strftime('%H:%M:%S'), "[data_processor] Get LD snp")
	final_result = {}
	for k, v in targets.items():
		final_result[k] = {}
		for snp_pos, snpID in v.items():
			if(snp_pos in final_result[k]): final_result[k][snp_pos].append(snpID)
			else: final_result[k][snp_pos] = [snpID]
		with open(input_directory+k+'.txt', 'r') as fp:
			for i in fp:
				pos, interacts = i.strip().split('\t')
				for snp_pos, snpID in v.items():
					if pos == snp_pos:
						for interact in interacts.split(','):
							ld_pos, count = interact.split(':')
							if int(count) > cutoff:
								if(ld_pos in final_result[k]): final_result[k][ld_pos].append(snpID)
								else: final_result[k][ld_pos] = [snpID]
	
	print(datetime.now().strftime('%H:%M:%S'), "[data_processor] Save SNP")
	for k, v in final_result.items(): final_result[k] = sorted(final_result[k].items(), key=lambda item:int(item[0]))
	with open(output_d+name+'_LDSnp.bed', 'w') as wp:
		wp.write("#"+'\t'.join(['chrID', 'start_pos', 'end_pos','tagSNP_rsID'])+'\n')
		for k, v in final_result.items():
			for pos, ID_list in v: wp.write('\t'.join([k, pos, str(int(pos)+1), ','.join(ID_list)])+'\n')
