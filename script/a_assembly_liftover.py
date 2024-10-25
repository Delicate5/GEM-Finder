from datetime import datetime
import os

output_d = '../data/GWAS_hg19_liftover/'

def liftover_input_format_maker(input_file = '../data/GWAS_Shrinked/association_shrinked.txt'):
	print(datetime.now().strftime('%H:%M:%S'), "[assembly_liftover] make bed file for use web liftover")
	with open(output_d+'association_shrinked_liftover_input.bed', 'w') as wp:
		with open(input_file, 'r') as fp:
			header = next(fp).strip().split('\t')
			for i in fp:
				line = i.strip().split('\t')
				wp.write('chr'+line[0]+' '+line[1]+' '+line[1]+' '+line[3]+'\n')
	print(datetime.now().strftime('%H:%M:%S'), "[assembly_liftover] liftover input prepared")

def liftover_converter(input_file = '../data/GWAS_Shrinked/association_shrinked.txt'):
	print(datetime.now().strftime('%H:%M:%S'), "[assembly_liftover] liftover shrinked association")
	with open(output_d+'association_shrinked_hg19.txt', 'w') as wp:
		with open(input_file, 'r') as fp1:
			wp.write(next(fp1))
			with open(output_d+'association_shrinked_liftover_output.bed', 'r') as fp2:
				for i in fp2:
					line = i.strip().split()
					oline = next(fp1).strip().split('\t')
					while(line[3] != oline[3]): oline = next(fp1).strip().split('\t')
					wp.write('\t'.join([oline[0], line[1]] + oline[2:])+'\n')
	print(datetime.now().strftime('%H:%M:%S'), "[assembly_liftover] save lifted data as "+output_d+'association_shrinked_hg19.txt')


