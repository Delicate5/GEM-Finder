#!/home/ajl2/anaconda3/bin/python



import os
import numpy



DIR=os.getcwd()
bin_path='/home/ajl1213/Projects/rasgrp3/data/hic/2_juicer/bin'
gtf='/home/ajl1213/genome.info/gencode/hg19.release38/GTF/gencode.v38.ptncoding.lncRNA.level_1_2.txt'
peak_file='/home/ajl1213/Projects/Endothelial/TEC_data/ChIP/Mapping/Peakcall/H3K27ac.mergepeaks.bed'
diffpeak_file='/home/ajl1213/Projects/Endothelial/TEC_data/ChIP/Processing/DiffPeaks/H3K27ac.DiffPeakLabel.txt'
rna_file='/home/ajl1213/Projects/Endothelial/TEC_data/RNA/Processing/CountMatrix/RNA.TPM.qq.stage.txt'


resolution=[20000, '20kb']

sample_list=[
'hESC',
'ME',
'EC4',
'EC12',
'EC24',
'EC48',
'nonEC12',
'nonEC48'
]

chr_list=[
'chr1','chr2','chr3','chr4', 'chr5','chr6','chr7','chr8', 'chr9','chr10',
'chr11','chr12','chr13','chr14','chr15','chr16','chr17', 'chr18','chr19','chr20',
'chr21','chr22','chrX'
]

os.system('mkdir '+DIR+'/bin')
os.system('mkdir '+DIR+'/anno_res')



def get_prom_bins():
        ##  
        prom_dict={}
        input1=open(gtf,'r')
        all_input1=input1.readlines()
        for line in all_input1[1:]:
                each=line.strip().split('\t')
                ensemble_id=each[0]
                gene=each[1]
                chrom=each[2]
                tss=int(each[3])
                tss_adj=int(tss/resolution[0])*resolution[0]
                key_id=chrom+':'+str(tss_adj)

                if key_id in prom_dict:
                        tmp=prom_dict[key_id]
                        tmp.append(gene)
                        prom_dict[key_id]=tmp
                else:
                        prom_dict[key_id]=[gene]
        input1.close()

        ## promoters in each bin
        output1=open(DIR+'/bin/'+resolution[1]+'.prom_list.allchr.bin','w')
        for chrom in chr_list:

                ##
                bin_list=[]

                input1=open(bin_path+'/'+chrom+'.'+resolution[1]+'.bin','r')
                all_input1=input1.readlines()
                for line in all_input1:
                        each=line.strip().split('\t')

                        bin_id=each[0]+':'+each[1]+'-'+each[2]
                        bin_list.append(bin_id)

                input1.close()

                ##
                for bin_id in bin_list:
                        key_id=bin_id.split('-')[0]

                        if key_id in prom_dict:
                                prom_list=';'.join(prom_dict[key_id])
                        else:
                                prom_list='none'

                        output1.write(bin_id+'\t'+str(prom_list)+'\n')

        output1.close()



def get_peak_bins():
        ##  
        peak_dict={}
        input1=open(peak_file,'r')
        all_input1=input1.readlines()
        for line in all_input1:
                each=line.strip().split('\t')

		chrom=each[0]
		pt1=each[1]
		pt2=each[2]
		
		mid_pt=int((int(pt1)+int(pt2))/2)
		mid_pt_adj=int(mid_pt/resolution[0])*resolution[0]
		peak=chrom+':'+pt1+'-'+pt2

		key_id=chrom+':'+str(mid_pt_adj)

                if key_id in peak_dict:
                        tmp=peak_dict[key_id]
                        tmp.append(peak)
                        peak_dict[key_id]=tmp
                else:
                        peak_dict[key_id]=[peak]
        input1.close()

        ## peaks in each bin
        output1=open(DIR+'/bin/'+resolution[1]+'.peak_list.allchr.bin','w')
        for chrom in chr_list:

                ##
                bin_list=[]

                input1=open(bin_path+'/'+chrom+'.'+resolution[1]+'.bin','r')
                all_input1=input1.readlines()
                for line in all_input1:
                        each=line.strip().split('\t')

                        bin_id=each[0]+':'+each[1]+'-'+each[2]
                        bin_list.append(bin_id)

                input1.close()

                ##
                for bin_id in bin_list:
                        key_id=bin_id.split('-')[0]

                        if key_id in peak_dict:
                                peak_list=';'.join(peak_dict[key_id])
                        else:
                                peak_list='none'

                        output1.write(bin_id+'\t'+str(peak_list)+'\n')

        output1.close()


def anno_inter():
	##
	peak_dict={}
	input1=open(DIR+'/bin/'+resolution[1]+'.peak_list.allchr.bin','r')
	all_input1=input1.readlines()
	for line in all_input1:
		each=line.strip().split('\t')
		chrom=each[0].split(':')[0]
		pt1=each[0].split(':')[1].split('-')[0]
		pt2=each[0].split(':')[1].split('-')[1]

		bin_id=chrom+'.'+pt1+'.'+pt2

		peak_list=each[1]

		peak_dict[bin_id]=peak_list
	input1.close()

	##
	prom_dict={}
	input1=open(DIR+'/bin/'+resolution[1]+'.prom_list.allchr.bin','r')
	all_input1=input1.readlines()
	for line in all_input1:
		each=line.strip().split('\t')
		chrom=each[0].split(':')[0]
		pt1=each[0].split(':')[1].split('-')[0]
		pt2=each[0].split(':')[1].split('-')[1]

		bin_id=chrom+'.'+pt1+'.'+pt2

		prom_list=each[1]

		prom_dict[bin_id]=prom_list
	input1.close()


	##
	for sample in sample_list:
		print(sample)

		input1=open(DIR+'/res_fithic/'+sample+'.'+resolution[1]+'.sig_inter.filtered.txt','r')
		all_input1=input1.readlines()

		output1=open(DIR+'/res_fithic/'+sample+'.'+resolution[1]+'.sig_inter.filtered.anno.txt','w')
		output1.write('inter_id\tinter_dir\tfrag1peak\tpeak_list\tfrag2prom\tprom_list\tdist\tinter_type\tinter_class\tqval\n')
		for line in all_input1[1:]:
			each=line.strip().split('\t')

			frag1=each[0]
			frag2=each[1]
			inter_id=frag1+'-'+frag2

			dist=each[2]
			qval=each[3]
			inter_type=each[4]

			# forward
			inter_dir='f'
			frag1peak=peak_dict[frag1]
			frag2prom=prom_dict[frag2]

			if (frag1peak!='none' and frag2prom!='none'):
				inter_class='peak_prom'
			else:
				if (frag1peak!='none' and frag2prom=='none'):
					inter_class='peak_none'
				else:
					inter_class='none_none'

			new_line=[inter_id, inter_dir, frag1, frag1peak, frag2, frag2prom, dist, inter_type, inter_class, qval]
			output1.write('\t'.join(new_line)+'\n')
	
			# reverse
			inter_dir='r'
			frag2peak=peak_dict[frag2]
			frag1prom=prom_dict[frag1]

			if (frag2peak!='none' and frag1prom!='none'):
				inter_class='peak_prom'
			else:
				if (frag2peak!='none' and frag1prom=='none'):
					inter_class='peak_none'
				else:
					inter_class='none_none'

			new_line=[inter_id, inter_dir, frag2, frag2peak, frag1, frag1prom, dist, inter_type, inter_class, qval]	
			output1.write('\t'.join(new_line)+'\n')

		input1.close()


	## anno merged siginter list
	input1=open(DIR+'/res_fithic/siginter.'+resolution[1]+'.merged.txt','r')
	all_input1=input1.readlines()

	output1=open(DIR+'/res_fithic/siginter.'+resolution[1]+'.merged.anno.txt','w')
	output1.write('inter_id\tinter_dir\tfrag1peak\tpeak_list\tfrag2prom\tprom_list\tdist\tinter_type\tinter_class\n')
	for line in all_input1[1:]:
		each=line.strip().split('\t')

		frag1=each[0]
		frag2=each[1]
		inter_id=frag1+'-'+frag2

		dist=each[2]
		inter_type=each[4]

		# forward
		inter_dir='f'
		frag1peak=peak_dict[frag1]
		frag2prom=prom_dict[frag2]

		if (frag1peak!='none' and frag2prom!='none'):
			inter_class='peak_prom'
		else:
			if (frag1peak!='none' and frag2prom=='none'):
				inter_class='peak_none'
			else:
				inter_class='none_none'

		new_line=[inter_id, inter_dir, frag1, frag1peak, frag2, frag2prom, dist, inter_type, inter_class]
		output1.write('\t'.join(new_line)+'\n')

		# reverse
		inter_dir='r'
		frag2peak=peak_dict[frag2]
		frag1prom=prom_dict[frag1]

		if (frag2peak!='none' and frag1prom!='none'):
			inter_class='peak_prom'
		else:
			if (frag2peak!='none' and frag1prom=='none'):
				inter_class='peak_none'
			else:
				inter_class='none_none'

		new_line=[inter_id, inter_dir, frag2, frag2peak, frag1, frag1prom, dist, inter_type, inter_class]	
		output1.write('\t'.join(new_line)+'\n')

	input1.close()


def attach_peak_label():
	#
	peak_dict={}
	input1=open(diffpeak_file,'r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
		each=line.strip().split('\t')

		peak=each[0]
		peak_label=each[2]

		peak_dict[peak]=peak_label
	input1.close()

	#
	output1=open(DIR+'/anno_res/sig_inter.peak_label.txt','w')
	output1.write('sample\tinter_id\tinter_type\tinter_class\tpeak\tpeak_label\n')
	for sample in sample_list:
		print(sample)
		input1=open(DIR+'/res_fithic/'+sample+'.'+resolution[1]+'.sig_inter.filtered.anno.txt','r')
		all_input1=input1.readlines()
		for line in all_input1[1:]:
			each=line.strip().split('\t')
			inter_id=each[0]
			peak_list=each[3]
			inter_type=each[7]
			inter_class=each[8]
			if inter_class in ['peak_prom','peak_none']:

				for peak in peak_list.split(';'):
					if peak_dict.has_key(peak):
						peak_label=peak_dict[peak]
						new_line=[sample, inter_id, inter_type, inter_class, peak, peak_label]
						output1.write('\t'.join(new_line)+'\n')
		input1.close()
	output1.close()

	#
	output1=open(DIR+'/anno_res/sig_inter.merged.peak_label.txt','w')
	output1.write('inter_id\tinter_type\tinter_class\tpeak\tpeak_label\n')

	input1=open(DIR+'/res_fithic/siginter.20kb.merged.anno.txt','r')
	all_input1=input1.readlines()
	for line in all_input1[1:]:
		each=line.strip().split('\t')
		inter_id=each[0]
		peak_list=each[3]
		inter_type=each[7]
		inter_class=each[8]
		if inter_class in ['peak_prom','peak_none']:

			for peak in peak_list.split(';'):
				if peak_dict.has_key(peak):
					peak_label=peak_dict[peak]
					new_line=[inter_id, inter_type, inter_class, peak, peak_label]
					output1.write('\t'.join(new_line)+'\n')
	input1.close()



def get_inter_count():

	#
	rna_dict={}
	input1=open(rna_file,'r')
	all_input1=input1.readlines()
	rna_samples=all_input1[0].strip().split('\t')[2:]
	for line in all_input1[1:]:
		each=line.strip().split('\t')
		gene=each[1]
		vals=each[2:]

		idx=0
		tmp_dict={}
		for val in vals:
			log2tpm=numpy.log2(float(val)+1)
			tmp_dict[rna_samples[idx]]=log2tpm
			idx+=1

		rna_dict[gene]=tmp_dict
	input1.close()

	#
	output1=open(DIR+'/anno_res/sig_inter_count.per_gene.txt','w')
	output1.write('sample\tgene\tprom_siginter\tlog2tpm\n')
	for sample in sample_list:
		input1=open(DIR+'/res_fithic/'+sample+'.'+resolution[1]+'.sig_inter.filtered.anno.txt','r')
		all_input1=input1.readlines()
		genecount_dict={}
		for line in all_input1[1:]:
			each=line.strip().split('\t')
			inter_id=each[0]
			prom_list=each[5]
			inter_class=each[8]
			if inter_class=='peak_prom':

				for gene in prom_list.split(';'):
					if genecount_dict.has_key(gene):
						genecount_dict[gene]+=1
					else:
						genecount_dict[gene]=1

		#
		for gene in rna_dict:

			if sample=='EC4':
				adj_sample='EC04'
			else:
				adj_sample=sample

			if genecount_dict.has_key(gene):
				genecount=genecount_dict[gene]
			else:
				genecount=0

			log2tpm=rna_dict[gene][adj_sample]

			new_line=[sample, gene, str(genecount), str(log2tpm)]
			output1.write('\t'.join(new_line)+'\n')

	
		input1.close()
	output1.close()




get_prom_bins()
get_peak_bins()
anno_inter()
attach_peak_label()
get_inter_count()


