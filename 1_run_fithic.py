#!/home/ajl1213/anaconda2/bin/python


import os
import gzip
import operator


DIR=os.getcwd()
fithic_path='/home/ajl1213/programs/fithic/fithic'
python_path='/home/ajl1213/anaconda2/bin'
bin_path='/home/ajl1213/Projects/rasgrp3/data/hic/2_juicer/bin'
bam_path='/home/ajl1213/Projects/rasgrp3/data/hic/1_mapping/SortedNodupBam'

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

chrom_list=[
'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6',  'chr7', 'chr8', 'chr9',
'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16','chr17', 'chr18', 'chr19', 'chr20',
'chr21', 'chr22', 'chrX'
]


min_dist=5000
max_dist=1000000
resolution='20kb'

min_qval=0.01

os.system('mkdir '+DIR+'/feature_vec')
os.system('mkdir '+DIR+'/res_fithic')


def run_fithic():
	ncpu=4

	for sample in sample_list:
		for chrom in chrom_list:
			print(sample, chrom)
			output1=open('run_fithic.'+sample+'.'+chrom+'.'+resolution+'.pbs','w')
			output1.write('#PBS -N fithic.'+sample+'.'+chrom+'.'+resolution+'\n')
			output1.write('#PBS -q workq\n')
			output1.write('#PBS -l nodes=1:ppn='+str(ncpu)+'\n')
			output1.write('#PBS -j oe\n')
			output1.write('\n')
			output1.write('# go workdir\n')
			output1.write('cd $PBS_O_WORKDIR\n')
			output1.write('\n')
			output1.write('# run command \n')
			output1.write('sleep 5\n')
			output1.write('\n')
			output1.write('echo -n \"I am on: \"\n')
			output1.write('hostname;\n')    
			output1.write('echo finding ssh-accessible nodes:\n')
			output1.write('echo -n \"running on: \"\n')
			output1.write('\n')
			output1.write('\n')
			output1.write('samtools_0.1.18 view '+bam_path+'/'+sample+'.sorted.nodup.bam | python bam2inter.py '+sample+' '+chrom+' '+str(min_dist)+' '+str(max_dist)+' '+resolution+' '+bin_path+' '+DIR+'/feature_vec\n')
			
			output1.write(python_path+'/python '+fithic_path+'/fithic.py -i '+DIR+'/feature_vec/intersfile.'+sample+'.'+chrom+'.'+resolution+'.gz -f '+DIR+'/feature_vec/fragsfile.'+sample+'.'+chrom+'.'+resolution+'.gz -o '+DIR+'/res_fithic/'+sample+'.'+chrom+'.'+resolution+' -U '+str(max_dist)+' -r 0 -p 2\n')
			output1.write('\n')
			output1.write('sleep 30\n')
			output1.write('exit 0')
			output1.close()

			os.system('qsub run_fithic.'+sample+'.'+chrom+'.'+resolution+'.pbs')


def get_siginter():

	for sample in sample_list:
		print(sample)
		output1=open(DIR+'/res_fithic/'+sample+'.'+resolution+'.sig_inter.txt','w')
		output1.write('frag1chr\tfrag1mid\tfrag2chr\tfrag2mid\tcontact_count\tp_value\tq_value\tbias1\tbias2\n')
		for chrom in chrom_list:
			print chrom
			input1=gzip.open(DIR+'/res_fithic/'+sample+'.'+chrom+'.'+resolution+'/FitHiC.spline_pass2.significances.txt.gz','r')
			all_input1=input1.readlines()
			for j in all_input1[1:]:
				each=j.split()
				if float(each[6]) < min_qval:
					output1.write(j)
		output1.close()


def match_ninter():

	if resolution=='40kb':
		bin_size=40000
	if resolution=='20kb':
		bin_size=20000
	if resolution=='10kb':
		bin_size=10000

	#
	inter_dict={}
	n_list=[]
	for sample in sample_list:
		print(sample)
		
		input1=open(DIR+'/res_fithic/'+sample+'.'+resolution+'.sig_inter.txt','r')
		all_input1=input1.readlines()

		tmp_dict={}
		idx=0
		for line in all_input1[1:]:
			each=line.strip().split('\t')

			frag1chr=each[0]
			frag1mid=int(each[1])
			frag1pt1=str(frag1mid - (bin_size/2))
			frag1pt2=str(frag1mid + (bin_size/2))

			frag2chr=each[2]
			frag2mid=int(each[3])
			frag2pt1=str(frag2mid - (bin_size/2))
			frag2pt2=str(frag2mid + (bin_size/2))

			frag1=frag1chr+'.'+frag1pt1+'.'+frag1pt2
			frag2=frag2chr+'.'+frag2pt1+'.'+frag2pt2

			inter_id=frag1+'-'+frag2
			qval=float(each[6])

			tmp_dict[inter_id]=qval
			idx+=1

		input1.close()

		inter_dict[sample]=tmp_dict
		n_list.append(idx)

	n_inter=min(n_list)

	#
	all_dict={}
	for sample in sample_list:
		print(sample)

		for key in sorted(inter_dict[sample].items(), key=operator.itemgetter(1))[:n_inter]:

			inter_id=key[0]
			frag1=inter_id.split('-')[0]
			frag2=inter_id.split('-')[1]

			dist=int(frag2.split('.')[1])-int(frag1.split('.')[1])

			qval=key[1]

			if all_dict.has_key(inter_id):
				tmp=all_dict[inter_id]
				tmp.append(sample)
				all_dict[inter_id]=tmp
			else:
				all_dict[inter_id]=[sample]
 

	## merged_list
        output1=open(DIR+'/res_fithic/siginter.'+resolution+'.merged.txt','w')
        output1.write('frag1\tfrag2\tdist\tnsample\tinter_type\n')
        for inter_id in all_dict:
                frag1=inter_id.split('-')[0]
                frag2=inter_id.split('-')[1]

                frag1pt1=int(frag1.split('.')[1])
                frag2pt1=int(frag2.split('.')[1])

                dist=frag2pt1-frag1pt1
		nsample=len(all_dict[inter_id])

		if nsample > 5:
			inter_type='static'
		else:
			inter_type='variable'

                new_line=[frag1, frag2, str(dist), str(nsample), inter_type]
                output1.write('\t'.join(new_line)+'\n')

        output1.close()

	## 
	for sample in sample_list:
		print(sample)
		output1=open(DIR+'/res_fithic/'+sample+'.'+resolution+'.sig_inter.filtered.txt','w')
		output1.write('frag1\tfrag2\tdist\tqval\tinter_type\n')
		for key in sorted(inter_dict[sample].items(), key=operator.itemgetter(1))[:n_inter]:
			inter_id=key[0]
			frag1=inter_id.split('-')[0]
			frag2=inter_id.split('-')[1]

			dist=int(frag2.split('.')[1])-int(frag1.split('.')[1])

			qval=key[1]
			nsample=len(all_dict[inter_id])

			if nsample > 5:
				inter_type='static'
			else:
				inter_type='variable'

			new_line=[frag1, frag2, str(dist), str(qval), inter_type]
			output1.write('\t'.join(new_line)+'\n')

		output1.close()




#run_fithic()
get_siginter()
match_ninter()





