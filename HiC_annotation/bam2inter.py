#!/usr/bin/python



import fileinput
import sys
import gzip


sample=sys.argv[1]
chrom=sys.argv[2]
min_dist=int(sys.argv[3])
max_dist=int(sys.argv[4])
resolution=sys.argv[5]
bin_path=sys.argv[6]
output_path=sys.argv[7]


if resolution=='5kb':
	bin_size=5000
if resolution=='10kb':
	bin_size=10000
if resolution=='20kb':
	bin_size=20000
if resolution=='40kb':
	bin_size=40000


inter_dict={}
cov_dict={}
for line in sys.stdin:
        # HWI-ST1113:549:HKWJFADXX:2:1115:20976:94471   163     chr22   16052163        41      64M36S  chrY    16541071        0       TCAAAAAACAAACAAACAAACAAACAAACAAACAAACTGTCAAAATCTGTACAGTATGTGAAGAAGCTTAGTTTTTTCGCTCTTTGCAATAAATCTTGCT    CCCFFFFFHHHHHJJJJJJJJJJJJJJJJJJJJJJJIGIJDGJJJJIJJJJJJIIIIJGHIIJJHHGHHFFCEFFDDDDDDDDDDDDEDDDDDEDDDDDC    NM:i:0  MD:Z:64 AS:i:64 XS:i:51 SA:Z:chrY,16540953,+,62S38M,27,0;
        [id1, flag1, pos1chr, loc_from1, mapq1, cigar1, pos2chr, loc_from2, dist, read1, read_qual1]=line.split('\t')[0:11]

	pos1=int(loc_from1)
	pos2=int(loc_from2)

	#
	if pos1chr==chrom:
		if pos2chr=='=':
			if abs(int(dist)) > min_dist:
				if abs(int(dist)) <= max_dist:

					pos1bin=int(loc_from1)/bin_size*bin_size
					pos2bin=int(loc_from2)/bin_size*bin_size

					#			 
					if pos1bin<pos2bin:
						inter_id=chrom+'.'+str(pos1bin)+'.'+str(pos2bin)
					else:
						inter_id=chrom+'.'+str(pos2bin)+'.'+str(pos1bin)

					if inter_dict.has_key(inter_id):
						inter_dict[inter_id]+=1
					else:
						inter_dict[inter_id]=1

					#
					if cov_dict.has_key(pos1bin):
						cov_dict[pos1bin]+=1
					else:
						cov_dict[pos1bin]=1

					if cov_dict.has_key(pos2bin):
						cov_dict[pos2bin]+=1
					else:
						cov_dict[pos2bin]=1
			
#
output1=gzip.open(output_path+'/intersfile.'+sample+'.'+chrom+'.'+resolution+'.gz','w')
for i in inter_dict:
	info=i.split('.')
	output1.write(info[0]+'\t'+str(int(info[1])+(bin_size/2))+'\t'+info[0]+'\t'+str(int(info[2])+(bin_size/2))+'\t'+str(inter_dict[i]/2)+'\n')
output1.close()	


#
bin_list=[]
input1=open(bin_path+'/'+chrom+'.'+resolution+'.bin','r')
all_input1=input1.readlines()
for line in all_input1:
	each=line.strip().split('\t')
	bin_id=each[0]+'.'+each[1]+'.'+each[2]
	bin_list.append(bin_id)
input1.close()


#
output1=gzip.open(output_path+'/fragsfile.'+sample+'.'+chrom+'.'+resolution+'.gz','w')
for bin_id in bin_list:

	pt1=int(bin_id.split('.')[1])

	mid=pt1+(bin_size/2)
	if cov_dict.has_key(pt1):
		coverage=cov_dict[pt1]
	else:
		coverage=0

	new_line=[chrom, bin_id, str(mid), str(coverage)]
	output1.write('\t'.join(new_line)+'\n')

output1.close()



		
