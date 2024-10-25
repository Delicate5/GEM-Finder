from datetime import datetime
import os, random

output_d 	= '../output/pvalue_maker/'
snp_d       = '../output/data_processor/'

def p_value_maker(chrlist, query, jobid, iterate, stage, peak_type):
    print(datetime.now().strftime('%H:%M:%S'), "["+jobid+"_pvalue_maker] calculate p-value")
    real_peak = '/home/pgs8070/oa13_final/bedfile/'+peak_type+'_'+stage+'.bed'

    real_count_dict = {}
    for name in query: real_count_dict[name] = int(os.popen('bedtools intersect -sorted -a "'+real_peak+'" -b "'+snp_d+name+'_LDSnp.bed" -wa -wb | wc -l').read().strip())

    results = {}
    alarm = max(iterate/10, 1)
    for name in query: results[name] = [0, 0]
    for i in range(iterate):
        if not (i % alarm): print(datetime.now().strftime('%H:%M:%S'), "["+jobid+"_pvalue_maker] "+str(i)+"th iteration...")
        os.system('bedtools shuffle -chrom -i "'+real_peak+'" -g /home/pgs8070/ref/rgt_config/hg19/chrom.sizes.hg19 | bedtools sort > tmp/'+peak_type+'_'+stage+'_'+jobid+'.bed')
        for name in query:
            rand_intersect = int(os.popen('bedtools intersect -sorted -a tmp/'+peak_type+'_'+stage+'_'+jobid+'.bed -b "'+snp_d+name+'_LDSnp.bed" -wa -wb | wc -l').read().strip())
            if real_count_dict[name] > rand_intersect: results[name][0] += 1
            results[name][1] += rand_intersect
    
    with open(output_d+'result_'+peak_type+'_'+stage+'_'+jobid+'.txt', 'w') as wp:
        wp.write('#trait_id\tmore_than_random_peak_while_'+str(iterate)+'\tsum_of_random_peak_while_'+str(iterate)+'\n')
        for k, v in results.items(): wp.write('\t'.join([k] + [str(j) for j in v])+'\n')

    os.system('rm tmp/'+peak_type+'_'+stage+'_'+jobid+'.bed')

def p_value_merger(query, njob, stage, peak_type):
    real_peak = '/home/pgs8070/oa13_final/bedfile/'+peak_type+'_'+stage+'.bed'

    results = {}
    whole_iter = 0
    for name in query: results[name] = [0, 0, 0]
    for name in query: results[name][0] = int(os.popen('bedtools intersect -sorted -a "'+real_peak+'" -b "'+snp_d+name+'_LDSnp.bed" -wa -wb | wc -l').read().strip())
    for i in range(int(njob)):
        with open(output_d+'result_'+peak_type+'_'+stage+'_'+str(i)+'.txt', 'r') as fp:
            whole_iter += int(next(fp).strip().split('_')[-1])
            for j in fp:
                try: name, count, npeak = j.strip().split('\t')
                except: name, npeak, count, npeak = j.strip().split('\t')
                results[name][1] += int(count)
                results[name][2] += int(npeak)
    with open(output_d+'results_'+peak_type+'_'+stage+'.txt', 'w') as wp:
        wp.write("#trait_id\treal_peak\tmore_than_random_peak_while_"+str(whole_iter)+"\tsum_of_random_peak_while_"+str(whole_iter)+"\n")
        for k, v in results.items(): wp.write('\t'.join([k] + [str(j) for j in v])+'\n')
