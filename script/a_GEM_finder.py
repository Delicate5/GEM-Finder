import matplotlib.pyplot as plt

def DEG_selector():
    DEG_list = {}
    for stage in stages:
        DEG_list[stage] = []
        with open('symbolic/DEGs/'+stage+'EC.valTable.txt', 'r') as fp:
            header = next(fp).strip().split('\t')
            for i in fp:
                line = i.strip().split('\t')
                if(abs(float(line[1])) > 1) : DEG_list[stage].append(line[0])

    DEG_symbols = {}
    for stage in stages: DEG_symbols[stage] = []
    with open('symbolic/DEGs/RNA.STable.txt', 'r') as fp:
        header = next(fp).strip().split('\t')
        for i in fp:
            line = i.strip().split('\t')
            for cstage, degs in DEG_list.items():
                if(line[0]) in degs: DEG_symbols[cstage].append(line[1])

    for stage, degs in DEG_symbols.items():
        with open('DEG/DEG_'+stage+'.txt', 'w') as wp: wp.write('#'+stage+'symbol\n'+'\n'.join(degs))

#DEG_contact_enhancer_finder. tried to put connected gene, hard at sorting. remove
def DEG_contact_enhancer_finder():
    DEG_list = {}
    for stage in stages:
        DEG_list[stage] = []
        with open('DEG/control_DEG_'+stage+'.txt', 'r') as fp:
            header = next(fp).strip().split('\t')
            for i in fp: DEG_list[stage].append(i.strip())

    for stage in stages:
        with open('bedfile/control_DEG_contact_enhancer_'+stage+'_unsorted.bed', 'w') as wp:
            wp.write('#chrID\tpos1\tpos2\n')
            with open('/home/pgs8070/oa13_final/symbolic/res_fithic/'+time_stages[stage]+'.20kb.sig_inter.filtered.anno.txt', 'r') as fp:
                header = next(fp).strip().split('\t')
                for i in fp:
                    line = i.strip().split('\t')
                    if line[8] == 'peak_prom':
                        for gene in line[5].split(';'):
                            if gene in DEG_list[stage]:
                                for enhancer in line[3].split(';'): wp.write(enhancer.replace(':', '\t').replace('-', '\t')+'\n')
