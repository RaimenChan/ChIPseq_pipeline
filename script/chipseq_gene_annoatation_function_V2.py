def type1_downstream_geneannotation(gff, Chr, summit,distance, type3): 
    ###基因重叠，一个peak可能有两个type3
    if(type3!=None):
        type3 = type3.split(', ')
    else:
        type3 = []
    
    # peak 的下游没有基因， min_start ='flag'; 有基因，min_start会变成数值
    min_start = 'flag'
    
    gff_chr = gff[gff['chr']==Chr]
    gff_chr_down = gff_chr[gff_chr['start']>=summit]
    gff_chr_down = gff_chr_down[gff_chr_down['strand'] == '+']
    gff_chr_down = gff_chr_down.sort_values('start', ascending=True)
    
    ##找出不在gene body上的最近的基因
    for i in range(len(gff_chr_down)):
        tem_gene = gff_chr_down.iloc[i]['ID']
        tem_gene = tem_gene.split(';')[0]
        tem_gene = tem_gene.split('-')[0]
        tem_gene = tem_gene.split('=')[1]
        if(tem_gene in type3):
            continue
        else:
            min_start = gff_chr_down.iloc[i]['start']
            break
    
    #peak的下游没有基因
    if(min_start=='flag'):
        return None, None
    
    d = min_start - summit
    if(d<=distance):
        gff_chr_down_target = gff_chr_down[gff_chr_down['start']==min_start]
        if(len(gff_chr_down_target)==1):
            target_ID_value = gff_chr_down_target.iloc[0]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('-')[0]
            target_ID = target_ID.split('=')[1]
            return target_ID, d
        else:
            print(Chr,'\t', summit, 'down stream abnormal')
            return 'type1_down_abnormal', None
    else:
        return None, None

    
def type1_upstream_geneannotation(gff, Chr, summit,distance, type3):
    if(type3 != None):
        type3 = type3.split(', ')
    else:
        type3 = []
    
    # peak 的上游没有基因, max_end = 'flag'; 有基因，max_end 会变成数值
    max_end = 'flag'
    
    gff_chr = gff[gff['chr']==Chr]
    gff_chr_up = gff_chr[gff_chr['end']<=summit]
    gff_chr_up = gff_chr_up[gff_chr_up['strand'] == '-']
    gff_chr_up = gff_chr_up.sort_values('end', ascending=False)
    
    ##找出不在gene body上的最近的基因
    for i in range(len(gff_chr_up)):
        tem_gene = gff_chr_up.iloc[i]['ID']
        tem_gene = tem_gene.split(';')[0]
        tem_gene = tem_gene.split('-')[0]
        tem_gene = tem_gene.split('=')[1]
        if(tem_gene in type3):
            continue
        else:
            max_end = gff_chr_up.iloc[i]['end']
            break
    
    # peak 的上游没有基因
    if(max_end == 'flag'):
        return None, None
    
    d = summit - max_end
    if(d<=distance):
        gff_chr_up_target = gff_chr_up[gff_chr_up['end']==max_end]
        if(len(gff_chr_up_target)==1):
            target_ID_value = gff_chr_up_target.iloc[0]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('-')[0]
            target_ID = target_ID.split('=')[1]
            return target_ID, d
        else:
            print(Chr,'\t', summit, 'up stream abnormal')
            return 'type1_down_abnormal', None
    else:
        return None, None



def type2_geneannotation(gff, Chr, summit):
    gff_chr = gff[gff['chr']==Chr]
    gff_chr = gff_chr[gff_chr['start']<=summit]
    gff_chr = gff_chr[gff_chr['end']>=summit]
    if(len(gff_chr)==1):
        target_ID_value = gff_chr.iloc[0]['ID']
        target_ID = target_ID_value.split('-')[0]
        target_ID = target_ID.split('=')[1]
        target_ID  = target_ID .split('_')[1]
        return target_ID
    
    elif(len(gff_chr)==0):
        return None
    
    #
    else:
        print('peak', '\t', Chr,'\t', summit, '\t', "type2 more than 1 gene(five_prime_UTR)")
        print(gff_chr)
        print()
        
        target_IDs = []
        for i in range(len(gff_chr)):
            target_ID_value = gff_chr.iloc[i]['ID']
            target_ID = target_ID_value.split('-')[0]
            target_ID = target_ID.split('=')[1]
            target_ID  = target_ID .split('_')[1]
            target_IDs.append(target_ID)
        
        gene_IDs = ''
        for i,ID in enumerate(target_IDs):
            gene_IDs = gene_IDs + ID
            if(i!=len(target_IDs)-1):
                gene_IDs = gene_IDs + ', '
        return gene_IDs



# this function will return genes also in type2, need to remove the type2  gene.
def type3_geneannotation(gff, Chr, summit):
    gff_chr = gff[gff['chr']==Chr]
    gff_chr = gff_chr[gff_chr['start']<=summit]
    gff_chr = gff_chr[gff_chr['end']>=summit]
    
    # peak 不绑定在gene body上
    if(len(gff_chr)==0):
        return None
    
    # peak 绑定在gene body上，正常情况
    elif(len(gff_chr)==1):
        target_ID_value = gff_chr.iloc[0]['ID']
        target_ID = target_ID_value.split(';')[0]
        target_ID = target_ID.split('=')[1]
        return target_ID
    
    # peak绑定在gene body上，但大于一个基因，异常情况
    else:
        print('peak', '\t', Chr,'\t', summit, 'type3 more than 1 gene')
        print(gff_chr)
        print()
        
        target_IDs = []
        for i in range(len(gff_chr)):
            target_ID_value = gff_chr.iloc[i]['ID']
            target_ID = target_ID_value.split(';')[0]
            target_ID = target_ID.split('=')[1]
            target_IDs.append(target_ID)
        
        gene_IDs = ''
        for i,ID in enumerate(target_IDs):
            gene_IDs = gene_IDs + ID
            if(i!=len(target_IDs)-1):
                gene_IDs = gene_IDs + ', '
        return gene_IDs


def remove_type3_in_type2(type2, type3):
    if(type2 == type3):
        return None
    else:
        return type3
    

def target_gene_list(type3, type2, type1_down, type1_up):
    target_genes =[]
    for genes in type3:
        target_genes = target_genes + genes.split(',')
    for genes in type2:
        target_genes = target_genes + genes.split(',')
    for genes in type1_down:
        target_genes = target_genes + genes.split(',')
    for genes in type1_up:
        target_genes = target_genes + genes.split(',')
    return target_genes


def two_type1_selection(d2, downg, down_d, upg, up_d):
    ##downg:down stream type1 gene
    ##upg: upstream type1 gene
    ##down_d, up_d: up/down distance to peak summit

    if(down_d!=None):
        down_d = int(down_d)
    if(up_d!=None):
        up_d = int(up_d)

    #只有上下游都有基因的时候再根据d2筛选,不然输入什么，输出什么
    if(downg != None and upg != None):
        #两边都小于d2
        if(down_d < d2 and up_d<d2):
            return downg, down_d, upg, up_d
        #只有down小于d2
        elif(down_d <d2 and up_d>=d2):
            return downg, down_d, None, None
        #只有up小于d2
        elif(down_d >= d2 and up_d<d2):
            return None, None, upg, up_d
        #都大于等于d2
        else:
            mind = min(down_d, up_d)
            if(down_d==mind and up_d==mind):
                return downg, down_d, upg, up_d
            elif(down_d==mind):
                return downg, down_d, None, None
            else:
                return None, None, upg, up_d
    else:
        return downg, down_d, upg, up_d
    