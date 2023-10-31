import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import argparse
import sys


##################################### arguments #####################################################
parser = argparse.ArgumentParser(description="script to calculate the correlation of two ChIPseq repeats\n",
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-peak_file", help = ".bed from chipr or other tab delimited file which the first three columns are chr, start, end")
parser.add_argument("-bdg_files", nargs='+', help = "bdg files, at least provide 2 bdg files. The bdg file can be generated from a mapping file using MACS2.")
parser.add_argument("-o", help = "output prefix", default='./output')
args = parser.parse_args()


#################################### read input files ###############################################
bdg_files = args.bdg_files
if(len(bdg_files)<2):
     sys.exit("Error: At least 2 bdg files are required.")

peak = pd.read_csv(args.peak_file, sep='\t', header=None)
peak = peak[[0,1,2]]


# bdg1 = pd.read_csv(args.bdg_file1, sep='\t', names=['chr', 'start', 'end', 'signal'])
# bdg2 = pd.read_csv(args.bdg_file2, sep='\t', names=['chr', 'start', 'end', 'signal'])

######################################################################################################
def get_bdg_info(bdg_file,Chr,start, end, start_time=time.time()):
    
    
    tem_bdg = bdg_file[bdg_file['chr']==Chr]
    
    tem_bdg_start = tem_bdg[tem_bdg['start']<=start]
    start_threshold = max(tem_bdg_start['start'])
    
    tem_bdg_end = tem_bdg[tem_bdg['end']>=end]
    end_threshold = min(tem_bdg_end['end'])
    
    tem_bdg = tem_bdg[tem_bdg['start']>=start_threshold]
    tem_bdg = tem_bdg[tem_bdg['end']<=end_threshold]
    
    scores = []
    for i in range(start,end+1):
        tem_bdg_i = tem_bdg[tem_bdg['start']<=i]
        tem_bdg_i = tem_bdg_i[tem_bdg_i['end']>=i]
        tem_bdg_i_score = tem_bdg_i.iloc[0,3]

        scores.append(tem_bdg_i_score)
    
    max_score = max(scores)
    sum_score = np.sum(scores)
    

    print(Chr, '\t', start, '\t', end, '\t', 'max_score=%.2f'%max_score, '\t', 'sum_score=%.2f'%sum_score, end="\t")    
    print('elapsed time:%.2fs'%(time.time()-start_time))
    return max_score, sum_score

####################################### main function ####################################################
peak['length'] = peak[2] - peak[1]


log = open("%s_bdg_log.txt"%args.o,'w+')
log.write("calculate correlation using %s"%args.peak_file)
log.write("\n")
for i, bdg_file in enumerate(bdg_files):
    bdg = pd.read_csv(bdg_file, sep='\t', names=['chr', 'start', 'end', 'signal'])
    ## get bdg info
    print("bdg%s = %s"%(i+1, bdg_file))
    log.write("bdg%s = %s"%(i+1, bdg_file))
    log.write("\n")
    peak[['bdg%s_max_signal'%(i+1), 'bdg%s_sum_signal'%(i+1)]] = None
    start_time = time.time()
    result1 = peak.apply(lambda X:get_bdg_info(bdg, X[0], X[1], X[2], start_time),axis=1)
    result1 = np.vstack([np.asarray(t) for t in result1])
    peak[['bdg%s_max_signal'%(i+1), 'bdg%s_sum_signal'%(i+1)]] = result1
    peak['bdg%s_signal_per_bp'%(i+1)] = peak['bdg%s_sum_signal'%(i+1)] /  peak['length']

log.close()



######################################### save result and visualization ###################################
## save data
peak.to_excel("%s_correlation_data.xlsx"%args.o, index=False)

## heatmap
index_feature_dict={
    4: 'max_signal',
    5: 'sum_signal',
    6: 'signal_per_bp'
}
for i in [4,5,6]:
    cols=[]
    for j in range(i,len(peak.iloc[0,:]),3):
        if(j<len(peak.iloc[0,:])):
            cols.append(j)

    corr_matrix = peak.iloc[:,cols].corr()
    fig_length = int(3*len(bdg_files)+1)
    plt.figure(figsize=(fig_length,fig_length))
    ax = sns.heatmap(corr_matrix, cmap='coolwarm', annot=True)
    ax.xaxis.tick_top()
    plt.savefig('%s_%s_correlation_heatmap.jpeg'%(args.o, index_feature_dict[i]), bbox_inches='tight', dpi=700)

