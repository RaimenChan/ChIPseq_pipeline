import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import time
import argparse


##################################### arguments #####################################################
parser = argparse.ArgumentParser(description="script to calculate the correlation of two ChIPseq repeats\n",
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-peak_file", help = ".bed from chipr or other tab delimited file which the first three columns are chr, start, end")
parser.add_argument("-bdg_file1", help = "Repeat1 bdg file. The bdg file can be generated from a mapping file using MACS2.")
parser.add_argument("-bdg_file2", help = "Repeat2 bdg file. The bdg file can be generated from a mapping file using MACS2.")
parser.add_argument("-o", help = "output prefix", default='./output')
args = parser.parse_args()


#################################### read input files ###############################################
peak = pd.read_csv(args.peak_file, sep='\t', header=None)
peak = peak[[0,1,2]]

bdg1 = pd.read_csv(args.bdg_file1, sep='\t', names=['chr', 'start', 'end', 'signal'])
bdg2 = pd.read_csv(args.bdg_file2, sep='\t', names=['chr', 'start', 'end', 'signal'])

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

## get bdg1 info
peak[['bdg1_max_signal', 'bdg1_sum_signal']] = None
start_time = time.time()
result1 = peak.apply(lambda X:get_bdg_info(bdg1, X[0], X[1], X[2], start_time),axis=1)
result1 = np.vstack([np.asarray(t) for t in result1])
peak[['bdg1_max_signal', 'bdg1_sum_signal']] = result1
peak['bdg1_signal_per_bp'] = peak['bdg1_sum_signal'] /  peak['length']

## get bdg2 info
peak[['bdg2_max_signal', 'bdg2_sum_signal']] = None
start_time = time.time()
result2 = peak.apply(lambda X:get_bdg_info(bdg2, X[0], X[1], X[2], start_time),axis=1)
result2 = np.vstack([np.asarray(t) for t in result2])
peak[['bdg2_max_signal', 'bdg2_sum_signal']] = result2
peak['bdg2_signal_per_bp'] = peak['bdg2_sum_signal'] /  peak['length']

corr_matrix = peak.iloc[:,4:].corr().iloc[:3,3:]
######################################### save result and visualization ###################################
## save data
peak.to_excel("%s_correlation_data.xlsx"%args.o, index=False)

## heatmap
plt.figure(figsize=(10,8))
ax = sns.heatmap(corr_matrix, cmap='coolwarm', annot=True)
ax.xaxis.tick_top()
plt.savefig('%s_correlation_heatmap.jpeg'%args.o, bbox_inches='tight', dpi=700)

## scatter plot
k=1
plt.figure(figsize=(25,25))
for bdg1_val in ['bdg1_max_signal', 'bdg1_sum_signal', 'bdg1_signal_per_bp']:
    for bdg2_val in ['bdg2_max_signal', 'bdg2_sum_signal', 'bdg2_signal_per_bp']:
        fig = plt.subplot(3,3,k)
        
        x_data = peak[bdg1_val]
        y_data = peak[bdg2_val]
        
        plt.scatter(x_data, y_data)
        
        x_pos = x_data.mean()
        y_pos = y_data.max() + (y_data.max() - y_data.min()) * 0.1
        
        plt.text(x_pos, y_pos, 'correlation of\n' + bdg1_val + ' and ' + bdg2_val, fontsize=15)
        k=k+1
        
plt.savefig('%s_correlation_scatter_plot.jpeg'%args.o, bbox_inches='tight', dpi=700)

print("end")

