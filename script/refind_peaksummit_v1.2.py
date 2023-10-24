import argparse
import pandas as pd
import numpy as np
import time

parser = argparse.ArgumentParser(description="Input file of this script is the output .bed from chipr and .bdg file.\n\
Aim of this script is to re-identify the peak summit.\n\
The first three columns of your input file should be chr, start, end.\n\
The summit index will be add to the last column.\n\
python package pandas, numpy and time are required.\n",\
formatter_class=argparse.RawTextHelpFormatter
)
parser.add_argument("-peak_file", help = ".bed from chipr or other tab delimited file which the first three columns are chr, start, end")
parser.add_argument("-bdg_file", help = "The bdg file can be generated from a mapping file using MACS2.")
parser.add_argument("-h_t", help = "high quality threshold. Peaks with summit bdg values greater than h_t will be classified as high quality. default=40", default='40')
parser.add_argument("-o", help = "output path", default='./output')
args = parser.parse_args()




peak = pd.read_csv(args.peak_file, delimiter='\t',  header=None)
bdg_cols = ['chr', 'start', 'end', 'score']
bdg = pd.read_table(args.bdg_file, names=bdg_cols)
h_t = int(args.h_t) #high quality threshold


# function for find peak summit
def summit(bdg_file,Chr,start, end, start_time=time.time()):
    
    
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
    max_index = scores.index(max_score)
    std = np.std(scores)
    

    print(Chr, '\t', start, '\t', end, '\t', max_index, end='\t')    
    print('elapsed time:%.2fs'%(time.time()-start_time))
    return max_index, max_score, std


# main function
peak[['summit', 'max_score', 'std']] = None
start_time = time.time()
result = peak.apply(lambda X:summit(bdg, X[0], X[1], X[2], start_time),axis=1)
result = np.vstack([np.asarray(t) for t in result])
peak[['summit', 'max_score', 'std']]  = result
peak['summit'] = peak['summit'].astype(int)
print('elapsed time:%.2fs'%(time.time()-start_time))



# save the result
peak.to_csv('%s_plus_summit_info.bed'%args.o, index=False, header=False, sep='\t')

narrowPeak = peak.iloc[:,:-2]
narrowPeak.to_csv('%s.narrowPeak'%args.o, index=False, header=False, sep='\t')

high_quality_peak = peak[peak['max_score']>h_t]
high_quality_narrowPeak = high_quality_peak.iloc[:,:-2]
high_quality_narrowPeak.to_csv('%s_high_quality_peaks_threshold_%d.narrowPeak'%(args.o, h_t), index=False, header=False, sep='\t')


print('end')
