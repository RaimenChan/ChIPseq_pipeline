import pandas as pd
from scipy.stats import hypergeom
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.cm as cm
import seaborn as sns
import os
import argparse

##################################### arguments #####################################################

parser = argparse.ArgumentParser(description="According to the hypergeometric distribution, calculate the enrichment probability of the target genes in various KEGG terms,\n\
 and select the KEGG terms with a p-value less than 0.05 for visualization.\n",
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-target_gene_file", help="target gene should save in a txt file format like \ngene1\ngene2\ngene3\n...")
parser.add_argument("-KEGG_relation", help = "in-home used file.")
parser.add_argument("-o", help="output", default='./output')

args = parser.parse_args()

#################################### read input files ###############################################
df =pd.read_csv(args.KEGG_relation)
with open(args.target_gene_file, 'r') as file:
    file_content = file.read()
Target_gene_list = file_content.split('\n')

KEGG_ID={}
group = df.groupby(['children3'])
for key,tem_df in group:
    term = key[0]
    Gene_IDs = tem_df['children4 (gene_ID)'].to_numpy().tolist()
    #同一条GO-ID记录可能出现多次,KEGG不清楚，也跑一遍
    Gene_IDs = list(set(Gene_IDs))
    KEGG_ID[term] = Gene_IDs


#################################### enrichment analysis ##############################################
results =[]

assoc = KEGG_ID
background_gene_ID = df['children4 (gene_ID)'].to_numpy().tolist()
background_gene_ID = list(set(background_gene_ID))
Study = list( 
        set(Target_gene_list).intersection( set(background_gene_ID) )
    )
for KEGG in KEGG_ID.keys():
    KEGG_result = {}
    Genes_in_KEGG = assoc[KEGG]
    Population = len(background_gene_ID)
    Term_size = len(Genes_in_KEGG)
    Intersection = list(
        set(Genes_in_KEGG).intersection( set(Study))
        )
    Study_size = len(Study)
    Intersection_size = len(Intersection)
    
    #Intersection_size - 1才是我们想要的pvalue
    p_value = hypergeom.sf(Intersection_size - 1, Population, Term_size, Study_size)
    
    KEGG_result['KEGG Term'] = KEGG
    KEGG_result['Category'] = 'KEGG'
    KEGG_result['Population'] = Population
    KEGG_result['Term_size'] = Term_size
    KEGG_result['Study_size'] = Study_size
    KEGG_result['Intersection_size'] = Intersection_size
    KEGG_result['p_value'] = p_value
    KEGG_result['KEGG_item'] = Genes_in_KEGG
    KEGG_result['Study_item'] = Study
    KEGG_result['Intersection'] = Intersection

    results.append(KEGG_result)

results_df = pd.DataFrame(results)
results_df.to_excel(args.o+'.xlsx', index=False)


################################ visualization #################################
results_df['Gene%']=results_df['Intersection_size']/results_df['Study_size']*100
results_df = results_df.sort_values('p_value')
df_005 = results_df[results_df['p_value']<0.05].copy()


##### if too many items in each category whose p < 0.05, show them in different picture, each picture show 20 items
for i in range(len(df_005)//20+1):
    if(20*i+20<=len(df_005)):
        df_005_i = df_005.iloc[20*i:20*i+20,:]
    else:
        df_005_i = df_005.iloc[20*i:len(df_005),:]
        
    if(len(df_005_i)==0):
        continue


    #### parameter used to draw the picture
    X = df_005_i['Gene%'].to_numpy().tolist()
    Y = df_005_i['KEGG Term'].to_numpy().tolist()
    P = df_005_i['p_value'].to_numpy().tolist()
    N = df_005_i['Intersection_size'].to_numpy().tolist()

    X = X[::-1]
    Y = Y[::-1]
    P = P[::-1]
    N = N[::-1]





    #### draw the picture
    plt.figure(figsize=(10,15))

    norm = plt.Normalize(0, 0.1)
    cmap = cm.get_cmap('Blues_r')  # Replace 'cool' with the desired colormap
    colors = cmap(norm(P)).tolist()

    plt.barh( Y, X, color=colors, height =0.8, alpha =1 )
    cbar = plt.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=norm), orientation='vertical', ax=plt.gca())
    cbar.ax.tick_params(labelsize=20)  # Increase the font size of the colorbar tick labels
    cbar.set_label('p value', fontsize=20, fontweight='bold')

    plt.xlabel('Gene[%]', fontsize=25)
    plt.xlim(0,max(X)+10)
    plt.yticks(fontsize=20, fontname='Times New Roman', fontweight='bold')
    plt.xticks(fontsize=20, fontname='Times New Roman')

    for i2, (x, y, p, n) in enumerate(zip(X, Y, P, N)):
        plt.text(x+0.3, i2, 'n=%d \np=%.2e'%(n, p), ha='left', va='center', font='Times New Roman', fontsize = 15)

    plt.savefig(args.o + '_KEGG_%d.jpeg'%i, dpi=700, bbox_inches='tight')
    #plt.show()

print("end")