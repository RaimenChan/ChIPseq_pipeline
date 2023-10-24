import pandas as pd
import time
from chipseq_gene_annoatation_function_V1 import type1_downstream_geneannotation, type1_upstream_geneannotation
from chipseq_gene_annoatation_function_V1 import type2_geneannotation, type3_geneannotation, remove_type3_in_type2
from chipseq_gene_annoatation_function_V1 import target_gene_list
import argparse

start_time = time.time()
print('\n\n')
print('start time:', time.ctime())


##################################### arguments #####################################################

parser = argparse.ArgumentParser(description="This is a in-home script developed by Raimen (ruiwenchen@um.edu.mo) for annotating target genes based on a peak file.\n\
 The script performs annotations for three types of target genes:\n\
    Type1 target genes: These are the nearest genes located upstream or downstream of the peak summit.\n\
    Type2 target genes: These are the genes where the peak summit falls within the 5'-UTR region.\n\
    Type3 target genes: These are the genes where the peak summit falls within the gene body region, excluding the type 2 target genes.\n\
    python package: pandas is required.\n\
    Before using this script, please ensure chipseq_gene_annoatation_function-VX.X.py is appear in the same folder.\n\
",\
formatter_class=argparse.RawTextHelpFormatter
)

parser.add_argument("-peak_file", help="Tab-delimited file containing peak information.\n"
                                        "If the peak_file is not in narrowPeak format, please use -chr_index, "
                                        "-start_index, -end_index, -summit_index to assign the column indices "
                                        "of chromosome, start, end, and summit."
                                        )

parser.add_argument("-gff", help = "gff_file")

parser.add_argument("-o", help = "output path, default=./output", default='./output')

parser.add_argument("-chr_index", help = "chromosome column index in peak_file, index starting from 0, default=0",\
 default=0)

parser.add_argument("-start_index", help = "start column index in peak_file, index starting from 0, default=1",\
 default=1)

parser.add_argument("-end_index", help = "end column index in peak_file, index starting from 0, default=2", default=2)

parser.add_argument("-summit_index", help = "summit column index in peak_file, index starting from 0, default=9",\
 default=9)

parser.add_argument("-type1_method", help = "use CDS or TSS to annotate type1 gene, defalut=CDS",\
 choices=['CDS', 'TSS'], default='CDS')

parser.add_argument("-distance", help = "max distance in type1 gene annotation, defalut=2000",  default=2000)

parser.add_argument("-ncRNA_gene", help="Consider ncRNA_gene or not. Set to '1' to include ncRNA_gene,\
 or '0' to exclude ncRNA_gene. Default is '0' (unconsider).", choices=['0', '1'], default='0')

args = parser.parse_args()

chr_index = int(args.chr_index)
start_index = int(args.start_index)
end_index = int(args.end_index)
summit_index = int(args.summit_index)
distance = int(args.distance)
type1_method = args.type1_method
ncRNA_gene = int(args.ncRNA_gene)

print("chr_index=%d, start_index=%d end_index=%d, summit_index=%d distance=%d, type1_method=%s"\
%(chr_index, start_index, end_index, summit_index, distance, type1_method))
print('\n\n')




#################################### read input files ###############################################

genome_gff_header = ['chr','source','Genome_attr','start','end','score','strand','phase','ID']
genome_gff = pd.read_csv(args.gff, sep='\t', comment='#', names=genome_gff_header)
CHIP_seq=pd.read_csv(args.peak_file, sep='\t', header = None)


CHIP_seq_output = CHIP_seq[[chr_index, start_index, end_index, summit_index]].copy()
CHIP_seq_output['abs_summit'] = CHIP_seq_output[start_index] + CHIP_seq_output[summit_index]
genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('protein_coding_gene','gene')
genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('pseudogene','gene')
if(ncRNA_gene == 1):
    genome_gff['Genome_attr'] = genome_gff['Genome_attr'].replace('ncRNA_gene','gene')

################## sub_genome_gff for annotating type1, type2 and type3 genes  #####################
#type1
if(type1_method == 'CDS'):
    genome_gff_type1 = genome_gff[genome_gff['Genome_attr']=='CDS'].copy()
elif(type1_method == 'TSS'):
    genome_gff_type1 = genome_gff[genome_gff['Genome_attr']=='gene'].copy()
#type2
genome_gff_fiveUTR = genome_gff[genome_gff['Genome_attr']=='five_prime_UTR'].copy()
#type3
genome_gff_gene = genome_gff[genome_gff['Genome_attr']=='gene'].copy()


###################################   main function   ################################################
print(time.ctime())
print("Annotating type3 target genes")
CHIP_seq_output['type3'] = CHIP_seq_output.apply(lambda X:
                                                      type3_geneannotation(genome_gff_gene, X[chr_index],
                                                      X['abs_summit']), axis=1)

print(time.ctime())
print("Annotating type2 target genes")
CHIP_seq_output['type2'] = CHIP_seq_output.apply(lambda X:
                                                      type2_geneannotation(genome_gff_fiveUTR, X[chr_index], 
                                                      X['abs_summit']), axis=1)

print(time.ctime())
print("Annotating type1 downstream target genes")
CHIP_seq_output['type1_down'] = CHIP_seq_output.apply(lambda X:
                                                      type1_downstream_geneannotation(genome_gff_type1,
                                                       X[chr_index], X['abs_summit'], distance, X['type3']), axis=1)

print(time.ctime())
print("Annotating type1 upstream target genes")
CHIP_seq_output['type1_up'] = CHIP_seq_output.apply(lambda X:
                                                      type1_upstream_geneannotation(genome_gff_type1, X[chr_index], 
                                                      X['abs_summit'], distance, X['type3']), axis=1)

print(time.ctime())
print("running function remove_type3_in_type2")
CHIP_seq_output['type3'] = CHIP_seq_output.apply(lambda X:
                                                      remove_type3_in_type2(X['type2'], X['type3']), axis=1)

target_genes = target_gene_list(CHIP_seq_output['type3'].dropna(),CHIP_seq_output['type2'].dropna(),
                                CHIP_seq_output['type1_down'].dropna(),CHIP_seq_output['type1_up'].dropna())
target_genes = list(set(target_genes))


######################################## output result #########################################################
CHIP_seq_output.to_excel('%s.xlsx'%args.o, index=False)

with open('%s.txt'%args.o, 'w+') as file:
    for gene in target_genes:
        file.write(str(gene) + '\n')

with open('%s_KEGG.txt'%args.o, 'w+') as file:
    for gene in target_genes:
        file.write(str(gene) +'.2' + '\n')

######################################## end ###################################################################
print('\n\n')
elapsed_time = time.time() - start_time
print('end time:', time.ctime())
print('running time = %.2fs'%elapsed_time)
print('end')