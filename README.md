# ChIPseq analysis pipeline

## Introduction to ChIPseq

## ChIPseq data analysis pipeline
The analysis pipeline had written in snakemake. After installing the required package, you can change the sample name and input name in config.yaml and change the parameters you want to use. Then run the command

```
snakemake -s TF_ChIPseq_snakemake_single_end.py
```



### Install the required package
```
conda create -n TF_ChIPseq_pipeline python=3.11
conda deactivate
conda activate TF_ChIPseq_pipeline

conda install -c bioconda -y snakemake=7.32.4
conda install -c bioconda fastaqc
conda install multiqc
conda install -c bioconda -y hisat2=2.2.1
conda install samtools=1.18
conda install bedtools=2.31.0

mkdir ~/biosofts
cd ~/biosofts
#bedGraphToBigWig 
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedGraphToBigWig
chmod +x bedGraphToBigWig
#bedsort
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64.v385/bedSort
chmod +x bedSort

pip install macs2
pip install snakemake


conda install pandas
conda install chip-r 
conda install openpyxl
conda install goatools
conda install seaborn
```



### QC
```
fastqc -o FastQC -t 4 sample_rep1.fastq.gz sample_rep2.fastq.gz input_rep1.fastq.gz input_rep2.fastq.gz
multiqc FastQC -f -o MultiQC MultiQC/logs/mutiqc.log 2>&1
```

### mapping
```
hisat2-build /genome_info/Specie.fasta /genome_info/Specie_hisat2_index 
## double end sequencing use hisat2 -X Specie_hisat2_index -1 sample_rep1_1.fastq.gz -2 sample_rep1_2.fastq.gz instead
hisat2 -X Specie_hisat2_index -u sample_rep1.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o sample_rep1_sorted.bam
hisat2 -X Specie_hisat2_index -u sample_rep2.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o sample_rep2_sorted.bam
hisat2 -X Specie_hisat2_index -u input_rep1.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o input_rep1_sorted.bam
hisat2 -X Specie_hisat2_index -u input_rep2.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o input_rep2_sorted.bam

```

###  generate.bdg and .bw file
bw file can be used in genome browser such as IGV.
```
macs2 pipeup -f BAM --extsize 200 -i sample_rep1_sorted.bam -sample_rep1_pipeup.bdg
macs2 pipeup -f BAM --extsize 200 -i sample_rep2_sorted.bam -sample_rep2_pipeup.bdg

samtools flagstat sample_rep1_sorted.bam sample_rep1_alignmentstat
samtools flagstat sample_rep2_sorted.bam sample_rep2_alignmentstat

#Change the sample_rep1 to sample_rep2, run again
mappedReads=$(awk 'NR==7 {print $1}' isample_rep1_alignment) #The first word in 7th row
scale=$(bc -l <<< "scale=3; 1000000 / $mappedReads")
echo " Normalizing sample_rep1_pipeup.bdg with factor $scale"
macs2 bdgopt -i sample_rep1_pipeup.bdg -m multiply -p $scale -o tem_tem_sample_rep1_normalized.bdg
sed -n '2,${{p}}' tem_tem_sample_rep1_normalized.bdg > tem_sample_rep1_normalized.bdg
bedSort tem_sample_rep1_normalized.bdg > sample_rep1_normalized.bdg

# generate fasta.fai file
samtools faidx Specie.fasta
# The first two columns of fasta.fai file is chromosome and chromosome size
cut -f1,2 Specie.fasta.fai > chrom.sizes
bedGraphToBigWig chrom.sizes bedGraphToBigWig 
```

### macs2 callpeak
```
macs2 callpeak -t sample_rep1_sorted.bam -c input_rep1_sorted.bam --name sample_rep1 --outdir /sample/rep1/macs2/withCtrl_narrow --nomodel -g 3e7 --extsize 200 -B --SPMR

macs2 callpeak -t sample_rep2_sorted.bam -c input_rep2_sorted.bam --name sample_rep2 --outdir /sample/rep2/macs2/withCtrl_narrow --nomodel -g 3e7 --extsize 200 -B --SPMR
```

### chipr select peaks
```
chipr -i /sample/rep1/macs2/withCtrl_narrow/sample_rep1_peaks.narrowpeak /sample/rep2/macs2/withCtrl_narrow/sample_rep2_peaks.narrowpeak -m 2 -o /sample/chipr/sample_chipr
```


### refind peak summit
result from chipr do not contain peak summit, refind peak summit according to the bdg file.
After chipr screening, there are still many peaks whose quality may not be very good. Select the peaks whose bdg signal of peak summit is greater than 40.
```
python refind_peaksummit_v1.2.py -peak_file /sample/chipr/sample_chipr_all.bed -bdg_file sample_rep1_normalized.bdg -h_t 40 -o sample_chipr_all
```


### peak annotation
run
`python chipseq_gene_annotation_V1.2.py -h`
 for detailed information
```
python chipseq_gene_annotation_V1.2.py -peak_file sample_chipr_all_high_quality_peaks_threshold_40.narrowPeak -gff Specie.gff -o sample_chipr_all_high_quality_peaks_threshold_40_target_gene -distance 2000 -type1_method CDS
```

### GO and KEGG analysis of target genes
ani_KEGG_info.csv is an internal use file downloaded from the KEGG website and organized into a custom format.
obo file can be download from Gene Ootology website. Here is the link http://purl.obolibrary.org/obo/go/go-basic.obo.
```
python GO_analysis_V1.py -target_gene_file sample_chipr_all_high_quality_peaks_threshold_40_target_gene.txt -gaf Specie.gaf -obo go-basic.obo -o sample_GO_result

python KEGG_analysis_V1.py -target_gene_file sample_chipr_all_high_quality_peaks_threshold_40_target_gene_KEGG.txt -KEGG_relation ani_KEGG_info.csv -o sample_KEGG_result
```

