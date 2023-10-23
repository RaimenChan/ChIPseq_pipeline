# ChIPseq_pipeline

## Introduction to ChIPseq

## ChIPseq data analysis pipeline
The analysis pipeline had written in snakemake. After installing the required package, you can change the sample name and input name in config.yaml and change the parameters you want to use. Then run the command

```
snakemake -s TF_ChIPseq_snakemake_single_end.py
```



### Install the required package





### QC
```
fastqc -o FastQC -t 4 fastq/sample_rep1.fastq.gz fastq/sample_rep2.fastq.gz fastq/input_rep1.fastq.gz fastq/input_rep2.fastq.gz
multiqc FastQC -f -o MultiQC MultiQC/logs/mutiqc.log 2>&1
```

### mapping
```
hisat2-build /genome_info/Specie.fasta /genome_info/Specie_hisat2_index 
## double end sequencing use hisat2 -X Specie_hisat2_index -1 fastq/sample_rep1_1.fastq.gz -2 fastq/sample_rep1_2.fastq.gz instead
hisat2 -X Specie_hisat2_index -u fastq/sample_rep1.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o sample_rep1_sorted.bam
hisat2 -X Specie_hisat2_index -u fastq/sample_rep2.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o sample_rep2_sorted.bam
hisat2 -X Specie_hisat2_index -u fastq/input_rep1.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o input_rep1_sorted.bam
hisat2 -X Specie_hisat2_index -u fastq/input_rep2.fastq.gz | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o input_rep2_sorted.bam

```

###  generate.bdg and .bw file
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
```
python refind_peaksummit_v1.2.py -peak_file /sample/chipr/sample_chipr_all.bed -bdg_file sample_rep1_normalized.bdg -h_t 40 -o sample/chipr/sample_chipr_all
```


### peak annotation


### GO and KEGG analysis of target genes

