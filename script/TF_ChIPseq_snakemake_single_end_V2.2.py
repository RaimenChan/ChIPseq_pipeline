# config.yaml
# FASTQ: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/fastq
# samples: ["CreAHA_Glucose", "CreAHA_Carbon_Free", "CreAHA_Proline"]
# Input: "Untag_control_for_HA"
# repeats: ["Repeat1", "Repeat2"]
# fqsuffix: fastq
# gff: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/genome_info/FungiDB-65_AnidulansFGSCA4.gff
# fasta: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/genome_info/FungiDB-65_AnidulansFGSCA4_Genome.fasta
# gfolder: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/genome_info
# output: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/output
# ncores : 4
# HT2index: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/genome_info/Anidulans_hisat2_index
# trim5: 8
# genomesize: 3e7
# extsize: 200
# bedGraphToBigWig: /home/ruiwenchen/Biosofts/bedGraphToBigWig
# bedSort: /home/ruiwenchen/Biosofts/bedSort
# m: 2 #chipr select peaks which appear in at least m repeats.
# h_t: 40 #high quality peak threshold, Peaks with summit bdg values greater than "h_t" will be classified as high quality.
# type1_method: CDS
# distance: 2000
# gaf: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/reference/FungiDB-65_AnidulansFGSCA4_Curated_GO.gaf
# obo: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/reference/go-basic.obo
# KEGG_relation: /home/ruiwenchen/CreA_test_ChIPseq_pipeline_V2/reference/ani_KEGG_info.csv

configfile: "./config_V2.yaml"

outputdir = config["output"]
FASTQdir = config["FASTQ"]
samples = config["samples"]
gff_file = config["gff"]
repeats = config["repeats"]
Input = config["Input"]

wildcard_constraints:
    sample= '|'.join([re.escape(x) for x in samples]),
    repeat = '|'.join([re.escape(x) for x in repeats])



rule all:
    input:
        ##QC
        expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", "{repeat}", "_fastqc.zip"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "FastQC", "".join([Input, "_", "{repeat}", "_fastqc.zip"])), repeat=repeats),
        os.path.join(outputdir, "MultiQC", "multiqc_report.html"),
        
        ## reference use in bedGraphToBigWig 
        config["fasta"] + '.fai',
        os.path.join(config["gfolder"], "chrom.sizes"),
        
        #macs2 result
        expand(os.path.join(outputdir, Input, "{repeat}", "hisat2", "".join([Input, "_", "{repeat}","_sorted.bam"])), repeat=repeats),
        expand(os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_sorted.bam",])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_pileup.bdg"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_alignmentstat",])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["tem_","{sample}", "_", "{repeat}", "_normalized.bdg"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_normalized.bdg"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_normalized.bw"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_peaks.narrowPeak"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_peaks.xls"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_summits.bed"])), sample=samples, repeat=repeats),
        
        ## chipr result
        expand(os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_all.bed"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_optimal.bed"])), sample=samples),

        ## inhome script result
        expand(os.path.join(outputdir,"{sample}","correlation","".join(["{sample}_correlation_data.xlsx"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all.narrowPeak"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all_plus_summit_info.bed"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), ".narrowPeak"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene.xlsx"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene.txt"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene_KEGG.txt"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result.xlsx"])), sample=samples),
        expand(os.path.join(outputdir, "{sample}", "KEGG","".join(["{sample}", "_KEGG_result.xlsx"])), sample=samples)





## ------------------------------------------------------------------------------------ ##
## Quality control
## ------------------------------------------------------------------------------------ ##
## FastQC, original reads
rule fastqc:
    input:
        os.path.join(FASTQdir, "".join(["{sample}_", "{repeat}", ".", str(config["fqsuffix"]), ".gz"]))
    output:
        os.path.join(outputdir, "FastQC", "".join(["{sample}_", "{repeat}", "_fastqc.zip"]))
    params:
        FastQC = lambda wildcards, output: os.path.dirname(output[0])  
    singularity:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads:
        config["ncores"]
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input}"

rule fastqc_Input:
    input:
        os.path.join(FASTQdir, "".join([Input, "_", "{repeat}", ".", str(config["fqsuffix"]), ".gz"]))
    output:
        os.path.join(outputdir, "FastQC", "".join([Input, "_", "{repeat}", "_fastqc.zip"]))
    params:
        FastQC = lambda wildcards, output: os.path.dirname(output[0]) 
    singularity:
        "docker://biocontainers/fastqc:v0.11.9_cv8"
    threads:
        config["ncores"]
    shell:
        "fastqc -o {params.FastQC} -t {threads} {input}"


## MultiQC
rule multiqc:
    input:
        expand(os.path.join(outputdir, "FastQC", "".join(["{sample}_", "{repeat}", "_fastqc.zip"])), sample=samples, repeat=repeats),
        expand(os.path.join(outputdir, "FastQC", "".join([Input, "_", "{repeat}", "_fastqc.zip"])), repeat=repeats)
    output:
        os.path.join(outputdir, "MultiQC", "multiqc_report.html")
    params:
        inputdirs = os.path.join(outputdir, "FastQC"),
        MultiQCdir = lambda wildcards, output: os.path.dirname(output[0])  ## dirname of first output
    log:
        os.path.join(outputdir, "logs", "multiqc.log")
    singularity:
        "docker://ewels/multiqc:1.7"
    shell:
        "multiqc {params.inputdirs} -f -o {params.MultiQCdir} > {log} 2>&1"


rule samtools_faidx:
    input:
        fasta = config["fasta"]
    output:
        fastafai = config["fasta"] + '.fai'
    shell:
        "samtools faidx {input.fasta}"

rule cutfastafai:
    input:
        config["fasta"] + '.fai'
    output:
        os.path.join(config["gfolder"], "chrom.sizes")
    shell:
        "cut -f1,2 {input} > {output}"

rule hisat2_Input:
    input:
        os.path.join(FASTQdir, "".join([Input, "_", "{repeat}", ".", str(config["fqsuffix"]), ".gz" ]))
    output:
        os.path.join(outputdir, Input, "{repeat}", "hisat2", "".join([Input, "_", "{repeat}","_sorted.bam"]))
    params:
        index = config["HT2index"],
        trim5 = config["trim5"]
    singularity:
        "docker://zlskidmore/hisat2:2.1.0"
    shell:
         "hisat2 --trim5 {params.trim5} -x {params.index} -U {input} | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o {output}"


rule hisat2:
    input:
        os.path.join(FASTQdir, "".join(["{sample}_", "{repeat}", ".", str(config["fqsuffix"]), ".gz" ]))
    output:
        os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_sorted.bam",]))
    params:
        index = config["HT2index"],
        trim5 = config["trim5"]
    singularity:
        "docker://zlskidmore/hisat2:2.1.0"
    shell:
         "hisat2 --trim5 {params.trim5} -x {params.index} -U {input}  | samtools view -@ 4 -bS - | samtools sort -@ 4  -O bam -o {output}"



rule pileup:
    input:
        os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_sorted.bam",]))
    output:
        os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_pileup.bdg"])) 
    shell:
        "macs2 pileup -f BAM --extsize 200 -i {input} -o {output}"

rule samtools_alignmentstat:
    input:
        os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_sorted.bam",]))
    output:
        os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_alignmentstat",]))
    shell:
        "samtools flagstat {input} > {output}"

rule bdgopt:
    input:
        pileup = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_pileup.bdg"])),
        alignment = os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_alignmentstat",]))
    output:
        tem1 = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["tem_","{sample}", "_", "{repeat}", "_normalized.bdg"])) 
    params:
        tem2 = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["tem_tem_","{sample}", "_", "{repeat}", "_normalized.bdg"])) 
    shell:
        """
        mappedReads=$(awk 'NR==7 {{print $1}}' {input.alignment})
        scale=$(bc -l <<< "scale=3; 1000000 / $mappedReads")
        echo " Normalizing {input.pileup} with factor $scale"
        macs2 bdgopt -i {input.pileup} -m multiply -p $scale -o {params.tem2}
        sed -n '2,${{p}}' {params.tem2} > {output.tem1}
        rm {params.tem2}
        """

rule bedGraphToBigWig:
    input:
        tem_bdg = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["tem_","{sample}", "_", "{repeat}", "_normalized.bdg"]))  ,
        chromsize =  os.path.join(config["gfolder"], "chrom.sizes")
    output:
        bdg = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_normalized.bdg"])) ,
        bw = os.path.join(outputdir,"{sample}", "{repeat}", "macs2", "".join(["{sample}", "_", "{repeat}", "_normalized.bw"]))
    params:
        bedSort = config["bedSort"],
        bedGraphToBigWig = config["bedGraphToBigWig"]
    shell:
        """
        {params.bedSort} {input.tem_bdg} {output.bdg}
        {params.bedGraphToBigWig} {output.bdg} {input.chromsize} {output.bw} 
        """

rule call_peak:
    input:
        s = os.path.join(outputdir, "{sample}", "{repeat}", "hisat2", "".join(["{sample}", "_", "{repeat}", "_sorted.bam",])),
        ctrl = os.path.join(outputdir, Input, "{repeat}", "hisat2", "".join([Input, "_", "{repeat}","_sorted.bam"]))
    output:
        os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_peaks.narrowPeak"])),
        os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_peaks.xls"])),
        os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow", "".join(["{sample}_", "{repeat}", "_summits.bed"]))
    params:
        name = "".join(["{sample}_", "{repeat}"]),
        outdir = os.path.join(outputdir, "{sample}", "{repeat}", "macs2", "withCtrl_narrow"),
        genomesize = config["genomesize"],
        extsize = config["extsize"]
    shell:
        "macs2 callpeak -t {input.s} -c {input.ctrl} --name {params.name} --outdir {params.outdir} --nomodel -g {params.genomesize} --extsize {params.extsize} -B --SPMR"


rule chipr:
    input:
        [os.path.join(outputdir, "{sample}", repeat, "macs2", "withCtrl_narrow", "{sample}_%s_peaks.narrowPeak"%repeat) for repeat in repeats]
    output:
        os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_all.bed"])),
        os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_optimal.bed"]))
    params:
        o =  os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr"])),
        m = config["m"]
    shell:
        "chipr -i {input} -m {params.m} -o {params.o}"

rule correlation:
    input:
        peak = os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_all.bed"])),
        bdg = [os.path.join(outputdir,"{sample}", repeat, "macs2", "".join(["{sample}", "_", repeat, "_normalized.bdg"])) for repeat in repeats]
    output:
        os.path.join(outputdir,"{sample}","correlation","".join(["{sample}_correlation_data.xlsx"]))
    params:
        o = os.path.join(outputdir,"{sample}","correlation","".join(["{sample}"]))
    shell:
        "python correlation_V2.py -peak_file {input.peak} -o {params.o} -bdg_files {input.bdg}"

rule refind_peak_summit:
    input:
        chipr_all = os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}","_chipr_all.bed"])),
        normalized_bdg = os.path.join(outputdir,"{sample}", repeats[0], "macs2", "".join(["{sample}", "_", repeats[0], "_normalized.bdg"]))
    output:
        os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all.narrowPeak"])),
        os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all_plus_summit_info.bed"])),
        os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), ".narrowPeak"]))
    params:
        h_t = config["h_t"],
        o = os.path.join(outputdir, "{sample}", "chipr", "".join(["{sample}", "_chipr_all"]))
    shell:
        "python refind_peaksummit_v1.2.py -peak_file {input.chipr_all} -bdg_file {input.normalized_bdg} -h_t {params.h_t} -o {params.o}"


rule peak_annotation:
    input:
        peak = os.path.join(outputdir, "{sample}", "chipr","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), ".narrowPeak"])),
        gff = gff_file
    output:
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene.xlsx"])),
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene.txt"])),
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene_KEGG.txt"]))
    params:
        o = os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene"])),
        type1_method = config["type1_method"],
        distance = config["distance"]
    shell:
        "python chipseq_gene_annotation_V1.2.py -peak_file {input.peak} -gff {input.gff} -o {params.o} -distance {params.distance} -type1_method {params.type1_method}"


rule GO_analysis:
    input:
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene.txt"]))
    output:
        os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result.xlsx"]))
    params:
        o = os.path.join(outputdir, "{sample}", "GO","".join(["{sample}", "_GO_result"])),
        gaf = config["gaf"],
        obo = config["obo"]
    shell:
        "python GO_analysis_V1.py -target_gene_file {input} -gaf {params.gaf} -obo {params.obo} -o {params.o}"

rule KEGG_analysis:
    input:
        os.path.join(outputdir, "{sample}", "target_gene","".join(["{sample}", "_chipr_all_high_quality_peaks_threshold_", str(config["h_t"]), "_target_gene_KEGG.txt"]))
    output:
        os.path.join(outputdir, "{sample}", "KEGG","".join(["{sample}", "_KEGG_result.xlsx"]))
    params:
        o = os.path.join(outputdir, "{sample}", "KEGG","".join(["{sample}", "_KEGG_result"])),
        KEGG_relation = config["KEGG_relation"]
    shell:
        "python KEGG_analysis_V1.py -target_gene_file {input}  -KEGG_relation {params.KEGG_relation} -o {params.o}"
    



