import itertools
import os
import sys
import gzip
import glob
import subprocess
from collections import defaultdict
from functools import partial
import fileinput
import pandas as pd
from Bio import SeqIO
import concurrent.futures
from scipy import stats
import numpy as np



configfile: "./sra_config.yaml" #"/storage/home/hcoda1/6/rridley3/scratch/subset_2/config.yaml"
localrules: all, get_sra_output

group_range = config['group_range']
group_list = ['G' + str(x).zfill(2) for x in list(range(int(group_range[0]),int(group_range[1])+1))]

out_name = config['out_dir']
out_name = os.getcwd() + '/' + out_name

data_file = config["data-file"]
df = pd.read_parquet(data_file).set_index("SampleID", drop=False).astype('str')
working_df = df.query("map_group in @group_list")

SAMPLES = list(working_df['SampleID'])
GROUPS = list(working_df['map_group'])
GROUPS_uniq = list(set(dict(zip(working_df.SampleID,working_df.map_group)).values()))
GROUP_IDS = list(set(working_df['map_group']))
temp_scratch_dir= str(config['scratch_dir']).strip() + '/smk_temp'
samples_in_group = defaultdict(list)
for item_s in list(zip(working_df.SampleID,working_df.map_group)):
    samples_in_group[item_s[1]].append(item_s[0])

map_groups = ['01_set','02_unbin'] #'04_unbin'

def get_threads(level='high'):
    return int(config['threads'][level])
def get_mem_cpu(level):
    return int(config['mem_per_cpu'][level])
def get_mem_mb(mem_level,thread_level='high'):
    return int(config['mem_per_cpu'][mem_level] * config['threads'][thread_level] * 1000)

def get_mem1_cpu(wildcards, attempt):
    level_dict =  { 1 : 'med', 2 : 'high', 3 : 'high_extra'}
    return get_mem_cpu(level_dict[attempt])
def get_mem1_mb(wildcards, attempt):
    level_dict =  { 1 : 'med', 2 : 'high', 3 : 'high_extra'}
    return get_mem_mb(mem_level=level_dict[attempt])

def get_mem2_cpu(wildcards, attempt):
    level_dict =  { 1 : 'high', 2 : 'high_extra', 3 : 'high_top'}
    return get_mem_cpu(level_dict[attempt])
def get_mem2_mb(wildcards, attempt):
    level_dict =  { 1 : 'high', 2 : 'high_extra', 3 : 'high_top'}
    return get_mem_mb(mem_level=level_dict[attempt])

def get_time_level(wildcards,attempt):
    time_level =  { 1 : '24:00:00', 2 : '36:00:00', 3 : '48:00:00'}
    return time_level[attempt]

rule all:
    input:
        #setup_sra = ancient(expand(f"{out_name}/touchfiles/{{group}}/{{sample}}_sra.touch",zip,sample=SAMPLES,
        #                            group=GROUPS)),
        #map_genomes = ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_{{map_group}}_map_done.touch",
        #                            group=GROUP_IDS,
        #                            map_group=map_groups)),
        #filter_files= ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_{{map_group}}_filter_done.touch",
        #                            group=GROUP_IDS,
        #                            map_group=map_groups)),
        #gff_abund = ancient(expand( f"{out_name}/touchfiles/{{group}}/parse_bedtools_{{map_group}}_done.touch",
        #                            group=GROUP_IDS,
        #                            map_group=map_groups)),
        map_genes = ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_set_genes_map_done.touch",
                                    group=GROUP_IDS)),
        tad_genes= ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_set_genes_tad_done.touch",
                                    group=GROUP_IDS)),
        map_unbin_genes= ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_unbin_genes_map_done.touch",
                                    group=GROUP_IDS)),
        unbin_tad_genes= ancient(expand(f"{out_name}/touchfiles/{{group}}/coverm_unbin_genes_tad_done.touch",
                                    group=GROUP_IDS)),
        #sample_geqs = ancient(expand(f"{out_name}/04_diversity/{{group}}/{{sample}}.geq.txt",zip,sample=SAMPLES,
        #                            group=GROUPS)),

    run:
        print(f"Run complete for groups {str(group_range)}")

rule setup_files:
    input:
        sra_file = lambda wildcards: working_df.loc[wildcards.sample,'sra_loc']
    output:
        setup_done = touch(f"{out_name}/touchfiles/{{group}}/{{sample}}_sra.touch"),
    params:
        bin = config["fasterq-dump"],
        out_dir = lambda wildcards: f"{out_name}/00_raw_reads/{wildcards.group}",
        run_acc = lambda wildcards: df.loc[wildcards.sample,'Run'],
        read1= lambda wildcards: f"{out_name}/00_raw_reads/{wildcards.group}/{wildcards.sample}_1.fq.gz",
        read2= lambda wildcards: f"{out_name}/00_raw_reads/{wildcards.group}/{wildcards.sample}_2.fq.gz",
        temp_dir = lambda  wildcards: f"{temp_scratch_dir}/sra_temp/{wildcards.group}"
    threads: get_threads('med')
    resources:
        mem_per_cpu = get_mem_cpu('med'),
        mem_mb = get_mem_mb(mem_level='med',thread_level='med')
    shell:"""
    mkdir -p {params.out_dir} {params.temp_dir}
    
    {params.bin} {input.sra_file} -O {params.out_dir}  --qual-defline '+' -e {threads} -t {params.temp_dir}
    pigz {params.out_dir}/{params.run_acc}*fastq
    mv {params.out_dir}/{params.run_acc}_1.fastq.gz {params.read1}
    mv {params.out_dir}/{params.run_acc}_2.fastq.gz {params.read2}
    rm {params.out_dir}/{params.run_acc}.fastq.gz || true  &> /dev/null
    
    """

rule fastp:
    input:
        setup_done = ancient(rules.setup_files.output.setup_done),
    output:
        trim1=f"{out_name}/01_trimmed/{{group}}/{{sample}}-t_1.fq.gz",
        trim2=f"{out_name}/01_trimmed/{{group}}/{{sample}}-t_2.fq.gz",
        html=f"{out_name}/reports/trim/{{group}}_{{sample}}-fastp.html",
        json=f"{out_name}/reports/trim/{{group}}_{{sample}}-fastp.json",
        trim_done= touch(f"{out_name}/touchfiles/{{group}}/{{sample}}_trim.touch")
    threads: get_threads('med')
    conda: "analysis_db"
    resources:
        time = '12:00:00',
        mem_per_cpu = get_mem_cpu('med'),
        mem_mb = get_mem_mb(mem_level='med',thread_level='med')
    log:
        f"{out_name}/reports/trim/{{group}}_{{sample}}-fastp.log"
    params:
        add_params = config['fastp']['add_params'],
        read1 = rules.setup_files.params.read1,
        read2 = rules.setup_files.params.read2
    shell:"""
    fastp -i {params.read1} -I {params.read2} -o {output.trim1} \
    -O {output.trim2} {params.add_params} -h {output.html} -j {output.json} &> {log}
    rm {params.read1} {params.read2}
    """

rule map_p1_set_genomes:
    input:
        read1 = lambda wildcards: ancient(expand(f"{out_name}/01_trimmed/{wildcards.group}/{{sample}}-t_1.fq.gz",
                                        sample=samples_in_group[wildcards.group])),
        read2 = lambda wildcards: ancient(expand(f"{out_name}/01_trimmed/{wildcards.group}/{{sample}}-t_2.fq.gz",
                                        sample=samples_in_group[wildcards.group])),
    output:
        out_cov = f"{out_name}/02_genome_analysis/00_coverage/01_set_tad80-{{group}}.tsv",
        bam_files = directory(f"{out_name}/02_genome_analysis/01_set/01_bam_files/{{group}}"),
        coverm_touch_done = touch(f"{out_name}/touchfiles/{{group}}/coverm_01_set_map_done.touch")
    log: f"{out_name}/reports/map/coverm_{{group}}_p1_set_genomes_map.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time = get_time_level,
        mem_per_cpu= get_mem2_cpu,
        mem_mb= get_mem2_mb,
    params:
        bins = config['bins']['01_set'],
        min_id = config['coverm']['min_id'],
        min_aln= config['coverm']['min_aln'],
        trim_min = config['coverm']['trim_min'],
        trim_max = config['coverm']['trim_max'],
        add_params = config['coverm']['add_params'],
        read_dir = lambda wildcards: f"{out_name}/01_trimmed/{wildcards.group}",
        tmpdir = lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/map"
    shell:"""
    mkdir -p {params.tmpdir} {output.bam_files}

    TMPDIR={params.tmpdir} coverm genome -m trimmed_mean -r {params.bins} \
    --output-file {output.out_cov} -p bwa-mem2 --trim-max {params.trim_max} --trim-min \
    {params.trim_min} --separator '~' --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    --bam-file-cache-directory {output.bam_files} -1 {params.read_dir}/*_1.fq.gz -2 {params.read_dir}/*_2.fq.gz | tee -a {log}
    """

rule filter_genomes:
    input:
        BAM_DIR=ancient(f"{out_name}/02_genome_analysis/{{map_group}}/01_bam_files/{{group}}"),
        coverm_touch_done=ancient(f"{out_name}/touchfiles/{{group}}/coverm_{{map_group}}_map_done.touch")
    output:
        MP_DIR=directory(f"{out_name}/02_genome_analysis/{{map_group}}/02_mapped/{{group}}"),
        UMP_DIR=directory(f"{out_name}/02_genome_analysis/{{map_group}}/03_unmapped/{{group}}"),
        FQ_DIR=directory(f"{out_name}/02_genome_analysis/{{map_group}}/04_ump_fq/{{group}}"),
        touch_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_{{map_group}}_filter_done.touch")
    log: f"{out_name}/reports/filter/coverm_{{map_group}}_{{group}}_filter.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        add_params=config['coverm']['add_params'],
        read_dir=f"{out_name}/01_trimmed/",
        tmpdir=lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.map_group}/{wildcards.group}/filter"

    shell: """        
        WORKDIR=$PWD
        mkdir -p {output.MP_DIR} {output.UMP_DIR} {output.FQ_DIR} {params.tmpdir}

        cd {input.BAM_DIR}
        for x in * ; do touch {output.MP_DIR}/$x  {output.UMP_DIR}/$x ; done
        cd $WORKDIR
        rename .gz .mp {output.MP_DIR}/*
        rename .gz .ump {output.UMP_DIR}*


        TMPDIR={params.tmpdir} coverm filter --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
        -t {threads}  -b {input.BAM_DIR}/* -o {output.MP_DIR}/* | tee -a {log}

        if [[ {wildcards.map_group} == '02_unbin' ]]
        then
        rm {input.BAM_DIR}/*
        rm {output.UMP_DIR}/*
        exit 0
        fi

        TMPDIR={params.tmpdir} coverm filter --inverse --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
        -t {threads} -b {input.BAM_DIR}/* -o {output.UMP_DIR}/* | tee -a {log}
        
        
        cd {output.UMP_DIR}

        for file in *.bam
        do
        OUT_FILE=`echo $file | cut -f2-3 -d. | sed 's|-t|-ump|g'`
        OUT_FILE2=`echo $OUT_FILE | sed 's|ump_1|ump_2|g'`
        TEMPFILE={params.tmpdir}/${{OUT_FILE}}_TEMPBAM
        samtools collate -O -@ {threads} $file $TEMPFILE | samtools fastq -1 {output.FQ_DIR}/${{OUT_FILE}}.gz -2 {output.FQ_DIR}/${{OUT_FILE2}}.gz \
        -s /dev/null -@ {threads} -
        done
        
        rm {input.BAM_DIR}/*
        rm {output.UMP_DIR}/*

    """


rule map_p2_unbin:
    input:
        in_dir = ancient(f"{out_name}/02_genome_analysis/01_set/04_ump_fq/{{group}}"),
        touch_done=ancient(f"{out_name}/touchfiles/{{group}}/coverm_01_set_filter_done.touch")
    output:
        out_cov= f"{out_name}/02_genome_analysis/00_coverage/02_unbin_tad80-{{group}}.tsv",
        bam_files=directory(f"{out_name}/02_genome_analysis/02_unbin/01_bam_files/{{group}}"),
        coverm_touch_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_02_unbin_map_done.touch")
    log: f"{out_name}/reports/map/coverm_{{group}}_p2_unbin_map.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        bins=config['bins']['02_unbin'],
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        trim_min=config['coverm']['trim_min'],
        trim_max=config['coverm']['trim_max'],
        add_params=config['coverm']['add_params'],
        tmpdir = lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/map"
    shell: """
    mkdir -p {params.tmpdir}
    TMPDIR={params.tmpdir} coverm genome -m trimmed_mean -r {params.bins} --single-genome \
    --output-file {output.out_cov} -p bwa-mem2 --trim-max {params.trim_max} --trim-min \
    {params.trim_min}  --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} -t {threads} \
    --bam-file-cache-directory {output.bam_files} -1 {input.in_dir}/*_1.fq.gz -2 {input.in_dir}/*_2.fq.gz | tee -a {log}
    
    rm {input.in_dir}/*
    """


rule bedtools_depth:
    input:
        bam_files = ancient(f"{out_name}/02_genome_analysis/{{map_group}}/02_mapped/{{group}}"),
        filter_done = ancient(f"{out_name}/touchfiles/{{group}}/coverm_{{map_group}}_filter_done.touch")
    output:
        done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/bedtools_depth_{{map_group}}_done.touch")
    log: f"{out_name}/reports/filter/bedtools_{{group}}_{{map_group}}_genomes_genes.log"
    conda: "analysis_db"
    threads: get_threads('low')
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        gff_file = lambda wildcards: config['gff_files'][wildcards.map_group],
        depth_out_dir = lambda wildcards: f"{out_name}/02_genome_analysis/{wildcards.map_group}/05_gene_counts/{wildcards.group}",
    shell: """
    mkdir -p {params.depth_out_dir}
    cd {input.bam_files}

    for x in *.bam 
    do
    FILE=`echo $x | cut -f2 -d. | cut -f1 -d-`
    echo "Parsing coverage of file $x " | tee -a {log}
    bedtools coverage -a {params.gff_file} -b $x  -hist -header > {params.depth_out_dir}/${{FILE}}_bedt.tsv | tee -a {log}
    done
    """


def process_sample(input_group):
    print(input_group,flush=True)
    tsv_file,abund_dir = input_group
    seq_cov_dict = defaultdict(dict)
    sample = '_'.join(tsv_file.split('_')[0:2])
    with open(tsv_file,'r') as file:
        for line in file:
            vals = line.strip().split()
            loc = vals[0]
            if line.startswith('#'):
                continue
            depth_val = int(vals[-4])
            depth_len = int(vals[-3])
            try:
                id_val = line.split('ID=')[1]
                id_val = id_val.split(';')[0]
            except:
                continue
            seq_cov_dict[f"{loc}|{id_val}"][depth_val] = depth_len

    with open(f'{abund_dir}/{sample}_gene_abund.tsv','w') as out:
        out.write(f'file~contig|ID\t{sample}-tad70\t{sample}-tad80\t{sample}-tad90\t{sample}-mean\n')
        for key, count_dict in seq_cov_dict.items():
            numbers_list = np.arange(0,0)
            if list(count_dict.keys()) == [0]:
                out.write(f'{key}\t0\t0\t0\t0\n')
            else:
                for depth, count in count_dict.items():
                    numbers_list = np.append(numbers_list,np.repeat(depth,count))
                tad70 = stats.trim_mean(numbers_list,0.15)
                tad80 = stats.trim_mean(numbers_list,0.10)
                tad90 = stats.trim_mean(numbers_list,0.05)
                mean_t = stats.tmean(numbers_list)
                out.write(f'{key}\t{tad70}\t{tad80}\t{tad90}\t{mean_t}\n')
    return f'Completed for {tsv_file}'



rule parse_bedtools:
    input:
        depth_done_touch = ancient(rules.bedtools_depth.output.done_group_touch)
    output:
        done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_{{map_group}}_done.touch")
    log: f"{out_name}/reports/filter/bedtools_{{group}}_{{map_group}}_parse_genes.log"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem1_cpu,
        mem_mb=get_mem1_mb,
    params:
        in_dir = rules.bedtools_depth.params.depth_out_dir,
        abund_dir = lambda wildcards: f"{out_name}/02_genome_analysis/{wildcards.map_group}/06_gene_abund/{wildcards.group}"
    run:
        if not sys.version_info >= (3, 7):
            raise OSError("python version must be 3.7 or higher for this script.")
        os.makedirs(params.abund_dir,exist_ok=True)

        os.chdir(params.in_dir)
        file_list = os.listdir()
        file_listr = [x for x in file_list if x.endswith('.tsv')]
        file_list = zip(file_listr,itertools.repeat(params.abund_dir))
        with concurrent.futures.ProcessPoolExecutor() as executor:
            results = executor.map(process_sample,file_list)
            for result in results:
                print(result)
        print("Completed parsing gff. zipping files...",flush=True)
        for file in file_listr:
            res = subprocess.run(f"lbzip2 {file}",shell=True)
            print(res,flush=True)
        os.chdir(params.abund_dir)
        file_list2 = os.listdir()
        for file in file_list2:
            res = subprocess.run(f"lbzip2 {file}",shell=True)
            print(res,flush=True)



rule map_genes:
    input:
        read1 = lambda wildcards: ancient(expand(f"{out_name}/01_trimmed/{wildcards.group}/{{sample}}-t_1.fq.gz",
                                        sample=samples_in_group[wildcards.group])),
        read2 = lambda wildcards: ancient(expand(f"{out_name}/01_trimmed/{wildcards.group}/{{sample}}-t_2.fq.gz",
                                        sample=samples_in_group[wildcards.group])),
    output:
        out_cov = f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_tad80-{{group}}.tsv",
        bam_files = directory(f"{out_name}/03_gene_analysis/01_set_genes/01_bam_files/{{group}}"),
        coverm_touch_done = touch(f"{out_name}/touchfiles/{{group}}/coverm_set_genes_map_done.touch"),

    log: f"{out_name}/reports/map/coverm_{{group}}_p1_set_genes_map.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        genes_ref=config['genes']['01_set'],
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        trim_min=config['coverm']['trim_min'],
        trim_max=config['coverm']['trim_max'],
        add_params=config['coverm']['add_params'],
        read_dir= lambda wildcards: f"{out_name}/01_trimmed/{wildcards.group}",
        tmpdir = lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/map_genes"
    shell: """
    mkdir -p {params.tmpdir}

    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean -r {params.genes_ref} \
    --output-file {output.out_cov} -p bwa-mem2 --trim-max {params.trim_max} --trim-min \
    {params.trim_min}  --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} -t {threads} \
    --bam-file-cache-directory {output.bam_files} -1 {params.read_dir}/*_1.fq.gz -2 {params.read_dir}/*_2.fq.gz | tee -a {log}
    """

rule filter_genes:
    input:
        BAM_DIR=ancient(rules.map_genes.output.bam_files),
        coverm_touch_done=ancient(rules.map_genes.output.bam_files)
    output:
        MP_DIR=directory(f"{out_name}/03_gene_analysis/01_set_genes/02_mapped/{{group}}"),
        UMP_DIR=directory(f"{out_name}/03_gene_analysis/01_set_genes/03_unmapped/{{group}}"),
        FQ_DIR=directory(f"{out_name}/03_gene_analysis/01_set_genes/04_ump_fq/{{group}}"),
        touch_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_set_genes_filter_done.touch")
    log: f"{out_name}/reports/map/coverm_{{group}}_p1_set_genes_filter.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        add_params=config['coverm']['add_params'],
        tmpdir=lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/filter_genes",

    shell: """
        WORKDIR=$PWD


        mkdir -p {output.MP_DIR} {output.UMP_DIR} {output.FQ_DIR} {params.tmpdir}

        cd {input.BAM_DIR}
        for x in * ; do touch {output.MP_DIR}/$x  {output.UMP_DIR}/$x ; done
        cd $WORKDIR
        rename .gz .mp {output.MP_DIR}/*
        rename .gz .ump {output.UMP_DIR}*

        TMPDIR={params.tmpdir} coverm filter --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
        -t {threads}  -b {input.BAM_DIR}/* -o {output.MP_DIR}/* | tee -a {log}

        TMPDIR={params.tmpdir} coverm filter --inverse --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
        -t {threads} -b {input.BAM_DIR}/* -o {output.UMP_DIR}/* | tee -a {log}
        
        

        cd {output.UMP_DIR}
        
    

        for file in *.bam
        do
        OUT_FILE=`echo $file | cut -f2-3 -d. | sed 's|-t|-ump|g'`
        OUT_FILE2=`echo $OUT_FILE | sed 's|ump_1|ump_2|g'`
        TEMPFILE={params.tmpdir}/${{OUT_FILE}}_TEMPBAM
        samtools collate -O -@ 24 $file $TEMPFILE | samtools fastq -1 {output.FQ_DIR}/${{OUT_FILE}}.gz \
        -2 {output.FQ_DIR}/${{OUT_FILE2}}.gz -s /dev/null -@ {threads} -
        done
        
        rm {input.BAM_DIR}/*
        rm {output.UMP_DIR}/*

    """

rule map_unbin_genes:
    input:
        read_dir = ancient(rules.filter_genes.output.FQ_DIR)
    output:
        out_cov = f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_tad80-{{group}}.tsv",
        bam_files = directory(f"{out_name}/03_gene_analysis/02_unbin_genes/01_bam_files/{{group}}"),
        coverm_touch_done = touch(f"{out_name}/touchfiles/{{group}}/coverm_unbin_genes_map_done.touch"),
    log: f"{out_name}/reports/map/coverm_{{group}}_p2_unbin_genes_map.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        genes_ref=config['genes']['02_unbin'],
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        trim_min=config['coverm']['trim_min'],
        trim_max=config['coverm']['trim_max'],
        add_params=config['coverm']['add_params'],
        tmpdir = lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/map_genes"
    shell: """
    mkdir -p {params.tmpdir}

    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean -r {params.genes_ref} \
    --output-file {output.out_cov} -p bwa-mem2 --trim-max {params.trim_max} --trim-min \
    {params.trim_min}  --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} -t {threads} \
    --bam-file-cache-directory {output.bam_files} -1 {input.read_dir}/*_1.fq.gz -2 {input.read_dir}/*_2.fq.gz | tee -a {log}
    
    
    rm {input.read_dir}/*
    """


rule get_tads_gene:
    input:
        bam_files = ancient(rules.filter_genes.output.MP_DIR),
        map_done = ancient(rules.filter_genes.output.touch_done)
    output:
        out_cov90 = f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_tad90-{{group}}.tsv",
        out_cov70 = f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_tad70-{{group}}.tsv",
        out_cov_mean= f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_mean-{{group}}.tsv",
        out_cov_count= f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_count-{{group}}.tsv",
        out_cov_other= f"{out_name}/03_gene_analysis/00_coverage/01_set_genes_oth-{{group}}.tsv",
        coverm_tads_done = touch(f"{out_name}/touchfiles/{{group}}/coverm_set_genes_tad_done.touch"),

    log: f"{out_name}/reports/map/coverm_{{group}}_p1_set_tad_genes.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        add_params=config['coverm']['add_params'],
        tmpdir = lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/tad_genes"
    shell: """
    mkdir -p {params.tmpdir}

    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean  \
    --output-file {output.out_cov90} --trim-max 95 --trim-min 5 \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    
    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean  \
    --output-file {output.out_cov70} --trim-max 85 --trim-min 15 \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    
    TMPDIR={params.tmpdir} coverm contig -m mean \
    --output-file {output.out_cov_mean}   \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    
    TMPDIR={params.tmpdir} coverm contig -m count \
    --output-file {output.out_cov_count}   \
     --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    
    TMPDIR={params.tmpdir} coverm contig -m length tpm rpkm \
    --output-file {output.out_cov_other}   --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    """

rule get_tads_unbin:
    input:
        bam_files=ancient(rules.map_unbin_genes.output.bam_files),
        map_done=ancient(rules.map_unbin_genes.output.coverm_touch_done)
    output:
        out_cov90=f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_tad90-{{group}}.tsv",
        out_cov70=f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_tad70-{{group}}.tsv",
        out_cov_mean=f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_mean-{{group}}.tsv",
        out_cov_count=f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_count-{{group}}.tsv",
        out_cov_other=f"{out_name}/03_gene_analysis/00_coverage/02_unbin_genes_oth-{{group}}.tsv",
        coverm_tads_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_unbin_genes_tad_done.touch"),

    log: f"{out_name}/reports/map/coverm_{{group}}_p2_unbin_tad_genes.log"
    conda: "coverm_bwa2"
    threads: get_threads()
    resources:
        time=get_time_level,
        mem_per_cpu=get_mem2_cpu,
        mem_mb=get_mem2_mb,
    params:
        min_id=config['coverm']['min_id'],
        min_aln=config['coverm']['min_aln'],
        add_params=config['coverm']['add_params'],
        tmpdir=lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/tad_unbin"
    shell: """
    mkdir -p {params.tmpdir}

    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean  \
    --output-file {output.out_cov90} --trim-max 95 --trim-min 5 \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}

    TMPDIR={params.tmpdir} coverm contig -m trimmed_mean  \
    --output-file {output.out_cov70} --trim-max 85 --trim-min 15 \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}

    TMPDIR={params.tmpdir} coverm contig -m mean \
    --output-file {output.out_cov_mean}   \
    --min-covered-fraction 0.1  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    
    TMPDIR={params.tmpdir} coverm contig -m count \
    --output-file {output.out_cov_count}   \
     --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}

    TMPDIR={params.tmpdir} coverm contig -m length rpkm \
    --output-file {output.out_cov_other}  --min-read-percent-identity {params.min_id} \
    --min-read-aligned-percent {params.min_aln} -t {threads} \
    -b {input.bam_files}/* | tee -a {log}
    """





rule microbe_census:
    conda: "microbe_census"
    input:
        read1 = ancient(rules.fastp.output.trim1),
        read2 = ancient(rules.fastp.output.trim2)
    output:
        out_file = f"{out_name}/04_diversity/{{group}}/{{sample}}.geq.txt"
    threads: get_threads()
    resources:
        time = '24:00:00',
        mem_per_cpu=get_mem_cpu('med_low'),
        mem_mb=get_mem_mb(mem_level='med_low'),
    params:
        n = config['microbe_census']['n'],
        q = config['microbe_census']['q']
    shell:
        """
        run_microbe_census.py -t {threads} -n {params.n} -q {params.q} {input.read1},{input.read2} {output.out_file} 
        """


## Single group rules

## GFF Gene abundances

# rule bedtools_depth_p1:
#     input:
#         bam_files = ancient(rules.filter_p1.output.MP_DIR)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/bedtools_abund_p1_done.touch")
#     log: f"{out_name}/reports/filter/bedtools_{{group}}_p1_set_genomes_genes.log"
#     conda: "analysis_db"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         gff_file = config['gff_files']['01_set'],
#         depth_out_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/05_gene_counts/{wildcards.group}",
#     shell: """
#     mkdir -p {params.depth_out_dir}
#     cd {input.bam_files}
#
#     for x in *.bam
#     do
#     FILE=`echo $x | cut -f2 -d. | cut -f1 -d-`
#     echo "Parsing coverage of file $x "
#     bedtools coverage -a {params.gff_file} -b $x  -hist -header > {params.depth_out_dir}/${{FILE}}_bedt.tsv
#     done
#     """
#
# rule parse_bedtools_p1:
#     input:
#         depth_done_touch = ancient(rules.bedtools_depth_p1.output.done_group_touch)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_p1_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_genes.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         abund_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/06_gene_abund/{wildcards.group}"
#     run:
#         if not sys.version_info >= (3, 7):
#             raise OSError("python version must be 3.7 or higher for this script.")
#         os.mkdirs(params.abund_dir,exist_ok=True)
#         main_dir = params.abund_dir
#         def process_sample(tsv_file):
#             seq_cov_dict = defaultdict(dict)
#             sample = '_'.join(tsv_file.split('_')[0:2])
#             with open(tsv_file,'r') as file:
#                 for line in file:
#                     vals = line.strip().split()
#                     loc = vals[0]
#                     depth_val = int(vals[-4])
#                     depth_len = int(vals[-3])
#                     try:
#                         id_val = line.split('ID=')[1]
#                         id_val = id_val.split(';')[0]
#                     except:
#                         continue
#                     seq_cov_dict[f"{loc}|{id_val}"][depth_val] = depth_len
#
#
#             with open(f'{main_dir}/{sample}_gene_abund.csv','w') as out:
#                 out.write(f'file~contig|ID,{sample}-tad70,{sample}-tad80,{sample}-tad90,{sample}-mean\n')
#                 for key, count_dict in seq_cov_dict.items():
#                     numbers_list = np.arange(0,0)
#                     if list(count_dict.keys()) == [0]:
#                         out.write(f'{key},0,0,0,0\n')
#                     else:
#                         for depth, count in count_dict.items():
#                             numbers_list = np.append(numbers_list,np.repeat(depth,count))
#                         tad70 = stats.trim_mean(numbers_list,0.15)
#                         tad80 = stats.trim_mean(numbers_list,0.10)
#                         tad90 = stats.trim_mean(numbers_list,0.05)
#                         mean_t = stats.tmean(numbers_list)
#                         out.write(f'{key},{tad70},{tad80},{tad90},{mean_t}\n')
#
#
#         os.chdir(f"{input.bam_files}]")
#         file_list = os.listdir()
#         file_list = [x for x in file_list if x.endswith('.tsv')]
#         with concurrent.futures.ProcessPoolExecutor() as executor:
#             executor.map(process_sample,file_list)
#
# rule bedtools_depth_p2:
#     input:
#         bam_files = ancient(rules.filter_p1.output.MP_DIR)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/bedtools_abund_p1_done.touch")
#     log: f"{out_name}/reports/filter/bedtools_{{group}}_p1_set_genomes_genes.log"
#     conda: "analysis_db"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         gff_file = config['gff_files']['01_set'],
#         depth_out_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/05_gene_counts/{wildcards.group}",
#     shell: """
#     mkdir -p {params.depth_out_dir}
#     cd {input.bam_files}
#
#     for x in *.bam
#     do
#     FILE=`echo $x | cut -f2 -d. | cut -f1 -d-`
#     echo "Parsing coverage of file $x "
#     bedtools coverage -a {params.gff_file} -b $x  -hist -header > {params.depth_out_dir}/${{FILE}}_bedt.tsv
#     done
#     """
#
# rule parse_bedtools_p2:
#     input:
#         depth_done_touch = ancient(rules.bedtools_depth_p1.output.done_group_touch)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_p1_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_genes.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         abund_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/06_gene_abund/{wildcards.group}"
#     run:
#         if not sys.version_info >= (3, 7):
#             raise OSError("python version must be 3.7 or higher for this script.")
#         os.mkdirs(params.abund_dir,exist_ok=True)
#         main_dir = params.abund_dir
#         def process_sample(tsv_file):
#             seq_cov_dict = defaultdict(dict)
#             sample = '_'.join(tsv_file.split('_')[0:2])
#             with open(tsv_file,'r') as file:
#                 for line in file:
#                     vals = line.strip().split()
#                     loc = vals[0]
#                     depth_val = int(vals[-4])
#                     depth_len = int(vals[-3])
#                     try:
#                         id_val = line.split('ID=')[1]
#                         id_val = id_val.split(';')[0]
#                     except:
#                         continue
#                     seq_cov_dict[f"{loc}|{id_val}"][depth_val] = depth_len
#
#
#             with open(f'{main_dir}/{sample}_gene_abund.csv','w') as out:
#                 out.write(f'file~contig|ID,{sample}-tad70,{sample}-tad80,{sample}-tad90,{sample}-mean\n')
#                 for key, count_dict in seq_cov_dict.items():
#                     numbers_list = np.arange(0,0)
#                     if list(count_dict.keys()) == [0]:
#                         out.write(f'{key},0,0,0,0\n')
#                     else:
#                         for depth, count in count_dict.items():
#                             numbers_list = np.append(numbers_list,np.repeat(depth,count))
#                         tad70 = stats.trim_mean(numbers_list,0.15)
#                         tad80 = stats.trim_mean(numbers_list,0.10)
#                         tad90 = stats.trim_mean(numbers_list,0.05)
#                         mean_t = stats.tmean(numbers_list)
#                         out.write(f'{key},{tad70},{tad80},{tad90},{mean_t}\n')
#
#
#         os.chdir(f"{input.bam_files}]")
#         file_list = os.listdir()
#         file_list = [x for x in file_list if x.endswith('.tsv')]
#         with concurrent.futures.ProcessPoolExecutor() as executor:
#             executor.map(process_sample,file_list)
#
# rule bedtools_depth_p3:
#     input:
#         bam_files = ancient(rules.filter_p1.output.MP_DIR)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/bedtools_abund_p1_done.touch")
#     log: f"{out_name}/reports/filter/bedtools_{{group}}_p1_set_genomes_genes.log"
#     conda: "analysis_db"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         gff_file = config['gff_files']['01_set'],
#         depth_out_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/05_gene_counts/{wildcards.group}",
#     shell: """
#     mkdir -p {params.depth_out_dir}
#     cd {input.bam_files}
#
#     for x in *.bam
#     do
#     FILE=`echo $x | cut -f2 -d. | cut -f1 -d-`
#     echo "Parsing coverage of file $x "
#     bedtools coverage -a {params.gff_file} -b $x  -hist -header > {params.depth_out_dir}/${{FILE}}_bedt.tsv
#     done
#     """
#
# rule parse_bedtools_p3:
#     input:
#         depth_done_touch = ancient(rules.bedtools_depth_p1.output.done_group_touch)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_p1_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_genes.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         abund_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/06_gene_abund/{wildcards.group}"
#     run:
#         if not sys.version_info >= (3, 7):
#             raise OSError("python version must be 3.7 or higher for this script.")
#         os.mkdirs(params.abund_dir,exist_ok=True)
#         main_dir = params.abund_dir
#         def process_sample(tsv_file):
#             seq_cov_dict = defaultdict(dict)
#             sample = '_'.join(tsv_file.split('_')[0:2])
#             with open(tsv_file,'r') as file:
#                 for line in file:
#                     vals = line.strip().split()
#                     loc = vals[0]
#                     depth_val = int(vals[-4])
#                     depth_len = int(vals[-3])
#                     try:
#                         id_val = line.split('ID=')[1]
#                         id_val = id_val.split(';')[0]
#                     except:
#                         continue
#                     seq_cov_dict[f"{loc}|{id_val}"][depth_val] = depth_len
#
#
#             with open(f'{main_dir}/{sample}_gene_abund.csv','w') as out:
#                 out.write(f'file~contig|ID,{sample}-tad70,{sample}-tad80,{sample}-tad90,{sample}-mean\n')
#                 for key, count_dict in seq_cov_dict.items():
#                     numbers_list = np.arange(0,0)
#                     if list(count_dict.keys()) == [0]:
#                         out.write(f'{key},0,0,0,0\n')
#                     else:
#                         for depth, count in count_dict.items():
#                             numbers_list = np.append(numbers_list,np.repeat(depth,count))
#                         tad70 = stats.trim_mean(numbers_list,0.15)
#                         tad80 = stats.trim_mean(numbers_list,0.10)
#                         tad90 = stats.trim_mean(numbers_list,0.05)
#                         mean_t = stats.tmean(numbers_list)
#                         out.write(f'{key},{tad70},{tad80},{tad90},{mean_t}\n')
#
#
#         os.chdir(f"{input.bam_files}]")
#         file_list = os.listdir()
#         file_list = [x for x in file_list if x.endswith('.tsv')]
#         with concurrent.futures.ProcessPoolExecutor() as executor:
#             executor.map(process_sample,file_list)
#
# rule bedtools_depth_p4:
#     input:
#         bam_files = ancient(rules.filter_p1.output.MP_DIR)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/bedtools_abund_p1_done.touch")
#     log: f"{out_name}/reports/filter/bedtools_{{group}}_p1_set_genomes_genes.log"
#     conda: "analysis_db"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         gff_file = config['gff_files']['01_set'],
#         depth_out_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/05_gene_counts/{wildcards.group}",
#     shell: """
#     mkdir -p {params.depth_out_dir}
#     cd {input.bam_files}
#
#     for x in *.bam
#     do
#     FILE=`echo $x | cut -f2 -d. | cut -f1 -d-`
#     echo "Parsing coverage of file $x "
#     bedtools coverage -a {params.gff_file} -b $x  -hist -header > {params.depth_out_dir}/${{FILE}}_bedt.tsv
#     done
#     """
#
# rule parse_bedtools_p4:
#     input:
#         depth_done_touch = ancient(rules.bedtools_depth_p1.output.done_group_touch)
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_p1_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_genes.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         abund_dir = lambda wildcards: f"{out_name}/02_genome_analysis/01_set_genomes/06_gene_abund/{wildcards.group}"
#     run:
#         if not sys.version_info >= (3, 7):
#             raise OSError("python version must be 3.7 or higher for this script.")
#         os.mkdirs(params.abund_dir,exist_ok=True)
#         main_dir = params.abund_dir
#         def process_sample(tsv_file):
#             seq_cov_dict = defaultdict(dict)
#             sample = '_'.join(tsv_file.split('_')[0:2])
#             with open(tsv_file,'r') as file:
#                 for line in file:
#                     vals = line.strip().split()
#                     loc = vals[0]
#                     depth_val = int(vals[-4])
#                     depth_len = int(vals[-3])
#                     try:
#                         id_val = line.split('ID=')[1]
#                         id_val = id_val.split(';')[0]
#                     except:
#                         continue
#                     seq_cov_dict[f"{loc}|{id_val}"][depth_val] = depth_len
#
#
#             with open(f'{main_dir}/{sample}_gene_abund.csv','w') as out:
#                 out.write(f'file~contig|ID,{sample}-tad70,{sample}-tad80,{sample}-tad90,{sample}-mean\n')
#                 for key, count_dict in seq_cov_dict.items():
#                     numbers_list = np.arange(0,0)
#                     if list(count_dict.keys()) == [0]:
#                         out.write(f'{key},0,0,0,0\n')
#                     else:
#                         for depth, count in count_dict.items():
#                             numbers_list = np.append(numbers_list,np.repeat(depth,count))
#                         tad70 = stats.trim_mean(numbers_list,0.15)
#                         tad80 = stats.trim_mean(numbers_list,0.10)
#                         tad90 = stats.trim_mean(numbers_list,0.05)
#                         mean_t = stats.tmean(numbers_list)
#                         out.write(f'{key},{tad70},{tad80},{tad90},{mean_t}\n')
#
#
#         os.chdir(f"{input.bam_files}]")
#         file_list = os.listdir()
#         file_list = [x for x in file_list if x.endswith('.tsv')]
#         with concurrent.futures.ProcessPoolExecutor() as executor:
#             executor.map(process_sample,file_list)

# rule join_tables_p1:
#     input:
#         bam_files = ancient(expand())
#     output:
#         done_group_touch = touch(f"{out_name}/touchfiles/{{group}}/parse_bedtools_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_genes.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         bam_gff_dir = directory(f"{out_name}/02_genome_analysis/01_set_genomes/05_gene_counts/{{group}}"),





rule cleanup_files:
    input:
        complete_genes = ancient(rules.map_unbin_genes.output.coverm_touch_done),
        complete_gffs = ancient(rules.parse_bedtools.output.done_group_touch),
        complete_genomes  = ancient(rules.map_p2_unbin.output.coverm_touch_done)
    output:
        cleanup_done_touch = touch(f"{out_name}/touchfiles/cleanup_{{map_group}}_done.touch")
    shell: """
    #find ./ -wholename '*/G{wildcards.map_group}/*.bam' -exec rm {{}} \;
    #find ./ -wholename '*/G{wildcards.map_group}/*.fq.gz' -exec rm {{}} \;
    
    
    
    """





# rule filter_p1:
#     input:
#         BAM_DIR=ancient(rules.map_p1_set_genomes.output.bam_files),
#         coverm_touch_done=ancient(rules.map_p1_set_genomes.output.coverm_touch_done)
#     output:
#         MP_DIR=directory(f"{out_name}/02_genome_analysis/01_set/02_mapped/{{group}}"),
#         UMP_DIR=directory(f"{out_name}/02_genome_analysis/01_set/03_unmapped/{{group}}"),
#         FQ_DIR=directory(f"{out_name}/02_genome_analysis/01_set/04_ump_fq/{{group}}"),
#         touch_done=touch(touch(f"{out_name}/touchfiles/{{group}}/coverm_set_filter_done.touch"))
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p1_set_genomes_filter.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         min_id=config['coverm']['min_id'],
#         min_aln=config['coverm']['min_aln'],
#         add_params=config['coverm']['add_params'],
#         read_dir=f"{out_name}/01_trimmed/",
#         tmpdir=lambda wildcards: f"{temp_scratch_dir}/coverm/{wildcards.group}/filter"
#
#     shell: """
#         WORKDIR=$PWD
#         mkdir -p {output.MP_DIR} {output.UMP_DIR} {output.FQ_DIR} {params.tmpdir}
#
#         cd {input.BAM_DIR}
#         for x in * ; do touch {output.MP_DIR}/$x  {output.UMP_DIR}/$x ; done
#         cd $WORKDIR
#         rename .gz .Smp {output.MP_DIR}/*
#         rename .gz .Sump {output.UMP_DIR}*
#
#
#         TMPDIR={params.tmpdir} coverm filter --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads}  -b {input.BAM_DIR}/* -o {output.MP_DIR}/* | tee -a {log}
#
#         TMPDIR={params.tmpdir} coverm filter --inverse --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads} -b {input.BAM_DIR}/* -o {output.UMP_DIR}/* | tee -a {log}
#
#         cd {output.UMP_DIR}
#
#         for file in *.bam
#         do
#         OUT_FILE=`echo $file | cut -f2-3 -d. | sed 's|-t|-ump|g'`
#         OUT_FILE2=`echo $OUT_FILE | sed 's|ump_1|ump_2|g'`
#         TEMPFILE={params.tmpdir}/${{OUT_FILE}}_TEMPBAM
#         samtools collate -O -@ {threads} $file $TEMPFILE | samtools fastq -1 {output.FQ_DIR}/${{OUT_FILE}}.gz -2 {output.FQ_DIR}/${{OUT_FILE2}}.gz \
#         -s /dev/null -@ {threads} -
#         done
#
#     """

# rule filter_p2:
#     input:
#         BAM_DIR=ancient(rules.map_p2_OD_genomes.output.bam_files),
#         coverm_touch_done=ancient(rules.map_p2_OD_genomes.output.coverm_touch_done)
#     output:
#         MP_DIR=directory(f"{out_name}/02_genome_analysis/02_OD/02_mapped/{{group}}"),
#         UMP_DIR=directory(f"{out_name}/02_genome_analysis/02_OD/03_unmapped/{{group}}"),
#         FQ_DIR=directory(f"{out_name}/02_genome_analysis/02_OD/04_ump_fq/{{group}}"),
#         touch_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_OD_filter_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p2_OD_genomes_filter.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         min_id=config['coverm']['min_id'],
#         min_aln=config['coverm']['min_aln'],
#         add_params=config['coverm']['add_params'],
#         tmpdir=f"{temp_scratch_dir}/OD_02/{group_id}",
#
#     shell: """
#         WORKDIR=$PWD
#         mkdir -p {output.MP_DIR} {output.UMP_DIR} {output.FQ_DIR} {params.tmpdir}
#
#         cd {input.BAM_DIR}
#         for x in * ; do touch {output.MP_DIR}/$x  {output.UMP_DIR}/$x ; done
#         cd $WORKDIR
#         rename .gz .Omp {output.MP_DIR}/*
#         rename .gz .Oump {output.UMP_DIR}*
#
#
#         TMPDIR={params.tmpdir} coverm filter --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads}  -b {input.BAM_DIR}/* -o {output.MP_DIR}/* | tee -a {log}
#
#         TMPDIR={params.tmpdir} coverm filter --inverse --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads} -b {input.BAM_DIR}/* -o {output.UMP_DIR}/* | tee -a {log}
#
#         cd {output.UMP_DIR}
#
#         for file in *.bam
#         do
#         OUT_FILE=`echo $file | cut -f2-3 -d. | sed 's|-t|-ump|g'`
#         OUT_FILE2=`echo $OUT_FILE | sed 's|ump_1|ump_2|g'`
#         TEMPFILE={params.tmpdir}/${{OUT_FILE}}_TEMPBAM
#         samtools collate -O -@ {threads} $file $TEMPFILE | samtools fastq -1 {output.FQ_DIR}/${{OUT_FILE}}.gz -2 {output.FQ_DIR}/${{OUT_FILE2}}.gz \
#         -s /dev/null -@ {threads} -
#         done
#     """

# rule filter_p3:
#     input:
#         BAM_DIR=ancient(rules.map_p3_cm2_genomes.output.bam_files),
#         coverm_touch_done=ancient(rules.map_p3_cm2_genomes.output.coverm_touch_done)
#     output:
#         MP_DIR=directory(f"{out_name}/02_genome_analysis/03_cm2/02_mapped/{{group}}"),
#         UMP_DIR=directory(f"{out_name}/02_genome_analysis/03_cm2/03_unmapped/{{group}}"),
#         FQ_DIR=directory(f"{out_name}/02_genome_analysis/03_cm2/04_ump_fq/{{group}}"),
#         touch_done=touch(f"{out_name}/touchfiles/{{group}}/coverm_cm2_filter_done.touch")
#     log: f"{out_name}/reports/filter/coverm_{{group}}_p3_cm2_genomes_filter.log"
#     conda: "coverm_bwa2"
#     threads: get_threads()
#     resources:
#         time='48:00:00',
#         mem_per_cpu=get_mem_cpu('med'),
#         mem_mb=get_mem_mb(mem_level='med'),
#     params:
#         min_id=config['coverm']['min_id'],
#         min_aln=config['coverm']['min_aln'],
#         add_params=config['coverm']['add_params'],
#         tmpdir=f"{temp_scratch_dir}/03_cm2/{group_id}",
#
#     shell: """
#         WORKDIR=$PWD
#         mkdir -p {output.MP_DIR} {output.UMP_DIR} {output.FQ_DIR} {params.tmpdir}
#
#         cd {input.BAM_DIR}
#         for x in * ; do touch {output.MP_DIR}/$x  {output.UMP_DIR}/$x ; done
#         cd $WORKDIR
#         rename .gz .Cmp {output.MP_DIR}/*
#         rename .gz .Cump {output.UMP_DIR}*
#
#
#         TMPDIR={params.tmpdir} coverm filter --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads}  -b {input.BAM_DIR}/* -o {output.MP_DIR}/* | tee -a {log}
#
#         TMPDIR={params.tmpdir} coverm filter --inverse --min-read-percent-identity {params.min_id} --min-read-aligned-percent {params.min_aln} \
#         -t {threads} -b {input.BAM_DIR}/* -o {output.UMP_DIR}/* | tee -a {log}
#
#         cd {output.UMP_DIR}
#
#         for file in *.bam
#         do
#         OUT_FILE=`echo $file | cut -f2-3 -d. | sed 's|-t|-ump|g'`
#         OUT_FILE2=`echo $OUT_FILE | sed 's|ump_1|ump_2|g'`
#         TEMPFILE={params.tmpdir}/${{OUT_FILE}}_TEMPBAM
#         samtools collate -O -@ {threads} $file $TEMPFILE | samtools fastq -1 {output.FQ_DIR}/${{OUT_FILE}}.gz -2 {output.FQ_DIR}/${{OUT_FILE2}}.gz \
#         -s /dev/null -@ {threads} -
#         done
#     """






rule get_sra_output:
    input:
        expand(f"{out_name}/{{group}}--{{sample}}_sra.touch",zip,sample=SAMPLES,group=GROUPS),
    output:
        "updated_metadata.csv",
        touch(f"{out_name}/download_complete.touch")
    run:
        with open(output[0], 'w') as out:
            for sample in SAMPLES:
                out_dir = os.path.dirname(working_df.loc[sample,'sra_file'])
                run_acc = working_df.loc[sample,'Run']
                #if df.loc[sample,'read_type'] == 'paired':
                f1 = str(glob.glob(f"{out_dir}/{run_acc}*.fastq.gz"))
                out.write(f"{f1},")
               # f2 = str(glob.glob(f"{out_dir}/{run_acc}*.fastq.gz"))
                #print(type(f1),f2)
                #df.loc[:,'File1'] = ' '
              #  df.loc[:,'File2'] = ' '
              #  df.loc[sample,'Read1'] = f1
              #  df.loc[sample, 'Read2'] = f2
              #  df.loc[sample,'sra_downloaded'] = True
            #elif df.loc[sample,'read_type'] == 'single':
             #   f1 = list(glob.glob(f"{out_dir}/{run_acc}*.fastq.gz"))
             #   df.loc[sample, 'Read1'] = f1
             #   df.loc[sample,'sra_downloaded'] = True
      #  df.to_csv(output[0],index=False)

