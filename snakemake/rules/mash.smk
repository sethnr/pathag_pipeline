import sys

ruleorder: mash_merge_calls > mash_call
localrules: mash_merge_calls

rule mash_index:
    input:
        fastas=expand(os.path.join(config['refsdir'],"{target}.fasta"),target=TARGETS),
    output:
        msh="results/mash/index_all.msh"
    log:
        err="logs/mash/all.mash.err"
    params:
        idx_script = os.path.join(os.getcwd(),"scripts/indexer.sh"),
        mashabs=os.path.join(os.getcwd(),"results/mash/index_all.msh"),
        genome_size = "11k",
        refdir=config["refsdir"],
        localfastas=expand("{target}.fasta",target=TARGETS)
    resources:
        partition="day",
        mem_mb="40G",
        cpus_per_task=4,
        runtime=300
    container: "docker://sethnr/pgcoe_processing:0.01"
    shell:"""
	cd {params.refdir}
        echo {params.idx_script} \
                -m {params.mashabs} -g {params.genome_size} \
                {input.fastas}
        {params.idx_script} \
                -m {params.mashabs} -g {params.genome_size} \
                {params.localfastas}
        """

rule mash_call:
    input:
        msh = "results/mash/index_all.msh",
        R1 = 'results/rawdata/{sample}_R1.fastq.gz',
        R2 = 'results/rawdata/{sample}_R2.fastq.gz',
    output:
        mashcalls=temporary("results/mash/{sample}_calls.txt"),
        mashout=temporary("results/mash/{sample}_mash.txt"),
    resources:
        partition="day",
        mem_mb="4G",
        cpus_per_task=1,
        runtime=30
    container: "docker://sethnr/pgcoe_processing:0.01"
    group: "mashcall"
    params:
        reads=10000, # compare top N reads to refs
        bloom=10,    # bloom filter kmers with < N coverage (seq errors)
        gsize="11k", # estimated genome size (for prob assignment)
        prob=1e-50,  # max mash prob to call
        dist=0.25,   # max mash dist to call
        masher = "scripts/masher.sh",
        prefix="results/mash/{sample}",
        refdir=config["refsdir"],
    log:
        "logs/getstrain_{sample}.log",
    shell:
        """
        {params.masher} -f {input.msh} -s {wildcards.sample}\
                -r {params.reads} -b {params.bloom} -g {params.gsize} \
                -d {params.dist} -p {params.prob} \
        -o {params.prefix} \
        {input.R1} {input.R2}
        """

checkpoint mash_merge_calls:
    input:
        calls=expand("results/mash/{sample}_calls.txt",sample=SAMPLES),
        mash=expand("results/mash/{sample}_mash.txt",sample=SAMPLES)
    output:
        mashout="results/mash/allmash.txt",
        mashcalls="results/mash/allcalls.txt",
    log:
        "logs/mash/all-calls.err"
    shell:
        """
        cat {input.calls} 1> {output.mashcalls} 2> {log}
        cat {input.mash} 1> {output.mashout} 2>> {log}
        """




def get_mash_targets(wildcards):
    #print("mash: requesting targets for "+wildcards.sample,file=sys.stderr)  
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        mashlines = f.read().splitlines()
    mytargets = [L.split()[1] for L in mashlines if L.split()[0]==wildcards.sample]
    #print("mash: returning "+":".join(mytargets)+" from "+wildcards.sample,file=sys.stderr)  
    return mytargets

def get_mash_samples(wildcards):
    #print("mash: requesting samples for "+wildcards.target,file=sys.stderr)  
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        mashlines = f.read().splitlines()
    mysamples = [L.split()[0] for L in mashlines if L.split()[1]==wildcards.target]
    #print("mash: returning "+";".join(mysamples)+" from "+wildcards.target,file=sys.stderr)  
    return mysamples

def get_valid_targets():
    with checkpoints.mash_merge_calls.get().output.mashcalls.open() as f:
        alltargets = [L.split()[1] for L in f.read().splitlines()]
	alltargets = [I for I in set(alltargets)]
    #print("mash: returning valid target "+";".join(alltargets),file=sys.stderr)  
    return alltargets

def get_mash_bams(wildcards):
    mytargets = get_mash_targets(wildcards)
    filenames = expand("results/bams/{sample}_{target}.bam",sample=wildcards.sample,target=mytargets)
    #print("mash: returning "+";".join(filenames)+" from "+wildcards.sample,file=sys.stderr)  
    return filenames


