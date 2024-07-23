
rule ivar_pipeline:
    input:
        consensus='results/ivar/{sample}_{target}_consensus.fa',
        ivariants='results/ivar/{sample}_{target}_ivariants.tsv'


rule ivar_primerclip:
    input:
        aligned = 'results/align/{sample}_{target}.bam',
        primers = os.path.join(config['ampsdir']+'{target}.bed'),
    output:
        trimmed_us=temporary('results/ivar/{sample}_{target}_itrim_us.bam'),
        trimmed='results/ivar/{sample}_{target}_itrim.bam',
	trimdex='results/ivar/{sample}_{target}_itrim.bam.bai'
    resources:
        mem_mb=8000,
        runtime=1440,
	cores=1,
    log:
        stderr="logs/ivar/{sample}_{target}_trim.err"
    shell:
        """
        ivar trim -i {input.aligned} -b {input.primers} -p {output.trimmed_us} -e 2>&1 > {log.stderr}
        samtools sort -@ {resources.cores} -o {output.trimmed} {output.trimmed_us} 2>&1 >> {log.stderr}
        samtools index {output.trimmed} 2>&1 >> {log.stderr}
        """


rule sam_pileup:
    input:
        tocall='results/ivar/{sample}_{target}_itrim.bam',
        indexed='results/ivar/{sample}_{target}_itrim.bam.bai'
    output:
        pileup=temporary('results/ivar/{sample}_{target}.pileup'),
    params:
        ref=os.path.join(config['refsdir'],"{target}.fasta"),
        threshold=0.2,
        maxdepth=10000,
    resources:
        mem_mb=8000,
        runtime=600,
    log:
        stderr="logs/ivar/pileup_{sample}_{target}.err"
    shell:
        """
        samtools mpileup -aa -A -Q 0 -d {params.maxdepth} -f {params.ref} \
                         -o {output.pileup} {input.tocall} 2>&1 > {log.stderr}
        """

rule ivar_consensus:
    input:
        pileup='results/ivar/{sample}_{target}.pileup'
    output:
        consensus='results/ivar/{sample}_{target}_consensus.fa'
    params:
        depth=20,
	    threshold=0.75,
        prefix='results/ivar/{sample}_{target}_consensus'
    resources:
        mem_mb=lambda wc, input: max(10 * input.size_mb, 4000),
        runtime=600,
    log:
        stderr="logs/ivar/consensus_{sample}_{target}.err",
    shell:
        """
        cat {input.pileup} | ivar consensus -t {params.threshold} \
           -m {params.depth} -p {params.prefix} -i {wildcards.sample} \
           2>&1 > {log.stderr}
        """

rule ivar_variants:
    input:
        pileup='results/ivar/{sample}_{target}.pileup'
    output:
        ivariants='results/ivar/{sample}_{target}_ivariants.tsv',
    params:
        ref=config['refsdir']+"/{target}.fa",
        threshold=0.2,
        depth=20,
        qual=20,
        prefix='results/ivar/{sample}_{target}_ivariants'
    resources:
        mem_mb=lambda wc, input: max(20 * input.size_mb, 4000),
        runtime=240,
    log:
        stderr="logs/ivar/{sample}_{target}_consensus.err"
    shell:
        """
        cat {input.pileup} | \
         ivar variants -q {params.qual} -r {params.ref} -t {params.threshold} \
           -m {params.depth} -p {params.prefix}  2>&1 >  {log.stderr}        
        """


