### STATS ###
rule depth:
    input:
        bams='results/align/{sample}_{target}.bam',
    output:
        depth='results/align/{sample}_{target}_depth.txt',
    resources:
        mem_mb=8000,
        runtime=120,
    group: "statsgroup"
    log:
        stderr="logs/depth/{sample}_{target}.err"
    container: "docker://sethnr/pgcoe_processing:0.01"
    shell:
        """
        samtools depth -a -H {input.bams} -o {output.depth} 2>&1 >  {log.stderr}
        """



rule depthist:
    input:
        depth='results/align/{sample}_{target}_depth.txt',
        subfactor = 'results/align/{sample}_{target}.substat'
    output:
        dwins=temporary('results/align/{sample}_{target}_depthwins.txt'),
        dhist=temporary('results/align/{sample}_{target}_depthhist.txt'),
    resources:
        mem_mb=8000,
        runtime=120,
    group: "statsgroup"
    params:
        maxdepth=0,
        minmapqual=60,
        minbasequal=13,
        winsize=10,
        prefix='results/align/{sample}_{target}'
    log:
        stderr="logs/depth/{sample}_{target}.err"
    container: "docker://sethnr/pgcoe_analysis:0.01"
    shell:
        """
        python scripts/get_depth_distribution.py -d {input.depth} \
            -s {wildcards.sample} -F `cat {input.subfactor} | cut -f 2,2` \
            -w {params.winsize} -o {params.prefix}
        """


rule alignstats:
    input: 
        subfactor = 'results/align/{sample}_{target}.substat',
        flagstats = 'results/align/{sample}_{target}.flagstats',
        dhist='results/align/{sample}_{target}_depthhist.txt',
    output:
        stats='results/align/{sample}_{target}_alignstats.txt',
    group: "statsgroup"
    params:
        mindepth=10
    run:
        sample = wildcards.sample
        target = wildcards.target

        #get subsampling factor
        with open(input.subfactor, "r") as f:
            l = f.read().split()
            subfact = l[1] 
        f.close()

        #get reads aligned
        reads = -1
        aligned = -1
        paired = -1
        with open(input.flagstats, "r") as f:
            for l in f:
                l = l.strip().split('\t')
                passreads = l[0]
                failreads = l[1]
                stat = l[2]
                if stat == 'total (QC-passed reads + QC-failed reads)': 
                    reads = int(passreads)
                elif stat == "mapped": 
                    aligned = int(passreads)
                elif stat == "properly paired": 
                    paired = int(passreads)
        f.close()

        #get coverage / depth
        goodcov = 0
        cov = 0
        gsize = 0
        dtotal = 0
        with open(input.dhist, "r") as f:
            for l in f:
                l = l.split("\t")
                depth = int(l[1])
                count = int(l[2])
                cdepth = int(l[3])
                dtotal += cdepth * count
                if cdepth >= params.mindepth:
                    gsize += count
                    cov += count
                    goodcov += count
                elif cdepth > 0:
                    gsize += count
                    cov += count
                elif cdepth == 0:
                    gsize += count
        f.close()
        meandepth = round(dtotal/gsize)


        #subprocess.run(["samtools","index",output.subsamp])
        f = open(output.stats, "w")
        print("\t".join(map(str,[sample, target, subfact,
	    reads, aligned, paired,
	    meandepth, goodcov, cov, gsize])),file=f)
        f.close()




rule depth_amplicons:
    input:
        bams='results/align/{sample}_{target}.bam',
        subfactor = 'results/align/{sample}_{target}.substat',
        depth='results/align/{sample}_{target}_depth.txt',
    output:
        depth=temporary('results/align/{sample}_{target}_ampdepth.txt'),
    resources:
        mem_mb=8000,
        runtime=60,
    group: "statsgroup"
    params:
        ampbed=os.path.join(config['ampsdir'],'{target}_amplicons.bed'),
    log:
        stderr="logs/depth/{sample}_{target}.err"
    container: "docker://sethnr/pgcoe_analysis:0.01"
    shell:
        """
        python scripts/get_depth_windows.py -d {input.depth} \
            -s {wildcards.sample} -F `cat {input.subfactor} | cut -f 2,2` \
            -b {params.ampbed} -o {output.depth}
        """


rule depth_genes:
    input:
        bams='results/align/{sample}_{target}.bam',
        subfactor = 'results/align/{sample}_{target}.substat',
        depth='results/align/{sample}_{target}_depth.txt',
    output:
        depth=temporary('results/align/{sample}_{target}_genedepth.txt'),
    resources:
        mem_mb=8000,
        runtime=60,
    group: "statsgroup"
    params:
        bed=os.path.join(config['ampsdir'],'{target}_genes.bed'),
    log:
        stderr="logs/depth/{sample}_{target}.err"
    container: "docker://sethnr/pgcoe_analysis:0.01"
    shell:
        """
        python scripts/get_depth_windows.py -d {input.depth} \
            -s {wildcards.sample} -F `cat {input.subfactor} | cut -f 2,2` \
            -b {params.bed} -o {output.depth}

        """





rule div_amplicons:
    input:
        vcf='results/bcftools/{target}_all.vcf.gz',
    output:
        div='results/summary/{target}_ampdiv.txt',
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        ampbed=os.path.join(config['ampsdir'],'{target}_amplicons.bed'),
    container: "docker://sethnr/pgcoe_analysis:0.01"
    group: "statsgroup"
    shell:
        """
        python scripts/get_refdist_distribution.py -v {input.vcf} \
            -b {params.ampbed} -o {output.div}
        """

rule div_genes:
    input:
        vcf='results/bcftools/{target}_all.vcf.gz',
    output:
        div='results/summary/{target}_genediv.txt',
    resources:
        mem_mb=8000,
        runtime=60,
    params:
        ampbed=os.path.join(config['ampsdir'],'{target}_genes.bed'),
    container: "docker://sethnr/pgcoe_analysis:0.01"
    group: "statsgroup"
    shell:
        """
        python scripts/get_refdist_distribution.py -v {input.vcf} \
            -b {params.ampbed} -o {output.div}
        """



rule catstats:
    input:
        dhists=lambda wildcards: expand('results/align/{sample}_{{target}}_depthhist.txt',sample=get_mash_samples(wildcards)),
        dwins=lambda wildcards: expand('results/align/{sample}_{{target}}_depthwins.txt',sample=get_mash_samples(wildcards)),
        alstats=lambda wildcards: expand('results/align/{sample}_{{target}}_alignstats.txt',sample=get_mash_samples(wildcards)),
        ampdepths=lambda wildcards: expand('results/align/{sample}_{{target}}_ampdepth.txt',sample=get_mash_samples(wildcards)),
        genedepths=lambda wildcards: expand('results/align/{sample}_{{target}}_genedepth.txt',sample=get_mash_samples(wildcards)),
        #prmdepths=lambda wildcards: expand('results/align/{sample}_{{target}}_prmdepth.txt',sample=get_mash_samples(wildcards)),        
    output:
        dhists=    'results/summary/{target}_depthhists.txt',
        dwins=     'results/summary/{target}_depthwins.txt',
        alstats=    'results/summary/{target}_alignstats.txt',
        ampdepths= 'results/summary/{target}_ampdepth.txt',
        genedepths='results/summary/{target}_genedepth.txt',
        #prmdepths='results/summary/{target}_prmdepth.txt',
    resources:
        mem_mb=8000,
        runtime=180,
    log:
        stderr="logs/depth/{target}_depths.err"
    shell:
        """
        cat results/align/*_{wildcards.target}_depthhist.txt > {output.dhists}
        cat results/align/*_{wildcards.target}_depthwins.txt > {output.dwins}
        cat results/align/*_{wildcards.target}_alignstats.txt > {output.alstats}
        cat results/align/*_{wildcards.target}_ampdepth.txt > {output.ampdepths}
        cat results/align/*_{wildcards.target}_genedepth.txt > {output.genedepths}
        """

rule runstats:
    input:
        dhists='results/summary/{target}_depthhists.txt',
        dwins='results/summary/{target}_depthwins.txt',
        dstats='results/summary/{target}_alignstats.txt',
        ampdepths='results/summary/{target}_ampdepth.txt',
        ampdiv='results/summary/{target}_ampdiv.txt',        
        genedepths='results/summary/{target}_genedepth.txt',
        genediv='results/summary/{target}_genediv.txt',                
    output:
        done=temporary('results/summary/{target}_stats'),
    shell:
        """
        touch {output.done}
        """

