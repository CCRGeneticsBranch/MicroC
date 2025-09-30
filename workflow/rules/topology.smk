# Convert to HiCExplorer native format

rule convert_hic_h5:
    """From Hsien-Chao's hicToh5:sh: converts .hic files to .cool and .h5 formats (10kb only)"""
    input:
        hic = join(RESULTSDIR,"hic","{sample}.hic")
    output:
        cool = join(RESULTSDIR,"h5","{sample}.10kb.cool"),
        h5 = join(RESULTSDIR,"h5","{sample}.10kb.h5")
    params:
        res = "10000"
    threads:
        getthreads("convert_hic_h5")
    envmodules:
        TOOLS["hicexplorer"]
    shell:
        """
        hicConvertFormat --matrices {input.hic} --outFileName {output.cool} --inputFormat hic --outputFormat cool --resolutions {params.res}
        hicConvertFormat --matrices {output.cool} --outFileName {output.h5} --inputFormat cool --outputFormat h5
        """


# AB Compartments

rule call_compartment_fanc:
    """Call AB compartments using FAN-C (see https://fan-c.readthedocs.io/en/latest/fanc-executable/fanc-analyse-hic/ab_compartments.html)"""
    input:
        hic = join(RESULTSDIR,"hic","{sample}.hic")
    output:
        ab = join(RESULTSDIR,"topology","compartment","{sample}_10kb.ab"),
        eigen = join(RESULTSDIR,"topology","compartment","{sample}_10kb_eigen.bed"),
        domain = join(RESULTSDIR,"topology","compartment","{sample}_10kb_AB.bed"),
        profile = join(RESULTSDIR,"topology","compartment","{sample}_10kb_profile.png")
    params:
        ref = config[genome]['ref'],
        res = "10kb"
    envmodules:
        TOOLS["fanc"]
    shell:
        """
        fanc compartments -f -v {output.eigen} -d {output.domain} -e {output.profile} -g {params.ref} \
        {input.hic}@{params.res} {output.ab}
        """


# Topologically Associated Domains

rule call_tads_arrowhead:
    """Call TADs using Juicer Arrowhead (https://github.com/aidenlab/juicer/wiki/Arrowhead)"""
    input:
        hic = join(RESULTSDIR, "hic", "{sample}.hic")
    output:
        tads_10kb = join(RESULTSDIR, "topology", "tads", "{sample}.10kb_juicer_tads.txt"),
        tads_5kb = join(RESULTSDIR, "topology", "tads", "{sample}.5kb_juicer_tads.txt")
    threads: 
        getthreads("call_tads_arrowhead")
    envmodules:
        TOOLS["juicer"]
    shell:
        """
        juicer_tools -Xmx96g arrowhead --threads {threads} -k KR -r 10000 {input.hic} {output.tads_10kb}
        juicer_tools -Xmx96g arrowhead --threads {threads} -k KR -r 5000 {input.hic} {output.tads_5kb}
        """

rule call_tads_hic:
    """From Hsien-Chao's findTAD.sh: identify TADs using hicexplorer"""
    input:
        h5 = join(RESULTSDIR,"h5","{sample}.10kb.h5")
    output:
        tad = join(RESULTSDIR,"topology","tads","{sample}.10kb_tad_score.bm")
    params:
        sample="{sample}"
    threads: 
        getthreads("call_tads_hic")
    envmodules:
        TOOLS["hicexplorer"]
    shell:
        """
        hicFindTADs -m {input.h5} --outPrefix {params.sample}.10kb -p {threads} --correctForMultipleTesting fdr
        """


# Chromatin Loop calling

rule call_loops_mustache:
    """Call chromatin loops using Mustache (https://github.com/ay-lab/mustache)"""
    input:
        hic = join(RESULTSDIR, "hic", "{sample}.hic")
    output:
        loops_10kb = join(RESULTSDIR, "topology", "loops", "{sample}_10kb_loops.tsv"),
        loops_5kb = join(RESULTSDIR, "topology", "loops", "{sample}_5kb_loops.tsv")
    threads: 
        getthreads("call_loops_mustache")
    envmodules:
        TOOLS["mustache"]
    shell:
        """
        mustache -p {threads} -f {input.hic} -r 10000 -o {output.loops_10kb}
        mustache -p {threads} -f {input.hic} -r 5000 -o {output.loops_5kb}
        """

rule call_loops_hic:
    """From Hsien-Chao's detectLoops.sh: identify chromatin loops using hicexplorer"""
    input:
        h5 = join(RESULTSDIR,"h5","{sample}.10kb.h5")
    output:
        bedpe = join(RESULTSDIR,"topology","loops","{sample}.10kb.bedpe")
    threads:
        getthreads("call_loops_hic")
    envmodules:
        TOOLS["hicexplorer"]
    shell:
        """
        hicDetectLoops -m {input.h5} -o {output.bedpe} --threads {threads}
        """