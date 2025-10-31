# Structural Variants
rule hic_breakfinder:
    """Detects structural variants using hic_breakfinder (https://github.com/dixonlab/hic_breakfinder) - currently not available on Biowulf"""
    input:
        bam = "{sample}.bam",
        exp_file_inter = "inter_expect_1Mb.hg38.txt",
        exp_file_intra = "intra_expect_100kb.hg38.txt"
    output:
        name = "hic_breakfinder/{sample}"
    threads:
        1
    shell:
        """
        hic_breakfinder --bam-file {input.bam} \
        --exp-file-inter {input.exp_file_inter} \
        --exp-file-intra {input.exp_file_intra} \
        --name {output.name} --min-1kb
        """

rule eaglec:
    """Detects structural variants using EagleC (https://github.com/XiaoTaoWang/EagleC)"""
    input:
        mcool = join(RESULTSDIR, "cool", "{sample}.1kb.mcool")
    output:
        sv = join(RESULTSDIR, "topology","sv","eaglec","{sample}.CNN_SVs.5K_combined.txt")
    params:
        sv = join(RESULTSDIR, "topology","sv","eaglec","{sample}")
    envmodules:
        TOOLS["eaglec"]
    shell:
        """
        predictSV --hic-5k {input.mcool}::/resolutions/5000 --hic-10k {input.mcool}::/resolutions/10000 --hic-50k {input.mcool}::/resolutions/50000 \
        -O {params.sv} -g hg38 --balance-type ICE
        """

rule call_cnv_neoloop:
    """Call CNVs using neoloopfinder"""
    input:
        mcool = join(RESULTSDIR, "cool", "{sample}.1kb.mcool")
    output:
        cnv_profile = join(RESULTSDIR, "topology","cnv", "{sample}_{res}.CNV-profile.bedGraph"),
        cnv_seg = join(RESULTSDIR, "topology","cnv", "{sample}_{res}.CNV-seg.bedGraph"),
        cnv_plot = join(RESULTSDIR, "topology","cnv", "{sample}_{res}.CNV.genome-wide.png"),
        log = join(RESULTSDIR, "topology","cnv", "{sample}_{res}.log")
    params:
        genome=config['genome'],
        binsize="{res}",
        ploidy=config['ploidy']
    threads:
        getthreads("call_cnv_neoloop")
    envmodules:
        TOOLS["neoloopfinder"]
    shell:
        """
        # Infer genome-wide CNV profile
        calculate-cnv -H {input.mcool}::/resolutions/{params.binsize} -g {params.genome} \
                      -e uniform --output {output.cnv_profile} --logFile {output.log}

        # Identify CNV segments
        segment-cnv --cnv-file {output.cnv_profile} --binsize {params.binsize} \
                    --ploidy {params.ploidy} --output {output.cnv_seg} --nproc {threads} \
                    --logFile {output.log}

        # Plot CNV
        plot-cnv --cnv-profile {output.cnv_profile} \
                 --cnv-segment {output.cnv_seg} \
                 --output-figure-name {output.cnv_plot} \
                 --dot-size 0.5 --dot-alpha 0.2 --line-width 1 --boundary-width 0.5 \
                 --label-size 7 --tick-label-size 6 --clean-mode

        # Correct CNV
        #correct-cnv -H {input.mcool}::/resolutions/{params.binsize} \
        #            --cnv-file {output.cnv_seg} --nproc {threads} -f
        """

rule reformat_eaglec:
    """Reformat EagleC output from 8 columns to 6 columns"""
    input:
        sv = join(RESULTSDIR, "topology","sv","eaglec","{sample}.CNN_SVs.5K_combined.txt")
    output:
        reformatted = join(RESULTSDIR, "topology","sv","eaglec","{sample}.CNN_SVs.5K_combined.reformat.txt")
    shell:
        """
        awk 'BEGIN {{OFS="\t"}} NR > 1 {{
            orientation = ($5 > $6 && $5 > $7 && $5 > $8) ? "++" : \
                         ($6 > $5 && $6 > $7 && $6 > $8) ? "+-" : \
                         ($7 > $5 && $7 > $6 && $7 > $8) ? "-+" : "--";
            type = ($1 != $3) ? "translocation" : \
                   (orientation == "+-") ? "deletion" : \
                   (orientation == "-+") ? "duplication" : "inversion";
            print $1, $3, orientation, $2, $4, type;
        }}' {input.sv} > {output.reformatted}
        """

rule assemble_sv:
    """Assemble structural variants for rearranged genomic regions."""
    input:
        sv = join(RESULTSDIR, "topology","sv","eaglec","{sample}.CNN_SVs.5K_combined.reformat.txt"),
        mcool = join(RESULTSDIR, "cool", "{sample}.1kb.mcool")
    output:
        assemblies = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}.assemblies.txt"),
        log = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}_assembly.log")
    params:
        sample = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}")
    threads:
        getthreads("assemble_sv")
    envmodules:
        TOOLS["neoloopfinder"]
    shell:
        """
        assemble-complexSVs -O {params.sample} -B {input.sv} --balance-type ICE --protocol insitu --nproc {threads} \
            -H {input.mcool}::/resolutions/25000 \
               {input.mcool}::/resolutions/10000 \
               {input.mcool}::/resolutions/5000 \
            --logFile {output.log} \
            -R 50000000 --minimum-size 50000
        """

rule call_neoloop:
    """Call neo-loops from assembled structural variants."""
    input:
        assemblies = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}.assemblies.txt"),
        mcool = join(RESULTSDIR, "cool", "{sample}.1kb.mcool")
    output:
        neoloops = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}.neo-loops.txt"),
        log = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}_neoloop.log")
    threads:
        getthreads("call_neoloop")
    envmodules:
        TOOLS["neoloopfinder"]
    shell:
        """
        neoloop-caller -O {output.neoloops} --assembly {input.assemblies} \
            --balance-type ICE --protocol insitu --prob 0.9 --nproc {threads} \
            -H {input.mcool}::/resolutions/25000 \
               {input.mcool}::/resolutions/10000 \
               {input.mcool}::/resolutions/5000 \
            --logFile {output.log} \
            -R 5000000
        """

rule call_neotad:
    """Call neo-TADs from assembled structural variants."""
    input:
        assemblies = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}.assemblies.txt"),
        mcool = join(RESULTSDIR, "cool", "{sample}.1kb.mcool")
    output:
        neotads = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}.neo-tads.txt"),
        log = join(RESULTSDIR, "topology","sv", "neoloopfinder","{sample}_neotad.log")
    threads:
        getthreads("call_neotad")
    envmodules:
        TOOLS["neoloopfinder"]
    shell:
        """
        neotad-caller -O {output.neotads} --assembly {input.assemblies} \
            --balance-type ICE --protocol insitu --nproc {threads} \
            -H {input.mcool}::/resolutions/5000 \
            --logFile {output.log} \
            -R 5000000
        """

# Copy Number Variation
#Recommended Software
#CNVkit https://cnvkit.readthedocs.io/en/stable/quickstart.html#install-cnvkit

#Example Command
#generate .cnr file
#cnvkit.py batch K562.bam -r FlatReference.cnn -p 8 -d K562

#segment into copy number calls
#cnvkit.py segment K562.cnr -o K562.cns
