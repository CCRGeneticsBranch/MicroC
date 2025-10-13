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

# Copy Number Variation
#Recommended Software
#CNVkit https://cnvkit.readthedocs.io/en/stable/quickstart.html#install-cnvkit

#Example Command
#generate .cnr file
#cnvkit.py batch K562.bam -r FlatReference.cnn -p 8 -d K562

#segment into copy number calls
#cnvkit.py segment K562.cnr -o K562.cns