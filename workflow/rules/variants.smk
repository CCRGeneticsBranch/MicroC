# Structural Variants
Recommended Software
hic_breakfinder https://github.com/dixonlab/hic_breakfinder

Example Command
./hic_breakfinder --bam-file K562.bam --exp-file-inter inter_expect_1Mb.hg38.txt --exp-file-intra intra_expect_100kb.hg38.txt --name hic_breakfinder/K562


# Copy Number Variation
Recommended Software
CNVkit https://cnvkit.readthedocs.io/en/stable/quickstart.html#install-cnvkit

Example Command
generate .cnr file
cnvkit.py batch K562.bam -r FlatReference.cnn -p 8 -d K562

segment into copy number calls
cnvkit.py segment K562.cnr -o K562.cns