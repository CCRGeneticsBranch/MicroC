rule call_compartment_hc:
    """From Hsien-Chao's run_ab_correct.sh: """
    input:
    output:
    envmodules:
        TOOLS["deeptools"],
        TOOLS["R"]
    shell:
        """
        multiBigwigSummary BED-file -b $bed -o RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/$bn.npz --outRawCounts RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/$bn.txt --BED RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/ev.txt
        Rscript /data/khanlab3/hsienchao/Yong/HiC/PFI90/ab_analysis/ab_correct.R RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/ev.txt RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/$bn.txt RH4_DMSO_24h_merged/compartment/$res/ev${ev_idx}/ev_corrected_$bn.txt        """
        """

##  (plot_ab.R)

ev_idx=$1
res=$2
bed=$3
#bed=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC.dd.25.RPM.bw
#bed=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_H3K4me1_201118_C_HV2Y5BGXG/Sample_RH4_DMSO_H3K4me1_201118_C_HV2Y5BGXG.25.RPM.bw
bn=$(basename $bed)
bn=`echo $bn | sed -r 's/.*(H3.*?)_.*/\1/'`

#generate ev_pick.txt manually
merge_ev.sh ab_analysis/ev_pick_500kb.txt 500kb
Rscript ab_analysis/plot_ab.R RH4_PFI90_24h_merged/compartment/$res/ev_corrected.txt RH4_DMSO_24h_merged/compartment/$res/ev_corrected.txt ab_analysis/ab_corrected_ev_mergd_500kb.pdf

#Rscript plot_ab.R RH4_PFI90_24h_merged/compartment/ev${ev_idx}/ev_corrected_$bn.txt RH4_DMSO_24h_merged/compartment/ev${ev_idx}/ev_corrected_$bn.txt ab_corrected_ev${ev_idx}_$bn.pdf
module load bedtools
bedtools intersect -wa -wb -a RH4_DMSO_24h_merged/compartment/$res/ev_corrected.txt -b RH4_PFI90_24h_merged/compartment/$res/ev_corrected.txt | cut -f1-5,10,11 > ab_analysis/ev_merged_$res.bed
awk '$4=="A" && $6=="A"' ab_analysis/ev_merged_$res.bed > ab_analysis/ev_merged_${res}_A.bed
awk '$4=="B" && $6=="B"' ab_analysis/ev_merged_$res.bed > ab_analysis/ev_merged_${res}_B.bed
awk '$4=="B" && $6=="A"' ab_analysis/ev_merged_$res.bed | awk '($7-$5) > 0.02' > ab_analysis/ev_merged_${res}_B_to_A.bed
awk '$4=="A" && $6=="B"' ab_analysis/ev_merged_$res.bed | awk '($5-$7) > 0.02' > ab_analysis/ev_merged_${res}_A_to_B.bed
#all genes, TPM >1
awk '$6=="protein_coding"' PFI90.expression.tpm.tsv | cut -f1,2,3,8,11,35 | grep -v seqname > PFI90.expression.tpm.coding.tsv
awk '$5>1 || $6>1' PFI90.expression.tpm.coding.tsv > PFI90.expression.tpm.coding.TPM1.tsv

bedtools intersect -wa -b ab_analysis/ev_merged_${res}_B_to_A.bed -a PFI90.expression.tpm.coding.TPM1.tsv > ab_analysis/ev_merged_${res}_B_to_A_coding_TPM1.tsv
bedtools intersect -wa -b ab_analysis/ev_merged_${res}_A_to_B.bed -a PFI90.expression.tpm.coding.TPM1.tsv > ab_analysis/ev_merged_${res}_A_to_B_coding_TPM1.tsv
cut -f4 ab_analysis/ev_merged_B_to_A_coding_TPM1.tsv | sort | uniq > ab_analysis/ev_merged_${res}_B_to_A_coding_TPM1.genes.list.tsv
cut -f4 ab_analysis/ev_merged_A_to_B_coding_TPM1.tsv | sort | uniq > ab_analysis/ev_merged_${res}_A_to_B_coding_TPM1.genes.list.tsv
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt ab_analysis/ev_merged_${res}_B_to_A_coding_TPM1.genes.list.tsv ab_analysis/ev_merged_${res}_B_to_A_coding_TPM1.NCI_GeneSet_v33.txt
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt ab_analysis/ev_merged_${res}_A_to_B_coding_TPM1.genes.list.tsv ab_analysis/ev_merged_${res}_A_to_B_coding_TPM1.NCI_GeneSet_v33.txt

p3f=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB/MACS_Out_p_1e-07/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB_peaks.narrowPeak.nobl.bed

bedtools intersect -wa -wb -a $p3f -b ab_analysis/ev_merged_${res}_A.bed > ab_analysis/ev_merged_${res}_P3F_A.bed
bedtools intersect -wa -wb -a $p3f -b ab_analysis/ev_merged_${res}_B.bed > ab_analysis/ev_merged_${res}_P3F_B.bed
bedtools intersect -wa -wb -a $p3f -b ab_analysis/ev_merged_${res}_B_to_A.bed > ab_analysis/ev_merged_${res}_P3F_B_to_A.bed
bedtools intersect -wa -wb -a $p3f -b ab_analysis/ev_merged_${res}_A_to_B.bed > ab_analysis/ev_merged_${res}_P3F_A_to_B.bed

cd ab_analysis
./runEDEN.sh ev_merged_${res}_P3F_A_to_B.bed ev_merged_${res}_P3F_B_to_A.bed
cd ..
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt ab_analysis/ev_P3F_A_to_B._TPM1_multi-genes.list.txt ab_analysis/ev_P3F_A_to_B_TPM1_multi-genes.NCI_GeneSet_v33.txt
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt ab_analysis/ev_P3F_B_to_A._TPM1_multi-genes.list.txt ab_analysis/ev_P3F_B_to_A_TPM1_multi-genes.NCI_GeneSet_v33.txt


#module load homer
#annotatePeaks.pl ev_P3F_B_to_A.bed hg19 -annStats ev_P3F_B_to_A.stat > ev_P3F_B_to_A.annotation.txt
#annotatePeaks.pl ev_P3F_A_to_B.bed hg19 -annStats ev_P3F_A_to_B.stat > ev_P3F_A_to_B.annotation.txt
#cut -f10,16 ev_P3F_A_to_B.annotation.txt > ev_P3F_A_to_B.annotation.genes.txt
#cut -f10,16 ev_P3F_B_to_A.annotation.txt > ev_P3F_B_to_A.annotation.genes.txt


# TADs
rule call_diff_tads:
    """From Hsien-Chao's diffTAD.sh: call differential TADs"""
    input:
        control_h5 = lambda wildcards: f"{{control_dir}}/{{control}}.5000.norm.h5",
        treatment_h5 = lambda wildcards: f"{{treatment_dir}}/{{treatment}}.5000.norm.h5",
        control_domains = lambda wildcards: f"{{control_dir}}/{{control}}.5000_domains.bed"
    output:
        diff_tad = lambda wildcards: f"{{out_dir}}/output_differential_tad"
    params:
        control_dir = "{control_dir}",
        treatment_dir = "{treatment_dir}",
        out_dir = "{out_dir}",
        control = "{control}",
        treatment = "{treatment}"
    threads: 
        getthreads("call_diff_tads")
    envmodules:
        TOOLS["hicexplorer"]
    shell:
        """
        hicDifferentialTAD -tm {input.treatment_h5} -cm {input.control_h5} -td {input.control_domains} -o {output.diff_tad}
        """

module load R/4.0.3
cutoff=0.05
Rscript ../calFDRdiffTAD.R $out_dir $cutoff

module load bedtools
sed 's/^chr//g' $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.bed > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.nochr.bed
sed 's/^chr//g' $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.bed > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.nochr.bed

tss_bed=/data/khanlab/projects/pipeline_production/khanlab_pipeline/ref/hg19.tss.coding.nochr.bed

pairToBed -a $control_dir/$control.merged.bedpe -b $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.nochr.bed > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.bedpe
pairToBed -a $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.bedpe -b $tss_bed | cut -f27 | sort | uniq > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.genes.txt
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.genes.txt $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.genes.NCI_GeneSet_v33.txt
pairToBed -a $control_dir/$control.merged.bedpe -b $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.nochr.bed > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.bedpe
pairToBed -a $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.bedpe -b $tss_bed | cut -f27 | sort | uniq > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.genes.txt
Rscript /data/khanlab/projects/hsienchao/scripts/fisherTestGeneSets.R /data/khanlab/projects/hsienchao/apps/msigdb_v7.5/NCI_GeneSet_v33.gmt $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.genes.txt $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.genes.NCI_GeneSet_v33.txt

grep -v geneset $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.genes.NCI_GeneSet_v33.txt | awk -F'\t' '$3 < 0.05' | cut -f8 | sort -n | perl -pe 's/\s/\n/g' | sort | uniq > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.genes.NCI_GeneSet_v33.sig_genes.txt
grep -v geneset $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.genes.NCI_GeneSet_v33.txt | awk -F'\t' '$3 < 0.05' | cut -f8 | sort -n | perl -pe 's/\s/\n/g' | sort | uniq > $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.genes.NCI_GeneSet_v33.sig_genes.txt

/data/khanlab/projects/hsienchao/Berkley/scripts/plotAQuA_APA_BEDPE.sh -1 RH4_DMSO_merged_all -2 RH4_PFI90_merged_all -l $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.dn.loops.bedpe -d . -o $out_dir/apa_plots_fdr$cutoff.dn/apa_plot_loops
/data/khanlab/projects/hsienchao/Berkley/scripts/plotAQuA_APA_BEDPE.sh -1 RH4_DMSO_merged_all -2 RH4_PFI90_merged_all -l $out_dir/output_differential_tad_rejected.diff_tad.fdr$cutoff.up.loops.bedpe -d . -o $out_dir/apa_plots_fdr$cutoff.up/apa_plot_loops

bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.up.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K4me3/H3K4me3_DB.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.up.H3K4me3_DB_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.dn.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K4me3/H3K4me3_DB.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.dn.H3K4me3_DB_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.up.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K9me2/H3K9me2_DB_up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.up.H3K9me2_DB_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.dn.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K9me2/H3K9me2_DB_up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.dn.H3K9me2_DB_up.bed

bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.up.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K4me3/H3K4me3_DB.10kb.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.up.H3K4me3_DB_10kb_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.dn.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K4me3/H3K4me3_DB.10kb.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.dn.H3K4me3_DB_10kb_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.up.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K9me2/H3K9me2_DB_10kb.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.up.H3K9me2_DB_10kb_up.bed
bedtools intersect -a output_differential_tad_rejected.diff_tad.fdr0.1.dn.bed -b /data/khanlab3/hsienchao/Yong/ChIPseq/diffBind_spikeIn/H3K9me2/H3K9me2_DB_10kb.up.bed -wa | uniq > output_differential_tad_rejected.diff_tad.fdr0.1.dn.H3K9me2_DB_10kb_up.bed

#python /data/khanlab/projects/hsienchao/scripts/bigwigConversion.py /data/khanlab/projects/hsienchao/ref/GRCh37_UCSC2ensembl.txt /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC.25.scaled.bw /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC.25.scaled.nochr.bw
#python /data/khanlab/projects/hsienchao/scripts/bigwigConversion.py /data/khanlab/projects/hsienchao/ref/GRCh37_UCSC2ensembl.txt /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB.dd.25.scaled.bw /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB.dd.25.scaled.nochr.bw
#python /data/khanlab/projects/hsienchao/scripts/bigwigConversion.py /data/khanlab/projects/hsienchao/ref/GRCh37_UCSC2ensembl.txt /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_PFI90_1uM_24h_H3K27ac_200125_C_H77JJBGXC/Sample_RH4_PFI90_1uM_24h_H3K27ac_200125_C_H77JJBGXC.25.scaled.bw /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_PFI90_1uM_24h_H3K27ac_200125_C_H77JJBGXC/Sample_RH4_PFI90_1uM_24h_H3K27ac_200125_C_H77JJBGXC.25.scaled.nochr.bw
#python /data/khanlab/projects/hsienchao/scripts/bigwigConversion.py /data/khanlab/projects/hsienchao/ref/GRCh37_UCSC2ensembl.txt /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_PFI90_1_24h_PAX3FOXO1_048_C_HNJTGBGXB/Sample_RH4_PFI90_1_24h_PAX3FOXO1_048_C_HNJTGBGXB.dd.25.scaled.bw /data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_PFI90_1_24h_PAX3FOXO1_048_C_HNJTGBGXB/Sample_RH4_PFI90_1_24h_PAX3FOXO1_048_C_HNJTGBGXB.dd.25.scaled.nochr.bw
#k27ac_bed=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC/MACS_Out_p_1e-07/Sample_RH4_DMSO_24h_H3K27ac_200125_C_H77JJBGXC_peaks.narrowPeak.nobl.bed
#P3F_bed=/data/khanlab/projects/ChIP_seq/DATA/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB/MACS_Out_p_1e-07/Sample_RH4_DMSO_24h_PAX3FOXO1_048_C_HNJTGBGXB_peaks.narrowPeak.nobl.bed


# Differential loop calling and APA
Recommended software: multiHiCcompare

## From Hsien-Chao (diff_loops.R, analysis.html)