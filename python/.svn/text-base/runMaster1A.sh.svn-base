#!/usr/bin/sh


# Generate report (from Anshul):
#python in2shape.py --path ~/meTRN --organism hs --infile akundaje/HumanFinalReport-TFs_SPP_pooled.tsv --folder extras --mode akundaje:report --name hs_modencode_DAT_xx


# Generate paralog files:
#python mapData.py --path ~/meTRN/ --mode paralogs --source extras --infile modencode.family.common.report.txt


# Make genomic windows:
#python in2shape.py --path ~/meTRN --organism ce --mode intervals --name ce_intervals_500_gx --parameters 500 --nuclear OFF
#python in2shape.py --path ~/meTRN --organism ce --mode intervals --name ce_intervals_1kb_gx --parameters 1000 --nuclear OFF


# Prepare worm expression annotations:
#python in2shape.py --path ~/meTRN --organism ce --infile waterston_EE_N2_EE-2.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:expression --name ce_expressed_RNA_ee
#python in2shape.py --path ~/meTRN --organism ce --infile waterston_LE_N2_LE-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:expression --name ce_expressed_RNA_le
#python in2shape.py --path ~/meTRN --organism ce --infile waterston_L1_N2_L1-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:expression --name ce_expressed_RNA_l1
#python in2shape.py --path ~/meTRN --organism ce --infile waterston_L2_N2_L2-4.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:expression --name ce_expressed_RNA_l2
#python in2shape.py --path ~/meTRN --organism ce --infile waterston_L3_N2_L3-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:expression --name ce_expressed_RNA_l3


# Merge worm expression files:
#cat in2shape_ce_expressed_RNA_ee.bed in2shape_ce_expressed_RNA_le.bed > in2shape_ce_expressed_RNA_em.bed
#cat in2shape_ce_expressed_RNA_ee.gff in2shape_ce_expressed_RNA_le.gff > in2shape_ce_expressed_RNA_em.gff
#cat in2shape_ce_expressed_RNA_ee.nh.bed in2shape_ce_expressed_RNA_le.nh.bed > in2shape_ce_expressed_RNA_em.nh.bed
#cat in2shape_ce_expressed_RNA_ee.nh.gff in2shape_ce_expressed_RNA_le.nh.gff > in2shape_ce_expressed_RNA_em.nh.gff


# Prepare worm feature annotation files:
#python in2shape.py --path ~/meTRN --organism ce --infile wormbase_ws220_annotations.gff --folder extras --mode wormbase:complete --name ce_wormbased_COM_gx
#python in2shape.py --path ~/meTRN --organism ce --infile wormbase_ws220_annotations.gff --folder extras --mode wormbase:gene --name ce_wormbased_GEN_gx
#python in2shape.py --path ~/meTRN --organism ce --infile wormbase_ws220_annotations.gff --folder extras --mode wormbase:transcript --name ce_wormbased_RNA_gx
#python in2shape.py --path ~/meTRN --organism ce --infile wormbase_ws220_annotations.gff --folder extras --mode wormbase:tes --name ce_wormbased_TES_gx
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_wbGene.bed --folder analysis --mode featureID:wormbase --name ce_wormbased_GEN_wb


# Prepare chromatin state files:
#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.worm_EE_final.bed --folder extras --mode import --name ce_modencode_HMM_ee
#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.worm_L3_final.bed --folder extras --mode import --name ce_modencode_HMM_l3

#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.fly_EL.bed --folder extras --mode import --name dm_modencode_HMM_le --cutChr ON
#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.fly_L3.bed --folder extras --mode import --name dm_modencode_HMM_pp --cutChr ON

#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.human_GM.bed --folder extras --mode import --name hs_modencode_HMM_gm
#python in2shape.py --path ~/meTRN --organism ce --infile iHMM/iHMM.M1K16.human_H1.bed --folder extras --mode import --name hs_modencode_HMM_h1


#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_worm-ee_dense.bed --folder extras --mode import --name ce_modencode_CMM_ee --cutChr ON
#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_worm-l3_dense.bed --folder extras --mode import --name ce_modencode_CMM_l3 --cutChr ON

#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_fly-el_dense.bed --folder extras --mode import --name dm_modencode_CMM_le --cutChr ON
#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_fly-l3_dense.bed --folder extras --mode import --name dm_modencode_CMM_pp --cutChr ON

#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_human-gm12878_dense.bed --folder extras --mode import --name hs_modencode_CMM_gm --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile cHMM/WFH_upwt_pval1e-3_n19_human-h1_dense.bed --folder extras --mode import --name hs_modencode_CMM_h1 --cutChr OFF

#python in2shape.py --path ~/meTRN --organism ce --infile xHMM/chromhmm.segway.gm12878.comb11.concord4.bed --folder extras --mode import --name hs_modencode_XMM_gm --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile xHMM/chromhmm.segway.h1hesc.comb11.concord4.bed --folder extras --mode import --name hs_modencode_XMM_h1 --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile xHMM/chromhmm.segway.hepg2.comb11.concord4.bed --folder extras --mode import --name hs_modencode_XMM_hg --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile xHMM/chromhmm.segway.helas3.comb11.concord4.bed --folder extras --mode import --name hs_modencode_XMM_hl --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile xHMM/chromhmm.segway.k562.comb11.concord4.bed --folder extras --mode import --name hs_modencode_XMM_k5 --cutChr OFF

#python in2shape.py --path ~/meTRN --organism ce --infile kHMM/gm12878.ChromHMM.bed --folder extras --mode import --name hs_modencode_KMM_gm --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile kHMM/h1hesc.ChromHMM.bed --folder extras --mode import --name hs_modencode_KMM_h1 --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile kHMM/hepg2.ChromHMM.bed --folder extras --mode import --name hs_modencode_KMM_hg --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile kHMM/helas3.ChromHMM.bed --folder extras --mode import --name hs_modencode_KMM_hl --cutChr OFF
#python in2shape.py --path ~/meTRN --organism ce --infile kHMM/k562.ChromHMM.bed --folder extras --mode import --name hs_modencode_KMM_k5 --cutChr OFF


# Prepare worm TSS annotation files (and adjust them):
#python in2shape.py --path ~/meTRN --organism ce --infile wormbase_ws220_annotations.gff --folder extras --mode wormbase:tss --name ce_wormbased_TSS_gx
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_wormbased_TSS_gx --folder annotations --mode slopbed --name ce_wormbased_TSS_gx --parameters "-l 5000 -r 500 -s"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_wormbased_TSS_gx --folder annotations --mode slopbed --name ce_wormbased_TSS_gx --parameters "-l 4000 -r 400 -s"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_wormbased_TSS_gx --folder annotations --mode slopbed --name ce_wormbased_TSS_gx --parameters "-l 3000 -r 300 -s"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_wormbased_TSS_gx --folder annotations --mode slopbed --name ce_wormbased_TSS_gx --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_wormbased_TSS_gx --folder annotations --mode slopbed --name ce_wormbased_TSS_gx --parameters "-l 1000 -r 100 -s"


# Prepare human TSS annotation files (and adjust them):
#python in2shape.py --path ~/meTRN --organism hs --infile Gencodev10_TSS_May2012.gff --folder extras --mode gff2bed --name hs_modencode_TSS_gx
#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_modencode_TSS_gx --folder annotations --mode slopbed --name hs_TSS_modencode_TSS_gx --parameters "-l 5000 -r 500 -s"


# Prepare fly TSS annotation files (and adjust them):
#python in2shape.py --path ~/meTRN --organism dm --infile Celniker_Drosophila_Annotation_20120616_1428.tss.gtf --folder extras --mode gff2bed --name dm_modencode_TSS_gx --cutChr ON
#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_modencode_TSS_gx --folder annotations --mode slopbed --name dm_modencode_TSS_gx --parameters "-l 2000 -r 200 -s"


# Import CapSeq and CIP-TAP TSSs from Gu et al. 2012 (Cell):
#python in2shape.py --path ~/meTRN --organism ce --infile gu_2012_table_s1b.txt --folder extras --mode gu2012:all --name ce_gu2012all_TSS_gx
#python in2shape.py --path ~/meTRN --organism ce --infile gu_2012_table_s1b.txt --folder extras --mode gu2012:pro --name ce_gu2012pro_TSS_gx

#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_gu2012all_TSS_gx --folder annotations --mode slopbed --name ce_gu2012all_TSS_gx --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_gu2012pro_TSS_gx --folder annotations --mode slopbed --name ce_gu2012pro_TSS_gx --parameters "-l 2000 -r 200 -s"


# Prepare species TSS annotation files (with feature names):
#python in2shape.py --path ~/meTRN --organism hs --infile Gencodev10_TSS_May2012.gff --folder extras --mode gencode:tss --name hs_gencode10_TSS_gx
#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_gencode10_TSS_gx --folder annotations --mode slopbed --name hs_gencode10_TSS_gx --parameters "-l 5000 -r 500 -s"

#python in2shape.py --path ~/meTRN --organism dm --infile Celniker_Drosophila_Annotation_20120616_1428.tss.gtf --folder extras --mode celniker:tss --name dm_celnikers_TSS_gx --cutChr ON
#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_celnikers_TSS_gx --folder annotations --mode slopbed --name dm_celnikers_TSS_gx --parameters "-l 2000 -r 200 -s"


# Prepare enhancer files:
#python in2shape.py --path ~/meTRN --organism ce --infile enhancers/worm_EE_CBP-based_enhancers.txt --folder extras --mode enhancers --name ce_modencode_DHS_ee --header ON --cutChr ON
#python in2shape.py --path ~/meTRN --organism ce --infile enhancers/worm_L3_CBP-based_enhancers.txt --folder extras --mode enhancers --name ce_modencode_DHS_l3 --header ON --cutChr ON

#python in2shape.py --path ~/meTRN --organism dm --infile enhancers/fly_LE_DHS-based_enhancers.txt --folder extras --mode enhancers --name dm_modencode_DHS_le --header ON --cutChr ON
#python in2shape.py --path ~/meTRN --organism dm --infile enhancers/fly_BG3_DHS-based_enhancers.txt --folder extras --mode enhancers --name dm_modencode_DHS_bg --header ON --cutChr ON
#python in2shape.py --path ~/meTRN --organism dm --infile enhancers/fly_Kc_DHS-based_enhancers.txt --folder extras --mode enhancers --name dm_modencode_DHS_kc --header ON --cutChr ON
#python in2shape.py --path ~/meTRN --organism dm --infile enhancers/fly_S2_DHS-based_enhancers.txt --folder extras --mode enhancers --name dm_modencode_DHS_s2 --header ON --cutChr ON

#python in2shape.py --path ~/meTRN --organism hs --infile enhancers/human_Gm12878_DHS-based_enhancers.txt --folder extras --mode enhancers --name hs_modencode_DHS_gm --header ON --cutChr OFF
#python in2shape.py --path ~/meTRN --organism hs --infile enhancers/human_H1_DHS-based_enhancers.txt --folder extras --mode enhancers --name hs_modencode_DHS_h1 --header ON --cutChr OFF
#python in2shape.py --path ~/meTRN --organism hs --infile enhancers/human_HeLa_DHS-based_enhancers.txt --folder extras --mode enhancers --name hs_modencode_DHS_hl --header ON --cutChr OFF
#python in2shape.py --path ~/meTRN --organism hs --infile enhancers/human_K562_DHS-based_enhancers.txt --folder extras --mode enhancers --name hs_modencode_DHS_k5 --header ON --cutChr OFF



# Prepare enhancer regions:
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_modencode_DHS_ee --folder annotations --mode slopbed --name ce_modencode_DHS_ee --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism ce --infile in2shape_ce_modencode_DHS_l3 --folder annotations --mode slopbed --name ce_modencode_DHS_l3 --parameters "-l 500 -r 500"

#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_modencode_DHS_le --folder annotations --mode slopbed --name dm_modencode_DHS_le --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_modencode_DHS_bg --folder annotations --mode slopbed --name dm_modencode_DHS_bg --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_modencode_DHS_kc --folder annotations --mode slopbed --name dm_modencode_DHS_kc --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism dm --infile in2shape_dm_modencode_DHS_s2 --folder annotations --mode slopbed --name dm_modencode_DHS_s2 --parameters "-l 500 -r 500"

#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_modencode_DHS_gm --folder annotations --mode slopbed --name hs_modencode_DHS_gm --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_modencode_DHS_h1 --folder annotations --mode slopbed --name hs_modencode_DHS_h1 --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_modencode_DHS_hl --folder annotations --mode slopbed --name hs_modencode_DHS_hl --parameters "-l 500 -r 500"
#python in2shape.py --path ~/meTRN --organism hs --infile in2shape_hs_modencode_DHS_k5 --folder annotations --mode slopbed --name hs_modencode_DHS_k5 --parameters "-l 500 -r 500"



# Generate TSS sites for the worm, as detected in the expression data:
#python in2shape.py --path ~/meTRN --infile waterston_EE_N2_EE-2.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:tss --name ce_modencode_TSS_ee
#python in2shape.py --path ~/meTRN --infile waterston_LE_N2_LE-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:tss --name ce_modencode_TSS_le
#python in2shape.py --path ~/meTRN --infile waterston_L1_N2_L1-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:tss --name ce_modencode_TSS_l1
#python in2shape.py --path ~/meTRN --infile waterston_L2_N2_L2-4.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:tss --name ce_modencode_TSS_l2
#python in2shape.py --path ~/meTRN --infile waterston_L3_N2_L3-1.integrated_transcripts.v1203f.ws220.gff3 --folder extras --mode waterston:tss --name ce_modencode_TSS_l3


# Generate promoter regions for the worm:
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_ee --folder annotations --mode slopbed --name ce_modencode_TSS_ee --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_le --folder annotations --mode slopbed --name ce_modencode_TSS_le --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l1 --folder annotations --mode slopbed --name ce_modencode_TSS_l1 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l2 --folder annotations --mode slopbed --name ce_modencode_TSS_l2 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l3 --folder annotations --mode slopbed --name ce_modencode_TSS_l3 --organism ce --parameters "-l 2000 -r 200 -s"

#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_ee --folder annotations --mode slopbed --name ce_modencode_TSS_ee --organism ce --parameters "-l 5000 -r 500 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_le --folder annotations --mode slopbed --name ce_modencode_TSS_le --organism ce --parameters "-l 5000 -r 500 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l1 --folder annotations --mode slopbed --name ce_modencode_TSS_l1 --organism ce --parameters "-l 5000 -r 500 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l2 --folder annotations --mode slopbed --name ce_modencode_TSS_l2 --organism ce --parameters "-l 5000 -r 500 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_modencode_TSS_l3 --folder annotations --mode slopbed --name ce_modencode_TSS_l3 --organism ce --parameters "-l 5000 -r 500 -s"


#cat in2shape_ce_modencode_TSS_ee.nh.bed in2shape_ce_modencode_TSS_le.nh.bed in2shape_ce_modencode_TSS_l1.nh.bed in2shape_ce_modencode_TSS_l2.nh.bed in2shape_ce_modencode_TSS_l3.nh.bed > in2shape_ce_modencode_TSS_gx.nh.bed
#mergeBed -i in2shape_ce_modencode_TSS_gx.nh.bed -nms -s > in2shape_ce_modencode_TSS_mx.nh.bed


#cat in2shape_ce_modencode_TSS_ee_slopbed_up2000_dn200.nh.bed in2shape_ce_modencode_TSS_le_slopbed_up2000_dn200.nh.bed in2shape_ce_modencode_TSS_l1_slopbed_up2000_dn200.nh.bed in2shape_ce_modencode_TSS_l2_slopbed_up2000_dn200.nh.bed in2shape_ce_modencode_TSS_l3_slopbed_up2000_dn200.nh.bed > in2shape_ce_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed
#mergeBed -i in2shape_ce_modencode_TSS_gx_slopbed_up2000_dn200.nh.bed -nms -s > in2shape_ce_modencode_TSS_mx_slopbed_up2000_dn200.nh.bed


#cat in2shape_ce_modencode_TSS_ee_slopbed_up5000_dn500.nh.bed in2shape_ce_modencode_TSS_le_slopbed_up5000_dn500.nh.bed in2shape_ce_modencode_TSS_l1_slopbed_up5000_dn500.nh.bed in2shape_ce_modencode_TSS_l2_slopbed_up5000_dn500.nh.bed in2shape_ce_modencode_TSS_l3_slopbed_up5000_dn500.nh.bed > in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed
#mergeBed -i in2shape_ce_modencode_TSS_gx_slopbed_up5000_dn500.nh.bed -nms -s > in2shape_ce_modencode_TSS_mx_slopbed_up5000_dn500.nh.bed


# Make TSS distance regions as (window distances from TSS!):
#python in2shape.py --path ~/meTRN --infile in2shape_ce_wormbased_TSS_gx.bed --folder annotations --mode windows --name 125kb --parameters 0:1000,1001:2000,2001:3000,3001:4000,4001:5000 --header ON
#python in2shape.py --path ~/meTRN --infile in2shape_dm_modencode_TSS_gx.bed --folder annotations --mode windows --name 125kb --parameters 0:1000,1001:2000,2001:3000,3001:4000,4001:5000 --header bed
#python in2shape.py --path ~/meTRN --infile in2shape_hs_modencode_TSS_gx.bed --folder annotations --mode windows --name 125kb --parameters 0:1000,1001:2000,2001:3000,3001:4000,4001:5000 --header bed

#python in2shape.py --path ~/meTRN --infile in2shape_ce_wormbased_TSS_gx.bed --folder annotations --mode windows --name 500bp --parameters 0:500,501:1000,1001:2000,2001:10000 --header ON
#python in2shape.py --path ~/meTRN --infile in2shape_dm_modencode_TSS_gx.bed --folder annotations --mode windows --name 500bp --parameters 0:500,501:1000,1001:2000,2001:10000 --header bed
#python in2shape.py --path ~/meTRN --infile in2shape_hs_modencode_TSS_gx.bed --folder annotations --mode windows --name 500bp --parameters 0:500,501:1000,1001:2000,2001:10000 --header bed


# Make combined regulatory regions (promoters and enhancers):
#python in2shape.py --path ~/meTRN --organism ce --folder annotations --mode regulatory.mix --infile in2shape_ce_wormbased_TSS_gx_windows_500bp.bed --parameters in2shape_ce_modencode_DHS_ee_slopbed_up500_dn500.bed --name ce_modencode_MIX_ee --target window
#python in2shape.py --path ~/meTRN --organism ce --folder annotations --mode regulatory.mix --infile in2shape_ce_wormbased_TSS_gx_windows_500bp.bed --parameters in2shape_ce_modencode_DHS_l3_slopbed_up500_dn500.bed --name ce_modencode_MIX_l3 --target window

#python in2shape.py --path ~/meTRN --organism dm --folder annotations --mode regulatory.mix --infile in2shape_dm_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_dm_modencode_DHS_le_slopbed_up500_dn500.bed --name dm_modencode_MIX_le --target window
#python in2shape.py --path ~/meTRN --organism dm --folder annotations --mode regulatory.mix --infile in2shape_dm_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_dm_modencode_DHS_bg_slopbed_up500_dn500.bed --name dm_modencode_MIX_bg --target window
#python in2shape.py --path ~/meTRN --organism dm --folder annotations --mode regulatory.mix --infile in2shape_dm_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_dm_modencode_DHS_kc_slopbed_up500_dn500.bed --name dm_modencode_MIX_kc --target window
#python in2shape.py --path ~/meTRN --organism dm --folder annotations --mode regulatory.mix --infile in2shape_dm_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_dm_modencode_DHS_s2_slopbed_up500_dn500.bed --name dm_modencode_MIX_s2 --target window

#python in2shape.py --path ~/meTRN --organism hs --folder annotations --mode regulatory.mix --infile in2shape_hs_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_hs_modencode_DHS_gm_slopbed_up500_dn500.bed --name hs_modencode_MIX_gm --target window
#python in2shape.py --path ~/meTRN --organism hs --folder annotations --mode regulatory.mix --infile in2shape_hs_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_hs_modencode_DHS_h1_slopbed_up500_dn500.bed --name hs_modencode_MIX_h1 --target window
#python in2shape.py --path ~/meTRN --organism hs --folder annotations --mode regulatory.mix --infile in2shape_hs_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_hs_modencode_DHS_hl_slopbed_up500_dn500.bed --name hs_modencode_MIX_hl --target window
#python in2shape.py --path ~/meTRN --organism hs --folder annotations --mode regulatory.mix --infile in2shape_hs_modencode_TSS_gx_windows_500bp.bed --parameters in2shape_hs_modencode_DHS_k5_slopbed_up500_dn500.bed --name hs_modencode_MIX_k5 --target window


# Generate promoter regions for expressed and repressed genes:
#python in2shape.py --path ~/meTRN --infile in2shape_ce_expressed_TSS_ee.bed --folder annotations --mode slopbed --name ce_expressed_TSS_ee --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_repressed_TSS_ee.bed --folder annotations --mode slopbed --name ce_repressed_TSS_ee --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_reference_TSS_ee.bed --folder annotations --mode slopbed --name ce_reference_TSS_ee --organism ce --parameters "-l 2000 -r 200 -s"

#python in2shape.py --path ~/meTRN --infile in2shape_ce_expressed_TSS_le.bed --folder annotations --mode slopbed --name ce_expressed_TSS_le --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_repressed_TSS_le.bed --folder annotations --mode slopbed --name ce_repressed_TSS_le --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_reference_TSS_le.bed --folder annotations --mode slopbed --name ce_reference_TSS_le --organism ce --parameters "-l 2000 -r 200 -s"

#python in2shape.py --path ~/meTRN --infile in2shape_ce_expressed_TSS_l1.bed --folder annotations --mode slopbed --name ce_expressed_TSS_l1 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_repressed_TSS_l1.bed --folder annotations --mode slopbed --name ce_repressed_TSS_l1 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_reference_TSS_l1.bed --folder annotations --mode slopbed --name ce_reference_TSS_l1 --organism ce --parameters "-l 2000 -r 200 -s"

#python in2shape.py --path ~/meTRN --infile in2shape_ce_expressed_TSS_l2.bed --folder annotations --mode slopbed --name ce_expressed_TSS_l2 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_repressed_TSS_l2.bed --folder annotations --mode slopbed --name ce_repressed_TSS_l2 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_reference_TSS_l2.bed --folder annotations --mode slopbed --name ce_reference_TSS_l2 --organism ce --parameters "-l 2000 -r 200 -s"

#python in2shape.py --path ~/meTRN --infile in2shape_ce_expressed_TSS_l3.bed --folder annotations --mode slopbed --name ce_expressed_TSS_l3 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_repressed_TSS_l3.bed --folder annotations --mode slopbed --name ce_repressed_TSS_l3 --organism ce --parameters "-l 2000 -r 200 -s"
#python in2shape.py --path ~/meTRN --infile in2shape_ce_reference_TSS_l3.bed --folder annotations --mode slopbed --name ce_reference_TSS_l3 --organism ce --parameters "-l 2000 -r 200 -s"