#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_cx --name basics
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name basics --tag GM-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name basics --tag H1-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name basics --tag HG-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name basics --tag HL-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name basics --tag K5-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.01 significance):
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name simP01 --cutoff 29 --tag GM-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name simP01 --cutoff 16 --tag H1-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name simP01 --cutoff 17 --tag HG-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name simP01 --cutoff 21 --tag HL-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name simP01 --cutoff 31 --tag K5-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.05 significance):
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name simP05 --cutoff 16 --tag GM-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name simP05 --cutoff 9 --tag H1-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name simP05 --cutoff 10 --tag HG-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name simP05 --cutoff 12 --tag HL-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name simP05 --cutoff 16 --tag K5-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (0.10 significance):
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name simP10 --cutoff 12 --tag GM-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name simP10 --cutoff 6 --tag H1-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name simP10 --cutoff 7 --tag HG-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name simP10 --cutoff 9 --tag HL-
python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name simP10 --cutoff 11 --tag K5-



# Test overlaps with HOT regions from Yip et al. 2012:
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_gm_simP01_cut29_conXX_compiled_fail.bed --name hs_selection_reg_gm_p01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_h1_simP01_cut16_conXX_compiled_fail.bed --name hs_selection_reg_h1_p01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hg_simP01_cut17_conXX_compiled_fail.bed --name hs_selection_reg_hg_p01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hl_simP01_cut21_conXX_compiled_fail.bed --name hs_selection_reg_hl_p01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_k5_simP01_cut31_conXX_compiled_fail.bed --name hs_selection_reg_k5_p01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_gm_simP01_cut29_conXX_compiled_pass.bed --name hs_selection_reg_gm_p01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_h1_simP01_cut16_conXX_compiled_pass.bed --name hs_selection_reg_h1_p01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hg_simP01_cut17_conXX_compiled_pass.bed --name hs_selection_reg_hg_p01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hl_simP01_cut21_conXX_compiled_pass.bed --name hs_selection_reg_hl_p01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_k5_simP01_cut31_conXX_compiled_pass.bed --name hs_selection_reg_k5_p01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_gm_simP05_cut16_conXX_compiled_fail.bed --name hs_selection_reg_gm_p05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_h1_simP05_cut09_conXX_compiled_fail.bed --name hs_selection_reg_h1_p05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hg_simP05_cut10_conXX_compiled_fail.bed --name hs_selection_reg_hg_p05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hl_simP05_cut12_conXX_compiled_fail.bed --name hs_selection_reg_hl_p05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_k5_simP05_cut16_conXX_compiled_fail.bed --name hs_selection_reg_k5_p05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_gm_simP05_cut16_conXX_compiled_pass.bed --name hs_selection_reg_gm_p05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_h1_simP05_cut09_conXX_compiled_pass.bed --name hs_selection_reg_h1_p05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hg_simP05_cut10_conXX_compiled_pass.bed --name hs_selection_reg_hg_p05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hl_simP05_cut12_conXX_compiled_pass.bed --name hs_selection_reg_hl_p05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_k5_simP05_cut16_conXX_compiled_pass.bed --name hs_selection_reg_k5_p05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_gm_simP10_cut12_conXX_compiled_fail.bed --name hs_selection_reg_gm_p10_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_h1_simP10_cut06_conXX_compiled_fail.bed --name hs_selection_reg_h1_p10_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_hg_simP10_cut07_conXX_compiled_fail.bed --name hs_selection_reg_hg_p10_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_hl_simP10_cut09_conXX_compiled_fail.bed --name hs_selection_reg_hl_p10_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_k5_simP10_cut11_conXX_compiled_fail.bed --name hs_selection_reg_k5_p10_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_gm_simP10_cut12_conXX_compiled_pass.bed --name hs_selection_reg_gm_p10_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_h1_simP10_cut06_conXX_compiled_pass.bed --name hs_selection_reg_h1_p10_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_hg_simP10_cut07_conXX_compiled_pass.bed --name hs_selection_reg_hg_p10_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_hl_simP10_cut09_conXX_compiled_pass.bed --name hs_selection_reg_hl_p10_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP10_mergeBed_hs_selection_reg_k5_simP10_cut11_conXX_compiled_pass.bed --name hs_selection_reg_k5_p10_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name simP01 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx

python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name simP05 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_simP01 --regions hs_selection_reg_gm_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_simP01 --regions hs_selection_reg_h1_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_simP01 --regions hs_selection_reg_hg_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_simP01 --regions hs_selection_reg_hl_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_simP01 --regions hs_selection_reg_k5_simP01 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_simP01 --regions hs_selection_reg_cx_simP01 --source overlap

python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_simP05 --regions hs_selection_reg_gm_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_simP05 --regions hs_selection_reg_h1_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_simP05 --regions hs_selection_reg_hg_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_simP05 --regions hs_selection_reg_hl_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_simP05 --regions hs_selection_reg_k5_simP05 --source analysis
python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_simP05 --regions hs_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets:
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xno --infile maphot_hs_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xct --infile maphot_hs_selection_reg_gm_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xct --infile maphot_hs_selection_reg_h1_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xct --infile maphot_hs_selection_reg_hg_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xct --infile maphot_hs_selection_reg_hl_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xct --infile maphot_hs_selection_reg_k5_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xno --infile maphot_hs_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xct --infile maphot_hs_selection_reg_gm_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xct --infile maphot_hs_selection_reg_h1_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xct --infile maphot_hs_selection_reg_hg_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xct --infile maphot_hs_selection_reg_hl_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xct --infile maphot_hs_selection_reg_k5_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xub --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_reg_gm_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_reg_h1_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_reg_hg_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_reg_hl_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_reg_k5_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/

rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_xct/
mkdir ~/meTRN/data/peaks/hs_selection_com_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_com_gm_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_com_h1_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_com_hg_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_com_hl_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
cp ~/meTRN/data/peaks/hs_selection_com_k5_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/


# Generate HOT-filtered peak sets:
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hno --infile maphot_hs_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hct --infile maphot_hs_selection_reg_gm_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hct --infile maphot_hs_selection_reg_h1_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hct --infile maphot_hs_selection_reg_hg_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hct --infile maphot_hs_selection_reg_hl_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hct --infile maphot_hs_selection_reg_k5_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hno --infile maphot_hs_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hct --infile maphot_hs_selection_reg_gm_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hct --infile maphot_hs_selection_reg_h1_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hct --infile maphot_hs_selection_reg_hg_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hct --infile maphot_hs_selection_reg_hl_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hct --infile maphot_hs_selection_reg_k5_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hub --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_reg_gm_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_reg_h1_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_reg_hg_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_reg_hl_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_reg_k5_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/

rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_hct/
mkdir ~/meTRN/data/peaks/hs_selection_com_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_com_gm_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_com_h1_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_com_hg_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_com_hl_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
cp ~/meTRN/data/peaks/hs_selection_com_k5_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF

#coverageBed -a ~/meTRN/data/peaks/mappeaks_hs_selection_com_k5_raw_compiled.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0258050 (0.9741951)
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_k5_simP05_rgb.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0199302
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_k5_simP05_hot.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0056409
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_k5_simP01_hot.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0028298

#coverageBed -a ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw_compiled.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0599497 (0.9400503)
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_simP05_any.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0105589 
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_simP05_all.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0006660
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_simP01_any.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0048583
#coverageBed -a ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_simP01_all.bed -b ~/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0000935

#RGB: 100*(0.0599497 - 0.0105589)
#HOT (any): 100*(0.0105589 - 0.0006660)
#HOT (all): 100*0.0006660
#XOT (any): 100*(0.0048583-0.0000935)
#XOT (all): 100*0.0000935



# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism c.elegans --source ce1
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_peaks_simple_gffkde2_hs_selection_reg_sx_simple_cut38_conXX_bw300_cs1_cp00001_pl30_peaks_extended_fail.bed --name hot_regions_simple
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/analysis/maphot_filter_select_sx_shared_peaks_basics_gffkde2_hs_selection_reg_sx_basics_cutXX_conXX_bw300_cs1_cp00001_pl30_peaks_extended_pass.bed --name hot_regions_filter
#python mapCircos.py --path ~/meTRN --mode import --organism c.elegans --source ce1 --infile ~/meTRN/data/hot/overlap/maphot_overlap_simple_hs_selection_reg_sx_shared --name hot_regions_filter

#circos -conf ./circos.conf 

#top