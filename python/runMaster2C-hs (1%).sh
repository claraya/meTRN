#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name basics --tag GM-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name basics --tag H1-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name basics --tag HG-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name basics --tag HL-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name basics --tag K5-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name occP01 --significance 01 --tag GM- --metric occupancy --cutoff 29
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name occP01 --significance 01 --tag H1- --metric occupancy --cutoff 16
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name occP01 --significance 01 --tag HG- --metric occupancy --cutoff 17
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name occP01 --significance 01 --tag HL- --metric occupancy --cutoff 21
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name occP01 --significance 01 --tag K5- --metric occupancy --cutoff 31


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name denP01 --significance 01 --tag GM- --metric density --cutoff 23.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name denP01 --significance 01 --tag H1- --metric density --cutoff 15.6
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name denP01 --significance 01 --tag HG- --metric density --cutoff 18.2
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name denP01 --significance 01 --tag HL- --metric density --cutoff 18.6
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name denP01 --significance 01 --tag K5- --metric density --cutoff 23.7


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.01 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name simP01 --significance 01 --tag GM- --metric combined --cutoff 29,23.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name simP01 --significance 01 --tag H1- --metric combined --cutoff 16,15.6
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name simP01 --significance 01 --tag HG- --metric combined --cutoff 17,18.2
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name simP01 --significance 01 --tag HL- --metric combined --cutoff 21,18.6
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name simP01 --significance 01 --tag K5- --metric combined --cutoff 31,23.7


# Test overlaps with HOT regions from Yip et al. 2012:
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_gm_occP01_sig01_conXX_occupan_fail.bed --name hs_selection_reg_gm_occP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_h1_occP01_sig01_conXX_occupan_fail.bed --name hs_selection_reg_h1_occP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_hg_occP01_sig01_conXX_occupan_fail.bed --name hs_selection_reg_hg_occP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_hl_occP01_sig01_conXX_occupan_fail.bed --name hs_selection_reg_hl_occP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_k5_occP01_sig01_conXX_occupan_fail.bed --name hs_selection_reg_k5_occP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_gm_occP01_sig01_conXX_occupan_pass.bed --name hs_selection_reg_gm_occP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_h1_occP01_sig01_conXX_occupan_pass.bed --name hs_selection_reg_h1_occP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_hg_occP01_sig01_conXX_occupan_pass.bed --name hs_selection_reg_hg_occP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_hl_occP01_sig01_conXX_occupan_pass.bed --name hs_selection_reg_hl_occP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP01_mergeBed_hs_selection_reg_k5_occP01_sig01_conXX_occupan_pass.bed --name hs_selection_reg_k5_occP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_gm_denP01_sig01_conXX_density_fail.bed --name hs_selection_reg_gm_denP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_h1_denP01_sig01_conXX_density_fail.bed --name hs_selection_reg_h1_denP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_hg_denP01_sig01_conXX_density_fail.bed --name hs_selection_reg_hg_denP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_hl_denP01_sig01_conXX_density_fail.bed --name hs_selection_reg_hl_denP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_k5_denP01_sig01_conXX_density_fail.bed --name hs_selection_reg_k5_denP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_gm_denP01_sig01_conXX_density_pass.bed --name hs_selection_reg_gm_denP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_h1_denP01_sig01_conXX_density_pass.bed --name hs_selection_reg_h1_denP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_hg_denP01_sig01_conXX_density_pass.bed --name hs_selection_reg_hg_denP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_hl_denP01_sig01_conXX_density_pass.bed --name hs_selection_reg_hl_denP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP01_mergeBed_hs_selection_reg_k5_denP01_sig01_conXX_density_pass.bed --name hs_selection_reg_k5_denP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_gm_simP01_sig01_conXX_combine_fail.bed --name hs_selection_reg_gm_simP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_h1_simP01_sig01_conXX_combine_fail.bed --name hs_selection_reg_h1_simP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hg_simP01_sig01_conXX_combine_fail.bed --name hs_selection_reg_hg_simP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hl_simP01_sig01_conXX_combine_fail.bed --name hs_selection_reg_hl_simP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_k5_simP01_sig01_conXX_combine_fail.bed --name hs_selection_reg_k5_simP01_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_gm_simP01_sig01_conXX_combine_pass.bed --name hs_selection_reg_gm_simP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_h1_simP01_sig01_conXX_combine_pass.bed --name hs_selection_reg_h1_simP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hg_simP01_sig01_conXX_combine_pass.bed --name hs_selection_reg_hg_simP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_hl_simP01_sig01_conXX_combine_pass.bed --name hs_selection_reg_hl_simP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP01_mergeBed_hs_selection_reg_k5_simP01_sig01_conXX_combine_pass.bed --name hs_selection_reg_k5_simP01_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name occP01 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name denP01 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name simP01 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_occP01 --regions hs_selection_reg_gm_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_occP01 --regions hs_selection_reg_h1_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_occP01 --regions hs_selection_reg_hg_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_occP01 --regions hs_selection_reg_hl_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_occP01 --regions hs_selection_reg_k5_occP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_occP01 --regions hs_selection_reg_cx_occP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_denP01 --regions hs_selection_reg_gm_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_denP01 --regions hs_selection_reg_h1_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_denP01 --regions hs_selection_reg_hg_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_denP01 --regions hs_selection_reg_hl_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_denP01 --regions hs_selection_reg_k5_denP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_denP01 --regions hs_selection_reg_cx_denP01 --source overlap

#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_simP01 --regions hs_selection_reg_gm_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_simP01 --regions hs_selection_reg_h1_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_simP01 --regions hs_selection_reg_hg_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_simP01 --regions hs_selection_reg_hl_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_simP01 --regions hs_selection_reg_k5_simP01 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_simP01 --regions hs_selection_reg_cx_simP01 --source overlap


# Generate HOT-filtered peak sets (occupancy P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xon --infile maphot_hs_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xot --infile maphot_hs_selection_reg_gm_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xot --infile maphot_hs_selection_reg_h1_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xot --infile maphot_hs_selection_reg_hg_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xot --infile maphot_hs_selection_reg_hl_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xot --infile maphot_hs_selection_reg_k5_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xon --infile maphot_hs_selection_reg_cx_occP01_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xot --infile maphot_hs_selection_reg_gm_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xot --infile maphot_hs_selection_reg_h1_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xot --infile maphot_hs_selection_reg_hg_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xot --infile maphot_hs_selection_reg_hl_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xot --infile maphot_hs_selection_reg_k5_occP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xob --infile maphot_hs_selection_reg_cx_occP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_xot/* ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_xot/* ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_xot/* ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_xot/* ~/meTRN/data/peaks/hs_selection_reg_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_xot/* ~/meTRN/data/peaks/hs_selection_reg_cx_xot/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_xot/* ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_xot/* ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_xot/* ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_xot/* ~/meTRN/data/peaks/hs_selection_com_cx_xot/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_xot/* ~/meTRN/data/peaks/hs_selection_com_cx_xot/


# Generate HOT-filtered peak sets (density P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xdn --infile maphot_hs_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xdt --infile maphot_hs_selection_reg_gm_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xdt --infile maphot_hs_selection_reg_h1_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xdt --infile maphot_hs_selection_reg_hg_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xdt --infile maphot_hs_selection_reg_hl_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xdt --infile maphot_hs_selection_reg_k5_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xdn --infile maphot_hs_selection_reg_cx_denP01_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xdt --infile maphot_hs_selection_reg_gm_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xdt --infile maphot_hs_selection_reg_h1_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xdt --infile maphot_hs_selection_reg_hg_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xdt --infile maphot_hs_selection_reg_hl_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xdt --infile maphot_hs_selection_reg_k5_denP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xdb --infile maphot_hs_selection_reg_cx_denP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_xdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_xdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_xdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_xdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_xdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_xdt/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_xdt/* ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_xdt/* ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_xdt/* ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_xdt/* ~/meTRN/data/peaks/hs_selection_com_cx_xdt/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_xdt/* ~/meTRN/data/peaks/hs_selection_com_cx_xdt/


# Generate HOT-filtered peak sets (combined P01 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xcn --infile maphot_hs_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xct --infile maphot_hs_selection_reg_gm_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xct --infile maphot_hs_selection_reg_h1_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xct --infile maphot_hs_selection_reg_hg_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xct --infile maphot_hs_selection_reg_hl_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xct --infile maphot_hs_selection_reg_k5_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xcn --infile maphot_hs_selection_reg_cx_simP01_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xct --infile maphot_hs_selection_reg_gm_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xct --infile maphot_hs_selection_reg_h1_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xct --infile maphot_hs_selection_reg_hg_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xct --infile maphot_hs_selection_reg_hl_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xct --infile maphot_hs_selection_reg_k5_simP01_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_xcb --infile maphot_hs_selection_reg_cx_simP01_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_xct/* ~/meTRN/data/peaks/hs_selection_reg_cx_xct/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_xct/* ~/meTRN/data/peaks/hs_selection_com_cx_xct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw_compiled.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0771868 (0.9228132)
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP01_any.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0082662
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP01_all.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0004812

#RGB: 100*(0.0771868 - 0.0082662)
#XOT (any): 100*(0.0082662 - 0.0004812)
#XOT (all): 100*0.0004812

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw*
#datasets: 340
#regions: 538665
#peaks: 4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_xob*
#datasets: 340
#regions: 538606
#peaks: 4346965
#percent: 4346965/4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_xot*
#datasets: 340
#regions: 541327
#peaks: 3379122
#percent: 3379122/4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_xon*
#datasets: 340
#regions: 523459
#peaks: 2900914
#percent: 2900914/4561908


#top
#bash "runMaster2C-hs (1%).sh"