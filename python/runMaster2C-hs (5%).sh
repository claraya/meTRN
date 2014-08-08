#!/usr/bin/sh


# Map HOT regions (fail) and non-HOT regions (pass); without any actual filtering cutoffs:
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_cx --name basics
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name basics --tag GM-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name basics --tag H1-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name basics --tag HG-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name basics --tag HL-
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name basics --tag K5-


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (occupancy 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name occP05 --significance 05 --tag GM- --metric occupancy --cutoff 16
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name occP05 --significance 05 --tag H1- --metric occupancy --cutoff 9
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name occP05 --significance 05 --tag HG- --metric occupancy --cutoff 10
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name occP05 --significance 05 --tag HL- --metric occupancy --cutoff 12
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name occP05 --significance 05 --tag K5- --metric occupancy --cutoff 16


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (density 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name denP05 --significance 05 --tag GM- --metric density --cutoff 17.7
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name denP05 --significance 05 --tag H1- --metric density --cutoff 12.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name denP05 --significance 05 --tag HG- --metric density --cutoff 14.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name denP05 --significance 05 --tag HL- --metric density --cutoff 14.4
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name denP05 --significance 05 --tag K5- --metric density --cutoff 18


# Map HOT regions (fail) and non-HOT regions (pass); applying simulated filtering cutoffs (combined 0.05 significance):
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_gm --name simP05 --significance 05 --tag GM- --metric combined --cutoff 16,17.7
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_h1 --name simP05 --significance 05 --tag H1- --metric combined --cutoff 9,12.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hg --name simP05 --significance 05 --tag HG- --metric combined --cutoff 10,14.5
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_hl --name simP05 --significance 05 --tag HL- --metric combined --cutoff 12,14.4
#python mapHOT.py --path ~/meTRN --organism hs --mode scan --target peaks --peaks hs_selection_reg_k5 --name simP05 --significance 05 --tag K5- --metric combined --cutoff 16,18


# Test overlaps with HOT regions from Yip et al. 2012:
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_gm_occP05_sig05_conXX_occupan_fail.bed --name hs_selection_reg_gm_occP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_h1_occP05_sig05_conXX_occupan_fail.bed --name hs_selection_reg_h1_occP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_hg_occP05_sig05_conXX_occupan_fail.bed --name hs_selection_reg_hg_occP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_hl_occP05_sig05_conXX_occupan_fail.bed --name hs_selection_reg_hl_occP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_k5_occP05_sig05_conXX_occupan_fail.bed --name hs_selection_reg_k5_occP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_gm_occP05_sig05_conXX_occupan_pass.bed --name hs_selection_reg_gm_occP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_h1_occP05_sig05_conXX_occupan_pass.bed --name hs_selection_reg_h1_occP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_hg_occP05_sig05_conXX_occupan_pass.bed --name hs_selection_reg_hg_occP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_hl_occP05_sig05_conXX_occupan_pass.bed --name hs_selection_reg_hl_occP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_occP05_mergeBed_hs_selection_reg_k5_occP05_sig05_conXX_occupan_pass.bed --name hs_selection_reg_k5_occP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_gm_denP05_sig05_conXX_density_fail.bed --name hs_selection_reg_gm_denP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_h1_denP05_sig05_conXX_density_fail.bed --name hs_selection_reg_h1_denP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_hg_denP05_sig05_conXX_density_fail.bed --name hs_selection_reg_hg_denP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_hl_denP05_sig05_conXX_density_fail.bed --name hs_selection_reg_hl_denP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_k5_denP05_sig05_conXX_density_fail.bed --name hs_selection_reg_k5_denP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_gm_denP05_sig05_conXX_density_pass.bed --name hs_selection_reg_gm_denP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_h1_denP05_sig05_conXX_density_pass.bed --name hs_selection_reg_h1_denP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_hg_denP05_sig05_conXX_density_pass.bed --name hs_selection_reg_hg_denP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_hl_denP05_sig05_conXX_density_pass.bed --name hs_selection_reg_hl_denP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_denP05_mergeBed_hs_selection_reg_k5_denP05_sig05_conXX_density_pass.bed --name hs_selection_reg_k5_denP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_gm_simP05_sig05_conXX_combine_fail.bed --name hs_selection_reg_gm_simP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_h1_simP05_sig05_conXX_combine_fail.bed --name hs_selection_reg_h1_simP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hg_simP05_sig05_conXX_combine_fail.bed --name hs_selection_reg_hg_simP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hl_simP05_sig05_conXX_combine_fail.bed --name hs_selection_reg_hl_simP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_k5_simP05_sig05_conXX_combine_fail.bed --name hs_selection_reg_k5_simP05_hot_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed

#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_gm_simP05_sig05_conXX_combine_pass.bed --name hs_selection_reg_gm_simP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Gm12878_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_h1_simP05_sig05_conXX_combine_pass.bed --name hs_selection_reg_h1_simP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_H1hesc_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hg_simP05_sig05_conXX_combine_pass.bed --name hs_selection_reg_hg_simP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Hepg2_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_hl_simP05_sig05_conXX_combine_pass.bed --name hs_selection_reg_hl_simP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_Helas3_merged.bed
#python mapHOT.py --path ~/meTRN --organism hs --mode compare --infile maphot_peaks_simP05_mergeBed_hs_selection_reg_k5_simP05_sig05_conXX_combine_pass.bed --name hs_selection_reg_k5_simP05_rgb_yip_total --source ~/meTRN/extras/yip_2012/HOT_K562_merged.bed


# Calculate overlap in HOT regions and filter regions against the ubiquitously HOT regions:
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name occP05 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name denP05 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx
#python mapHOT.py --path ~/meTRN --organism hs --mode overlap --name simP05 --target GM:hs_selection_reg_gm,H1:hs_selection_reg_h1,HG:hs_selection_reg_hg,HL:hs_selection_reg_hl,K5:hs_selection_reg_k5 --overlap hs_selection_reg_cx


# Store HOT and RGB region files for each context and the global context file:
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_occP05 --regions hs_selection_reg_gm_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_occP05 --regions hs_selection_reg_h1_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_occP05 --regions hs_selection_reg_hg_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_occP05 --regions hs_selection_reg_hl_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_occP05 --regions hs_selection_reg_k5_occP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_occP05 --regions hs_selection_reg_cx_occP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_denP05 --regions hs_selection_reg_gm_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_denP05 --regions hs_selection_reg_h1_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_denP05 --regions hs_selection_reg_hg_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_denP05 --regions hs_selection_reg_hl_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_denP05 --regions hs_selection_reg_k5_denP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_denP05 --regions hs_selection_reg_cx_denP05 --source overlap

#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_gm_simP05 --regions hs_selection_reg_gm_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_h1_simP05 --regions hs_selection_reg_h1_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hg_simP05 --regions hs_selection_reg_hg_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_hl_simP05 --regions hs_selection_reg_hl_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_k5_simP05 --regions hs_selection_reg_k5_simP05 --source analysis
#python mapHOT.py --path ~/meTRN --organism hs --mode regions --name hs_selection_reg_cx_simP05 --regions hs_selection_reg_cx_simP05 --source overlap


# Generate HOT-filtered peak sets (occupancy P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hon --infile maphot_hs_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hot --infile maphot_hs_selection_reg_gm_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hot --infile maphot_hs_selection_reg_h1_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hot --infile maphot_hs_selection_reg_hg_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hot --infile maphot_hs_selection_reg_hl_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hot --infile maphot_hs_selection_reg_k5_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hon --infile maphot_hs_selection_reg_cx_occP05_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hot --infile maphot_hs_selection_reg_gm_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hot --infile maphot_hs_selection_reg_h1_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hot --infile maphot_hs_selection_reg_hg_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hot --infile maphot_hs_selection_reg_hl_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hot --infile maphot_hs_selection_reg_k5_occP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hob --infile maphot_hs_selection_reg_cx_occP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_hot/* ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_hot/* ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_hot/* ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_hot/* ~/meTRN/data/peaks/hs_selection_reg_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_hot/* ~/meTRN/data/peaks/hs_selection_reg_cx_hot/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_hot/* ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_hot/* ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_hot/* ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_hot/* ~/meTRN/data/peaks/hs_selection_com_cx_hot/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_hot/* ~/meTRN/data/peaks/hs_selection_com_cx_hot/


# Generate HOT-filtered peak sets (density P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hdn --infile maphot_hs_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hdt --infile maphot_hs_selection_reg_gm_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hdt --infile maphot_hs_selection_reg_h1_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hdt --infile maphot_hs_selection_reg_hg_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hdt --infile maphot_hs_selection_reg_hl_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hdt --infile maphot_hs_selection_reg_k5_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hdn --infile maphot_hs_selection_reg_cx_denP05_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hdt --infile maphot_hs_selection_reg_gm_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hdt --infile maphot_hs_selection_reg_h1_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hdt --infile maphot_hs_selection_reg_hg_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hdt --infile maphot_hs_selection_reg_hl_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hdt --infile maphot_hs_selection_reg_k5_denP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hdb --infile maphot_hs_selection_reg_cx_denP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_hdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_hdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_hdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_hdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_hdt/* ~/meTRN/data/peaks/hs_selection_reg_cx_hdt/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_hdt/* ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_hdt/* ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_hdt/* ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_hdt/* ~/meTRN/data/peaks/hs_selection_com_cx_hdt/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_hdt/* ~/meTRN/data/peaks/hs_selection_com_cx_hdt/


# Generate HOT-filtered peak sets (combined P05 cutoff):
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hcn --infile maphot_hs_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hct --infile maphot_hs_selection_reg_gm_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hct --infile maphot_hs_selection_reg_h1_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hct --infile maphot_hs_selection_reg_hg_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hct --infile maphot_hs_selection_reg_hl_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hct --infile maphot_hs_selection_reg_k5_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_cx_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_gm_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_h1_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hg_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_hl_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_reg_k5_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_reg_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hcn --infile maphot_hs_selection_reg_cx_simP05_any.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hct --infile maphot_hs_selection_reg_gm_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hct --infile maphot_hs_selection_reg_h1_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hct --infile maphot_hs_selection_reg_hg_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hct --infile maphot_hs_selection_reg_hl_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hct --infile maphot_hs_selection_reg_k5_simP05_hot.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_cx_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_cx_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_gm_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_gm_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_h1_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_h1_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hg_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hg_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_hl_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_hl_raw
#python mapHOT.py --path ~/meTRN --organism hs --mode filter:remove --peaks hs_selection_com_k5_hcb --infile maphot_hs_selection_reg_cx_simP05_all.bed --source ~/meTRN/data/peaks/hs_selection_com_k5_raw

#rm -rf ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#mkdir ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_reg_gm_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_reg_h1_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_reg_hg_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_reg_hl_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_reg_k5_hct/* ~/meTRN/data/peaks/hs_selection_reg_cx_hct/

#rm -rf ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#mkdir ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_com_gm_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_com_h1_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_com_hg_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_com_hl_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/
#cp ~/meTRN/data/peaks/hs_selection_com_k5_hct/* ~/meTRN/data/peaks/hs_selection_com_cx_hct/


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF


#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw_compiled.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0771868 (0.9228132)
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP05_any.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0158601 
#coverageBed -a /Volumes/HD1/Users/claraya/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP05_all.bed -b /Volumes/HD1/Users/claraya/meTRN/input/ucsc_hg19_nuclear_sizes.bed -hist
#0.0017175

#RGB: 100*(0.0771868 - 0.0158601)
#HOT (any): 100*(0.0158601 - 0.0017175)
#HOT (all): 100*0.0017175

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_raw*
#datasets: 340
#regions: 538665
#peaks: 4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_hob*
#datasets: 340
#regions: 537994
#peaks: 3868719
#percent: 3868719/4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_hot*
#datasets: 340
#regions: 536450
#peaks: 2485836
#percent: 2485836/4561908

#wc -l ~/meTRN/data/peaks/mappeaks_hs_selection_com_cx_hon*
#datasets: 340
#regions: 498521
#peaks: 2084994
#percent: 2084994/4561908


#top
#bash "runMaster2C-hs (5%).sh"