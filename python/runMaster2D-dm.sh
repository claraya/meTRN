#!/usr/bin/sh


# Generate merged region files:
#mergeBed -i ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP01_any.bed -nms > maphot_dm_selection_reg_cx_occP01_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_denP01_any.bed -nms > maphot_dm_selection_reg_cx_denP01_mrg.bed

#mergeBed -i ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP05_any.bed -nms > maphot_dm_selection_reg_cx_occP05_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_denP05_any.bed -nms > maphot_dm_selection_reg_cx_denP05_mrg.bed


# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism dm
#python mapCircos.py --path ~/meTRN --mode import --organism dm --infile ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_cx_occP05_region_shared.bed --max 20 --hi 0.90 --lo 0.98 --color dgrey --fillcolor dgrey
#python mapCircos.py --path ~/meTRN --mode extend --organism dm --name maphot_overlap_dm_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_ee_occP05_region_subset.bed --track ee --max 20 --hi 0.80 --lo 0.88 --color contextorange --fillcolor contextorange
#python mapCircos.py --path ~/meTRN --mode extend --organism dm --name maphot_overlap_dm_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_le_occP05_region_subset.bed --track le --max 20 --hi 0.70 --lo 0.78 --color contextmandarin --fillcolor contextmandarin
#python mapCircos.py --path ~/meTRN --mode extend --organism dm --name maphot_overlap_dm_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_pp_occP05_region_subset.bed --track pp --max 20 --hi 0.60 --lo 0.68 --color contextblue --fillcolor contextblue
#python mapCircos.py --path ~/meTRN --mode launch --organism dm --name maphot_overlap_dm_selection_reg_cx_occP05_region_shared


# Generate HOT-region Circos plots:
python mapCircos.py --path ~/meTRN --mode karyotype --organism dm
python mapCircos.py --path ~/meTRN --mode import --organism dm --infile ~/meTRN/data/hot/regions/maphot_dm_selection_reg_cx_occP05_any.bed --max 20 --hi 0.70 --lo 0.98 --color dmColor --fillcolor dmColor
python mapCircos.py --path ~/meTRN --mode launch --organism dm --name maphot_dm_selection_reg_cx_occP05_any


#top
#bash runMaster2D-dm.sh