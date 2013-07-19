#!/usr/bin/sh


# Generate merged region files:
#mergeBed -i ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP01_any.bed -nms > maphot_ce_selection_reg_cx_occP01_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_denP01_any.bed -nms > maphot_ce_selection_reg_cx_denP01_mrg.bed

#mergeBed -i ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP05_any.bed -nms > maphot_ce_selection_reg_cx_occP05_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_denP05_any.bed -nms > maphot_ce_selection_reg_cx_denP05_mrg.bed


# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism ce
#python mapCircos.py --path ~/meTRN --mode import --organism ce --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP05_region_shared.bed --max 50 --hi 0.90 --lo 0.98 --color dgrey --fillcolor dgrey
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_ex_occP05_region_subset.bed --track ex --max 50 --hi 0.80 --lo 0.88 --color contextyellow --fillcolor contextyellow
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_l1_occP05_region_subset.bed --track l1 --max 50 --hi 0.70 --lo 0.78 --color contextgreen --fillcolor contextgreen
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_l2_occP05_region_subset.bed --track l2 --max 50 --hi 0.60 --lo 0.68 --color contextaqua --fillcolor contextaqua
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_l3_occP05_region_subset.bed --track l3 --max 50 --hi 0.50 --lo 0.58 --color contextsky --fillcolor contextsky
#python mapCircos.py --path ~/meTRN --mode extend --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_l4_occP05_region_subset.bed --track l4 --max 50 --hi 0.40 --lo 0.48 --color contextblue --fillcolor contextblue
#python mapCircos.py --path ~/meTRN --mode launch --organism ce --name maphot_overlap_ce_selection_reg_cx_occP05_region_shared


# Generate HOT-region Circos plots:
python mapCircos.py --path ~/meTRN --mode karyotype --organism ce
python mapCircos.py --path ~/meTRN --mode import --organism ce --infile ~/meTRN/data/hot/regions/maphot_ce_selection_reg_cx_occP05_any.bed --max 50 --hi 0.70 --lo 0.98 --color ceColor --fillcolor ceColor
python mapCircos.py --path ~/meTRN --mode launch --organism ce --name maphot_ce_selection_reg_cx_occP05_any


#top
#bash runMaster2D-ce.sh