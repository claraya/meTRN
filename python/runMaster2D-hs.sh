#!/usr/bin/sh


# Generate merged region files:
#mergeBed -i ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP01_any.bed -nms > maphot_hs_selection_reg_cx_occP01_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_denP01_any.bed -nms > maphot_hs_selection_reg_cx_denP01_mrg.bed

#mergeBed -i ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP05_any.bed -nms > maphot_hs_selection_reg_cx_occP05_mrg.bed
#mergeBed -i ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_denP05_any.bed -nms > maphot_hs_selection_reg_cx_denP05_mrg.bed


# Generate HOT-region Circos plots:
#python mapCircos.py --path ~/meTRN --mode karyotype --organism hs
#python mapCircos.py --path ~/meTRN --mode import --organism hs --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_cx_occP05_region_shared.bed --max 100 --hi 0.90 --lo 0.98 --color dgrey --fillcolor dgrey
#python mapCircos.py --path ~/meTRN --mode extend --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_h1_occP05_region_subset.bed --track h1 --max 100 --hi 0.80 --lo 0.88 --color contextyellow --fillcolor contextyellow
#python mapCircos.py --path ~/meTRN --mode extend --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_gm_occP05_region_subset.bed --track gm --max 100 --hi 0.70 --lo 0.78 --color contextgreen --fillcolor contextgreen
#python mapCircos.py --path ~/meTRN --mode extend --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_hg_occP05_region_subset.bed --track hg --max 100 --hi 0.60 --lo 0.68 --color contextaqua --fillcolor contextaqua
#python mapCircos.py --path ~/meTRN --mode extend --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_hl_occP05_region_subset.bed --track hl --max 100 --hi 0.50 --lo 0.58 --color contextsky --fillcolor contextsky
#python mapCircos.py --path ~/meTRN --mode extend --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared --infile ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_k5_occP05_region_subset.bed --track k5 --max 100 --hi 0.40 --lo 0.48 --color contextblue --fillcolor contextblue
#python mapCircos.py --path ~/meTRN --mode launch --organism hs --name maphot_overlap_hs_selection_reg_cx_occP05_region_shared


# Generate HOT-region Circos plots:
python mapCircos.py --path ~/meTRN --mode karyotype --organism hs
python mapCircos.py --path ~/meTRN --mode import --organism hs --infile ~/meTRN/data/hot/regions/maphot_hs_selection_reg_cx_occP05_any.bed --max 100 --hi 0.70 --lo 0.98 --color hsColor --fillcolor hsColor
python mapCircos.py --path ~/meTRN --mode launch --organism hs --name maphot_hs_selection_reg_cx_occP05_any


#top
#bash runMaster2D-hs.sh