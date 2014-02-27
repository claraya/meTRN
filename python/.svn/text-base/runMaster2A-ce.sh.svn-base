#!/usr/bin/sh


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF
#python mapPeaks.py --path ~/meTRN --mode core --overwrite OFF

# Transfer density files from raw peaks to density path:
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_cx_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_cx_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_ex_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_ex_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_l1_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_l1_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_l2_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_l2_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_l3_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_l3_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_ce_selection_reg_l4_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_ce_selection_reg_l4_compiled.bed

# Perform kernel density estimation of TFBSs (advance):
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/ce_selection_reg_ex/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_ce_selection_reg_ex_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/ce_selection_reg_l1/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_ce_selection_reg_l1_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/ce_selection_reg_l2/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_ce_selection_reg_l2_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/ce_selection_reg_l3/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_ce_selection_reg_l3_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/ce_selection_reg_l4/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_ce_selection_reg_l4_bw300_cs1_cp00001_pl30