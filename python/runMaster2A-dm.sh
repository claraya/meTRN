#!/usr/bin/sh


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF

# Transfer density files from raw peaks to density path:
cp ~/meTRN/data/peaks/mappeaks_dm_selection_reg_cx_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_dm_selection_reg_cx_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_dm_selection_reg_ee_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_dm_selection_reg_ee_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_dm_selection_reg_le_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_dm_selection_reg_le_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_dm_selection_reg_pp_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_dm_selection_reg_pp_compiled.bed

# Perform kernel density estimation of TFBSs (advance):
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/dm_selection_reg_ee/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_dm_selection_reg_ee_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/dm_selection_reg_le/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_dm_selection_reg_le_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/dm_selection_reg_pp/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_dm_selection_reg_pp_bw300_cs1_cp00001_pl30
