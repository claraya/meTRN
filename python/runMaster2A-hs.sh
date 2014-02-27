#!/usr/bin/sh


# Generate complete, collapsed, density and report files:
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF

# Transfer density files from raw peaks to density path:
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_cx_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_cx_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_gm_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_gm_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_h1_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_h1_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_hg_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_hg_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_hl_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_hl_compiled.bed
cp ~/meTRN/data/peaks/mappeaks_hs_selection_reg_k5_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_selection_reg_k5_compiled.bed

#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_cx_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_cx_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_gm_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_gm_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_h1_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_h1_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_hg_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_hg_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_hl_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_hl_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_standards_reg_k5_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_standards_reg_k5_compiled.bed

#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_cx_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_cx_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_gm_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_gm_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_h1_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_h1_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_hg_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_hg_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_hl_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_hl_compiled.bed
#cp ~/meTRN/data/peaks/mappeaks_hs_redundant_reg_k5_raw_compiled.bed ~/meTRN/data/hot/density/mergeBed_hs_redundant_reg_k5_compiled.bed


# Perform kernel density estimation of TFBSs (advance):
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/hs_standards_reg_gm/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_hs_standards_reg_gm_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/hs_standards_reg_h1/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_hs_standards_reg_h1_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/hs_standards_reg_hg/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_hs_standards_reg_hg_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/hs_standards_reg_hl/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_hs_standards_reg_hl_bw300_cs1_cp00001_pl30
#perl ~/meTRN/scripts/GFFKDE2.pl ~/meTRN/idr/peaks/hs_standards_reg_k5/ 300 0.1 0.00001 30 ~/meTRN/data/hot/density/gffkde2_hs_standards_reg_k5_bw300_cs1_cp00001_pl30