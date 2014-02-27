#!/usr/bin/sh


# Prepare inputs for HOT region analyses:
#mkdir ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/
#mkdir ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/
#mkdir ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/

#mkdir ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/input/
#mkdir ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/input/
#mkdir ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP05_region_subset.bed ~/meTRN/data/go/ce_HOTregion_occ_cx_p05/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP05_region_subset.bed ~/meTRN/data/go/dm_HOTregion_occ_cx_p05/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP05_region_subset.bed ~/meTRN/data/go/hs_HOTregion_occ_cx_p05/input/


# Prepare inputs for XOT region analyses:
#mkdir ~/meTRN/data/go/ce_XOTregion_occ_cx_p01/
#mkdir ~/meTRN/data/go/dm_XOTregion_occ_cx_p01/
#mkdir ~/meTRN/data/go/hs_XOTregion_occ_cx_p01/

#mkdir ~/meTRN/data/go/ce_XOTregion_occ_cx_p01/input/
#mkdir ~/meTRN/data/go/dm_XOTregion_occ_cx_p01/input/
#mkdir ~/meTRN/data/go/hs_XOTregion_occ_cx_p01/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/go/ce_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/go/ce_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP01_region_subset.bed ~/meTRN/data/go/ce_XOTregion_occ_cx_p01/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/go/dm_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/go/dm_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP01_region_subset.bed ~/meTRN/data/go/dm_XOTregion_occ_cx_p01/input/

#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/go/hs_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/go/hs_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP01_region_subset.bed ~/meTRN/data/go/hs_XOTregion_occ_cx_p01/input/


# Launch HOT-region GO analysis:
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_hot-hs.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_hot-ce.r
#r --slave < ~/Desktop/Dropbox/meTRN/Code/R/mapGO_hot-dm.r


# Parse GO analysis from HOT regions:
python mapGO.py --path ~/meTRN --organism ce --peaks ce_HOTregion_occ_cx_p05 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
python mapGO.py --path ~/meTRN --organism dm --peaks dm_HOTregion_occ_cx_p05 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
python mapGO.py --path ~/meTRN --organism hs --peaks hs_HOTregion_occ_cx_p05 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05


# Parse GO analysis from XOT regions:
python mapGO.py --path ~/meTRN --organism ce --peaks ce_XOTregion_occ_cx_p01 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
python mapGO.py --path ~/meTRN --organism dm --peaks dm_XOTregion_occ_cx_p01 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05
python mapGO.py --path ~/meTRN --organism hs --peaks hs_XOTregion_occ_cx_p01 --analysis p5e-1 --target hot.regions --minPvalue 0.5 --hitCount 1 --hitPvalue 0.05


#top