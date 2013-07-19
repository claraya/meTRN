#!/usr/bin/sh


# Determine overlap with chromatin states as a function of occupancy:
#python mapBinding.py --path ~/meTRN --organism dm --mode map:overlap --peaks dm_selection_reg_le_raw --infile in2shape_dm_modencode_HMM_le.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism dm --mode map:overlap --peaks dm_selection_reg_pp_raw --infile in2shape_dm_modencode_HMM_pp.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap

#python mapBinding.py --path ~/meTRN --organism dm --mode map:overlap --peaks dm_selection_reg_le_raw --infile in2shape_dm_modencode_CMM_le.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature
#python mapBinding.py --path ~/meTRN --organism dm --mode map:overlap --peaks dm_selection_reg_pp_raw --infile in2shape_dm_modencode_CMM_pp.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature


# Determine overlap with chromatin states for each individual factor:
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_selection_reg_ee_raw --infile in2shape_dm_modencode_HMM_le.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_selection_reg_le_raw --infile in2shape_dm_modencode_HMM_le.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_selection_reg_pp_raw --infile in2shape_dm_modencode_HMM_pp.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap


# Determine overlap with promoter/enhancer distances for each individual factor:
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_selection_reg_ee_raw --infile in2shape_dm_modencode_MIX_le.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_selection_reg_le_raw --infile in2shape_dm_modencode_MIX_le.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_extension_reg_kc_raw --infile in2shape_dm_modencode_MIX_kc.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism dm --mode map:dataset --peaks dm_extension_reg_s2_raw --infile in2shape_dm_modencode_MIX_s2.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500


# Prepare inputs for HOT region analyses (specific):
rm -rf ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/
mkdir ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/
mkdir ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/input/
rm -rf ~/meTRN/data/binding/dm_HOTregion_occ_cx_p05/input/maphot_overlap_dm_selection_reg_cx_occP05_region_unique.bed


# Prepare inputs for XOT region analyses (specific):
rm -rf ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/
mkdir ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/
mkdir ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_dm_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/input/
rm -rf ~/meTRN/data/binding/dm_XOTregion_occ_cx_p01/input/maphot_overlap_dm_selection_reg_cx_occP01_region_unique.bed


# Determine overlap with chromatin states for each individual factor:
python mapBinding.py --path ~/meTRN --organism dm --mode map:regions --peaks dm_HOTregion_occ_cx_p05 --infile in2shape_dm_modencode_HMM_le.bed --name genomic.states.le --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_dm_selection_reg_:,occP05_region_:,_:- --order cx-shared,ee-unique,le-unique,pp-unique
python mapBinding.py --path ~/meTRN --organism dm --mode map:regions --peaks dm_XOTregion_occ_cx_p01 --infile in2shape_dm_modencode_HMM_le.bed --name genomic.states.le --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_dm_selection_reg_:,occP01_region_:,_:- --order cx-shared,ee-unique,le-unique,pp-unique

python mapBinding.py --path ~/meTRN --organism dm --mode map:regions --peaks dm_HOTregion_occ_cx_p05 --infile in2shape_dm_modencode_HMM_pp.bed --name genomic.states.pp --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_dm_selection_reg_:,occP05_region_:,_:- --order cx-shared,ee-unique,le-unique,pp-unique
python mapBinding.py --path ~/meTRN --organism dm --mode map:regions --peaks dm_XOTregion_occ_cx_p01 --infile in2shape_dm_modencode_HMM_pp.bed --name genomic.states.pp --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_dm_selection_reg_:,occP01_region_:,_:- --order cx-shared,ee-unique,le-unique,pp-unique


# Determine overlap with promoter/enhancer distances for HOT regions:
#python mapBinding.py --path ~/meTRN --organism dm --mode map:regions --peaks dm_HOTregion_occ_cx_p05 --infile in2shape_dm_modencode_MIX_le.bed --name enhancer.overlaps.le --queries auto --headerDict auto --target window --label factor --rename maphot_overlap_dm_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer


#top
#bash runMaster2E-dm.sh


#cp ~/meTRN/data/binding/ce_selection_reg_ex_raw/enhancer.overlaps/summary/mapbinding_ce_selection_reg_ex_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/
#cp ~/meTRN/data/binding/ce_selection_reg_l3_raw/enhancer.overlaps/summary/mapbinding_ce_selection_reg_l3_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/

#cp ~/meTRN/data/binding/dm_selection_reg_le_raw/enhancer.overlaps/summary/mapbinding_dm_selection_reg_le_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/

#cp ~/meTRN/data/binding/hs_selection_reg_gm_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_gm_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/
#cp ~/meTRN/data/binding/hs_selection_reg_h1_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/
#cp ~/meTRN/data/binding/hs_selection_reg_hl_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_hl_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/
#cp ~/meTRN/data/binding/hs_selection_reg_k5_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_k5_raw_compiled_dataset_summary_values ~/Desktop/Dropbox/
