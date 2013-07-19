#!/usr/bin/sh


# Determine overlap with chromatin states as a function of occupancy:
#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_selection_reg_ex_raw --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_selection_reg_l3_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap

#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_selection_reg_ex_raw --infile in2shape_ce_modencode_CMM_ee.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature
#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_selection_reg_l3_raw --infile in2shape_ce_modencode_CMM_l3.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature


# Determine overlap with chromatin states for each individual factor (selection):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_ex_raw --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l1_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l2_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l3_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l4_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap


# Determine overlap with chromatin states for each individual factor (extension):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_ex_raw --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l1_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l2_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l3_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l4_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap


# Determine overlap with promoter/enhancer distances for each individual factor (selection):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_ex_raw --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l1_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l2_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l3_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_selection_reg_l4_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500


# Determine overlap with promoter/enhancer distances for each individual factor (extension):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_ex_raw --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l1_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l2_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l3_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_extension_reg_l4_raw --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500


#python mapBinding.py --path ~/meTRN --organism ce --mode map:dataset --peaks ce_reporting_com_xx_raw --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps --headerDict auto --target window --label dataset --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500


# Determine correlation between promoter/enhancer distances and number of peaks for each individual factor:
#python mapBinding.py --path ~/meTRN --organism ce --mode map:peaks --peaks ce_selection_reg_ex_raw --name enhancer.overlaps
#python mapBinding.py --path ~/meTRN --organism ce --mode map:peaks --peaks ce_selection_reg_l1_raw --name enhancer.overlaps
#python mapBinding.py --path ~/meTRN --organism ce --mode map:peaks --peaks ce_selection_reg_l2_raw --name enhancer.overlaps
#python mapBinding.py --path ~/meTRN --organism ce --mode map:peaks --peaks ce_selection_reg_l3_raw --name enhancer.overlaps
#python mapBinding.py --path ~/meTRN --organism ce --mode map:peaks --peaks ce_selection_reg_l4_raw --name enhancer.overlaps


# Prepare inputs for HOT region analyses (specific):
#rm -rf ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/
#mkdir ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/
#mkdir ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/input/
#rm -rf ~/meTRN/data/binding/ce_HOTregion_occ_cx_p05/input/maphot_overlap_ce_selection_reg_cx_occP05_region_unique.bed


# Prepare inputs for XOT region analyses (specific):
#rm -rf ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/
#mkdir ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/
#mkdir ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/input/
#rm -rf ~/meTRN/data/binding/ce_XOTregion_occ_cx_p01/input/maphot_overlap_ce_selection_reg_cx_occP01_region_unique.bed


# Determine overlap with chromatin states for HOT regions (specific):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_occ_cx_p05 --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states.ex --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --order cx-shared,ex-unique,l1-unique,l2-unique,l3-unique,l4-unique
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_XOTregion_occ_cx_p01 --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states.ex --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP01_region_:,_:- --order cx-shared,ex-unique,l1-unique,l2-unique,l3-unique,l4-unique
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_occ_cx_p05 --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states.l3 --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --order cx-shared,ex-unique,l1-unique,l2-unique,l3-unique,l4-unique
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_XOTregion_occ_cx_p01 --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states.l3 --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP01_region_:,_:- --order cx-shared,ex-unique,l1-unique,l2-unique,l3-unique,l4-unique


# Determine overlap with promoter/enhancer distances for HOT regions (specific):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_occ_cx_p05 --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps.ee --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_occ_cx_p05 --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps.l3 --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer


# Determine overlap with promoter (only) distances for HOT regions (specific):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_occ_cx_p05 --infile in2shape_ce_modencode_MIX_ee.bed --name promoter.overlaps --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize 0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000


# Prepare inputs for HOT region analyses (subset):
#rm -rf ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/
#mkdir ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/
#mkdir ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/input/
#cp ~/meTRN/data/hot/overlap/maphot_overlap_ce_selection_reg_**_occP05_region_subset.bed ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/input/
#rm -rf ~/meTRN/data/binding/ce_HOTregion_sub_cx_p05/input/maphot_overlap_ce_selection_reg_cx_occP05_region_subset.bed


# Determine overlap with chromatin states for HOT regions (subset):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_sub_cx_p05 --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states.ex --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --order cx-shared,ex-subset,l1-subset,l2-subset,l3-subset,l4-subset
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_sub_cx_p05 --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states.l3 --queries auto --target feature --label factor --exclude 17_Unmap --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --order cx-shared,ex-subset,l1-subset,l2-subset,l3-subset,l4-subset


# Determine overlap with promoter/enhancer distances for HOT regions (subset):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_sub_cx_p05 --infile in2shape_ce_modencode_MIX_ee.bed --name enhancer.overlaps.ee --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_sub_cx_p05 --infile in2shape_ce_modencode_MIX_l3.bed --name enhancer.overlaps.l3 --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer


# Determine overlap with promoter (only) distances for HOT regions (subset):
#python mapBinding.py --path ~/meTRN --organism ce --mode map:regions --peaks ce_HOTregion_sub_cx_p05 --infile in2shape_ce_modencode_MIX_ee.bed --name promoter.overlaps --headerDict auto --target window --label factor --rename maphot_overlap_ce_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize 0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000


# Examine temperature changes in stage-specific (unique) HOT regions:
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_ex_occP05_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l1_occP05_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l2_occP05_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l3_occP05_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l4_occP05_region_unique.bed --contexts order.extended --source overlap


# Examine temperature changes in stage-specific (unique) XOT regions:
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_ex_occP01_region_unique.bed --contexts order.extended --source overlap 
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l1_occP01_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l2_occP01_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l3_occP01_region_unique.bed --contexts order.extended --source overlap
#python mapHOT.py --path ~/meTRN --organism ce --mode temperature --peaks ce_selection_reg_cx --target dataset --infile maphot_overlap_ce_selection_reg_l4_occP01_region_unique.bed --contexts order.extended --source overlap


# Generate the set of peaks for factors that were assayed in embryonic and larval L3 stages:
#python mapPeaks.py --path ~/meTRN --mode comparison --rename selection --peaks ce_selection_reg_xx_raw --nametag comparing --contexts ex,l3
#python mapPeaks.py --path ~/meTRN --mode build --overwrite OFF
		

# Determine overlap with chromatin states as a function of occupancy in the context-comparative peak sets:
#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_comparing_reg_ex_raw --infile in2shape_ce_modencode_HMM_ee.bed --name genomic.states --start 1 --stop 12 --queries auto --target feature
#python mapBinding.py --path ~/meTRN --organism ce --mode map:overlap --peaks ce_comparing_reg_l3_raw --infile in2shape_ce_modencode_HMM_l3.bed --name genomic.states --start 1 --stop 12 --queries auto --target feature
		

#top
#bash runMaster2E-ce.sh