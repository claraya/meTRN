#!/usr/bin/sh


# Determine overlap with chromatin states as a function of occupancy (NOTE: These crash on memory issues!):
#python mapBinding.py --path ~/meTRN --organism hs --mode map:overlap --peaks hs_selection_reg_gm_raw --infile in2shape_hs_modencode_HMM_gm.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism hs --mode map:overlap --peaks hs_selection_reg_h1_raw --infile in2shape_hs_modencode_HMM_h1.bed --name genomic.states --start 1 --stop 40 --queries auto --target feature --exclude 17_Unmap

#python mapBinding.py --path ~/meTRN --organism hs --mode map:overlap --peaks hs_selection_reg_gm_raw --infile in2shape_hs_modencode_CMM_gm.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature
#python mapBinding.py --path ~/meTRN --organism hs --mode map:overlap --peaks hs_selection_reg_h1_raw --infile in2shape_hs_modencode_CMM_h1.bed --name chrmHMM.states --start 1 --stop 40 --queries auto --target feature


# Determine overlap with chromatin states for each individual factor:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_gm_raw --infile in2shape_hs_modencode_HMM_gm.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_h1_raw --infile in2shape_hs_modencode_HMM_h1.bed --name genomic.states --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap


# Determine overlap with promoter/enhancer distances for each individual factor:
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_gm_raw --infile in2shape_hs_modencode_MIX_gm.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_h1_raw --infile in2shape_hs_modencode_MIX_h1.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_hl_raw --infile in2shape_hs_modencode_MIX_hl.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500
#python mapBinding.py --path ~/meTRN --organism hs --mode map:dataset --peaks hs_selection_reg_k5_raw --infile in2shape_hs_modencode_MIX_k5.bed --name enhancer.overlaps --queries auto --headerDict auto --target window --label factor --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer --reference 0:500



# Prepare inputs for HOT region analyses (specific):
rm -rf ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/
mkdir ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/
mkdir ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_cx_occP05_region_shared.bed ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP05_region_unique.bed ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/input/
rm -rf ~/meTRN/data/binding/hs_HOTregion_occ_cx_p05/input/maphot_overlap_hs_selection_reg_cx_occP05_region_unique.bed


# Prepare inputs for XOT region analyses (specific):
rm -rf ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/
mkdir ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/
mkdir ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_cx_occP01_region_shared.bed ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/input/
cp ~/meTRN/data/hot/overlap/maphot_overlap_hs_selection_reg_**_occP01_region_unique.bed ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/input/
rm -rf ~/meTRN/data/binding/hs_XOTregion_occ_cx_p01/input/maphot_overlap_hs_selection_reg_cx_occP01_region_unique.bed


# Determine overlap with chromatin states for each individual factor:
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_HMM_gm.bed --name genomic.states.gm --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:-
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_XOTregion_occ_cx_p01 --infile in2shape_hs_modencode_HMM_gm.bed --name genomic.states.gm --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap --rename maphot_overlap_hs_selection_reg_:,occP01_region_:,_:-

python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_HMM_h1.bed --name genomic.states.h1 --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:-
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_XOTregion_occ_cx_p01 --infile in2shape_hs_modencode_HMM_h1.bed --name genomic.states.h1 --queries auto --target feature --label factor --reference 1_Pro --exclude 17_Unmap --rename maphot_overlap_hs_selection_reg_:,occP01_region_:,_:-


# Determine overlap with promoter/enhancer distances for HOT regions:
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_MIX_gm.bed --name enhancer.overlaps.gm --queries auto --headerDict auto --target window --label factor --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_MIX_h1.bed --name enhancer.overlaps.h1 --queries auto --headerDict auto --target window --label factor --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_MIX_hl.bed --name enhancer.overlaps.hl --queries auto --headerDict auto --target window --label factor --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer
python mapBinding.py --path ~/meTRN --organism hs --mode map:regions --peaks hs_HOTregion_occ_cx_p05 --infile in2shape_hs_modencode_MIX_k5.bed --name enhancer.overlaps.k5 --queries auto --headerDict auto --target window --label factor --rename maphot_overlap_hs_selection_reg_:,occP05_region_:,_:- --others ON --ids feature --prioritize enhancer,0:500,501:1000,1001:2000,2001:10000 --fraction 0.1 --elsewhere ON --queries 0:500,501:1000,1001:2000,2001:10000,enhancer


#top
#bash runMaster2E-hs.sh