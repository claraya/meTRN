#!/usr/bin/sh


# Compare ortholog binding in chromatin regions:
python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,ce --orthology groups --source binding --A hs_selection_reg_h1_raw/genomic.states/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_normal --B ce_extension_reg_ex_raw/genomic.states/summary/mapbinding_ce_extension_reg_ex_raw_compiled_dataset_summary_normal --name genomic.states --label rebuild --target factor --indexes iHMM

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,ce --orthology groups --source binding --A hs_selection_reg_gm_raw/genomic.states/summary/mapbinding_hs_selection_reg_gm_raw_compiled_dataset_summary_normal --B ce_extension_reg_l3_raw/genomic.states/summary/mapbinding_ce_extension_reg_l3_raw_compiled_dataset_summary_normal --name genomic.states --label rebuild --target factor --indexes iHMM

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,dm --orthology groups --source binding --A hs_selection_reg_h1_raw/genomic.states/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_normal --B dm_selection_reg_le_raw/genomic.states/summary/mapbinding_dm_selection_reg_le_raw_compiled_dataset_summary_normal --name genomic.states --label rebuild --target factor --indexes iHMM

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,dm --orthology groups --source binding --A hs_selection_reg_gm_raw/genomic.states/summary/mapbinding_hs_selection_reg_gm_raw_compiled_dataset_summary_normal --B dm_selection_reg_pp_raw/genomic.states/summary/mapbinding_dm_selection_reg_pp_raw_compiled_dataset_summary_normal --name genomic.states --label rebuild --target factor --indexes iHMM


# Compare ortholog binding in promoter/enhancer regions:
python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,ce --orthology groups --source binding --A hs_selection_reg_h1_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_normal --B ce_extension_reg_ex_raw/enhancer.overlaps/summary/mapbinding_ce_extension_reg_ex_raw_compiled_dataset_summary_normal --name enhancer.overlaps --label rebuild --target factor --indexes 125EN

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,ce --orthology groups --source binding --A hs_selection_reg_gm_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_gm_raw_compiled_dataset_summary_normal --B ce_extension_reg_l3_raw/enhancer.overlaps/summary/mapbinding_ce_extension_reg_l3_raw_compiled_dataset_summary_normal --name enhancer.overlaps --label rebuild --target factor --indexes 125EN

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,dm --orthology groups --source binding --A hs_selection_reg_h1_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_normal --B dm_selection_reg_le_raw/enhancer.overlaps/summary/mapbinding_dm_selection_reg_le_raw_compiled_dataset_summary_normal --name enhancer.overlaps --label rebuild --target factor --indexes 125EN

python mapHybrid.py --path ~/meTRN --mode merge.binding --organism hs --species hs,dm --orthology groups --source binding --A hs_selection_reg_h1_raw/enhancer.overlaps/summary/mapbinding_hs_selection_reg_h1_raw_compiled_dataset_summary_normal --B dm_selection_reg_ee_raw/enhancer.overlaps/summary/mapbinding_dm_selection_reg_ee_raw_compiled_dataset_summary_normal --name enhancer.overlaps --label rebuild --target factor --indexes 125EN


# Generate combined pairwse version of the binding matrixes (from Alan's pairwise file):
#python mapBinding.py --path ~/meTRN --mode pairwise --organism hs --source extras --infile boyle_ortholog_context_matrix_matched.mat --peaks xx_orthologs_reg_xx_raw --name enhancer.overlaps


### Convert ortholog binding in promoter/enhancer regions to our format (retired):
###python mapBinding.py --path ~/meTRN --mode convert --organism hs --source extras --infile boyle_ortholog_context_list.txt --peaks xx_orthologs_reg_xx_raw --name enhancer.overlaps
###python mapBinding.py --path ~/meTRN --mode aggregate --organism hs --source hs_selection_reg_gm_raw,hs_selection_reg_h1_raw,hs_selection_reg_hl_raw,hs_selection_reg_k5_raw,ce_extension_reg_ex_raw,ce_extension_reg_l3_raw,dm_selection_reg_le_raw --peaks xx_orthologs_reg_xx_raw --name enhancer.overlaps


#top
#bash runMaster2F.sh