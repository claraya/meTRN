#!/usr/bin/sh


# Filter cell annotations:
#python mapCells.py --path ~/meTRN --organism ce --mode filter --source extras --infile waterston_cell_original/pedigree.csv --target waterston_cell_original/filter.csv --name waterston_cell_pedigree.csv
#python mapCells.py --path ~/meTRN --organism ce --mode filter --source extras --infile waterston_cell_original/tissues.csv --target waterston_cell_original/filter.csv --name waterston_cell_tissues-v1.csv


# Simplify cell annotations (converts RHW updated mapping and tissue annotations to the V1 format used elsewhere):
#python mapCells.py --path ~/meTRN --organism ce --mode simply --source extras --infile waterston_cell_mapping-rhw.csv --name waterston_cell_mapping-v2.csv
#python mapCells.py --path ~/meTRN --organism ce --mode simply --source extras --infile waterston_cell_tissues-rhw.csv --name waterston_cell_tissues-v2.csv


# Fill-in cell annotations from Murray:
###python mapCells.py --path ~/meTRN --organism ce --mode fillin --source extras --infile murray_cell_tissues-v1.csv --tissues murray_cell_tissues-v2.csv


# Import cellular expression data (Re-done for new annotations!):
#python mapCells.py --path ~/meTRN --organism ce --mode import --infile waterston_avgExpression.csv --pedigree waterston_cell_pedigree.csv --name waterston --measure avg.expression --ascendants 5 --tissues murray_cell_tissues-v2.csv --target AMA-1,LSY-2,LIN-13,MEP-1,CES-1,CEH-26,NHR-2,CEH-39,MAB-5,HLH-1,UNC-39,F23F12.9,NHR-28 --exclude LIR-2
#python mapCells.py --path ~/meTRN --organism ce --mode import --infile waterston_avgExpression.csv --pedigree waterston_cell_pedigree.csv --name waterston --measure max.expression --ascendants 5 --tissues murray_cell_tissues-v2.csv --target AMA-1,LSY-2,LIN-13,MEP-1,CES-1,CEH-26,NHR-2,CEH-39,MAB-5,HLH-1,UNC-39,F23F12.9,NHR-28 --exclude LIR-2


# Inherit cellular expression down through ancestry (using leaf cells from Murray et al. 2012; re-done for new annotations!):
#python mapCells.py --path ~/meTRN --organism ce --mode inherit --infile waterston_avgExpression.csv --pedigree waterston_cell_pedigree.csv --name waterston --measure avg.expression --mapping murray_2012_SD1_per_gene.txt --inherit last
#python mapCells.py --path ~/meTRN --organism ce --mode inherit --infile waterston_avgExpression.csv --pedigree waterston_cell_pedigree.csv --name waterston --measure max.expression --mapping murray_2012_SD1_per_gene.txt --inherit last


# Rank and correct cellular expression estimates (Re-done for new annotations!):
#python mapCells.py --path ~/meTRN --organism ce --mode correct --infile waterston_avgExpression.csv --pedigree waterston_cell_pedigree.csv --name waterston --measure avg.expression --limit 1500 --parameters 3


# Identify the fractional expression cutoff that minimizes lineage-distance correlation and maximizes cellular-expression correlation (f=0.X):
#python mapCells.py --path ~/meTRN --organism ce --mode cell.distance --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --overwrite ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.distance --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_inleafs --overwrite ON


# Identify the numbers of assayed, started, tracked, focused, inherited and maximal (other form of inherited) cells per time point in development:
#python mapCells.py --path ~/meTRN --organism ce --mode cell.times --times waterston_cell_times.csv --name waterston --measure avg.expression
#python mapCells.py --path ~/meTRN --organism ce --mode cell.cubism --times waterston_cell_times.csv --name waterston --measure avg.expression


# Check data overlap between cellular-expression and ChIP-seq experiments:
#python mapCells.py --path ~/meTRN --organism ce --mode check.status --peaks ce_selection_reg_cx_raw --name waterston --measure avg.expression --ascendants OFF
#python mapCells.py --path ~/meTRN --organism ce --mode check.status --peaks ce_selection_reg_ex_raw --name waterston --measure avg.expression --ascendants OFF

#python mapCells.py --path ~/meTRN --organism ce --mode check.status --peaks ce_selection_com_cx_raw --name waterston --measure avg.expression --ascendants OFF
#python mapCells.py --path ~/meTRN --organism ce --mode check.status --peaks ce_selection_com_ex_raw --name waterston --measure avg.expression --ascendants OFF


# Collect cells expressed per gene:
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_assayed --collection cetrn_assayed_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_started --collection cetrn_started_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.0 --fraction 0.0
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.0 --fraction 0.0

#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_assayed --collection cetrn_assayed_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_started --collection cetrn_started_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.1 --fraction 0.1

#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_assayed --collection cetrn_assayed_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_started --collection cetrn_started_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.2 --fraction 0.2
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.2 --fraction 0.2


# Collect genes expressed per cell:
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_assayed --collection cetrn_assayed_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_started --collection cetrn_started_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.1 --fraction 0.1
#python mapCells.py --path ~/meTRN --organism ce --mode gene.collection --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.1 --fraction 0.1


# Generate gene and cell expression reports:
#python mapCells.py --path ~/meTRN --organism ce --mode reports --collection cetrn_assayed_f0.1,cetrn_started_f0.1,cetrn_tracked_f0.1,cetrn_focused_f0.1,cetrn_inherit_f0.1,cetrn_inleafs_f0.1,cetrn_maximal_f0.1,cetrn_mxleafs_f0.1


# Transfer cells expressed per gene:
#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times tracked --collection cetrn_tracked_f0.0 --name cetrn_time --nametag _f0.0 --fraction 0.0 --start 50 --stop 250 --total 350
#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times tracked --collection cetrn_tracked_f0.1 --name cetrn_time --nametag _f0.1 --fraction 0.1 --start 50 --stop 250 --total 350

#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times focused --collection cetrn_focused_f0.0 --name cetrn_fime --nametag _f0.0 --fraction 0.0 --start 50 --stop 250 --total 350
#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times focused --collection cetrn_focused_f0.1 --name cetrn_fime --nametag _f0.1 --fraction 0.1 --start 50 --stop 250 --total 350

#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times inherit --collection cetrn_inherit_f0.0 --name cetrn_iime --nametag _f0.0 --fraction 0.0 --start 50 --stop 250 --total 350
#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times inherit --collection cetrn_inherit_f0.1 --name cetrn_iime --nametag _f0.1 --fraction 0.1 --start 50 --stop 250 --total 350

#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times maximal --collection cetrn_maximal_f0.0 --name cetrn_mime --nametag _f0.0 --fraction 0.0 --start 50 --stop 250 --total 350
#python mapCells.py --path ~/meTRN --organism ce --mode cell.transfer --times maximal --collection cetrn_maximal_f0.1 --name cetrn_mime --nametag _f0.1 --fraction 0.1 --start 50 --stop 250 --total 350


# Build cellular expression overlap matrix (across all cells tracked and focused):
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.2

#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.1 --extend ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_focused_f0.2 --extend ON

#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.1 --extend ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_inherit_f0.2 --extend ON

#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --extend ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.2 --extend ON

#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.1 --extend ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_maximal_f0.2 --extend ON

#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.1 --extend ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.overlap --expression mapcells_avgExp_waterston_expression_mxleafs --collection cetrn_mxleafs_f0.2 --extend ON


# Build cellular expression overlap matrix (across developmental time-points):
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_time --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --chunks 5 --threads 4 --name OFF
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_time --nametag _f0.1 --start 244 --stop 244 --total 350 --step 1 --chunks 1 --threads 4 --name OFF

#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_fime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --chunks 5 --threads 4 --extend ON --name OFF
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_focused --collection cetrn_fime --nametag _f0.1 --start 244 --stop 244 --total 350 --step 1 --chunks 1 --threads 4 --extend ON --name OFF

#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_iime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --chunks 5 --threads 4 --extend ON --name OFF
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_inherit --collection cetrn_iime --nametag _f0.1 --start 244 --stop 244 --total 350 --step 1 --chunks 1 --threads 4 --extend ON --name OFF

#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_mime --nametag _f0.1 --start 50 --stop 250 --total 350 --step 5 --chunks 5 --threads 4 --extend ON --name OFF
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.overlap --expression mapcells_avgExp_waterston_expression_maximal --collection cetrn_mime --nametag _f0.1 --start 244 --stop 244 --total 350 --step 1 --chunks 1 --threads 4 --extend ON --name OFF


# Build cellular expression SOM matrixes:
#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1 --name binary --fraction 0.1 --technique binary
#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1 --name signal --fraction 0.1 --technique fraction
#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_tracked --collection cetrn_tracked_f0.1 --name normal --fraction 0.1 --technique normal

#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --name binary --fraction 0.1 --technique binary
#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --name signal --fraction 0.1 --technique fraction
#python mapCells.py --path ~/meTRN --organism ce --mode cell.matrix --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1 --name normal --fraction 0.1 --technique normal


# Map lineage and tissue enrichments for factors (Re-done for new annotations!):
#python mapCells.py --path ~/meTRN --organism ce --mode test.tissues --infile mapcells_avgExp_waterston_expression_tissues --expression mapcells_avgExp_waterston_expression_inleafs --collection cetrn_inleafs_f0.1


# Build in silico cellular peaks (across all cells tracked):
#python mapCells.py --path ~/meTRN --organism ce --mode cell.peaks --peaks ce_selection_reg_ex_raw --expression mapcells_avgExp_waterston_expression_focused --name ce_xpfocused_reg_ex_all --fraction 0.1 --overwrite ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.peaks --peaks ce_selection_reg_ex_raw --expression mapcells_avgExp_waterston_expression_tracked --name ce_xptracked_reg_ex_all --fraction 0.1 --overwrite ON
#python mapCells.py --path ~/meTRN --organism ce --mode cell.peaks --peaks ce_selection_reg_ex_raw --expression mapcells_avgExp_waterston_expression_inherit --name ce_xpinherit_reg_ex_all --fraction 0.1 --overwrite ON


# Build in silico cellular peaks (across developmental time-points 170-250):
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.peaks --peaks ce_selection_reg_ex_xot --expression mapcells_avgExp_waterston_expression_focused --name ce_xpfocused_reg_ex_ --fraction 0.1 --times tracked --start 170 --stop 250 --total 350 --threads 4 --chunks 1
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.peaks --peaks ce_selection_reg_ex_xot --expression mapcells_avgExp_waterston_expression_tracked --name ce_xptracked_reg_ex_ --fraction 0.1 --times tracked --start 170 --stop 250 --total 350 --threads 4 --chunks 1
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.peaks --peaks ce_selection_reg_ex_xot --expression mapcells_avgExp_waterston_expression_inherit --name ce_xpinherit_reg_ex_ --fraction 0.1 --times tracked --start 170 --stop 250 --total 350 --threads 4 --chunks 1


# Generate cellular-resolution TSS annotation files (for all cells tracked, focused, and inherited):
#python mapCells.py --path ~/meTRN --organism ce --mode cell.annotation --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --collection cetrn_focused_f0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.annotation --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --collection cetrn_tracked_f0.1
#python mapCells.py --path ~/meTRN --organism ce --mode cell.annotation --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --collection cetrn_inherit_f0.1


# Generate cellular-resolution TSS annotation files (for across developmental time-points):
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.annotation --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --times tracked --name cetrn_time --nametag _f0.1 --start 170 --stop 250 --total 350 --step 10 --chunks 10 --threads 4
#python mapCells.py --path ~/meTRN --organism ce --mode master:cell.annotation --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --times focused --name cetrn_fime --nametag _f0.1 --start 170 --stop 250 --total 350 --step 10 --chunks 10 --threads 4


#top
#bash runMaster6A.sh
