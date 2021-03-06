#!/usr/bin/sh


# Generate cellular-resolution matrix files:
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_focused --peaks ce_selection_com_ex_raw --contexts order.condensed --name focus --target factor --include AMA-1 --expression 0.1 --technique binary
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_focused --peaks ce_selection_com_ex_raw --contexts order.condensed --name focus --target organism --include AMA-1 --expression 0.1 --technique binary

#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_inleafs --peaks ce_selection_com_ex_raw --contexts order.condensed --name leafs --target factor --include AMA-1 --expression 0.1 --technique binary
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_inleafs --peaks ce_selection_com_ex_raw --contexts order.condensed --name leafs --target organism --include AMA-1 --expression 0.1 --technique binary

#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_focused --peaks ce_selection_com_ex_xot --contexts order.condensed --name focus --target factor --include AMA-1 --expression 0.1 --technique binary
#python mapNeurons.py --path ~/meTRN --organism ce --mode matrix.cells --infile mapcells_avgExp_waterston_expression_focused --peaks ce_selection_com_ex_xot --contexts order.condensed --name focus --target organism --include AMA-1 --expression 0.1 --technique binary


# Generate cellular-resolution signature reports:
#python mapNeurons.py --path ~/meTRN --organism ce --mode signatures --infile mapneurons_matrix_cells_focus_xp1_FAC_ce_selection_com_ex_raw.txt --peaks ce_selection_com_ex_raw --technique binary
#python mapNeurons.py --path ~/meTRN --organism ce --mode signatures --infile mapneurons_matrix_cells_focus_xp1_ORG_ce_selection_com_ex_raw.txt --peaks ce_selection_com_ex_raw --technique binary

#python mapNeurons.py --path ~/meTRN --organism ce --mode signatures --infile mapneurons_matrix_cells_focus_xp1_FAC_ce_selection_com_ex_xot.txt --peaks ce_selection_com_ex_xot --technique binary
#python mapNeurons.py --path ~/meTRN --organism ce --mode signatures --infile mapneurons_matrix_cells_focus_xp1_ORG_ce_selection_com_ex_xot.txt --peaks ce_selection_com_ex_xot --technique binary


# Generate BED files from SOM regions:
#python mapNeurons.py --path ~/meTRN --organism ce --mode convert.regions --peaks ce_selection_com_ex_raw --neurons cell.ex.som


# Generate neuron reports from SOM codes and regions:
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.summary --peaks ce_selection_com_ex_raw --neurons cell.ex.som --fast ON


# Assign features to neurons:
#python mapNeurons.py --path ~/meTRN --organism ce --mode neuron.mapping --peaks ce_selection_com_ex_raw --neurons cell.ex.som --source annotations --infile in2shape_ce_wormbased_TSS_gx_slopbed_up2000_dn200.nh.bed --name promoter_regions --query cetrn_focused_f0.1 --parameters cetrn_tracked_f0.1 --overwrite OFF


# Collect cells from neurons:
#python mapCells.py --path ~/meTRN --organism ce --mode cell.collection --peaks ce_selection_com_ex_raw --technique binary --neurons cell.ex.som --collection cetrn_cellular_raw


# Check composition of neurons (Re-done for new annotations!):
#python mapCells.py --path ~/meTRN --organism ce --mode test.composition --infile mapcells_avgExp_waterston_expression_tissues --expression mapcells_avgExp_waterston_expression_focused --peaks ce_selection_com_ex_raw --neurons cell.ex.som --name cetrn_cellular_raw



# Examine lineage enrichments in cellular-resolution neurons:
#python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --name waterston.focused --method builder --lineages tracked --descendants OFF --ascendants OFF --limit 100000 --collection cetrn_cellular_raw
#python mapCells.py --path ~/meTRN --organism ce --mode test.comparison --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --query cetrn_focused_f0.1 --target cetrn_cellular_raw --name cetrn_cellular_raw
#python mapCells.py --path ~/meTRN --organism ce --mode test.regions --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --peaks ce_selection_com_ex_raw --neurons cell.ex.som --query cetrn_focused_f0.1 --target cetrn_cellular_raw --name cetrn_cellular_raw --up 2000 --dn 200
#python mapCells.py --path ~/meTRN --organism ce --mode test.fdr --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --peaks ce_selection_com_ex_raw --neurons cell.ex.som --query cetrn_focused_f0.1 --target cetrn_cellular_raw --name cetrn_cellular_raw --up 2000 --dn 200




# UNDER DEVELOPMENT:
#python mapCells.py --path ~/meTRN --organism ce --mode test.similarity --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --peaks ce_selection_com_ex_raw --neurons cell.ex.som --query cetrn_focused_f0.1 --target cetrn_cellular_raw --name cetrn_cellular_raw --reference in2shape_ce_wormbased_TSS_gx.bed


# UNDER THE HOOD:
# Build lineage enrichment trees (.json) for visualization:
#python mapCells.py --path ~/meTRN --organism ce --mode tree.build --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages complete --ascendants OFF --name complete
#python mapCells.py --path ~/meTRN --organism ce --mode tree.build --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --ascendants 5 --name tracked
#python mapCells.py --path ~/meTRN --organism ce --mode tree.color --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --lineages tracked --ascendants 5 --name tracked --infile /Users/claraya/ceTRN/data/cells/compare/tracked_vs_cells/hyper/mapcells_comparison_focused_vs_cells_focused_neurons.txt


#top