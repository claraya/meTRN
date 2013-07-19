#!/usr/bin/sh


# Note:	These files were copied from the ceTRN run! 
#		The commands here were modified so as to record
#		the structure of the runs.
# Note: The focused and tracked runs were performed!


# Build lineage cells as descendants from nodes (descender):
python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method descender --lineages tracked --descendants OFF --ascendants OFF
python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method descender --lineages complete --descendants OFF --ascendants OFF

python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --name waterston.focused --method descender --lineages tracked --descendants OFF --ascendants OFF
python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --name waterston.focused --method descender --lineages complete --descendants OFF --ascendants OFF

python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_assayed --name waterston.assayed --method descender --lineages tracked --descendants OFF --ascendants OFF
python mapCells.py --path ~/meTRN --organism ce --mode build.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_assayed --name waterston.assayed --method descender --lineages complete --descendants OFF --ascendants OFF


# Test lineage cells as descendants from nodes (descender):
python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method descender --lineages tracked --descendants OFF --ascendants OFF --collection cetrn_tracked_f0.1
python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_tracked --name waterston.tracked --method descender --lineages complete --descendants OFF --ascendants OFF --collection cetrn_tracked_f0.1

python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --name waterston.focused --method descender --lineages tracked --descendants OFF --ascendants OFF --collection cetrn_focused_f0.1
python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_focused --name waterston.focused --method descender --lineages complete --descendants OFF --ascendants OFF --collection cetrn_focused_f0.1

python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_assayed --name waterston.assayed --method descender --lineages tracked --descendants OFF --ascendants OFF --collection cetrn_assayed_f0.1
python mapCells.py --path ~/meTRN --organism ce --mode test.lineages --pedigree waterston_cell_pedigree.csv --expression mapcells_avgExp_waterston_expression_assayed --name waterston.assayed --method descender --lineages complete --descendants OFF --ascendants OFF --collection cetrn_assayed_f0.1


#top
#bash runMaster6B.sh